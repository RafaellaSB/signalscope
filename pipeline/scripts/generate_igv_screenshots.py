#!/usr/bin/env python3
"""
generate_igv_screenshots.py — SignalScope IGV screenshot generator

Reads SNV, indel, and SV shortlists, renders one IGV pileup PNG per variant
using a SINGLE IGV batch session (start IGV once, load BAM once, take all
screenshots, exit), and writes a manifest TSV:

    key                              b64_png
    snv:1_12345_A_G                  <base64 data>
    indel:1_12345_ACGT_A             <base64 data>
    sv:16_89623987_DEL               <base64 data>

Key formats match load_igv_manifest() in generate_report_html.py:
    snv:   snv:{chrom}_{pos}_{ref}_{alt}
    indel: indel:{chrom}_{pos}_{ref}_{alt}
    sv:    sv:{chrom}_{start}_{svtype}

RENDERING STRATEGY (in order of preference):
    1. IGV desktop (IGV_Linux via xvfb-run) — SINGLE BATCH SESSION, all
       variants in one IGV invocation. Best quality, strand-coloured reads.
    2. samtools tview + ImageMagick convert — per-variant fallback for any
       variants the IGV batch failed to render.
    3. Placeholder PNG — blank grey image with variant label (always works).

Why single-batch: the BAM may be 9+ GB. Launching a fresh IGV per variant
re-loads/indexes the BAM each time and blows past any reasonable timeout.
One batch loads it once and snapshots cheaply.

USAGE
    python3 scripts/generate_igv_screenshots.py \\
        --snv_shortlist   ... \\
        --indel_shortlist ... \\
        --sv_shortlist    ... \\
        --bam             samples/barcode06/alignments/barcode06.primary.bam \\
        --ref_fasta       /path/to/GRCh38.fa \\
        --out_dir         samples/barcode06/report/igv \\
        --sample          barcode06 \\
        --igv_sh          igv \\
        --flank           200
"""
import argparse
import base64
import csv
import glob
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Generate IGV screenshots for SignalScope report")
    p.add_argument("--snv_shortlist",   required=True)
    p.add_argument("--indel_shortlist", required=True)
    p.add_argument("--sv_shortlist",    required=True)
    p.add_argument("--bam",             required=True)
    p.add_argument("--ref_fasta",       required=True)
    p.add_argument("--igv_genome",      default="",
                   help="IGV genome id or path to a local .json genome definition. "
                        "Defaults to hosted hg38 if empty (requires network).")
    p.add_argument("--out_dir",         required=True)
    p.add_argument("--sample",          required=True)
    p.add_argument("--igv_sh",          default="",
                   help="Path to igv.sh (IGV desktop launcher). Auto-detected if not set.")
    p.add_argument("--flank",           type=int, default=200,
                   help="bp either side of variant centre for IGV window (default 200)")
    p.add_argument("--height",          type=int, default=300,
                   help="PNG height in pixels — used by tview fallback only (default 300)")
    p.add_argument("--width",           type=int, default=900,
                   help="PNG width in pixels — used by tview fallback only (default 900)")
    p.add_argument("--igv_timeout",     type=int, default=0,
                   help="IGV batch timeout in seconds. 0 = auto-scale "
                        "(max(600, 60 + 8 * n_variants)). Override for tuning.")
    p.add_argument("--igv_memory",      default="4g",
                   help="Java heap for IGV (-Xmx flag). Default 4g.")
    p.add_argument("--preflight",       default="warn",
                   choices=["warn", "fail", "clean", "skip"],
                   help="What to do if orphan IGV/Xvfb processes or stale X "
                        "lockfiles are detected before launch. "
                        "'warn' (default) prints a warning and continues. "
                        "'fail' aborts with an error. "
                        "'clean' removes stale lockfiles and kills orphans "
                        "owned by this user, then continues. "
                        "'skip' bypasses the check entirely.")
    return p.parse_args()


# =============================================================================
# I/O HELPERS
# =============================================================================

def opt_file(path):
    """Return path if non-empty file exists, else None."""
    return path if path and os.path.isfile(path) and os.path.getsize(path) > 0 else None


def read_tsv(path):
    """Read TSV with header; return list of row dicts. Returns [] if absent."""
    if not opt_file(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def b64_encode_file(path):
    """Return base64 string of file contents."""
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


# =============================================================================
# VARIANT COLLECTION
# =============================================================================

def collect_variants(snv_path, indel_path, sv_path):
    """
    Returns list of dicts:
        key    — manifest key string
        chrom  — chromosome (no 'chr' prefix, matching GRCh38 ref)
        start  — 1-based genomic start
        end    — 1-based genomic end
        label  — human-readable variant label for alt-text
    """
    variants = []

    # ── SNVs ─────────────────────────────────────────────────────────────────
    for r in read_tsv(snv_path):
        variant_str = r.get("variant", "")
        gene        = r.get("gene", "")
        try:
            parts = variant_str.split("_")
            chrom = parts[0]
            pos   = int(parts[1])
            alleles = parts[2].split("/")
            ref = alleles[0]; alt = alleles[1] if len(alleles) > 1 else "."
            key = f"snv:{chrom}_{pos}_{ref}_{alt}"
            variants.append({
                "key":   key,
                "chrom": chrom,
                "start": pos,
                "end":   pos,
                "label": f"{gene} {variant_str}",
            })
        except (IndexError, ValueError) as e:
            print(f"  [IGV] WARNING: skipping malformed SNV variant '{variant_str}': {e}",
                  file=sys.stderr)

    # ── Indels ────────────────────────────────────────────────────────────────
    for r in read_tsv(indel_path):
        chrom = r.get("chrom", r.get("CHROM", ""))
        gene  = r.get("gene", "")
        try:
            pos  = int(r.get("pos", r.get("POS", 0)))
            ref  = r.get("ref", r.get("REF", "-"))
            alt  = r.get("alt", r.get("ALT", "-"))
            key  = f"indel:{chrom}_{pos}_{ref}_{alt}"
            end  = pos + max(len(ref), len(alt))
            variants.append({
                "key":   key,
                "chrom": chrom,
                "start": pos,
                "end":   end,
                "label": f"{gene} {chrom}:{pos} {ref}>{alt}",
            })
        except (ValueError, TypeError) as e:
            print(f"  [IGV] WARNING: skipping malformed indel row: {e}", file=sys.stderr)

    # ── SVs ───────────────────────────────────────────────────────────────────
    for r in read_tsv(sv_path):
        chrom   = r.get("chrom", "")
        gene    = r.get("gene", "")
        sv_type = r.get("sv_type", "")
        try:
            start = int(r.get("start", 0))
            end   = int(r.get("end", start + 1000))
            key   = f"sv:{chrom}_{start}_{sv_type}"
            variants.append({
                "key":   key,
                "chrom": chrom,
                "start": start,
                "end":   end,
                "label": f"{gene} {sv_type} {chrom}:{start}-{end}",
            })
        except (ValueError, TypeError) as e:
            print(f"  [IGV] WARNING: skipping malformed SV row: {e}", file=sys.stderr)

    print(f"  [IGV] Collected {len(variants)} variants for screenshots", file=sys.stderr)
    return variants


# =============================================================================
# BACKEND DETECTION
# =============================================================================

def _find_igv_sh(igv_sh_arg):
    """Locate the IGV launcher.

    Resolution order:
      1. The explicit --igv_sh argument, if it is an existing file or a
         command resolvable on PATH.
      2. An executable named "igv" or "igv.sh" on PATH (installed via the
         conda environment).
    Returns the path/command to use, or None if IGV cannot be found.
    """
    if igv_sh_arg:
        if os.path.isfile(igv_sh_arg):
            return igv_sh_arg
        found = shutil.which(igv_sh_arg)
        if found:
            return found
    for name in ("igv", "igv.sh"):
        found = shutil.which(name)
        if found:
            return found
    return None


def _check_xvfb():
    """Return True if xvfb-run is available."""
    return shutil.which("xvfb-run") is not None


def _check_imagemagick():
    """Return True if ImageMagick 'convert' command is available."""
    return shutil.which("convert") is not None


def _check_samtools():
    """Return True if samtools is available."""
    return shutil.which("samtools") is not None


def _preflight_environment():
    """
    Refuse to launch if the environment is dirty.

    A previous IGV run that didn't clean up its process group leaves:
      - orphan java/igv.sh/Xvfb processes (will compete for X displays)
      - stale /tmp/.X*-lock files (xvfb-run wastes time walking past them)

    Both correlate with the hangs we've seen on this box. Better to refuse
    fast with an actionable message than to launch into a flaky run.
    """
    issues = []

    # Stale X lockfiles from dead servers.
    locks = sorted(glob.glob("/tmp/.X*-lock"))
    stale_locks = []
    for lock in locks:
        try:
            with open(lock) as fh:
                pid = int(fh.read().strip())
        except (OSError, ValueError):
            continue
        # If the PID isn't a live process, the lock is stale.
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            stale_locks.append(lock)
        except PermissionError:
            pass   # exists, owned by another user — leave it alone
    if stale_locks:
        issues.append(
            f"{len(stale_locks)} stale X lockfile(s): "
            + " ".join(stale_locks[:5])
            + (f" ... ({len(stale_locks)-5} more)" if len(stale_locks) > 5 else "")
        )

    # Orphan IGV / Xvfb processes from this user.
    #
    # Approach: list all matching processes via pgrep with a deliberately
    # liberal regex, then filter by:
    #   (a) excluding our own pid and parent pid (in case our command line
    #       contains tokens like '/path/to/igv.sh' as an argument);
    #   (b) keeping only entries whose *executable* (first whitespace token
    #       of the command line) looks like a real IGV/Xvfb invocation:
    #         - basename is 'Xvfb'
    #         - basename is 'igv.sh' or 'sh' (when igv.sh runs as /bin/sh igv.sh)
    #         - the cmdline contains 'org.broad.igv' (the Java entry point)
    # This double check stops the preflight matching `python3 … --igv_sh
    # /path/to/igv.sh …`, which would otherwise be a false positive.
    my_pid = os.getpid()
    my_ppid = os.getppid()
    try:
        out = subprocess.check_output(
            ["pgrep", "-u", str(os.getuid()), "-af",
             r"Xvfb|igv\.sh|org\.broad\.igv"],
            text=True,
        )
        all_matches = [l for l in out.strip().split("\n") if l]
    except subprocess.CalledProcessError:
        all_matches = []   # pgrep returns 1 when no matches
    except FileNotFoundError:
        all_matches = []   # no pgrep — skip this check

    orphan_lines = []
    for line in all_matches:
        # pgrep -af output: "<pid> <command line>"
        parts = line.split(None, 1)
        if len(parts) < 2:
            continue
        try:
            pid = int(parts[0])
        except ValueError:
            continue
        if pid in (my_pid, my_ppid):
            continue
        cmdline = parts[1]
        argv = cmdline.split()
        if not argv:
            continue
        # First token is the executable (or interpreter, e.g. /bin/sh for
        # shell scripts). Classify the process by looking at its argv:
        first_base = os.path.basename(argv[0])
        # Real orphan if any of:
        #   - the binary itself is Xvfb
        #   - the binary is a shell whose script argument basename is igv.sh
        #   - the cmdline contains the org.broad.igv java entry point
        is_xvfb = first_base == "Xvfb"
        is_igv_sh = first_base == "igv.sh"
        if not is_igv_sh and first_base in ("sh", "bash") and len(argv) > 1:
            is_igv_sh = os.path.basename(argv[1]) == "igv.sh"
        is_igv_java = "org.broad.igv" in cmdline
        if is_xvfb or is_igv_sh or is_igv_java:
            orphan_lines.append(line)

    if orphan_lines:
        issues.append(
            f"{len(orphan_lines)} orphan IGV/Xvfb process(es) still running:"
        )
        for line in orphan_lines[:5]:
            issues.append(f"    {line}")
        if len(orphan_lines) > 5:
            issues.append(f"    ... ({len(orphan_lines)-5} more)")

    return issues, stale_locks


# =============================================================================
# OUTPUT PATH HELPER
# =============================================================================

def _out_png_path(out_dir, sample, key):
    """Map a variant key to its on-disk PNG path."""
    type_prefix = key.split(":")[0]   # snv / indel / sv
    coords      = key.split(":")[1]   # e.g. 19_35850409_T_A
    # Sanitize: replace filesystem-unsafe characters (/ in multi-allelic notation, etc.)
    coords_safe = coords.replace("/", "-").replace("\\", "-").replace(":", "-")
    return os.path.join(out_dir, f"{sample}.{type_prefix}.{coords_safe}.igv.png")


# =============================================================================
# BACKEND 1: IGV DESKTOP — SINGLE BATCH SESSION
# =============================================================================

def render_igv_batch(variants, bam, ref_fasta, out_dir, sample,
                     igv_sh, flank, timeout, igv_memory, igv_genome="hg38"):
    """
    Launch IGV ONCE with a batch script that loads the BAM and snapshots
    every variant. Returns the set of variant keys that produced a valid
    PNG (> 5000 bytes).

    Batch script structure (validated against gold-standard images):
        new
        genome hg38
        load <recovered.bam>
        snapshotDirectory <out_dir>
        maxPanelHeight 500
        # ── per variant ──
        goto <chrom>:<centre±50bp for SNV/indel, ±200bp for SV>
        colorBy READ_STRAND
        shadeBases QUALITY
        sort BASE <chrom>:<centre>
        expanded
        snapshot <filename.png>
        # ── end per variant ──
        exit
    """
    if not _check_xvfb():
        print("  [IGV] xvfb-run not available; skipping IGV desktop backend", file=sys.stderr)
        return set()

    if not variants:
        return set()

    out_dir_abs = os.path.abspath(out_dir)
    bam_abs     = os.path.abspath(bam)

    # Per-type flank in bp around the variant centre:
    #   SNV / indel — narrow (~50bp) so the Sequence track shows
    #                  reference base letters at the variant position.
    #   SV          — wider (~200bp) so the structural event is visible.
    # The --flank CLI arg sets the SV flank (kept for backwards compatibility);
    # SNV/indel always uses the narrow value.
    sv_flank_bp  = flank
    snv_flank_bp = 50

    # Build the batch file. One session, one BAM load, N snapshots.
    # Recipe matches the validated gold-standard IGV screenshots:
    #   genome hg38         — IGV's hosted hg38 (brings RefSeq gene annotation
    #                         and the Sequence track that's missing if we load
    #                         a bare local FASTA).
    #   colorBy READ_STRAND — strand-colour reads (the enum is READ_STRAND,
    #                         not STRAND; STRAND throws IllegalArgumentException).
    #   shadeBases QUALITY  — fade low-quality bases so high-confidence calls
    #                         stand out.
    #   sort BASE {pos}     — group reads carrying the variant allele to the
    #                         top of the pile.
    #   expanded            — give each read its own row; squish/collapse hides
    #                         the read separation needed to count alleles.
    #   maxPanelHeight 500  — matches the panel height used for the gold images.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".batch",
                                     delete=False, dir=out_dir_abs) as fh:
        batch_file = fh.name
        fh.write("new\n")
        fh.write(f"genome {igv_genome}\n")
        fh.write(f"load {bam_abs}\n")
        fh.write(f"snapshotDirectory {out_dir_abs}\n")
        fh.write("maxPanelHeight 500\n")

        for v in variants:
            chrom    = v["chrom"]
            centre   = v["start"] + max(0, (v["end"] - v["start"]) // 2)
            var_type = v["key"].split(":", 1)[0]   # snv / indel / sv
            this_flank = sv_flank_bp if var_type == "sv" else snv_flank_bp
            region   = f"{chrom}:{max(1, centre - this_flank)}-{centre + this_flank}"
            png_name = os.path.basename(_out_png_path(out_dir, sample, v["key"]))

            fh.write(f"goto {region}\n")
            fh.write("colorBy READ_STRAND\n")
            fh.write("shadeBases QUALITY\n")
            fh.write(f"sort BASE {chrom}:{centre}\n")
            fh.write("expanded\n")
            fh.write(f"snapshot {png_name}\n")

        fh.write("exit\n")

    print(f"  [IGV] launching single batch session for {len(variants)} variants "
          f"(timeout={timeout}s, heap={igv_memory})", file=sys.stderr)
    print(f"  [IGV] batch file: {batch_file}", file=sys.stderr)

    # Pass -Xmx to IGV via the env var that igv.sh respects.
    env = os.environ.copy()
    # Seed a private IGV HOME with the window-bounds pref so headless snapshots
    # render off-screen on any host (a fresh ~/igv without this hangs the batch).
    # Bounds 1150x800 fix the snapshot canvas to the validated dimensions.
    _igv_home = os.path.join(os.path.abspath(out_dir), ".igv_home")
    os.makedirs(os.path.join(_igv_home, "igv"), exist_ok=True)
    with open(os.path.join(_igv_home, "igv", "prefs.properties"), "w") as _pf:
        _pf.write("IGV.Bounds=0,0,1150,800\n")
    env["HOME"] = _igv_home
    # IGV_Linux's igv.sh reads IGV_JAVA_ARGS for extra JVM flags.
    existing = env.get("IGV_JAVA_ARGS", "")
    env["IGV_JAVA_ARGS"] = f"-Xmx{igv_memory} {existing}".strip()

    # Expected output paths — used to detect "all done" without waiting for
    # IGV to exit cleanly (it often doesn't on Linux/xvfb).
    expected_paths = [
        _out_png_path(out_dir, sample, v["key"]) for v in variants
    ]

    # Snapshot pre-existing file mtimes so we only count *newly written* PNGs.
    # (Re-runs leave old PNGs from previous attempts in place.)
    pre_existing_mtimes = {}
    for p in expected_paths:
        if os.path.isfile(p):
            try:
                pre_existing_mtimes[p] = os.path.getmtime(p)
            except OSError:
                pass

    def count_fresh_snapshots():
        """Return number of expected PNGs that look freshly written."""
        n = 0
        for p in expected_paths:
            if not os.path.isfile(p) or os.path.getsize(p) <= 5000:
                continue
            mt = os.path.getmtime(p)
            if p in pre_existing_mtimes and mt <= pre_existing_mtimes[p]:
                continue   # stale leftover from a previous run
            n += 1
        return n

    def kill_process_group(p, reason):
        """
        Kill the entire process group rooted at p.

        This is critical: xvfb-run is a shell wrapper that fork-execs IGV
        which fork-execs Java. Killing only the parent leaves Java alive,
        holding the X display lock. Putting the launch into its own session
        (start_new_session=True) lets us SIGTERM/SIGKILL the whole group
        with os.killpg.
        """
        if p.poll() is not None:
            return   # already gone
        try:
            pgid = os.getpgid(p.pid)
        except ProcessLookupError:
            return
        print(f"  [IGV] {reason} — killing process group {pgid}",
              file=sys.stderr)
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            return
        # Give it up to 8s to die cleanly, then SIGKILL the survivors.
        for _ in range(16):
            if p.poll() is not None:
                return
            time.sleep(0.5)
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            return
        try:
            p.wait(timeout=5)
        except subprocess.TimeoutExpired:
            pass

    # start_new_session=True puts xvfb-run + igv.sh + java in a fresh
    # process group whose pgid == proc.pid. That's what makes os.killpg
    # below kill the whole tree, not just the shell wrapper.
    proc = subprocess.Popen(
        ["xvfb-run", "--auto-servernum",
         igv_sh, "-b", batch_file],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, env=env, start_new_session=True,
    )

    poll_interval = 2.0          # seconds between checks
    settle_seconds = 3.0         # wait after the last PNG appears, so the
                                 # final snapshot has time to fully flush
    deadline = time.monotonic() + timeout
    all_done_at = None
    last_progress_n = -1
    last_progress_t = time.monotonic()

    try:
        while True:
            # IGV exited on its own — great, we're done.
            if proc.poll() is not None:
                break

            n_done = count_fresh_snapshots()

            # Per-snapshot progress logging (matches the old log format
            # so behaviour is easy to compare against the gold-standard run).
            if n_done > last_progress_n:
                for v in variants[last_progress_n + 1 if last_progress_n >= 0 else 0 : n_done]:
                    path = _out_png_path(out_dir, sample, v["key"])
                    print(f"  [IGV] OK: {os.path.basename(path)}",
                          file=sys.stderr)
                last_progress_n = n_done
                last_progress_t = time.monotonic()

            if n_done >= len(expected_paths):
                # All PNGs are on disk. Give IGV a brief moment to flush the
                # last write, then kill it ourselves rather than waiting on
                # `exit` which often hangs under xvfb.
                if all_done_at is None:
                    all_done_at = time.monotonic()
                    print(f"  [IGV] all {n_done} snapshots written — "
                          f"shutting IGV down", file=sys.stderr)
                elif time.monotonic() - all_done_at >= settle_seconds:
                    kill_process_group(proc, "clean shutdown after success")
                    break

            if time.monotonic() > deadline:
                kill_process_group(
                    proc,
                    f"IGV batch timed out after {timeout}s "
                    f"({n_done}/{len(expected_paths)} snapshots written)",
                )
                break

            time.sleep(poll_interval)

        # Drain remaining stderr for diagnostics if IGV returned non-zero.
        # Process group is already dead by this point.
        try:
            _, stderr_tail = proc.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            kill_process_group(proc, "stderr drain timeout")
            try:
                _, stderr_tail = proc.communicate(timeout=2)
            except subprocess.TimeoutExpired:
                stderr_tail = ""
        rc = proc.returncode
        if rc not in (0, None) and rc != -15 and rc != -9:
            # -15 = SIGTERM (us), -9 = SIGKILL (us) — both expected when we
            # shut IGV down after success. Only complain about other rcs.
            print(f"  [IGV] IGV batch returned rc={rc}", file=sys.stderr)
            if stderr_tail:
                print(f"  [IGV] stderr tail: {stderr_tail[-400:]}", file=sys.stderr)
    except Exception as e:
        print(f"  [IGV] IGV batch raised: {e}", file=sys.stderr)
        kill_process_group(proc, "exception in poll loop")
    finally:
        # Keep the batch file around on failure for debugging; remove on success.
        # Comment out the unlink if you want to inspect it.
        try:
            os.unlink(batch_file)
        except OSError:
            pass

    # Audit which PNGs actually landed and are non-trivial.
    succeeded = set()
    for v in variants:
        path = _out_png_path(out_dir, sample, v["key"])
        if os.path.isfile(path) and os.path.getsize(path) > 5000:
            succeeded.add(v["key"])

    print(f"  [IGV] batch produced {len(succeeded)}/{len(variants)} valid screenshots",
          file=sys.stderr)
    return succeeded


# =============================================================================
# BACKEND 2: samtools tview fallback (per-variant)
# =============================================================================

def render_samtools_tview(variant, bam, ref_fasta, out_png, flank, height, width):
    """
    Per-variant fallback when IGV didn't produce a usable PNG for this variant.
    Not strand-coloured but functional.
    """
    if not (_check_samtools() and _check_imagemagick()):
        return False

    chrom  = variant["chrom"]
    centre = variant["start"] + max(0, (variant["end"] - variant["start"]) // 2)
    region = f"{chrom}:{max(1, centre - flank)}-{centre + flank}"

    with tempfile.TemporaryDirectory() as tmpdir:
        txt_out = os.path.join(tmpdir, "tview.txt")
        try:
            result = subprocess.run(
                ["samtools", "tview", "-d", "T", "-p", region, bam, ref_fasta],
                capture_output=True, text=True, timeout=30
            )
            if result.returncode != 0:
                print(f"  [IGV] samtools tview failed for {region}: {result.stderr[:200]}",
                      file=sys.stderr)
                return False
            with open(txt_out, "w") as fh:
                fh.write(result.stdout)
        except subprocess.TimeoutExpired:
            print(f"  [IGV] samtools tview timed out for {region}", file=sys.stderr)
            return False

        # Convert text → PNG with ImageMagick.
        # Strip coordinate ruler line (pure digits) and cap at 15 lines
        # to stay within ImageMagick canvas policy limits.
        try:
            with open(txt_out, "r") as fh:
                txt_content = fh.read()
            lines_tv = [l for l in txt_content.splitlines()
                        if not l.strip().replace(" ", "").isdigit()]
            txt_content = "\n".join(lines_tv[:15])
            subprocess.run([
                "convert",
                "-size", f"{width}x{height}",
                "-background", "#1a1a2e",
                "-fill", "#c5d8f0",
                "-font", "Courier",
                "-pointsize", "11",
                f"caption:{txt_content}",
                "-trim",
                out_png
            ], check=True, capture_output=True, timeout=30, text=True)
            return os.path.isfile(out_png) and os.path.getsize(out_png) > 0
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            print(f"  [IGV] ImageMagick convert failed: {e}", file=sys.stderr)
            return False


# =============================================================================
# BACKEND 3: placeholder PNG (always succeeds)
# =============================================================================

def render_placeholder(variant, out_png, width, height):
    """
    Minimal placeholder PNG with variant label.
    Always works — no external tools required.
    """
    label = variant.get("label", variant.get("key", "unknown"))

    try:
        from PIL import Image, ImageDraw
        img  = Image.new("RGB", (width, height), color="#1a1a2e")
        draw = ImageDraw.Draw(img)
        draw.text((10, height // 2 - 10), f"No IGV screenshot\n{label}", fill="#c5d8f0")
        img.save(out_png, "PNG")
        return True
    except ImportError:
        pass

    # Pure stdlib fallback — minimal valid PNG
    import struct, zlib

    def png_chunk(chunk_type, data):
        c = chunk_type + data
        return struct.pack(">I", len(data)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)

    w = h = 1
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    raw  = b"\x00\x80\x80\x80"
    idat = zlib.compress(raw)
    png_data = (
        b"\x89PNG\r\n\x1a\n"
        + png_chunk(b"IHDR", ihdr)
        + png_chunk(b"IDAT", idat)
        + png_chunk(b"IEND", b"")
    )
    with open(out_png, "wb") as fh:
        fh.write(png_data)
    return True


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    manifest_path = os.path.join(args.out_dir, f"{args.sample}.igv_manifest.tsv")

    # Detect available backends.
    igv_sh          = _find_igv_sh(args.igv_sh)
    has_xvfb        = _check_xvfb()
    has_igv_desktop = igv_sh is not None and has_xvfb
    has_imagemagick = _check_imagemagick()
    has_samtools    = _check_samtools()

    print(f"[generate_igv_screenshots] sample={args.sample}", file=sys.stderr)
    print(f"  backends: igv-desktop={has_igv_desktop}  "
          f"samtools+imagemagick={has_samtools and has_imagemagick}", file=sys.stderr)
    if has_igv_desktop:
        print(f"  igv_sh:   {igv_sh}", file=sys.stderr)
    print(f"  bam:      {args.bam}", file=sys.stderr)

    variants = collect_variants(args.snv_shortlist, args.indel_shortlist, args.sv_shortlist)

    # ── Pre-flight environment check ─────────────────────────────────────────
    # Stale X locks and orphan IGV/Xvfb processes correlate with hangs.
    # Detect early; behaviour controlled by --preflight.
    if args.preflight != "skip" and has_igv_desktop and variants:
        issues, stale_locks = _preflight_environment()
        if issues:
            for msg in issues:
                print(f"  [preflight] {msg}", file=sys.stderr)
            if args.preflight == "fail":
                print("  [preflight] aborting (--preflight=fail). "
                      "Clean up with: "
                      "pkill -KILL -f 'org.broad.igv|igv.sh|Xvfb' && "
                      "rm -f /tmp/.X*-lock",
                      file=sys.stderr)
                sys.exit(2)
            elif args.preflight == "clean":
                # Kill our own orphans, then remove stale locks we proved are stale.
                # Same tight regex as the detection step, for the same reason
                # (don't match the current python invocation's own args).
                pkill_pattern = r"^Xvfb\b|/igv\.sh\b|org\.broad\.igv\b"
                for sig in (signal.SIGTERM, signal.SIGKILL):
                    subprocess.run(
                        ["pkill", f"-{sig.value}", "-u", str(os.getuid()),
                         "-f", pkill_pattern],
                        check=False,
                    )
                    time.sleep(2)
                for lock in stale_locks:
                    try:
                        os.unlink(lock)
                    except OSError:
                        pass
                print(f"  [preflight] cleaned {len(stale_locks)} stale lock(s) "
                      f"and killed orphan processes", file=sys.stderr)
            # 'warn' falls through.

    # Compute the IGV batch timeout.
    if args.igv_timeout > 0:
        igv_timeout = args.igv_timeout
    else:
        # 60s for IGV startup + BAM index load, plus ~8s per snapshot, floored at 600s.
        igv_timeout = max(600, 60 + 8 * len(variants))

    # ── Backend 1: one IGV batch for all variants ────────────────────────────
    igv_succeeded = set()
    if has_igv_desktop and variants:
        igv_succeeded = render_igv_batch(
            variants, args.bam, args.ref_fasta, args.out_dir, args.sample,
            igv_sh, args.flank, igv_timeout, args.igv_memory,
            igv_genome=(args.igv_genome or "hg38"),
        )

    # ── Per-variant fallback pass for anything IGV missed ────────────────────
    rows = []
    for i, v in enumerate(variants):
        key     = v["key"]
        out_png = _out_png_path(args.out_dir, args.sample, key)

        print(f"  [{i+1}/{len(variants)}] {key} ...", file=sys.stderr)

        if key in igv_succeeded:
            print(f"    rendered via IGV desktop (batch)", file=sys.stderr)
            rendered = True
        else:
            rendered = False
            # Backend 2: samtools tview + ImageMagick
            if has_samtools and has_imagemagick:
                rendered = render_samtools_tview(v, args.bam, args.ref_fasta,
                                                 out_png, args.flank,
                                                 args.height, args.width)
                if rendered:
                    print(f"    rendered via samtools tview (fallback)", file=sys.stderr)

            # Backend 3: placeholder PNG
            if not rendered:
                rendered = render_placeholder(v, out_png, args.width, args.height)
                print(f"    rendered placeholder (no rendering backend available)",
                      file=sys.stderr)

        if rendered and os.path.isfile(out_png):
            b64 = b64_encode_file(out_png)
            rows.append({"key": key, "b64_png": b64})
        else:
            print(f"  [IGV] WARNING: failed to generate screenshot for {key}",
                  file=sys.stderr)
            rows.append({"key": key, "b64_png": ""})

    # Write manifest TSV
    with open(manifest_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["key", "b64_png"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # Summary
    n_igv  = len(igv_succeeded)
    n_total = len(variants)
    print(f"  [IGV] summary: {n_igv}/{n_total} via IGV desktop, "
          f"{n_total - n_igv} via fallback/placeholder", file=sys.stderr)
    print(f"  Written: {manifest_path} ({len(rows)} variants)", file=sys.stderr)


if __name__ == "__main__":
    main()
