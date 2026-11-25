#!/usr/bin/env python3
"""
check_model_match.py - does the Clair3 model match the basecaller model that made the reads?

Reads the basecall_model tag(s) from the BAM header (@RG DS:basecall_model=...) and compares
them with the Clair3 model name. Compared, when both names carry them:
  accuracy  fast / hac / sup
  version   v4.3.0 in the reads -> v430 in the model name
  chemistry r10.4.1 -> r1041, e8.2 -> e82, 400bps
A model name that does not carry a field (for example a folder renamed without the chemistry)
is not compared on that field; the check says so.

Exit 0 = match (or nothing to compare and --on_unknown pass); 1 = mismatch; writes a report.
With --model auto the model is CHOSEN: among the installed models (--model_dir, names starting with
--prefix) the one whose accuracy, version and chemistry equal the basecaller model of the reads;
the report then has a line 'model: <name>' that the Clair3 rule reads. When the BAM header has no
basecall_model tag, --basecall_model (from the config) is used and the report says so; with
neither, or with reads from several basecaller models, the run stops instead of guessing.
  python3 check_model_match.py --bam x.bam --model clair3_pytorch_hac_v430 --out check.txt \
      [--samtools samtools] [--on_mismatch stop|warn] [--on_unknown pass|stop]
"""
import argparse
import os
import re
import subprocess
import sys

ACC = ("fast", "hac", "sup")
R10_DEFAULT = "r1041_e82_400bps"


def read_fields(name):
    """{'acc', 'ver', 'chem'} found in a basecaller or Clair3 model name (lower case)."""
    n = name.lower()
    out = {}
    for a in ACC:
        if re.search(r"(^|[_@\-])%s([_@\-]|$)" % a, n):
            out["acc"] = a
            break
    m = re.search(r"@v(\d+)\.(\d+)\.(\d+)", n)              # basecaller: ..._hac@v4.3.0
    if m:
        out["ver"] = "".join(m.groups())
    else:
        m = re.findall(r"(?:^|_)v(\d{3})(?:_|$)", n)          # Clair3: ..._hac_v430
        if m:
            out["ver"] = m[-1]
    chem = []
    m = re.search(r"(?:^|_)r(10|9)\.?(\d)\.?(\d)", n)        # r10.4.1 / r1041 / r9.4.1 / r941
    if m:
        chem.append("r" + "".join(m.groups()))
    if re.search(r"e8\.?2", n):
        chem.append("e82")
    m = re.search(r"(\d{3})bps", n)
    if m:
        chem.append(m.group(1) + "bps")
    if chem:
        out["chem"] = "_".join(chem)
    return out


def basecall_models(header_text):
    return sorted(set(re.findall(r"basecall_model=([^\s\t]+)", header_text)))


def compare(read_models, clair3_model):
    """(status, lines); status 'match', 'mismatch' or 'unknown'."""
    lines = ["Clair3 model: %s" % clair3_model]
    if not read_models:
        return "unknown", lines + ["BAM header has no basecall_model tag: the match cannot be checked"]
    mf = read_fields(clair3_model)
    lines.append("Clair3 model fields: %s" % (", ".join("%s=%s" % kv for kv in sorted(mf.items())) or "none"))
    status = "match"
    for rm in read_models:
        rf = read_fields(rm)
        lines.append("reads basecalled with %s (%s)" % (rm, ", ".join("%s=%s" % kv for kv in sorted(rf.items()))))
        for k in ("acc", "ver", "chem"):
            if k not in mf:
                lines.append("  %s: not in the Clair3 model name, not compared" % k)
            elif k not in rf:
                lines.append("  %s: not in the basecaller model name, not compared" % k)
            elif mf[k] != rf[k]:
                lines.append("  %s: MISMATCH (reads %s, Clair3 model %s)" % (k, rf[k], mf[k]))
                status = "mismatch"
            else:
                lines.append("  %s: match (%s)" % (k, rf[k]))
    if len(read_models) > 1:
        lines.append("note: the BAM holds reads from %d basecaller models" % len(read_models))
    return status, lines


def select_model(read_models, available, prefix="clair3_pytorch_"):
    """Choose the Clair3 model for the basecaller model of the reads.
    Candidates: names in `available` starting with `prefix` whose accuracy and version equal those
    of the reads (and chemistry, when the model name carries it). Returns (name or None, lines)."""
    lines = []
    if not read_models:
        return None, ["no basecaller model known: the BAM has no basecall_model tag and none is given in the "
                      "config (clair3.basecall_model); set clair3.model or clair3.basecall_model"]
    if len(read_models) > 1:
        return None, ["the reads come from %d basecaller models (%s): one Clair3 model cannot fit all; "
                      "set clair3.model explicitly" % (len(read_models), ", ".join(read_models))]
    rf = read_fields(read_models[0])
    if "acc" not in rf or "ver" not in rf:
        return None, ["cannot read accuracy and version from '%s'" % read_models[0]]
    cands = []
    for m in sorted(available):
        if not m.startswith(prefix):
            continue
        mf = read_fields(m)
        if mf.get("acc") != rf["acc"] or mf.get("ver") != rf["ver"]:
            continue
        # a model name without chemistry is taken as R10.4.1 (the chemistry this pipeline and its
        # DeepVariant ONT_R104 model are built for); it is never chosen for other chemistries
        if "chem" in rf and mf.get("chem", R10_DEFAULT) != rf["chem"]:
            continue
        cands.append(m)
    lines.append("reads basecalled with %s (acc=%s, ver=%s%s)" % (read_models[0], rf["acc"], rf["ver"],
                 ", chem=%s" % rf["chem"] if "chem" in rf else ""))
    lines.append("Clair3 models available with prefix '%s': %s" % (
        prefix, ", ".join(m for m in sorted(available) if m.startswith(prefix)) or "none"))
    if not cands:
        lines.append("no installed Clair3 model matches; install the %s v%s model for this chemistry in "
                     "clair3.model_dir, or set clair3.model" % (rf["acc"], rf["ver"]))
        return None, lines
    if len(cands) > 1:
        lines.append("several models match (%s): the first is used" % ", ".join(cands))
    return cands[0], lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--model", required=True, help="Clair3 model name or folder, or 'auto'")
    ap.add_argument("--model_dir", default="", help="folder of installed Clair3 models (for --model auto)")
    ap.add_argument("--prefix", default="clair3_pytorch_", help="model names usable by this Clair3 version")
    ap.add_argument("--basecall_model", default="",
                    help="basecaller model given by the user, used when the BAM header has none")
    ap.add_argument("--out", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--on_mismatch", choices=["stop", "warn"], default="stop")
    ap.add_argument("--on_unknown", choices=["pass", "stop"], default="pass")
    a = ap.parse_args()
    hdr = subprocess.run([a.samtools, "view", "-H", a.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         universal_newlines=True, check=True).stdout
    reads = basecall_models(hdr)
    source = "BAM header"
    if not reads and a.basecall_model:
        reads, source = [a.basecall_model], "config clair3.basecall_model (not in the BAM header)"
    if a.model == "auto":
        avail = [d for d in os.listdir(a.model_dir)] if a.model_dir and os.path.isdir(a.model_dir) else []
        chosen, lines = select_model(reads, avail, a.prefix)
        lines.insert(0, "basecaller model from: %s" % (source if reads else "nowhere (no tag in the BAM header, none in the config)"))
        if chosen is None:
            lines.insert(0, "model check: NO MODEL SELECTED")
            out_lines = lines
        else:
            out_lines = ["model check: OK (selected automatically)", "model: %s" % chosen] + lines
        with open(a.out, "w") as fh:
            fh.write("\n".join(out_lines) + "\n")
        print("\n".join(out_lines))
        if chosen is None:
            print("[STOP] no Clair3 model could be selected (see %s)" % a.out, file=sys.stderr)
            sys.exit(1)
        return
    status, lines = compare(reads, a.model.rstrip("/").split("/")[-1])
    verdict = {"match": "OK", "unknown": "UNKNOWN", "mismatch": "MISMATCH"}[status]
    lines = ["model check: %s" % verdict, "model: %s" % a.model.rstrip("/").split("/")[-1],
             "basecaller model from: %s" % (source if reads else "nowhere (no tag, none in the config)")] + lines
    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if (status == "mismatch" and a.on_mismatch == "stop") or (status == "unknown" and a.on_unknown == "stop"):
        print("[STOP] the Clair3 model does not match the basecaller model (set clair3.model: auto or the right "
              "model, or clair3.model_check: warn to continue anyway)", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
