#!/usr/bin/env python3
"""
build_sequencing_performance.py

Builds {sample}.sequencing_performance.tsv by combining:
  - Whole-dataset metrics from build_run_summary (run_summary.tsv)
  - Targeted-region metrics computed from NanoPlot on the BED-filtered BAM

Output TSV columns:
    metric  whole_dataset  targeted_regions

Usage:
    python3 build_sequencing_performance.py \
        --run_summary    samples/{sample}/report/data/{sample}.run_summary.tsv \
        --targeted_bam   samples/{sample}/qc/targeted/{sample}.targeted.bam \
        --targeted_nanostats samples/{sample}/qc/targeted/nanoplot/{sample}_NanoStats.txt \
        --out_tsv        samples/{sample}/report/data/{sample}.sequencing_performance.tsv \
        --sample         barcode06
"""

import argparse
import gzip
import os
import re
import subprocess
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run_summary",        required=True,
                   help="Whole-run {sample}.run_summary.tsv from build_run_summary")
    p.add_argument("--targeted_bam",       required=True,
                   help="BED-filtered primary BAM")
    p.add_argument("--targeted_nanostats", required=True,
                   help="NanoStats.txt from NanoPlot run on targeted BAM")
    p.add_argument("--out_tsv",            required=True)
    p.add_argument("--sample",             required=True)
    p.add_argument("--seq_summary",        default="",
                   help="MinKNOW sequencing_summary file (for duration extraction)")
    p.add_argument("--barcode",            default="",
                   help="Barcode name to filter seq_summary (default: same as sample)")
    return p.parse_args()


def read_run_summary(path):
    """Read {sample}.run_summary.tsv → dict of metric: value."""
    d = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.strip().split("\t", 1)
            if len(parts) == 2:
                d[parts[0].strip()] = parts[1].strip()
    return d


def count_reads_and_yield(bam_path, samtools="samtools"):
    """
    Count reads and total bases in a BAM using samtools flagstat + idxstats.
    Returns (n_reads, total_bases).
    """
    # Read count from flagstat (total reads, including unmapped)
    try:
        result = subprocess.run(
            [samtools, "flagstat", bam_path],
            capture_output=True, text=True, check=True
        )
        n_reads = 0
        for line in result.stdout.splitlines():
            m = re.match(r"^(\d+) \+ \d+ in total", line)
            if m:
                n_reads = int(m.group(1))
                break
    except subprocess.CalledProcessError as e:
        print(f"  [ERROR] samtools flagstat failed: {e}", file=sys.stderr)
        sys.exit(1)

    # Total bases from samtools stats
    try:
        result = subprocess.run(
            [samtools, "stats", bam_path],
            capture_output=True, text=True, check=True
        )
        total_bases = 0
        for line in result.stdout.splitlines():
            if line.startswith("SN\ttotal length:"):
                total_bases = int(line.split("\t")[2])
                break
    except subprocess.CalledProcessError as e:
        print(f"  [ERROR] samtools stats failed: {e}", file=sys.stderr)
        sys.exit(1)

    return n_reads, total_bases


def parse_nanostats(path):
    """
    Parse NanoStats.txt for mean read length, N50, mean quality.
    Fails hard if any are missing.
    """
    patterns = {
        "mean_length":  re.compile(r"Mean\s+read\s+length[:\s]+([\d,.]+)", re.IGNORECASE),
        "n50":          re.compile(r"Read(?:\s+length)?\s+N50[:\s]+([\d,.]+)", re.IGNORECASE),
        "mean_quality": re.compile(r"Mean\s+read\s+quality[:\s]+([\d.]+)", re.IGNORECASE),
    }
    results = {}
    with open(path) as fh:
        for line in fh:
            for key, pat in list(patterns.items()):
                if key in results:
                    continue
                m = pat.search(line)
                if m:
                    results[key] = float(m.group(1).replace(",", ""))

    missing = [k for k in patterns if k not in results]
    if missing:
        print(f"  [ERROR] NanoStats missing: {', '.join(missing)}", file=sys.stderr)
        print(f"  File: {path}", file=sys.stderr)
        sys.exit(1)

    return results



def compute_duration(seq_summary_path, barcode):
    """Extract sequencing duration in hours from MinKNOW sequencing_summary."""
    if not seq_summary_path or not os.path.isfile(seq_summary_path):
        return None
    try:
        with open(seq_summary_path) as fh:
            header = fh.readline().strip().split("\t")
        start_idx = header.index("start_time")
        dur_idx   = header.index("duration")
        bc_idx    = header.index("barcode_arrangement")

        max_end = 0.0
        min_start = float("inf")
        count = 0
        with open(seq_summary_path) as fh:
            fh.readline()  # skip header
            for line in fh:
                cols = line.strip().split("\t")
                if cols[bc_idx] != barcode:
                    continue
                st = float(cols[start_idx])
                du = float(cols[dur_idx])
                if st < min_start:
                    min_start = st
                end = st + du
                if end > max_end:
                    max_end = end
                count += 1
        if count == 0 or min_start == float("inf"):
            return None
        duration_hours = round((max_end - min_start) / 3600, 1)
        print(f"  Duration: {duration_hours} hours ({count} reads for {barcode})", file=sys.stderr)
        return duration_hours
    except Exception as e:
        print(f"  [WARN] Could not compute duration: {e}", file=sys.stderr)
        return None

def main():
    args = parse_args()
    print(f"[build_sequencing_performance] sample={args.sample}", file=sys.stderr)

    # --- Whole-dataset metrics from run_summary.tsv --------------------------
    whole = read_run_summary(args.run_summary)
    w_reads   = whole.get("Number_of_reads",   "")
    w_yield   = whole.get("Total_yield_Gb",     "")
    w_length  = whole.get("Mean_read_length_bp","")
    w_n50     = whole.get("Read_N50_bp",         "")
    w_quality = whole.get("Mean_read_quality",  "")
    print(f"  Whole-run reads: {w_reads}", file=sys.stderr)

    # --- Targeted metrics from targeted BAM + NanoStats ----------------------
    t_reads, t_bases = count_reads_and_yield(args.targeted_bam)
    t_yield = round(t_bases / 1e9, 4)
    print(f"  Targeted reads: {t_reads}  yield: {t_yield} Gb", file=sys.stderr)

    nano = parse_nanostats(args.targeted_nanostats)
    t_length  = nano["mean_length"]
    t_n50     = int(nano["n50"])
    t_quality = nano["mean_quality"]
    print(f"  Targeted length: {t_length} bp  N50: {t_n50} bp  quality: {t_quality}",
          file=sys.stderr)

    # --- Write output TSV ----------------------------------------------------
    # --- Duration from sequencing_summary -----------------------------------
    barcode = args.barcode if args.barcode else args.sample
    duration = compute_duration(args.seq_summary, barcode)

    rows = [
        ("Number_of_reads",    w_reads,   t_reads),
        ("Total_yield_Gb",     w_yield,   t_yield),
        ("Mean_read_length_bp", w_length, t_length),
        ("Read_N50_bp",        w_n50,     t_n50),
        ("Mean_read_quality",  w_quality, t_quality),
    ]
    if duration is not None:
        rows.append(("Sequencing_duration_hours", duration, "—"))

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("metric\twhole_dataset\ttargeted_regions\n")
        for metric, whole_val, targ_val in rows:
            fh.write(f"{metric}\t{whole_val}\t{targ_val}\n")

    print(f"  Written: {args.out_tsv}", file=sys.stderr)
    print(f"[build_sequencing_performance] complete", file=sys.stderr)


if __name__ == "__main__":
    main()
