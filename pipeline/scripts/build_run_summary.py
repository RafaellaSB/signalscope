#!/usr/bin/env python3
"""
build_run_summary.py

Combines per-barcode run-level metrics from the sequencing summary file
with read-quality metrics from NanoStats.txt into two output tables:

  {sample}.run_summary.tsv         - machine-readable key/value table
  {sample}.run_summary_display.tsv - human-readable display version

Run-level metrics (from sequencing_summary_*.txt, column barcode_arrangement):
  Number_of_reads, Total_yield_Gb, Run_duration_hours

NanoPlot-derived metrics (from NanoStats.txt):
  Mean_read_length_bp, Read_N50_bp, Mean_read_quality

Usage:
    python3 build_run_summary.py \
        --seq_summary    run_metadata/sequencing_summary_*.txt \
        --nanostats      samples/{sample}/qc/nanoplot/{sample}_NanoStats.txt \
        --out_summary    samples/{sample}/report/data/{sample}.run_summary.tsv \
        --out_display    samples/{sample}/report/data/{sample}.run_summary_display.tsv \
        --sample         barcode06 \
        --barcode_col    25
"""

import argparse
import os
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seq_summary",  required=True)
    p.add_argument("--nanostats",    required=True)
    p.add_argument("--out_summary",  required=True)
    p.add_argument("--out_display",  required=True)
    p.add_argument("--sample",       required=True)
    p.add_argument("--barcode_col",  type=int, default=25,
                   help="1-based column index of barcode_arrangement (default: 25)")
    return p.parse_args()


def parse_seq_summary(path, sample, barcode_col_1based):
    col_idx = barcode_col_1based - 1
    n_reads = 0
    total_bases = 0
    max_start_time = 0.0
    start_time_idx = None
    sequence_length_idx = None

    print(f"  Scanning sequencing summary for barcode: {sample}", file=sys.stderr)

    with open(path) as fh:
        for i, line in enumerate(fh):
            cols = line.rstrip("\n").split("\t")
            if i == 0:
                # Validate that the configured column is actually barcode_arrangement
                if len(cols) > col_idx:
                    actual = cols[col_idx].strip()
                    if actual != "barcode_arrangement":
                        print(
                            f"  [ERROR] Column {barcode_col_1based} is '{actual}', "
                            f"expected 'barcode_arrangement'.",
                            file=sys.stderr
                        )
                        print(
                            f"  Check barcode_col in config.yaml or the file header.",
                            file=sys.stderr
                        )
                        sys.exit(1)
                try:
                    start_time_idx = cols.index("start_time")
                except ValueError:
                    start_time_idx = None
                for name in ("sequence_length_template", "sequence_length"):
                    if name in cols:
                        sequence_length_idx = cols.index(name)
                        break
                continue

            if len(cols) <= col_idx or cols[col_idx].strip() != sample:
                continue

            n_reads += 1
            if sequence_length_idx is not None and len(cols) > sequence_length_idx:
                try:
                    total_bases += int(cols[sequence_length_idx])
                except ValueError:
                    pass
            if start_time_idx is not None and len(cols) > start_time_idx:
                try:
                    t = float(cols[start_time_idx])
                    if t > max_start_time:
                        max_start_time = t
                except ValueError:
                    pass

    if n_reads == 0:
        print(f"  [ERROR] No rows found for barcode '{sample}' in column {barcode_col_1based}",
              file=sys.stderr)
        print(f"  File: {path}", file=sys.stderr)
        sys.exit(1)

    yield_gb = round(total_bases / 1e9, 2)
    duration_h = round(max_start_time / 3600, 2)
    print(f"  Reads found    : {n_reads:,}", file=sys.stderr)
    print(f"  Total yield    : {yield_gb} Gb", file=sys.stderr)
    print(f"  Run duration   : {duration_h} h", file=sys.stderr)
    return n_reads, yield_gb, duration_h


def parse_nanostats(path):
    results = {}
    patterns = {
        "mean_length":  re.compile(r"Mean\s+read\s+length[:\s]+([\d,.]+)", re.IGNORECASE),
        "n50":          re.compile(r"Read(?:\s+length)?\s+N50[:\s]+([\d,.]+)", re.IGNORECASE),
        "mean_quality": re.compile(r"Mean\s+read\s+quality[:\s]+([\d.]+)", re.IGNORECASE),
        "n_reads":      re.compile(r"Number\s+of\s+reads[:\s]+([\d,.]+)", re.IGNORECASE),
        "total_bases":  re.compile(r"Total\s+bases[:\s]+([\d,.]+)", re.IGNORECASE),
    }
    with open(path) as fh:
        for line in fh:
            for key, pat in list(patterns.items()):
                if key in results:
                    continue
                m = pat.search(line)
                if m:
                    results[key] = float(m.group(1).replace(",", ""))

    required = ["mean_length", "n50", "mean_quality"]
    missing = [k for k in required if k not in results]
    if missing:
        print(f"  [ERROR] Missing NanoStats metrics: {', '.join(missing)}", file=sys.stderr)
        print(f"  File: {path}", file=sys.stderr)
        sys.exit(1)

    print(f"  Mean length    : {results['mean_length']} bp", file=sys.stderr)
    print(f"  N50            : {results['n50']} bp", file=sys.stderr)
    print(f"  Mean quality   : {results['mean_quality']}", file=sys.stderr)
    return results


def write_summary(path, n_reads, yield_gb, duration_h, nano):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"Number_of_reads\t{n_reads}\n")
        fh.write(f"Total_yield_Gb\t{yield_gb}\n")
        fh.write(f"Run_duration_hours\t{duration_h}\n")
        fh.write(f"Mean_read_length_bp\t{nano['mean_length']}\n")
        fh.write(f"Read_N50_bp\t{int(nano['n50'])}\n")
        fh.write(f"Mean_read_quality\t{nano['mean_quality']}\n")
    print(f"  Written: {path}", file=sys.stderr)


def write_display(path, n_reads, yield_gb, duration_h, nano):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("Metric\tValue\n")
        _nr = f"{n_reads:,}" if isinstance(n_reads, int) else str(n_reads)
        fh.write(f"Number of reads\t{_nr}\n")
        fh.write(f"Total yield (Gb)\t{yield_gb}\n")
        fh.write(f"Run duration (hours)\t{duration_h}\n")
        fh.write(f"Mean read length (bp)\t{nano['mean_length']}\n")
        fh.write(f"Read N50 (bp)\t{int(nano['n50'])}\n")
        fh.write(f"Mean read quality\t{nano['mean_quality']}\n")
    print(f"  Written: {path}", file=sys.stderr)


def seq_summary_is_real(path):
    """A real MinKNOW sequencing_summary has many data rows and no 'dummy' placeholder."""
    if not path or not os.path.isfile(path):
        return False
    try:
        with open(path) as fh:
            first = fh.readline()          # header
            n = 0
            for line in fh:
                if "dummy" in line.lower():
                    return False
                n += 1
                if n > 5:                  # more than a stub → treat as real
                    return True
        return n > 1
    except Exception:
        return False


def main():
    args = parse_args()
    print(f"[build_run_summary] sample={args.sample}", file=sys.stderr)

    # Whole-dataset read count, yield, and QC all come from NanoStats (real reads),
    # NOT the sequencing_summary (which may be absent or a stub for public/BAM data).
    nano = parse_nanostats(args.nanostats)

    if "n_reads" in nano:
        n_reads = int(nano["n_reads"])
    else:
        n_reads = "NA"
    if "total_bases" in nano:
        yield_gb = round(nano["total_bases"] / 1e9, 2)
    else:
        yield_gb = "NA"
    print(f"  Whole-dataset reads (NanoStats): {n_reads}  yield: {yield_gb} Gb", file=sys.stderr)

    # Duration only from a REAL sequencing_summary; else N/A (a BAM/public dataset
    # has no run-duration information).
    duration_h = "NA"
    if seq_summary_is_real(args.seq_summary):
        try:
            _, _, duration_h = parse_seq_summary(
                args.seq_summary, args.sample, args.barcode_col
            )
        except SystemExit:
            duration_h = "NA"
    else:
        print("  [info] No real sequencing_summary — duration reported as NA.", file=sys.stderr)

    write_summary(args.out_summary, n_reads, yield_gb, duration_h, nano)
    write_display(args.out_display, n_reads, yield_gb, duration_h, nano)
    print(f"[build_run_summary] complete", file=sys.stderr)


if __name__ == "__main__":
    main()
