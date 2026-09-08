#!/usr/bin/env python3
"""
parse_nanostats.py

Extracts basecall quality threshold percentages from NanoPlot's NanoStats.txt
and writes a clean TSV for downstream plotting.

Output format (tab-separated):
    threshold   percent_reads
    Q10         99.7
    Q15         78.2
    Q20         30.5
    Q25         4.5
    Q30         0.4

Usage:
    python3 parse_nanostats.py \
        --nanostats  samples/{sample}/qc/nanoplot/NanoStats.txt \
        --out_tsv    samples/{sample}/report/data/{sample}.qc_thresholds.tsv \
        --sample     barcode04
"""

import argparse
import os
import re
import sys


THRESHOLDS = ["Q10", "Q15", "Q20", "Q25", "Q30"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nanostats", required=True,
                   help="Path to NanoStats.txt")
    p.add_argument("--out_tsv",   required=True,
                   help="Output TSV path")
    p.add_argument("--sample",    required=True,
                   help="Sample identifier (for log messages)")
    return p.parse_args()


def extract_thresholds(nanostats_path):
    """
    Parse NanoStats.txt for Q10-Q30 percentage values.

    NanoStats lines of interest look like:
        >Q10:   17457842 (99.7%) 9977.9Mb
        >Q15:   13705034 (78.2%) 7838.1Mb
        >Q20:   5329286 (30.5%) 2975.6Mb
        >Q25:   791157 (4.5%) 452.5Mb
        >Q30:   74228 (0.4%) 42.5Mb

    The leading '>' may or may not be present depending on NanoPlot version.
    We match on the threshold label and capture the parenthesised percentage.
    """
    results = {}

    with open(nanostats_path) as fh:
        for line in fh:
            line = line.strip()
            for q in THRESHOLDS:
                if re.match(rf"^>?{q}[:\s]", line, re.IGNORECASE):
                    m = re.search(r"\((\d+(?:\.\d+)?)\s*%\)", line)
                    if m:
                        results[q] = float(m.group(1))
                        break

    missing = [q for q in THRESHOLDS if q not in results]
    if missing:
        print(f"  [ERROR] Could not find thresholds: {', '.join(missing)}",
              file=sys.stderr)
        print(f"  Check NanoStats.txt format: {nanostats_path}", file=sys.stderr)
        sys.exit(1)

    return results


def main():
    args = parse_args()
    print(f"[parse_nanostats] sample={args.sample}", file=sys.stderr)
    print(f"  input : {args.nanostats}", file=sys.stderr)

    results = extract_thresholds(args.nanostats)

    for q, pct in results.items():
        print(f"  {q}: {pct}%", file=sys.stderr)

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("threshold\tpercent_reads\n")
        for q in THRESHOLDS:
            fh.write(f"{q}\t{results[q]}\n")

    print(f"  output: {args.out_tsv}", file=sys.stderr)
    print(f"[parse_nanostats] complete", file=sys.stderr)


if __name__ == "__main__":
    main()
