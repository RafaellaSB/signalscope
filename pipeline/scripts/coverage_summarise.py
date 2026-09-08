#!/usr/bin/env python3
"""
coverage_summarise.py — SignalScope

Reads coverage.final.tsv and writes:
  {sample}.coverage_summary.tsv  — summary metrics (mean, median, counts)
  {sample}.genes_below_10x.tsv   — genes with mean depth < 10x
  {sample}.genes_below_20x.tsv   — genes with mean depth < 20x

Input format (from coverage_build_tables.py):
    gene    mean_depth
    AGXT    11.74
    AHI1    16.47

Usage:
    python3 coverage_summarise.py \
        --final_tsv   {sample}.coverage.final.tsv \
        --out_summary {sample}.coverage_summary.tsv \
        --out_below10 {sample}.genes_below_10x.tsv \
        --out_below20 {sample}.genes_below_20x.tsv \
        --sample      barcode06
"""
import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--final_tsv",   required=True)
    p.add_argument("--out_summary", required=True)
    p.add_argument("--out_below10", required=True)
    p.add_argument("--out_below20", required=True)
    p.add_argument("--sample",      required=True)
    return p.parse_args()


def main():
    args = parse_args()
    print(f"[coverage_summarise] sample={args.sample}", file=sys.stderr)

    rows = []
    with open(args.final_tsv) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if i == 0 or not line:
                continue  # skip header
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0].strip()
            try:
                depth = float(parts[1])
            except ValueError:
                print(f"  [WARN] Could not parse depth for gene '{gene}': {parts[1]}",
                      file=sys.stderr)
                continue
            rows.append((gene, depth))

    if not rows:
        print(f"  [ERROR] No data rows found in {args.final_tsv}", file=sys.stderr)
        sys.exit(1)

    print(f"  Genes loaded: {len(rows)}", file=sys.stderr)

    depths     = [d for _, d in rows]
    n          = len(depths)
    mean_depth = sum(depths) / n
    sorted_d   = sorted(depths)
    # Median: middle value for odd n, average of two middle for even n
    if n % 2 == 1:
        median = sorted_d[n // 2]
    else:
        median = (sorted_d[n // 2 - 1] + sorted_d[n // 2]) / 2.0

    below10 = [(g, d) for g, d in rows if d < 10.0]
    below20 = [(g, d) for g, d in rows if d < 20.0]

    print(f"  Mean depth  : {mean_depth:.2f}x", file=sys.stderr)
    print(f"  Median depth: {median:.2f}x", file=sys.stderr)
    print(f"  Below 10x   : {len(below10)} genes", file=sys.stderr)
    print(f"  Below 20x   : {len(below20)} genes", file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_summary)), exist_ok=True)

    with open(args.out_summary, "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"n_genes\t{n}\n")
        fh.write(f"mean_depth\t{mean_depth:.2f}\n")
        fh.write(f"median_depth\t{median:.2f}\n")
        fh.write(f"genes_below_10x\t{len(below10)}\n")
        fh.write(f"genes_below_20x\t{len(below20)}\n")
    print(f"  Written: {args.out_summary}", file=sys.stderr)

    with open(args.out_below10, "w") as fh:
        fh.write("gene\tmean_depth\n")
        for g, d in sorted(below10, key=lambda x: x[1]):
            fh.write(f"{g}\t{d:.2f}\n")
    print(f"  Written: {args.out_below10}", file=sys.stderr)

    with open(args.out_below20, "w") as fh:
        fh.write("gene\tmean_depth\n")
        for g, d in sorted(below20, key=lambda x: x[1]):
            fh.write(f"{g}\t{d:.2f}\n")
    print(f"  Written: {args.out_below20}", file=sys.stderr)

    print(f"[coverage_summarise] complete", file=sys.stderr)


if __name__ == "__main__":
    main()
