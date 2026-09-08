#!/usr/bin/env python3
"""
generate_shared_target_flags.py

Detects semicolon-joined gene labels in the mosdepth regions BED
(e.g. COL4A3;COL4A4) and writes a flags TSV used by the coverage plot
and HTML report.

Output: {sample}.shared_targets.tsv
Format:
    gene    shared_target
    COL4A3  COL4A3/COL4A4
    COL4A4  COL4A3/COL4A4

Rule: depth is NOT divided or modified. Annotate with * in plot only.

Usage:
    python3 generate_shared_target_flags.py \
        --regions_bed samples/{sample}/coverage/{sample}.regions.bed.gz \
        --out_flags   samples/{sample}/report/data/{sample}.shared_targets.tsv \
        --sample      barcode06
"""

import argparse
import gzip
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions_bed", required=True)
    p.add_argument("--out_flags",   required=True)
    p.add_argument("--sample",      required=True)
    return p.parse_args()


def main():
    args = parse_args()
    print(f"[generate_shared_target_flags] sample={args.sample}", file=sys.stderr)

    shared_pairs = {}
    opener = gzip.open if args.regions_bed.endswith(".gz") else open
    with opener(args.regions_bed, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            label = cols[3]
            # Accept either ";" (legacy) or "_" (new naming, e.g. COL4A4_COL4A3) as separator
            sep = None
            if ";" in label:
                sep = ";"
            elif "_" in label:
                # Only treat "_" as a separator if it splits into >= 2 valid gene-like tokens
                parts = [p.strip() for p in label.split("_") if p.strip()]
                if len(parts) >= 2:
                    sep = "_"
            if sep:
                genes = [g.strip() for g in label.split(sep) if g.strip()]
                display = "/".join(sorted(genes))
                for g in genes:
                    shared_pairs[g] = display

    os.makedirs(os.path.dirname(args.out_flags), exist_ok=True)
    with open(args.out_flags, "w") as fh:
        fh.write("gene\tshared_target\n")
        for gene, label in sorted(shared_pairs.items()):
            fh.write(f"{gene}\t{label}\n")

    if shared_pairs:
        print(f"  Shared genes: {', '.join(sorted(shared_pairs.keys()))}", file=sys.stderr)
    else:
        print(f"  No shared-target genes detected", file=sys.stderr)
    print(f"  Written: {args.out_flags}", file=sys.stderr)
    print(f"[generate_shared_target_flags] complete", file=sys.stderr)


if __name__ == "__main__":
    main()
