#!/usr/bin/env python3
"""
coverage_build_tables.py — SignalScope

Reads the mosdepth regions BED.gz (one row per gene, mean depth pre-computed)
and the gene intervals TSV (canonical gene list), then writes two TSVs:

  {sample}.coverage_by_gene.tsv  — all genes with depth (intermediate)
  {sample}.coverage.final.tsv    — same content, used by R plot

Mosdepth BED format (no header):
    chrom  start  end  gene  mean_depth
    1      5862811  6002473  NPHP4  12.94

Gene intervals TSV format (no header):
    chrom  start  end  gene  strand
    2      240858824  240880502  AGXT  +

The canonical gene list comes from gene_intervals_ref.
Any gene in the canonical list but missing from mosdepth output gets depth=0.0
and a warning — this should never happen in a correctly run pipeline.

Semicolon-joined genes (e.g. COL4A3;COL4A4) in the BED are split and each
gene gets the same depth value (depth is NOT divided).

Usage:
    python3 coverage_build_tables.py \
        --regions_bed  {sample}.regions.bed.gz \
        --gene_ref     kidgen33.gene_intervals.tsv \
        --out_by_gene  {sample}.coverage_by_gene.tsv \
        --out_final    {sample}.coverage.final.tsv \
        --sample       barcode06
"""
import argparse
import gzip
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regions_bed",  required=True)
    p.add_argument("--gene_ref",     required=True)
    p.add_argument("--out_by_gene",  required=True)
    p.add_argument("--out_final",    required=True)
    p.add_argument("--sample",       required=True)
    return p.parse_args()


def load_canonical_genes(path):
    """
    Load ordered canonical gene list from gene_intervals.tsv.
    Format: chrom  start  end  gene  strand  (no header, tab-separated)
    Returns list of gene names in file order.
    """
    genes = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            gene = parts[3].strip()
            if gene and gene not in seen:
                genes.append(gene)
                seen.add(gene)
    return genes


def load_mosdepth_depths(path):
    """
    Load per-gene mean depth from mosdepth regions BED.gz.
    Format: chrom  start  end  gene  mean_depth  (no header)
    Handles semicolon-joined genes: COL4A3;COL4A4 → both get same depth.
    Returns dict: gene -> mean_depth (float)
    """
    depths = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            name  = parts[3].strip()
            try:
                depth = float(parts[4])
            except ValueError:
                continue
            # Split semicolon-joined genes — each gets same depth, no division
            for gene in name.split(";"):
                gene = gene.strip()
                if gene:
                    depths[gene] = depth
    return depths


def main():
    args = parse_args()
    print(f"[coverage_build_tables] sample={args.sample}", file=sys.stderr)
    print(f"  regions_bed: {args.regions_bed}", file=sys.stderr)
    print(f"  gene_ref   : {args.gene_ref}", file=sys.stderr)

    canonical = load_canonical_genes(args.gene_ref)
    print(f"  Canonical genes: {len(canonical)}", file=sys.stderr)

    depths = load_mosdepth_depths(args.regions_bed)
    print(f"  Mosdepth genes found: {len(depths)}", file=sys.stderr)

    # Warn on any missing genes
    missing = [g for g in canonical if g not in depths]
    if missing:
        print(f"  [WARN] Genes in canonical list but missing from mosdepth: "
              f"{', '.join(missing)}", file=sys.stderr)
        print(f"  These will be written with depth=0.0", file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_by_gene)), exist_ok=True)

    # Write both outputs — same content, canonical gene order
    for outpath in (args.out_by_gene, args.out_final):
        with open(outpath, "w") as fh:
            fh.write("gene\tmean_depth\n")
            for gene in canonical:
                depth = depths.get(gene, 0.0)
                fh.write(f"{gene}\t{depth:.2f}\n")
        print(f"  Written: {outpath}", file=sys.stderr)

    print(f"[coverage_build_tables] complete", file=sys.stderr)


if __name__ == "__main__":
    main()
