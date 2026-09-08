#!/usr/bin/env python3
"""
build_gene_features_from_exons.py — SignalScope v3

Builds a detailed gene_features.tsv from:
  1. Panel BED file (gene coordinates with flanking)
  2. Exon boundaries TSV (from build_exon_boundaries.py / Ensembl)

Output format matches the CKD pipeline's kidgen33.gene_features.tsv:
    chrom  start  end  gene  feature_type

Feature types: Flanking Up, Flanking Down, Exon/CDS, Intron

Usage:
  python3 build_gene_features_from_exons.py \
      --bed panel.bed \
      --exon_boundaries exon_boundaries.tsv \
      --flank 5000 \
      --out gene_features.tsv
"""

import argparse
import sys
from collections import defaultdict


def read_bed(path):
    """Read panel BED: chrom, start, end, gene."""
    genes = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, gene = parts[0], int(parts[1]), int(parts[2]), parts[3]
            # Keep widest span if gene appears multiple times
            if gene in genes:
                old = genes[gene]
                genes[gene] = (chrom, min(old[1], start), max(old[2], end), gene)
            else:
                genes[gene] = (chrom, start, end, gene)
    return genes


def read_exon_boundaries(path):
    """Read exon boundaries: gene, chrom, exon_start, exon_end."""
    exons = defaultdict(list)
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            gene = parts[0]
            chrom = parts[1]
            start = int(parts[2])
            end = int(parts[3])
            exons[gene].append((chrom, start, end))
    return exons


def normalize_chrom(chrom, target_has_chr):
    """Ensure chrom matches target convention."""
    if target_has_chr and not chrom.startswith("chr"):
        return "chr" + chrom
    elif not target_has_chr and chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def build_features(gene_name, gene_chrom, gene_start, gene_end, exon_list, flank):
    """
    Build feature intervals for one gene:
    Flanking Up → [Exon/CDS → Intron → Exon/CDS → ...] → Flanking Down
    """
    features = []

    # Determine chr convention from gene_chrom
    has_chr = gene_chrom.startswith("chr")

    # Sort exons by start position
    sorted_exons = sorted(exon_list, key=lambda x: x[1])

    # Normalize exon chroms to match gene_chrom
    sorted_exons = [(normalize_chrom(c, has_chr), s, e) for c, s, e in sorted_exons]

    # Gene body boundaries (excluding flanking)
    body_start = gene_start + flank
    body_end = gene_end - flank

    # If BED doesn't have flanking built in, use gene boundaries directly
    if body_start >= body_end:
        body_start = gene_start
        body_end = gene_end

    # Flanking Up
    if gene_start < body_start:
        features.append((gene_chrom, gene_start, body_start, gene_name, "Flanking Up"))

    # Build exon/intron intervals within gene body
    # Clip exons to gene body boundaries
    clipped_exons = []
    for chrom, es, ee in sorted_exons:
        cs = max(es, body_start)
        ce = min(ee, body_end)
        if cs < ce:
            clipped_exons.append((cs, ce))

    if clipped_exons:
        # Before first exon: intron
        if clipped_exons[0][0] > body_start:
            features.append((gene_chrom, body_start, clipped_exons[0][0], gene_name, "Intron"))

        for i, (es, ee) in enumerate(clipped_exons):
            # Exon
            features.append((gene_chrom, es, ee, gene_name, "Exon/CDS"))

            # Intron between this exon and the next
            if i + 1 < len(clipped_exons):
                next_start = clipped_exons[i + 1][0]
                if ee < next_start:
                    features.append((gene_chrom, ee, next_start, gene_name, "Intron"))

        # After last exon: intron
        if clipped_exons[-1][1] < body_end:
            features.append((gene_chrom, clipped_exons[-1][1], body_end, gene_name, "Intron"))
    else:
        # No exons mapped — entire gene body is one block
        features.append((gene_chrom, body_start, body_end, gene_name, "Gene Body"))

    # Flanking Down
    if body_end < gene_end:
        features.append((gene_chrom, body_end, gene_end, gene_name, "Flanking Down"))

    return features


def main():
    p = argparse.ArgumentParser(description="Build detailed gene_features from exon boundaries")
    p.add_argument("--bed", required=True, help="Panel BED file (4-col: chrom, start, end, gene)")
    p.add_argument("--exon_boundaries", required=True, help="Exon boundaries TSV from build_exon_boundaries.py")
    p.add_argument("--flank", type=int, default=5000, help="Flanking region size in bp (default: 5000)")
    p.add_argument("--out", required=True, help="Output gene_features.tsv")
    args = p.parse_args()

    # Read inputs
    genes = read_bed(args.bed)
    exons = read_exon_boundaries(args.exon_boundaries)

    print(f"Panel: {len(genes)} genes from {args.bed}", file=sys.stderr)
    print(f"Exon boundaries: {sum(len(v) for v in exons.values())} exons for {len(exons)} genes", file=sys.stderr)

    # Detect chr convention from BED
    sample_chrom = list(genes.values())[0][0] if genes else "chr1"
    bed_has_chr = sample_chrom.startswith("chr")
    print(f"BED chr convention: {'chr-prefixed' if bed_has_chr else 'no-chr'}", file=sys.stderr)

    all_features = []
    genes_with_exons = 0
    genes_without_exons = 0

    for gene_name, (chrom, start, end, _) in sorted(genes.items()):
        gene_exons = exons.get(gene_name, [])

        if gene_exons:
            genes_with_exons += 1
            features = build_features(gene_name, chrom, start, end, gene_exons, args.flank)
        else:
            genes_without_exons += 1
            # Fallback: coarse Gene Body + Flanking
            body_start = start + args.flank
            body_end = end - args.flank
            if body_start >= body_end:
                body_start = start
                body_end = end
            features = []
            if start < body_start:
                features.append((chrom, start, body_start, gene_name, "Flanking Up"))
            features.append((chrom, body_start, body_end, gene_name, "Gene Body"))
            if body_end < end:
                features.append((chrom, body_end, end, gene_name, "Flanking Down"))

        all_features.extend(features)

    # Write output
    with open(args.out, "w") as f:
        f.write("chrom\tstart\tend\tgene\tfeature_type\n")
        for chrom, start, end, gene, ftype in all_features:
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\t{ftype}\n")

    print(f"\nWritten: {args.out}", file=sys.stderr)
    print(f"  Total features: {len(all_features)}", file=sys.stderr)
    print(f"  Genes with exon detail: {genes_with_exons}", file=sys.stderr)
    print(f"  Genes with coarse fallback: {genes_without_exons}", file=sys.stderr)

    # Summary of feature types
    from collections import Counter
    counts = Counter(f[4] for f in all_features)
    for ftype, count in sorted(counts.items()):
        print(f"  {ftype}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
