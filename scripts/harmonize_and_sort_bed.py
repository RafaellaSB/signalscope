#!/usr/bin/env python3
"""
harmonize_and_sort_bed.py

Harmonise a panel BED's contig naming to match the reference, then sort it
into reference contig order. Prevents two failure modes:

  1. chr/no-chr naming mismatch between BED and reference. If the BED uses
     'chr1' but the reference uses '1' (or vice versa), a naive contig-order
     sort silently drops EVERY interval -> empty BED -> callers get no regions
     -> blank report. We detect the mismatch and rewrite the BED to match the
     reference naming instead of silently dropping intervals.

  2. Unsorted BED causing Clair3 to miss contigs. After harmonising, we sort
     by reference contig order + start coordinate.

Fails LOUDLY (non-zero exit + clear message) if harmonisation cannot reconcile
the BED to the reference, rather than emitting an empty/partial BED.

Usage:
    harmonize_and_sort_bed.py --bed panel.bed --fai ref.fa.fai --out sorted.bed
"""
import argparse
import sys
import os


def read_fai_contigs(fai_path):
    """Return (ordered_list, set) of contig names from a .fai, in file order."""
    contigs = []
    with open(fai_path) as fh:
        for line in fh:
            if line.strip():
                contigs.append(line.split("\t")[0])
    return contigs, set(contigs)


def detect_chr_style(names):
    """Return 'chr' if the majority of autosome-like names start with 'chr',
    else 'nochr'. Looks at chr1/1 .. chr22/22, X, Y, M/MT."""
    chr_prefixed = sum(1 for n in names if n.startswith("chr"))
    return "chr" if chr_prefixed > len(names) / 2 else "nochr"


def harmonise_name(bed_contig, ref_style, ref_contig_set):
    """Rewrite a single BED contig name to the reference naming style.
    Returns the harmonised name if it exists in the reference, else None."""
    # Already matches the reference as-is?
    if bed_contig in ref_contig_set:
        return bed_contig

    # Mitochondrial special cases: chrM <-> chrMT <-> MT <-> M
    mito_variants = {"chrM", "chrMT", "MT", "M"}
    if bed_contig in mito_variants:
        for cand in mito_variants:
            if cand in ref_contig_set:
                return cand

    # Toggle chr prefix to match reference style
    if ref_style == "chr":
        cand = bed_contig if bed_contig.startswith("chr") else "chr" + bed_contig
    else:  # ref is nochr
        cand = bed_contig[3:] if bed_contig.startswith("chr") else bed_contig

    return cand if cand in ref_contig_set else None


def main():
    p = argparse.ArgumentParser(description="Harmonise + sort a panel BED to a reference .fai")
    p.add_argument("--bed", required=True, help="Input panel BED")
    p.add_argument("--fai", required=True, help="Reference .fai (samtools faidx index)")
    p.add_argument("--out", required=True, help="Output harmonised+sorted BED")
    args = p.parse_args()

    if not os.path.isfile(args.bed):
        sys.exit(f"[FATAL] BED not found: {args.bed}")
    if not os.path.isfile(args.fai):
        sys.exit(f"[FATAL] reference .fai not found: {args.fai}")

    ref_order, ref_set = read_fai_contigs(args.fai)
    ref_rank = {name: i for i, name in enumerate(ref_order)}
    ref_style = detect_chr_style(ref_order)

    # Read BED, detect its style, harmonise each interval
    rows = []            # (rank, chrom, full_line_fields)
    bed_contigs_seen = set()
    dropped = []         # contigs we could not reconcile
    n_in = 0
    with open(args.bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            n_in += 1
            fields = line.rstrip("\n").split("\t")
            bed_chrom = fields[0]
            bed_contigs_seen.add(bed_chrom)
            new_chrom = harmonise_name(bed_chrom, ref_style, ref_set)
            if new_chrom is None:
                dropped.append(bed_chrom)
                continue
            fields[0] = new_chrom
            rows.append((ref_rank[new_chrom], int(fields[1]), fields))

    bed_style = detect_chr_style(list(bed_contigs_seen))

    # Report what we detected
    print(f"[INFO] reference naming style: {ref_style}; BED naming style: {bed_style}", file=sys.stderr)
    if ref_style != bed_style:
        print(f"[INFO] naming mismatch detected -> harmonised BED contigs to '{ref_style}' style", file=sys.stderr)

    # GUARD 1: unreconcilable contigs
    unique_dropped = sorted(set(dropped))
    if unique_dropped:
        sys.exit(
            f"[FATAL] {len(unique_dropped)} BED contig(s) do not exist in the reference "
            f"even after chr/no-chr harmonisation: {', '.join(unique_dropped[:10])}"
            f"{' ...' if len(unique_dropped) > 10 else ''}\n"
            f"        This usually means the BED was built for a DIFFERENT assembly "
            f"(e.g. GRCh38 vs T2T) than the reference. Check that your BED and reference "
            f"are the same genome build."
        )

    # GUARD 2: empty output (should be caught by GUARD 1, but belt-and-braces)
    if not rows:
        sys.exit(
            f"[FATAL] harmonisation produced an EMPTY BED (input had {n_in} intervals). "
            f"Refusing to continue — an empty BED would make variant callers return nothing. "
            f"Check BED/reference contig naming."
        )

    # Sort by reference contig order, then start coordinate
    rows.sort(key=lambda r: (r[0], r[1]))

    with open(args.out, "w") as out:
        for _, _, fields in rows:
            out.write("\t".join(fields) + "\n")

    print(f"[INFO] harmonised+sorted BED written: {args.out} "
          f"({len(rows)}/{n_in} intervals kept)", file=sys.stderr)


if __name__ == "__main__":
    main()
