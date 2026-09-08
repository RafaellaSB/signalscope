#!/usr/bin/env python3
"""
pseudogene_recovery.py  (SignalScope v3)

Generalised pseudogene read recovery via competitive alignment.

Recovers genuine reads from secondary alignments by competitively aligning
candidate reads against the real gene and its pseudogene(s). Reads scoring
higher for the real gene are confirmed and merged into the primary BAM.

Supported genes (shipped with pipeline):
  - PKD1:  6 pseudogenes (PKD1P1-P6), chr16 — validated on HG002
  - PMS2:  1 pseudogene (PMS2CL), chr7 — exons 11-15, 98% identity

Architecture:
  For each pseudogene-affected gene in the user's panel:
    1. Extract reads aligned to pseudogene region(s) from the all-reads BAM
    2. Competitively align these candidates against:
       a. Real gene reference (extracted from genome FASTA)
       b. Pseudogene reference(s) (extracted from genome FASTA)
    3. Compare alignment scores — reads scoring better for real gene are confirmed
    4. Merge confirmed reads into primary BAM

Usage:
  python3 pseudogene_recovery.py \\
      --bam_all   sample.bam \\
      --bam_pri   sample.primary.bam \\
      --ref       GRCh38.fa \\
      --genes     PKD1,PMS2 \\
      --out_dir   recovery/ \\
      --sample    barcode01 \\
      --threads   8

Output:
  {sample}.recovered.bam      — primary BAM with confirmed reads added
  {sample}.recovery_summary.tsv — per-gene recovery statistics
"""

import argparse
import os
import subprocess
import sys


# =============================================================================
# PSEUDOGENE DATABASE
#
# Each entry defines:
#   real_region:    GRCh38 coordinates of the real gene
#   pseudo_regions: list of GRCh38 coordinates for each pseudogene copy
#   description:    human-readable note for the report
#
# To add a new gene:
#   1. Add an entry to PSEUDOGENE_DB below
#   2. No other code changes needed — the recovery loop handles it generically
# =============================================================================

PSEUDOGENE_DB = {
    "PKD1": {
        "real_region": "16:2088708-2145898",
        "pseudo_regions": [
            "16:14886639-14928823",   # PKD1P1
            "16:15092009-15134193",   # PKD1P2
            "16:15297380-15339564",   # PKD1P3
            "16:15502750-15544934",   # PKD1P4
            "16:15708121-15750305",   # PKD1P5
            "16:15913491-15955675",   # PKD1P6
        ],
        "description": "6 pseudogenes (PKD1P1-P6) on chr16, ~97% identity",
    },
    "PMS2": {
        "real_region": "7:5970925-6009019",
        "pseudo_regions": [
            "7:6012870-6048756",      # PMS2CL
        ],
        "description": "1 pseudogene (PMS2CL) on chr7, 98% identity to exons 11-15",
    },
    # ── Future genes (add coordinates when validated) ──
    # "CHEK2": {
    #     "real_region": "22:28687743-28742422",
    #     "pseudo_regions": [
    #         "22:...",  # CHEK2P2 and processed pseudogenes
    #     ],
    #     "description": "Multiple processed pseudogenes affecting exons 10-14",
    # },
    # "SMN1": {
    #     "real_region": "5:70924941-70953015",
    #     "pseudo_regions": [
    #         "5:70049523-70077595",  # SMN2 (inverted duplicate)
    #     ],
    #     "description": "SMN2 inverted duplicate on chr5",
    # },
    # "CYP21A2": {
    #     "real_region": "6:32038265-32041670",
    #     "pseudo_regions": [
    #         "6:32005983-32009388",  # CYP21A1P
    #     ],
    #     "description": "CYP21A1P pseudogene on chr6, within HLA region",
    # },
}

# Genes where recovery is known to be ineffective (documented, not attempted)
RECOVERY_SKIP = {
    # Genes documented as recovery-ineffective can be listed here.
    # MUC1 removed (dropped from analysis). Leave empty unless intentionally skipping a gene.
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def run(cmd, log=None):
    """Execute a shell command, exit on failure."""
    if log:
        print(f"  CMD: {cmd}", file=log)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    return result.stdout.strip()


def detect_chr_naming(path, samtools, is_bam):
    """Return True if the BAM (@SQ) or FASTA (first header) uses 'chr' contig names."""
    if is_bam:
        out = subprocess.run(f"{samtools} view -H {path}", shell=True,
                             capture_output=True, text=True).stdout
        return any(l.startswith("@SQ") and "SN:chr" in l for l in out.splitlines())
    else:
        out = subprocess.run(f"head -1 {path}", shell=True,
                             capture_output=True, text=True).stdout
        return out.startswith(">chr")

def norm_region(region, uses_chr):
    """Normalise '7:100-200' or 'chr7:100-200' to match chr/no-chr naming."""
    chrom, coords = region.split(":", 1)
    chrom = chrom[3:] if chrom.startswith("chr") else chrom
    return (("chr" + chrom) if uses_chr else chrom) + ":" + coords

def get_depth(bam, region, samtools):
    """Calculate mean depth over a region."""
    out = subprocess.run(
        f"{samtools} depth -r {region} {bam} | "
        f"awk '{{sum+=$3; n++}} END{{if(n>0) print sum/n; else print 0}}'",
        shell=True, capture_output=True, text=True
    )
    try:
        return float(out.stdout.strip())
    except (ValueError, AttributeError):
        return 0.0


def count_fastq_reads(fastq_path):
    """Count reads in a FASTQ file."""
    out = subprocess.run(
        f"grep -c '^@' {fastq_path}",
        shell=True, capture_output=True, text=True
    )
    return int(out.stdout.strip()) if out.returncode == 0 else 0


def extract_region_reads(bam, regions, out_fastq, samtools):
    """Extract reads from specific regions as FASTQ."""
    regions_str = " ".join(f"'{r}'" for r in regions)
    run(f"{samtools} view -h {bam} {regions_str} | {samtools} fastq - > {out_fastq}")


# =============================================================================
# COMPETITIVE ALIGNMENT
# =============================================================================

def competitive_align(gene_name, fastq, ref, real_region, pseudo_regions,
                      out_dir, minimap2, samtools, threads):
    """
    Competitively align candidate reads against real gene vs pseudogene(s).

    Returns:
        confirmed_reads: list of read names scoring better for real gene
        n_evaluated_real: number of reads that aligned to real gene
        n_evaluated_pseudo: number of reads that aligned to pseudogene(s)
    """
    prefix = gene_name.lower()
    real_fa   = os.path.join(out_dir, f"{prefix}_real.fa")
    pseudo_fa = os.path.join(out_dir, f"{prefix}_pseudogenes.fa")
    sam_real  = os.path.join(out_dir, f"vs_{prefix}_real.sam")
    sam_pseudo = os.path.join(out_dir, f"vs_{prefix}_pseudogenes.sam")

    # Extract reference sequences for real gene and pseudogene(s)
    run(f"{samtools} faidx {ref} '{real_region}' > {real_fa}")

    pseudo_args = " ".join(f"'{r}'" for r in pseudo_regions)
    run(f"{samtools} faidx {ref} {pseudo_args} > {pseudo_fa}")

    # Align candidate reads against real gene
    run(f"{minimap2} -ax map-ont --secondary=no -t {threads} "
        f"{real_fa} {fastq} > {sam_real} 2>/dev/null")

    # Align candidate reads against pseudogene(s)
    run(f"{minimap2} -ax map-ont --secondary=no -t {threads} "
        f"{pseudo_fa} {fastq} > {sam_pseudo} 2>/dev/null")

    # Parse alignment scores
    real_scores = _parse_alignment_scores(sam_real)
    pseudo_scores = _parse_alignment_scores(sam_pseudo)

    # Compare: reads scoring strictly better for real gene are confirmed
    confirmed = []
    for read_name in real_scores:
        score_real   = real_scores.get(read_name, -1)
        score_pseudo = pseudo_scores.get(read_name, -1)
        if score_real > score_pseudo:
            confirmed.append(read_name)

    return confirmed, len(real_scores), len(pseudo_scores)


def _parse_alignment_scores(sam_path):
    """Parse AS:i: alignment scores from a SAM file."""
    scores = {}
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            flag = int(parts[1])
            if flag & 4:  # unmapped
                continue
            name = parts[0]
            for tag in parts[11:]:
                if tag.startswith("AS:i:"):
                    scores[name] = int(tag[5:])
                    break
    return scores


# =============================================================================
# READ EXTRACTION AND BAM MERGING
# =============================================================================

def extract_confirmed_reads(bam_all, read_names, out_bam, samtools, tmp_dir):
    """Extract confirmed reads by name from the all-reads BAM."""
    if not read_names:
        # Create empty BAM with header only
        run(f"{samtools} view -H {bam_all} | {samtools} view -b > {out_bam}")
        run(f"{samtools} index {out_bam}")
        return 0

    names_file = os.path.join(tmp_dir, "confirmed_reads.txt")
    with open(names_file, "w") as f:
        for name in read_names:
            f.write(name + "\n")

    run(f"{samtools} view -h -N {names_file} {bam_all} | "
        f"{samtools} view -b | {samtools} sort > {out_bam}")
    run(f"{samtools} index {out_bam}")
    return len(read_names)


# =============================================================================
# AUTO-DETECTION
# =============================================================================

def detect_pseudogenes_in_panel(bed_path):
    """
    Check which genes in the user's panel BED have known pseudogene issues.
    Returns list of gene names that need recovery.

    No internet required — checks against the bundled PSEUDOGENE_DB.
    """
    panel_genes = set()
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                panel_genes.add(parts[3].upper())

    detected = []
    for gene in PSEUDOGENE_DB:
        if gene.upper() in panel_genes:
            detected.append(gene)

    skipped = []
    for gene in RECOVERY_SKIP:
        if gene.upper() in panel_genes:
            skipped.append(gene)

    return detected, skipped


# =============================================================================
# MAIN
# =============================================================================

def main():
    p = argparse.ArgumentParser(
        description="Pseudogene read recovery via competitive alignment"
    )
    p.add_argument("--bam_all",   required=True, help="All-reads BAM (primary + secondary)")
    p.add_argument("--bam_pri",   required=True, help="Primary-filtered BAM")
    p.add_argument("--ref",       required=True, help="Reference genome FASTA (GRCh38)")
    p.add_argument("--genes",     required=False, default="",
                   help="Comma-separated list of genes to recover. "
                        "If empty, auto-detects from --panel_bed")
    p.add_argument("--panel_bed", required=False, default="",
                   help="Panel BED file for auto-detection of pseudogene genes")
    p.add_argument("--out_dir",   required=True, help="Output directory for recovery files")
    p.add_argument("--sample",    required=True, help="Sample name")
    p.add_argument("--threads",   type=int, default=8)
    p.add_argument("--samtools",  default="samtools")
    p.add_argument("--minimap2",  default="minimap2")
    args = p.parse_args()

    tmp_dir = os.path.join(args.out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    sam = args.samtools
    mm2 = args.minimap2

    # ── Determine which genes to recover ──────────────────────────────────────
    if args.genes:
        target_genes = [g.strip().upper() for g in args.genes.split(",")]
    elif args.panel_bed:
        target_genes, skipped_genes = detect_pseudogenes_in_panel(args.panel_bed)
        if skipped_genes:
            for g in skipped_genes:
                info = RECOVERY_SKIP[g]
                print(f"  [{g}] Skipped: {info['reason']}", file=sys.stderr)
    else:
        target_genes = list(PSEUDOGENE_DB.keys())  # recover all known genes

    # Filter to only genes in the database
    valid_genes = [g for g in target_genes if g in PSEUDOGENE_DB]
    unknown = [g for g in target_genes if g not in PSEUDOGENE_DB]
    if unknown:
        print(f"  WARNING: No pseudogene data for: {', '.join(unknown)}. Skipping.",
              file=sys.stderr)

    print(f"[pseudogene_recovery] sample={args.sample}", file=sys.stderr)
    print(f"  Genes to recover: {', '.join(valid_genes) if valid_genes else 'none'}",
          file=sys.stderr)

    # ── Per-gene recovery loop ────────────────────────────────────────────────
    all_confirmed = []
    results = []

    bam_uses_chr = detect_chr_naming(args.bam_all, sam, is_bam=True)
    ref_uses_chr = detect_chr_naming(args.ref, sam, is_bam=False)
    print(f"  [naming] bam_uses_chr={bam_uses_chr} ref_uses_chr={ref_uses_chr}", file=sys.stderr)
    for gene in valid_genes:
        db = PSEUDOGENE_DB[gene]
        real_region = db["real_region"]
        pseudo_regions = db["pseudo_regions"]
        real_region_bam    = norm_region(real_region, bam_uses_chr)
        pseudo_regions_bam = [norm_region(r, bam_uses_chr) for r in pseudo_regions]
        real_region_ref    = norm_region(real_region, ref_uses_chr)
        pseudo_regions_ref = [norm_region(r, ref_uses_chr) for r in pseudo_regions]

        print(f"\n  [{gene}] {db['description']}", file=sys.stderr)

        # Measure baseline depth
        depth_before = get_depth(args.bam_pri, real_region_bam, sam)
        print(f"  [{gene}] Baseline depth: {depth_before:.1f}x", file=sys.stderr)

        # Extract candidate reads from pseudogene regions
        fastq = os.path.join(tmp_dir, f"{gene.lower()}_candidates.fastq")
        extract_region_reads(args.bam_all, pseudo_regions_bam, fastq, sam)
        n_candidates = count_fastq_reads(fastq)
        print(f"  [{gene}] Candidate reads from pseudogene regions: {n_candidates}",
              file=sys.stderr)

        confirmed = []
        n_evaluated = 0
        if n_candidates > 0:
            print(f"  [{gene}] Running competitive alignment...", file=sys.stderr)
            confirmed, n_evaluated, _ = competitive_align(
                gene, fastq, args.ref, real_region_ref, pseudo_regions_ref,
                tmp_dir, mm2, sam, args.threads
            )
            print(f"  [{gene}] Confirmed genuine reads: {len(confirmed)}", file=sys.stderr)

        all_confirmed.extend(confirmed)

        results.append({
            "gene": gene,
            "candidates": n_candidates,
            "evaluated": n_evaluated,
            "confirmed": len(confirmed),
            "depth_before": depth_before,
            "method": "competitive_alignment",
            "description": db["description"],
        })

    # ── Add skip-documented genes to results (e.g. MUC1) ─────────────────────
    for gene, info in RECOVERY_SKIP.items():
        # Only include if gene is in the panel (or always if no BED provided)
        include = True
        if args.panel_bed:
            panel_genes = set()
            with open(args.panel_bed) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        panel_genes.add(parts[3].upper())
            include = gene.upper() in panel_genes

        if include:
            depth = get_depth(args.bam_pri, norm_region(info["region"], bam_uses_chr), sam)
            results.append({
                "gene": gene,
                "candidates": 0,
                "evaluated": 0,
                "confirmed": 0,
                "depth_before": depth,
                "method": "not_attempted",
                "description": info["reason"],
            })

    # ── Merge confirmed reads into primary BAM ────────────────────────────────
    all_confirmed = list(set(all_confirmed))  # deduplicate
    print(f"\n  Total confirmed reads to merge: {len(all_confirmed)}", file=sys.stderr)

    confirmed_bam = os.path.join(tmp_dir, "confirmed_reads.bam")
    extract_confirmed_reads(args.bam_all, all_confirmed, confirmed_bam, sam, tmp_dir)

    recovered_bam = os.path.join(args.out_dir, f"{args.sample}.recovered.bam")
    if all_confirmed:
        run(f"{sam} merge -f {recovered_bam} {args.bam_pri} {confirmed_bam}")
    else:
        run(f"cp {args.bam_pri} {recovered_bam}")
    run(f"{sam} index {recovered_bam}")
    print(f"  Written: {recovered_bam}", file=sys.stderr)

    # ── Measure post-recovery depths and write summary ────────────────────────
    tsv_path = os.path.join(args.out_dir, f"{args.sample}.recovery_summary.tsv")
    with open(tsv_path, "w") as f:
        f.write("gene\tcandidates\tevaluated\tconfirmed\tdepth_before\t"
                "depth_after\tmethod\tdescription\n")

        for r in results:
            if r["method"] == "competitive_alignment":
                depth_after = get_depth(recovered_bam, norm_region(PSEUDOGENE_DB[r["gene"]]["real_region"], bam_uses_chr), sam)
            else:
                depth_after = r["depth_before"]  # unchanged for skipped genes

            f.write(f"{r['gene']}\t{r['candidates']}\t{r['evaluated']}\t"
                    f"{r['confirmed']}\t{r['depth_before']:.2f}\t{depth_after:.2f}\t"
                    f"{r['method']}\t{r['description']}\n")

    print(f"  Written: {tsv_path}", file=sys.stderr)
    print(f"[pseudogene_recovery] complete — {len(valid_genes)} gene(s) processed",
          file=sys.stderr)


if __name__ == "__main__":
    main()
