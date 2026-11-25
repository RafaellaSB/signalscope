#!/usr/bin/env python3
"""check_clinical_methylation.py — Screen curated clinical methylation loci.

Reads the genome-wide modkit bedMethyl pileup and a curated BED of clinically
relevant methylation loci. For each locus, computes the mean methylation
fraction and flags abnormal values.

Output: TSV with one row per clinical locus.

Columns:
  gene, region, chrom, start, end, clinical_context,
  cpg_count, mean_methylation_pct, expected_methylation,
  status (NORMAL / ABNORMAL / LOW_COVERAGE / NO_DATA),
  detail

Usage:
    python check_clinical_methylation.py \
        --bedmethyl coverage/sample.cpg.combined.bed.gz \
        --clinical_bed pipeline/assets/clinical_methylation_loci.bed \
        --output report/data/sample.clinical_methylation.tsv \
        --min_coverage 5

Part of SignalScope v2.0-wgs methylation module.
"""

import argparse
import logging
import sys

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)

# Thresholds for flagging abnormal methylation
# "low" expected = promoter normally unmethylated.  Flag if mean > 25%.
# "intermediate" expected = imprinting, normally ~50%. Flag if < 20% or > 80%.
# "high" expected = normally methylated. Flag if mean < 50%.
THRESHOLDS = {
    "low":          {"abnormal_above": 25.0, "abnormal_below": None},
    "intermediate": {"abnormal_above": 80.0, "abnormal_below": 20.0},
    "high":         {"abnormal_above": None,  "abnormal_below": 50.0},
}

MIN_CPGS_FOR_CALL = 3   # need at least this many CpGs to make a call


def parse_args():
    p = argparse.ArgumentParser(description="Screen clinical methylation loci.")
    p.add_argument("--bedmethyl", required=True,
                   help="Bgzipped + tabix-indexed bedMethyl (.bed.gz).")
    p.add_argument("--clinical_bed", required=True,
                   help="Curated clinical methylation loci BED (with header).")
    p.add_argument("--output", required=True, help="Output TSV.")
    p.add_argument("--min_coverage", type=int, default=5,
                   help="Minimum per-CpG read coverage (default: 5).")
    return p.parse_args()


def load_clinical_loci(path):
    """Read the curated clinical BED file.
    Expected columns: chrom, start, end, gene, region, expected_methylation, clinical_context
    """
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            rows.append({
                "chrom": parts[0],
                "start": int(parts[1]),
                "end":   int(parts[2]),
                "gene":  parts[3],
                "region": parts[4],
                "expected_methylation": parts[5],
                "clinical_context": parts[6],
            })
    log.info("Loaded %d clinical methylation loci.", len(rows))
    return rows


def query_region(bedmethyl_path, chrom, start, end, min_cov):
    """Return list of pct_modified values for CpGs in region."""
    import pysam

    pcts = []
    try:
        tbx = pysam.TabixFile(bedmethyl_path)
        query_chrom = chrom
        if query_chrom not in tbx.contigs:
            alt = chrom.replace("chr", "") if chrom.startswith("chr") else f"chr{chrom}"
            if alt in tbx.contigs:
                query_chrom = alt
            else:
                return pcts
        for row in tbx.fetch(query_chrom, max(0, start), end):
            fields = row.split("\t")
            if len(fields) < 11:
                continue
            cov = int(fields[9])
            if cov < min_cov:
                continue
            pcts.append(float(fields[10]))
    except Exception as e:
        log.warning("Tabix query failed for %s:%d-%d: %s", chrom, start, end, e)
    return pcts


def evaluate_locus(pcts, expected):
    """Determine status for a clinical locus.
    Returns (status, detail) tuple.
    """
    if not pcts:
        return "NO_DATA", "No CpG sites found with sufficient coverage"

    if len(pcts) < MIN_CPGS_FOR_CALL:
        mean_pct = sum(pcts) / len(pcts)
        return "LOW_COVERAGE", f"Only {len(pcts)} CpG(s) (need {MIN_CPGS_FOR_CALL}+). Mean {mean_pct:.1f}%"

    mean_pct = sum(pcts) / len(pcts)
    thresh = THRESHOLDS.get(expected, THRESHOLDS["low"])

    above = thresh.get("abnormal_above")
    below = thresh.get("abnormal_below")

    if above is not None and mean_pct > above:
        if expected == "low":
            return "ABNORMAL", f"Hypermethylated: {mean_pct:.1f}% (expected <{above:.0f}% for unmethylated promoter)"
        else:
            return "ABNORMAL", f"Hypermethylated: {mean_pct:.1f}% (expected <{above:.0f}%)"

    if below is not None and mean_pct < below:
        if expected == "intermediate":
            return "ABNORMAL", f"Hypomethylated: {mean_pct:.1f}% (expected >{below:.0f}% for imprinted region)"
        else:
            return "ABNORMAL", f"Hypomethylated: {mean_pct:.1f}% (expected >{below:.0f}%)"

    return "NORMAL", f"Mean {mean_pct:.1f}% (within expected range for {expected} methylation)"


def main():
    args = parse_args()

    loci = load_clinical_loci(args.clinical_bed)

    results = []
    abnormal_count = 0

    for locus in loci:
        pcts = query_region(args.bedmethyl, locus["chrom"], locus["start"],
                            locus["end"], args.min_coverage)

        status, detail = evaluate_locus(pcts, locus["expected_methylation"])

        mean_pct = round(sum(pcts) / len(pcts), 1) if pcts else float("nan")

        results.append({
            "gene": locus["gene"],
            "region": locus["region"],
            "chrom": locus["chrom"],
            "start": locus["start"],
            "end": locus["end"],
            "expected_methylation": locus["expected_methylation"],
            "clinical_context": locus["clinical_context"],
            "cpg_count": len(pcts),
            "mean_methylation_pct": mean_pct,
            "status": status,
            "detail": detail,
        })

        if status == "ABNORMAL":
            abnormal_count += 1
            log.warning("ABNORMAL: %s %s — %s", locus["gene"], locus["region"], detail)

    df = pd.DataFrame(results)
    df.to_csv(args.output, sep="\t", index=False)

    log.info("Clinical methylation screening complete:")
    log.info("  Loci screened:  %d", len(results))
    log.info("  NORMAL:         %d", sum(1 for r in results if r["status"] == "NORMAL"))
    log.info("  ABNORMAL:       %d", abnormal_count)
    log.info("  LOW_COVERAGE:   %d", sum(1 for r in results if r["status"] == "LOW_COVERAGE"))
    log.info("  NO_DATA:        %d", sum(1 for r in results if r["status"] == "NO_DATA"))
    log.info("  Output: %s", args.output)


if __name__ == "__main__":
    main()
