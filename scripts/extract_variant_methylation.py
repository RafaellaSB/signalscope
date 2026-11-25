#!/usr/bin/env python3
"""extract_variant_methylation.py — Add local methylation context to shortlisted variants.

For each shortlisted SNV/indel, extracts CpG methylation in a ±500 bp window
from the modkit bedMethyl pileup and adds summary columns to the shortlist TSV.

Columns added:
  meth_cpg_count     — number of CpGs in window with sufficient coverage
  meth_mean_pct      — mean methylation % across those CpGs
  meth_min_pct       — lowest single-CpG methylation %
  meth_max_pct       — highest single-CpG methylation %
  meth_low_cpgs      — count of CpGs with methylation < 20% (hypomethylated)
  meth_high_cpgs     — count of CpGs with methylation > 80% (hypermethylated)

Usage:
    python extract_variant_methylation.py \
        --shortlist results/shortlist_snv.tsv \
        --bedmethyl coverage/sample.cpg.combined.bed.gz \
        --output results/shortlist_snv_meth.tsv \
        --window 500 \
        --min_coverage 5

Part of SignalScope v2.0-wgs methylation module.
"""

import argparse
import gzip
import logging
import sys

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description="Add methylation context to shortlisted variants.")
    p.add_argument("--shortlist", required=True, help="Input shortlist TSV.")
    p.add_argument("--bedmethyl", required=True,
                   help="Bgzipped + tabix-indexed bedMethyl from modkit pileup (.bed.gz).")
    p.add_argument("--output", required=True, help="Output annotated shortlist TSV.")
    p.add_argument("--window", type=int, default=500,
                   help="Flanking window (bp) around each variant (default: 500).")
    p.add_argument("--min_coverage", type=int, default=5,
                   help="Minimum read coverage to include a CpG site (default: 5).")
    return p.parse_args()


def load_bedmethyl_tabix(path, chrom, start, end, min_cov=5):
    """Query a tabix-indexed bedMethyl for CpG sites in a region.

    Returns list of (position, pct_modified, coverage) tuples.
    bedMethyl columns (modkit --cpg --combine-strands --combine-mods):
      0: chrom, 1: start, 2: end, 3: mod_code, 4: score,
      5: strand, 6: thickStart, 7: thickEnd, 8: color,
      9: N_valid_cov, 10: pct_modified
    """
    import pysam

    results = []
    try:
        tbx = pysam.TabixFile(path)
        # Normalise chrom (try with and without 'chr' prefix)
        query_chrom = chrom
        if query_chrom not in tbx.contigs:
            alt = chrom.replace("chr", "") if chrom.startswith("chr") else f"chr{chrom}"
            if alt in tbx.contigs:
                query_chrom = alt
            else:
                return results

        for row in tbx.fetch(query_chrom, max(0, start), end):
            fields = row.split("\t")
            if len(fields) < 11:
                continue
            cov = int(fields[9])
            if cov < min_cov:
                continue
            pos = int(fields[1])
            pct = float(fields[10])
            results.append((pos, pct, cov))
    except Exception as e:
        log.warning("Tabix query failed for %s:%d-%d: %s", chrom, start, end, e)
    return results


def annotate_row(row, bedmethyl_path, window, min_cov):
    """Compute methylation summary for one variant row."""
    import math

    chrom = str(row.get("CHROM", row.get("chrom", "")))
    try:
        pos = int(row.get("POS", row.get("pos", 0)))
    except (ValueError, TypeError):
        return _empty_meth()

    if not chrom or pos == 0:
        return _empty_meth()

    start = pos - window
    end = pos + window

    cpgs = load_bedmethyl_tabix(bedmethyl_path, chrom, start, end, min_cov)

    if not cpgs:
        return _empty_meth()

    pcts = [pct for _, pct, _ in cpgs]
    return {
        "meth_cpg_count": len(pcts),
        "meth_mean_pct": round(sum(pcts) / len(pcts), 1),
        "meth_min_pct": round(min(pcts), 1),
        "meth_max_pct": round(max(pcts), 1),
        "meth_low_cpgs": sum(1 for p in pcts if p < 20),
        "meth_high_cpgs": sum(1 for p in pcts if p > 80),
    }


def _empty_meth():
    return {
        "meth_cpg_count": 0,
        "meth_mean_pct": float("nan"),
        "meth_min_pct": float("nan"),
        "meth_max_pct": float("nan"),
        "meth_low_cpgs": 0,
        "meth_high_cpgs": 0,
    }


def main():
    args = parse_args()

    log.info("Reading shortlist: %s", args.shortlist)
    try:
        df = pd.read_csv(args.shortlist, sep="\t", dtype=str)
    except FileNotFoundError:
        log.warning("Shortlist not found: %s. Writing empty output.", args.shortlist)
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        return

    if df.empty:
        log.info("Empty shortlist. Writing as-is.")
        df.to_csv(args.output, sep="\t", index=False)
        return

    log.info("Annotating %d variants with methylation context (window=±%d bp, min_cov=%d)...",
             len(df), args.window, args.min_coverage)

    meth_records = []
    for _, row in df.iterrows():
        meth_records.append(annotate_row(row, args.bedmethyl, args.window, args.min_coverage))

    meth_df = pd.DataFrame(meth_records)
    for col in meth_df.columns:
        df[col] = meth_df[col].values

    scored = sum(1 for m in meth_records if m["meth_cpg_count"] > 0)
    log.info("Methylation context added: %d/%d variants had CpGs in window.", scored, len(df))

    df.to_csv(args.output, sep="\t", index=False)
    log.info("Wrote: %s", args.output)


if __name__ == "__main__":
    main()
