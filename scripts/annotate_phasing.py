#!/usr/bin/env python3
"""annotate_phasing.py — Add phasing info and flag compound heterozygotes.

Reads the WhatsHap-phased VCF and shortlisted variant TSVs. For each
shortlisted variant, adds:
  phase_gt         — phased genotype (0|1 or 1|0), or unphased (0/1)
  phase_set        — PS tag value (phase block ID)
  phase_block_size — number of phased variants in this block

Then detects compound het candidates: two heterozygous shortlisted variants
in the same gene, in the same phase set, on different haplotypes (trans).

Outputs:
  - Annotated shortlist TSVs (in-place, adds phase columns)
  - Compound het report TSV

Usage:
    python annotate_phasing.py \
        --phased_vcf results/sample.phased.vcf.gz \
        --shortlist_snv report/snv/sample.snv.shortlist.tsv \
        --shortlist_indel report/indel/sample.indel.shortlist.tsv \
        --output_compound_hets report/data/sample.compound_hets.tsv

Part of SignalScope v2.0-wgs phasing module.
"""

import argparse
import gzip
import logging
import sys
from collections import defaultdict

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description="Add phasing info and flag compound hets.")
    p.add_argument("--phased_vcf", required=True,
                   help="WhatsHap phased VCF (.vcf or .vcf.gz).")
    p.add_argument("--shortlist_snv", required=True, help="SNV shortlist TSV.")
    p.add_argument("--shortlist_indel", required=True, help="Indel shortlist TSV.")
    p.add_argument("--output_snv", required=True, help="Output annotated SNV shortlist.")
    p.add_argument("--output_indel", required=True, help="Output annotated indel shortlist.")
    p.add_argument("--output_compound_hets", required=True,
                   help="Output compound het candidates TSV.")
    return p.parse_args()


def parse_phased_vcf(vcf_path):
    """Parse a WhatsHap-phased VCF and extract phasing info per variant.

    Returns dict keyed by (chrom, pos, ref, alt) ->
        {"gt": "0|1", "phase_set": "12345", "phased": True}

    Also returns phase_block_sizes: dict of phase_set -> count of variants.
    """
    variants = {}
    phase_block_sizes = defaultdict(int)

    opener = gzip.open if vcf_path.endswith(".gz") else open

    try:
        with opener(vcf_path, "rt") as f:
            sample_idx = None
            format_fields = []

            for line in f:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    cols = line.strip().split("\t")
                    # Sample is the last column (single-sample VCF)
                    sample_idx = len(cols) - 1
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alts = parts[4].split(",")
                fmt = parts[8].split(":")
                sample_data = parts[sample_idx].split(":")

                # Build format -> value mapping
                fmt_dict = {}
                for i, field in enumerate(fmt):
                    if i < len(sample_data):
                        fmt_dict[field] = sample_data[i]

                gt = fmt_dict.get("GT", "./.")
                ps = fmt_dict.get("PS", ".")
                is_phased = "|" in gt

                if is_phased and ps != ".":
                    ps_key = f"{chrom}:{ps}"
                    phase_block_sizes[ps_key] += 1

                for alt in alts:
                    key = (chrom, pos, ref, alt)
                    variants[key] = {
                        "gt": gt,
                        "phase_set": ps if is_phased else ".",
                        "phased": is_phased,
                    }
    except Exception as e:
        log.error("Failed to parse phased VCF %s: %s", vcf_path, e)
        return {}, {}

    log.info("Parsed %d variants from phased VCF, %d phase blocks.",
             len(variants), len(phase_block_sizes))
    return variants, phase_block_sizes


def normalise_chrom(chrom):
    """Normalise chromosome name (strip 'chr' prefix for matching)."""
    return chrom.replace("chr", "") if chrom.startswith("chr") else chrom


def lookup_variant(variants, chrom, pos, ref, alt):
    """Look up a variant in the phased VCF dict, trying chr/no-chr variants."""
    pos = int(pos)
    key = (chrom, pos, ref, alt)
    if key in variants:
        return variants[key]
    # Try with/without chr prefix
    alt_chrom = f"chr{chrom}" if not chrom.startswith("chr") else chrom.replace("chr", "")
    alt_key = (alt_chrom, pos, ref, alt)
    if alt_key in variants:
        return variants[alt_key]
    return None


def annotate_shortlist(df, variants, phase_block_sizes):
    """Add phase columns to a shortlist DataFrame."""
    phase_gts = []
    phase_sets = []
    block_sizes = []

    for _, row in df.iterrows():
        chrom = str(row.get("CHROM", row.get("chrom", "")))
        pos = row.get("POS", row.get("pos", ""))
        ref = str(row.get("REF", row.get("ref", "")))
        alt = str(row.get("ALT", row.get("alt", "")))

        info = lookup_variant(variants, chrom, pos, ref, alt)

        if info:
            phase_gts.append(info["gt"])
            ps = info["phase_set"]
            phase_sets.append(ps)
            if info["phased"] and ps != ".":
                ps_key = f"{normalise_chrom(chrom)}:{ps}"
                # Try both normalised and original
                bs = phase_block_sizes.get(ps_key,
                     phase_block_sizes.get(f"chr{normalise_chrom(chrom)}:{ps}",
                     phase_block_sizes.get(f"{chrom}:{ps}", 0)))
                block_sizes.append(bs)
            else:
                block_sizes.append(0)
        else:
            phase_gts.append(".")
            phase_sets.append(".")
            block_sizes.append(0)

    df["phase_gt"] = phase_gts
    df["phase_set"] = phase_sets
    df["phase_block_size"] = block_sizes

    phased_count = sum(1 for g in phase_gts if "|" in g)
    log.info("Phasing annotated: %d/%d variants phased.", phased_count, len(df))
    return df


def detect_compound_hets(snv_df, indel_df):
    """Detect compound het candidates across SNV and indel shortlists.

    A compound het is two heterozygous variants in the same gene, in the
    same phase set, on DIFFERENT haplotypes (one 0|1, the other 1|0 = trans).

    Returns a DataFrame of compound het pairs.
    """
    # Combine all phased shortlisted variants
    all_variants = []

    for source, df in [("snv", snv_df), ("indel", indel_df)]:
        if df.empty:
            continue
        for _, row in df.iterrows():
            gt = row.get("phase_gt", ".")
            ps = row.get("phase_set", ".")
            gene = row.get("gene", "")

            if not gene or "|" not in gt or ps == ".":
                continue

            chrom = str(row.get("CHROM", row.get("chrom", "")))

            all_variants.append({
                "gene": gene,
                "chrom": chrom,
                "pos": row.get("POS", row.get("pos", "")),
                "ref": row.get("REF", row.get("ref", "")),
                "alt": row.get("ALT", row.get("alt", "")),
                "variant_key": row.get("variant", row.get("variant_key", "")),
                "consequence": row.get("consequence", ""),
                "impact": row.get("impact", ""),
                "clinvar": row.get("clinvar", ""),
                "phase_gt": gt,
                "phase_set": ps,
                "type": source,
            })

    if not all_variants:
        log.info("No phased variants to check for compound hets.")
        return pd.DataFrame()

    # Group by gene + phase_set
    groups = defaultdict(list)
    for v in all_variants:
        key = (v["gene"], f"{v['chrom']}:{v['phase_set']}")
        groups[key].append(v)

    compound_hets = []
    for (gene, ps_key), variants in groups.items():
        if len(variants) < 2:
            continue

        # Split by haplotype
        hap1 = [v for v in variants if v["phase_gt"] in ("0|1", "0|2", "0|3")]
        hap2 = [v for v in variants if v["phase_gt"] in ("1|0", "2|0", "3|0")]

        if hap1 and hap2:
            # Trans configuration: variants on different haplotypes
            for v1 in hap1:
                for v2 in hap2:
                    compound_hets.append({
                        "gene": gene,
                        "variant_1": v1["variant_key"],
                        "type_1": v1["type"],
                        "consequence_1": v1["consequence"],
                        "phase_gt_1": v1["phase_gt"],
                        "variant_2": v2["variant_key"],
                        "type_2": v2["type"],
                        "consequence_2": v2["consequence"],
                        "phase_gt_2": v2["phase_gt"],
                        "phase_set": ps_key,
                        "configuration": "trans (compound het candidate)",
                    })

    if compound_hets:
        log.warning("COMPOUND HET CANDIDATES DETECTED: %d pair(s)", len(compound_hets))
        for ch in compound_hets:
            log.warning("  %s: %s (%s) + %s (%s)",
                        ch["gene"], ch["variant_1"], ch["phase_gt_1"],
                        ch["variant_2"], ch["phase_gt_2"])

    return pd.DataFrame(compound_hets)


def main():
    args = parse_args()

    # Parse phased VCF
    log.info("Parsing phased VCF: %s", args.phased_vcf)
    variants, phase_block_sizes = parse_phased_vcf(args.phased_vcf)

    # Annotate SNV shortlist
    log.info("Reading SNV shortlist: %s", args.shortlist_snv)
    try:
        snv_df = pd.read_csv(args.shortlist_snv, sep="\t", dtype=str)
    except FileNotFoundError:
        log.warning("SNV shortlist not found. Creating empty output.")
        snv_df = pd.DataFrame()

    if not snv_df.empty:
        snv_df = annotate_shortlist(snv_df, variants, phase_block_sizes)
    snv_df.to_csv(args.output_snv, sep="\t", index=False)

    # Annotate indel shortlist
    log.info("Reading indel shortlist: %s", args.shortlist_indel)
    try:
        indel_df = pd.read_csv(args.shortlist_indel, sep="\t", dtype=str)
    except FileNotFoundError:
        log.warning("Indel shortlist not found. Creating empty output.")
        indel_df = pd.DataFrame()

    if not indel_df.empty:
        indel_df = annotate_shortlist(indel_df, variants, phase_block_sizes)
    indel_df.to_csv(args.output_indel, sep="\t", index=False)

    # Detect compound hets
    log.info("Checking for compound het candidates...")
    ch_df = detect_compound_hets(snv_df, indel_df)
    ch_df.to_csv(args.output_compound_hets, sep="\t", index=False)
    log.info("Wrote compound het report: %s (%d pairs)", args.output_compound_hets, len(ch_df))

    log.info("Phasing annotation complete.")


if __name__ == "__main__":
    main()
