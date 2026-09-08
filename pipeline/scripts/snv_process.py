#!/usr/bin/env python3
"""
snv_process.py — SignalScope SNV processing (updated)

Changes vs original:
  - ClinVar rescue: pathogenic/likely_pathogenic bypass MAX_AF filter
  - MAX_AF missing defaults to 0.0 (not 1.0) — no longer drops unannotated variants
  - Location annotation from gene_features.tsv (exon/intron/flanking)
  - Intronic/flanking SNVs only shortlisted if consequence is clinically strong
  - location column added to shortlist output

Usage:
    python3 snv_process.py \
        --clair3_vcf       ... \
        --deepvariant_vcf  ... \
        --panel_bed        ... \
        --vep_tsv          ... \
        --gene_context_ref ... \
        --gene_features    ... \
        --out_dir          ... \
        --sample           barcode06 \
        --max_af           0.01
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

RELEVANCE_COL = "renal_relevance"  # overridden in main()

def _detect_relevance_col(context_path):
    """Auto-detect the relevance column name from gene_context.tsv header."""
    try:
        with open(context_path) as f:
            headers = f.readline().strip().split("\t")
            for h in headers:
                if "relevance" in h.lower():
                    return h
    except:
        pass
    return "renal_relevance"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clair3_vcf",       required=True)
    p.add_argument("--deepvariant_vcf", required=True)
    p.add_argument("--panel_bed",        required=True)
    p.add_argument("--vep_tsv",          required=True)
    p.add_argument("--deepvariant_vep",  required=False, default="")
    p.add_argument("--gene_context_ref", required=True)
    p.add_argument("--gene_features",    required=False, default="",
                   help="kidgen33.gene_features.tsv (exon/intron/flanking) — optional")
    p.add_argument("--out_dir",          required=True)
    p.add_argument("--sample",           required=True)
    p.add_argument("--max_af",           type=float, default=0.01)
    return p.parse_args()


# ---------------------------------------------------------------------------
# VCF reading helpers
# ---------------------------------------------------------------------------

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_snps_from_vcf(path):
    snps = []
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, pos, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            alt = alt.split(",")[0]
            vcf_filter = cols[6] if len(cols) > 6 else "."
            if vcf_filter not in ("PASS", "."):
                continue
            if len(ref) == 1 and len(alt) == 1 and ref != "." and alt != ".":
                snps.append((chrom, pos, ref, alt))
    return snps


def read_genotype_info_from_vcf(path):
    info = {}
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos, _, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            qual       = cols[5]
            vcf_filter = cols[6]
            fmt        = cols[8].split(":")
            sample     = cols[9].split(":")
            alt        = alt.split(",")[0]
            key        = (chrom, pos, ref, alt)
            fmt_idx    = {f: i for i, f in enumerate(fmt)}
            gt = sample[fmt_idx["GT"]] if "GT" in fmt_idx and len(sample) > fmt_idx["GT"] else ""
            dp = sample[fmt_idx["DP"]] if "DP" in fmt_idx and len(sample) > fmt_idx["DP"] else ""
            af = sample[fmt_idx["AF"]] if "AF" in fmt_idx and len(sample) > fmt_idx["AF"] else ""
            info[key] = {"gt": gt, "dp": dp, "af": af, "qual": qual, "filter": vcf_filter}
    return info


def read_panel_regions(bed_path):
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 3:
                regions.append((cols[0], int(cols[1]), int(cols[2])))
    return regions


def in_panel(chrom, pos, regions):
    p = int(pos)
    for r_chrom, r_start, r_end in regions:
        if chrom == r_chrom and r_start <= p < r_end:
            return True
    return False


# ---------------------------------------------------------------------------
# Gene features — location lookup (exon/intron/flanking)
# ---------------------------------------------------------------------------

# New 6-feature scheme (strand-aware): Exon/CDS, UTR5, UTR3, Intron, Flanking UP, Flanking Down
FEATURE_PRIORITY = {
    "Exon/CDS":      6,
    "UTR5":          5,
    "UTR3":          4,
    "Intron":        3,
    "Flanking UP":   2,
    "Flanking Down": 1,
}
# Variants in these locations are always retained (clinically relevant)
LOCATION_ALWAYS_KEEP = {"Exon/CDS", "UTR5"}


def load_gene_features(path):
    """Load kidgen33.gene_features.tsv into a list of (chrom,start,end,gene,ftype)."""
    features = []
    if not path or not os.path.isfile(path):
        return features
    with open(path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if i == 0 and parts[0] == "chrom":
                continue
            if len(parts) < 5:
                continue
            features.append((parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]))
    return features


def lookup_location(chrom, pos, features):
    """
    Return (gene, location) for the best-priority feature overlapping (chrom, pos).
    location is one of: exon, intron, flanking.
    Returns ("", "") if no feature overlaps.
    """
    pos_int = int(pos)
    best_priority = -1
    best_gene     = ""
    best_location = ""

    for fchrom, fstart, fend, gene, ftype in features:
        if chrom != fchrom:
            continue
        if not (fstart <= pos_int < fend):
            continue
        pri   = FEATURE_PRIORITY.get(ftype, 0)
        label = ftype  # use full feature_type name (Exon/CDS, UTR5, UTR3, Intron, Flanking UP, Flanking Down)
        if pri > best_priority:
            best_priority = pri
            best_gene     = gene
            best_location = label

    return best_gene, best_location


# ---------------------------------------------------------------------------
# Caller comparison
# ---------------------------------------------------------------------------

def make_key(chrom, pos, ref, alt):
    return f"{chrom}_{pos}_{ref}/{alt}"


def compare_callers(clair3_snps, dv_snps, panel_regions):
    clair3_keys   = {make_key(*s): s for s in clair3_snps}
    dv_panel = [s for s in dv_snps if in_panel(s[0], s[1], panel_regions)]
    dv_keys = {make_key(*s): s for s in dv_panel}
    all_keys = set(clair3_keys) | set(dv_keys)

    rows = []
    for key in sorted(all_keys):
        in_c = key in clair3_keys
        in_dv = key in dv_keys
        rec  = clair3_keys.get(key) or dv_keys.get(key)
        chrom, pos, ref, alt = rec
        group   = "shared" if (in_c and in_dv) else ("clair3_only" if in_c else "deepvariant_only")
        support = 2 if (in_c and in_dv) else 1
        rows.append({
            "variant_key":      key,
            "CHROM":            chrom,
            "POS":              pos,
            "REF":              ref,
            "ALT":              alt,
            "caller_clair3":    "1" if in_c else "0",
            "caller_deepvariant": "1" if in_dv else "0",
            "caller_support":   str(support),
            "comparison_group": group,
        })
    return rows


# ---------------------------------------------------------------------------
# VEP TSV parsing
# ---------------------------------------------------------------------------

def parse_vep_tsv(path):
    records = defaultdict(list)
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                continue
            if header is None:
                continue
            cols = line.strip().split("\t")
            if len(cols) < len(header):
                cols += [""] * (len(header) - len(cols))
            row = dict(zip(header, cols))
            uv  = row.get("Uploaded_variation", row.get("#Uploaded_variation", ""))
            m   = re.match(r"^(.+)_(\d+)_(.+)/(.+)$", uv)
            if not m:
                continue
            key = (m.group(1), m.group(2), m.group(3), m.group(4))
            records[key].append(row)

    # Compute spliceai_max across ALL transcript rows (not just canonical)
    spliceai_max_by_key = {}
    def _spliceai_max_from_rows(rows):
        best_score = None
        for r in rows:
            raw = r.get("SpliceAI_pred", "")
            if not raw or raw in (".", "-"):
                continue
            # Format: SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
            parts = raw.split("|")
            if len(parts) < 5:
                continue
            try:
                scores = [float(parts[i]) for i in (1, 2, 3, 4)]
                m = max(scores)
                if best_score is None or m > best_score:
                    best_score = m
            except (ValueError, IndexError):
                continue
        return f"{best_score:.4f}" if best_score is not None else ""
    for key, rows in records.items():
        spliceai_max_by_key[key] = _spliceai_max_from_rows(rows)

    best = {}
    for key, rows in records.items():
        chosen = None
        for r in rows:
            if r.get("MANE_SELECT", "-") not in ("-", ""):
                chosen = r
                break
        if chosen is None:
            for r in rows:
                if r.get("CANONICAL", "-") == "YES":
                    chosen = r
                    break
        if chosen is None and rows:
            chosen = rows[0]
        if chosen:
            chosen["spliceai_max"] = spliceai_max_by_key.get(key, "")
            best[key] = chosen
    return best


KEEP_VEP_COLS = [
    "SYMBOL", "Feature", "Consequence", "IMPACT", "BIOTYPE",
    "CANONICAL", "MANE_SELECT", "Existing_variation",
    "HGVSc", "HGVSp",
    "SIFT", "PolyPhen", "REVEL", "EXON", "INTRON",
    "gnomADe_AF", "gnomADg_AF", "MAX_AF", "MAX_AF_POPS", "CLIN_SIG",
    "FILTER",
    "SpliceAI_pred",
]


def join_vep(merged_rows, vep_best):
    annotated    = []
    unannotated  = 0
    for row in merged_rows:
        key = (row["CHROM"], row["POS"], row["REF"], row["ALT"])
        vep = vep_best.get(key, {})
        if not vep:
            unannotated += 1
        new_row = dict(row)
        for col in KEEP_VEP_COLS:
            new_row[col] = vep.get(col, "")
        # spliceai_max is a derived field attached to the vep row — copy explicitly
        new_row["spliceai_max"] = vep.get("spliceai_max", "")
        annotated.append(new_row)
    if unannotated:
        print(f"  [WARN] {unannotated} variants had no VEP annotation", file=sys.stderr)
    return annotated


# ---------------------------------------------------------------------------
# ClinVar normalisation
# ---------------------------------------------------------------------------

CLINVAR_MAP = {
    "pathogenic":             "pathogenic",
    "likely_pathogenic":      "likely_pathogenic",
    "uncertain_significance": "VUS",
    "benign":                 "benign",
    "likely_benign":          "benign",
}

CLINVAR_RESCUE = {"pathogenic", "likely_pathogenic"}

# Panel genes — only variants in these genes are shortlisted
# Panel genes — loaded dynamically from BED at runtime (see main())
# Populated by _load_panel_genes() before shortlisting
PANEL_GENES = set()

def _load_panel_genes(bed_path):
    """Load gene names from column 4 of the panel BED file."""
    global PANEL_GENES
    genes = set()
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and not line.startswith("#"):
                genes.add(parts[3].strip())
    PANEL_GENES = genes
    print(f"  Panel genes loaded from BED: {len(PANEL_GENES)} genes", __import__("sys").stderr)

# Consequences strong enough to shortlist intronic/flanking variants
STRONG_CONSEQUENCES = {
    "splice_donor_variant",
    "splice_acceptor_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "stop_gained",
    "stop_lost",
    "frameshift_variant",
    "start_lost",
    "transcript_ablation",
}



# =============================================================================
# =============================================================================
# ACMG-lite classification
# Simplified rule-based ACMG classification using available VEP annotations.
# Evidence criteria and thresholds follow ClinGen SVI recommendations:
#   - PP3/BP4 REVEL thresholds: Pejaver et al. 2022 (PMID: 36413997)
#   - SpliceAI thresholds: Walker et al. 2023 (PMID: 37352859)
#   - Point system: Tavtigian et al. 2020 (PMID: 32720330)
# Not a substitute for full InterVar/ACMG classification but provides
# evidence-based pathogenicity assessment for clinical prioritisation.
# =============================================================================
PVS1_CONSEQUENCES = {
    "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
    "splice_donor_variant", "splice_acceptor_variant",
    "transcript_ablation"
}

def acmg_lite(r):
    """
    Returns (classification, evidence_codes) tuple.
    classification: Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign

    Scoring follows the Tavtigian et al. 2020 Bayesian point system:
      Supporting = ±1, Moderate = ±2, Strong = ±4, Very Strong = ±8

    PP3/BP4 use REVEL as the sole calibrated missense predictor per ClinGen
    SVI (Pejaver et al. 2022). SIFT and PolyPhen are displayed in the report
    but excluded from scoring (they did not meet PP3 Supporting thresholds in
    the ClinGen calibration). SpliceAI is used for PP3 on non-PVS1 variants
    (Walker et al. 2023).
    """
    consequences = set(r.get("Consequence", "").split(","))
    impact       = r.get("IMPACT", "")

    try:
        max_af = float(r.get("MAX_AF", "") or 0.0)
    except:
        max_af = 0.0

    try:
        revel = float(r.get("REVEL", "") or -1)
    except:
        revel = -1

    try:
        spliceai = float(r.get("spliceai_max", "") or -1)
    except:
        spliceai = -1

    evidence = []
    path_score   = 0  # positive = pathogenic evidence
    benign_score = 0  # positive = benign evidence

    # ── Pathogenic criteria ──────────────────────────────────────────────────
    # PVS1: null variant in gene where LOF is disease mechanism (Very Strong, +8)
    pvs1_applied = False
    if consequences & PVS1_CONSEQUENCES:
        evidence.append("PVS1")
        path_score += 8
        pvs1_applied = True

    # PM2: absent/rare in population databases (Moderate, +2)
    if max_af < 0.001:
        evidence.append("PM2")
        path_score += 2

    # PP3: computational evidence supports deleterious effect
    # REVEL thresholds calibrated by ClinGen SVI (Pejaver et al. 2022, Table 2)
    if revel >= 0.932:
        evidence.append("PP3_REVEL")
        path_score += 4   # Strong
    elif revel >= 0.773:
        evidence.append("PP3_REVEL")
        path_score += 2   # Moderate
    elif revel >= 0.644:
        evidence.append("PP3_REVEL")
        path_score += 1   # Supporting

    # SpliceAI PP3: only when PVS1 is NOT applied, to avoid double-counting
    # splice evidence (Walker et al. 2023, ClinGen SVI Splicing Subgroup)
    if not pvs1_applied and spliceai >= 0.2:
        evidence.append("PP3_SpliceAI")
        path_score += 1   # Supporting

    # ── Benign criteria ──────────────────────────────────────────────────────
    # BA1: allele frequency >5% in population databases (Standalone Benign, +8)
    if max_af > 0.05:
        evidence.append("BA1")
        benign_score += 8

    # BP4: computational evidence supports benign effect
    # REVEL thresholds calibrated by ClinGen SVI (Pejaver et al. 2022, Table 2)
    if revel >= 0 and revel <= 0.183:
        evidence.append("BP4_REVEL")
        benign_score += 2   # Moderate
    elif revel >= 0 and revel <= 0.290:
        evidence.append("BP4_REVEL")
        benign_score += 1   # Supporting

    # BP7: synonymous variant with no predicted splice impact
    if "synonymous_variant" in consequences and "splice" not in " ".join(consequences):
        evidence.append("BP7")
        benign_score += 1

    # ── Classification ───────────────────────────────────────────────────────
    net = path_score - benign_score

    if benign_score >= 8:
        classification = "Benign"
    elif net >= 10:
        classification = "Pathogenic"
    elif net >= 6:
        classification = "Likely Pathogenic"
    elif benign_score >= 3:
        classification = "Likely Benign"
    elif net >= 2:
        classification = "VUS"
    else:
        classification = "VUS"

    # NOTE: ACMG-lite is intentionally INDEPENDENT of ClinVar. ClinVar's aggregate
    # classification is reported separately (its own column) as a distinct opinion.
    # Per ClinGen SVI guidance (retirement of PP5/BP6), a database verdict is not
    # used as ACMG evidence here; this classification rests on primary evidence
    # (consequence, population frequency, computational predictors) only.

    return classification, "|".join(evidence)
def normalize_clinvar(raw):
    if not raw or raw in ("-", "."):
        return "no_data"
    r = raw.lower().replace(" ", "_")
    # Detect conflicting submissions: both pathogenic-side AND benign-side present.
    # (e.g. "benign,pathogenic" — submitters disagree; not a clean pathogenic call.)
    has_path   = ("pathogenic" in r)          # covers pathogenic + likely_pathogenic
    has_benign = ("benign" in r)              # covers benign + likely_benign
    if has_path and has_benign:
        return "conflicting"
    # Otherwise map by priority: pathogenic > likely_pathogenic > VUS > benign.
    if "likely_pathogenic" in r:
        return "likely_pathogenic"
    if "pathogenic" in r:
        return "pathogenic"
    if "uncertain_significance" in r:
        return "VUS"
    if "benign" in r:
        return "benign"
    return "no_data"


# ---------------------------------------------------------------------------
# AF helper — FIXED: default 0.0 not 1.0
# ---------------------------------------------------------------------------

def _af(r):
    """
    Return MAX_AF as float.
    Missing/empty AF defaults to 0.0 (rare assumed) not 1.0.
    This was the root cause of pathogenic variants being dropped.
    """
    v = r.get("MAX_AF", "")
    try:
        return float(v)
    except (ValueError, TypeError):
        return 0.0


def is_clinvar_rescue(r):
    """True if variant is pathogenic/likely_pathogenic — bypasses AF filter."""
    return normalize_clinvar(r.get("CLIN_SIG", "")) in CLINVAR_RESCUE


def passes_location_filter(r, location):
    """
    Variants in Exon/CDS or UTR5 are always retained.
    Variants in UTR3, Intron, Flanking UP, Flanking Down are retained only
    if their consequence is clinically strong.
    Variants with empty location (not in any feature) pass (defensive).
    """
    if location == "" or location in LOCATION_ALWAYS_KEEP:
        return True
    consequences = set(r.get("Consequence", "").split(","))
    return bool(consequences & STRONG_CONSEQUENCES)


# ---------------------------------------------------------------------------
# Write helpers
# ---------------------------------------------------------------------------

def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    print(f"  Written: {path} ({len(rows)} rows)", file=sys.stderr)


def write_kv_tsv(path, pairs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("metric\tcount\n")
        for k, v in pairs:
            fh.write(f"{k}\t{v}\n")
    print(f"  Written: {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    _load_panel_genes(args.panel_bed)
    global RELEVANCE_COL
    RELEVANCE_COL = _detect_relevance_col(args.gene_context_ref)
    print(f"  Relevance column: {RELEVANCE_COL}", __import__("sys").stderr)
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"[snv_process] sample={args.sample}", file=sys.stderr)

    # --- 1. Load gene features (optional) ------------------------------------
    features = load_gene_features(args.gene_features)
    if features:
        print(f"  Gene features loaded: {len(features)} intervals", file=sys.stderr)
    else:
        print("  Gene features: not provided — location column will be empty",
              file=sys.stderr)

    # --- 2. Read SNPs + genotype info ----------------------------------------
    print("  Reading Clair3 SNPs...", file=sys.stderr)
    clair3_snps = read_snps_from_vcf(args.clair3_vcf)
    print(f"  Clair3 SNPs: {len(clair3_snps)}", file=sys.stderr)

    geno_info = read_genotype_info_from_vcf(args.clair3_vcf)
    print(f"  Genotype records read: {len(geno_info)}", file=sys.stderr)

    print("  Reading DeepVariant SNPs...", file=sys.stderr)
    dv_snps = read_snps_from_vcf(args.deepvariant_vcf)
    print(f"  DeepVariant SNPs (raw): {len(dv_snps)}", file=sys.stderr)

    panel_regions = read_panel_regions(args.panel_bed)

    # --- 3. Compare callers --------------------------------------------------
    merged_rows = compare_callers(clair3_snps, dv_snps, panel_regions)
    print(f"  Merged SNV set: {len(merged_rows)}", file=sys.stderr)

    merged_cols = [
        "variant_key", "CHROM", "POS", "REF", "ALT",
        "caller_clair3", "caller_deepvariant", "caller_support", "comparison_group",
    ]
    write_tsv(
        os.path.join(args.out_dir, f"{args.sample}.snv.merged.tsv"),
        merged_rows, merged_cols
    )

    # --- 4. Parse VEP and join -----------------------------------------------
    print("  Parsing VEP TSV...", file=sys.stderr)
    vep_best = parse_vep_tsv(args.vep_tsv)
    if args.deepvariant_vep:
        _dv_vep = parse_vep_tsv(args.deepvariant_vep)
        _added = 0
        for _k, _v in _dv_vep.items():
            if _k not in vep_best:
                vep_best[_k] = _v
                _added += 1
        print(f"  + DeepVariant VEP: {len(_dv_vep)} keys, {_added} new DV-only", file=sys.stderr)
    print(f"  VEP annotated variants: {len(vep_best)}", file=sys.stderr)

    annotated = join_vep(merged_rows, vep_best)
    _unann = sum(1 for _r in annotated if not _r.get("IMPACT"))
    if merged_rows and _unann / len(merged_rows) > 0.05:
        sys.exit(f"[FATAL] {_unann}/{len(merged_rows)} SNVs unannotated (>5%) - check VEP inputs")

    ann_cols = merged_cols + KEEP_VEP_COLS
    write_tsv(
        os.path.join(args.out_dir, f"{args.sample}.snv.annotated.tsv"),
        annotated, ann_cols
    )

    # --- 5. Summary counts ---------------------------------------------------
    clair3_total   = sum(1 for r in merged_rows if r["caller_clair3"]   == "1")
    dv_total = sum(1 for r in merged_rows if r["caller_deepvariant"] == "1")
    shared         = sum(1 for r in merged_rows if r["comparison_group"] == "shared")
    clair3_only    = sum(1 for r in merged_rows if r["comparison_group"] == "clair3_only")
    dv_only  = sum(1 for r in merged_rows if r["comparison_group"] == "deepvariant_only")

    impact_filtered = [r for r in annotated if r.get("IMPACT", "") in ("HIGH", "MODERATE")]
    impact_rare     = impact_filtered  # AF filter removed — all HIGH/MODERATE shown

    write_kv_tsv(
        os.path.join(args.out_dir, f"{args.sample}.snv.summary.tsv"),
        [
            ("clair3_total",    clair3_total),
            ("deepvariant_total", dv_total),
            ("shared",          shared),
            ("clair3_only",     clair3_only),
            ("deepvariant_only",  dv_only),
            ("impact_high_moderate", len(impact_filtered)),
        ]
    )

    # --- 6. ClinVar summary --------------------------------------------------
    clinvar_counts = defaultdict(int)
    for r in annotated:
        norm = normalize_clinvar(r.get("CLIN_SIG", ""))
        clinvar_counts[norm] += 1

    order = ["pathogenic", "likely_pathogenic", "VUS", "conflicting", "benign", "no_data"]
    write_tsv(
        os.path.join(args.out_dir, f"{args.sample}.snv.clinvar_summary.normalized.tsv"),
        [{"clinvar_class": k, "count": clinvar_counts.get(k, 0)} for k in order],
        ["clinvar_class", "count"]
    )

    # --- 7. Shortlist — FIXED filtering logic --------------------------------
    shortlist_display = []
    rescued           = 0
    location_filtered = 0

    for r in annotated:
        impact   = r.get("IMPACT", "")
        clinvar  = normalize_clinvar(r.get("CLIN_SIG", ""))
        af_val   = _af(r)
        rescue   = is_clinvar_rescue(r)

        # Must be in panel genes
        if r.get("SYMBOL", "") not in PANEL_GENES:
            continue
        # Must be HIGH or MODERATE impact
        if impact not in ("HIGH", "MODERATE"):
            continue

        # Compute ACMG-lite early so the shortlist decision can use it.
        acmg_class, acmg_ev = acmg_lite(r)
        r["acmg_class"]    = acmg_class
        r["acmg_evidence"] = acmg_ev
        # Exclude ONLY when BOTH ClinVar AND ACMG-lite agree the variant is benign.
        # If they disagree (e.g. ClinVar pathogenic but ACMG benign by high frequency),
        # keep the variant so a human can adjudicate the discordance.
        clinvar_benign = clinvar in ("benign", "likely_benign")
        acmg_benign    = acmg_class in ("Benign", "Likely Benign")
        if clinvar_benign and acmg_benign:
            continue

        # Population AF filter — exclude common variants (>1%) unless ClinVar P/LP
        if af_val > 0.01 and not rescue:
            continue

        # ClinVar rescue tracking
        if rescue:
            rescued += 1

        # Location lookup
        _, location = lookup_location(r["CHROM"], r["POS"], features)

        # Intronic/flanking: only keep if strong consequence
        if not passes_location_filter(r, location):
            location_filtered += 1
            continue

        # Genotype lookup
        vkey = (r.get("CHROM",""), r.get("POS",""), r.get("REF",""), r.get("ALT",""))
        geno = geno_info.get(vkey, {})

        exon   = r.get("EXON", "")
        intron = r.get("INTRON", "")
        exon_intron = (f"exon {exon}"     if exon   and exon   not in ("-", "", ".")
                  else f"intron {intron}" if intron and intron not in ("-", "", ".")
                  else "")

        shortlist_display.append({
            "gene":              r.get("SYMBOL", ""),
            "variant":           r.get("variant_key", ""),
            "varsome_url":       (f"https://varsome.com/variant/hg38/{r.get('CHROM','').replace('chr','')}-{r.get('POS','')}-{r.get('REF','')}-{r.get('ALT','')}"
                                  if len(r.get('REF','')) <= 50 and len(r.get('ALT','')) <= 50
                                  else ""),
            "novel_candidate":   "yes" if (normalize_clinvar(r.get("CLIN_SIG","")) == "no_data" and r.get("IMPACT","") in ("HIGH","MODERATE")) else "no",
            "consequence":       r.get("Consequence", ""),
            "impact":            impact,
            "clinvar":           clinvar,
            "location":          location,
            "caller_support":    r.get("caller_support", ""),
            "comparison_group":  r.get("comparison_group", ""),
            "gt":                geno.get("gt", ""),
            "depth":             geno.get("dp", ""),
            "af":                geno.get("af", ""),
            "qual":              geno.get("qual", ""),
            "filter":            geno.get("filter", ""),
            "max_af":            r.get("MAX_AF", ""),
            "gnomADe_AF":        r.get("gnomADe_AF", ""),
            "gnomADg_AF":        r.get("gnomADg_AF", ""),
            "HGVSc":             r.get("HGVSc", ""),
            "HGVSp":             r.get("HGVSp", ""),
            "exon_intron":       exon_intron,
            "SIFT":              r.get("SIFT", ""),
            "PolyPhen":          r.get("PolyPhen", ""),
            "REVEL":             r.get("REVEL", ""),
            "spliceai_max":      r.get("spliceai_max", ""),
            "acmg_class":        r.get("acmg_class", ""),
            "acmg_evidence":     r.get("acmg_evidence", ""),
            "Existing_variation": r.get("Existing_variation", ""),
            "Feature":           r.get("Feature", ""),
            "BIOTYPE":           r.get("BIOTYPE", ""),
        })

    print(f"  Shortlist variants: {len(shortlist_display)}", file=sys.stderr)
    if rescued:
        print(f"  ClinVar rescued (pathogenic/LP bypassed AF filter): {rescued}",
              file=sys.stderr)
    if location_filtered:
        print(f"  Location filtered (intronic/flanking, weak consequence): {location_filtered}",
              file=sys.stderr)

    _shortlist_cols = [
        "gene", "variant", "varsome_url", "novel_candidate", "consequence", "impact", "clinvar", "location",
        "caller_support", "comparison_group",
        "gt", "depth", "af", "qual", "filter",
        "max_af", "gnomADe_AF", "gnomADg_AF",
        "HGVSc", "HGVSp", "exon_intron",
        "SIFT", "PolyPhen", "REVEL", "spliceai_max", "acmg_class", "acmg_evidence", "Existing_variation",
        "Feature", "BIOTYPE",
    ]

    write_tsv(
        os.path.join(args.out_dir, f"{args.sample}.snv.shortlist.tsv"),
        shortlist_display, _shortlist_cols
    )

    # --- 8. Join gene context ------------------------------------------------
    CONTEXT_FILE_COLS = ["gene", "disease_summary", RELEVANCE_COL,
                         "inheritance", "phenotype_group", "omim_gene_id"]
    gene_context = {}
    try:
        with open(args.gene_context_ref) as fh:
            for i, line in enumerate(fh):
                if i == 0:
                    continue
                cols = line.strip().split("\t")
                if not cols or not cols[0].strip():
                    continue
                padded = cols + [""] * len(CONTEXT_FILE_COLS)
                rec    = dict(zip(CONTEXT_FILE_COLS, padded))
                gene   = rec.get("gene", "").strip()
                if gene:
                    gene_context[gene] = rec
        print(f"  Gene context loaded: {len(gene_context)} genes", file=sys.stderr)
    except FileNotFoundError:
        print(f"  [WARN] gene_context_ref not found: {args.gene_context_ref}",
              file=sys.stderr)

    context_cols = ["disease_summary", RELEVANCE_COL, "inheritance", "phenotype_group", "omim_gene_id"]
    final_cols   = _shortlist_cols + context_cols

    final_rows = []
    for r in shortlist_display:
        gene = r.get("gene", "")
        ctx  = gene_context.get(gene, {})
        row  = dict(r)
        for col in context_cols:
            row[col] = ctx.get(col, "")
        final_rows.append(row)

    write_tsv(
        os.path.join(args.out_dir,
                     f"{args.sample}.snv.shortlist.with_gene_context.header.tsv"),
        final_rows, final_cols
    )

    print(f"[snv_process] complete — {len(shortlist_display)} shortlisted variants",
          file=sys.stderr)


if __name__ == "__main__":
    main()
