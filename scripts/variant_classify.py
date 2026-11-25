#!/usr/bin/env python3
"""
variant_classify.py - SignalScope shared ACMG-lite classification and tiering.

Imported by snv_process.py and indel_process.py (same scripts/ folder).
Python 3.8 compatible. No third-party dependencies.

Design rules (Tavtigian et al. 2020 points):
  1. The ACMG-lite label comes ONLY from the point total:
        P >= 10 | LP 6..9 | VUS 0..5 | LB -1..-6 | B <= -7
     ClinVar is reported separately and never changes the label.
  2. BA1 (non-bottleneck popmax AF > 5%) is stand-alone Benign, except for
     variants on the ClinGen SVI BA1 exception list (assets/ba1_exceptions.tsv).
  3. In-silico evidence is counted ONCE per variant:
        missense            -> REVEL (Pejaver et al. 2022 calibrated bands)
        non-missense/splice -> SpliceAI (Walker et al. 2023)
     SIFT and PolyPhen are display-only. Genes with a ClinGen VCEP (optional
     list) are capped at Supporting (+1 / -1), because most VCEPs calibrate
     their own tool and do not allow more than Supporting.
  4. PVS1 only for null variants in genes where loss of function is an
     established disease mechanism (ClinGen HI 3/30, or recessive panel genes).
     Simplified, conservative Abou Tayoun 2018 tree (PMID 30192042):
        nonsense / frameshift, NMD predicted          PVS1 (+8)
        nonsense / frameshift, NMD escape             PVS1_Moderate (+2)  [review]
            (VEP NMD plugin: last exon or last 50 bp of the penultimate exon;
             last exon is also used when the NMD column is absent)
        canonical splice +-1/2                        PVS1_Strong (+4)    [review]
        start lost                                    PVS1_Moderate (+2)
        whole-transcript deletion                     PVS1 (+8)
     Gene-disease validity not Definitive/Strong (ClinGen) -> down 2 levels
     (+8 -> +2, +4 -> +1, +2 -> 0). Without a ClinGen validity file, PanelApp
     green is used as the proxy for an established gene-disease relationship.
     [review] means the full tree (truncated region, exon skipping frame,
     alternative transcripts) must be checked by the curator.
  5. PS1 (+4) and PM5 (+2) are applied when the caller passes
     row["codon_evidence"] = "PS1" or "PM5": a DIFFERENT nucleotide change at the
     same codon is ClinVar P/LP (rescue-level review) and gives the same amino
     acid (PS1) or another missense (PM5). The lookup lives in
     classify_small_variants.py. The variant's own ClinVar record would be PP5,
     which ClinGen SVI retired; it is not scored.
  6. ClinVar may only RESCUE a variant past the AF filter when it is a
     germline Pathogenic / Likely pathogenic, allele-matched record from the
     ClinVar VCF with >= 2 review stars (or 1 star and no conflict).
     The rsID-level CLIN_SIG column from the VEP cache is display-only.

Every criterion that fires is recorded with its points in `acmg_evidence`
(e.g. "PVS1(+8)|PM2(+2)") so the label is always reproducible by hand.
"""

import os
import re
import sys

# ---------------------------------------------------------------------------
# Defaults (override by passing a dict to Classifier(cfg=...))
# ---------------------------------------------------------------------------
DEFAULTS = {
    "ba1_af": 0.05,              # stand-alone benign
    "bs1_enabled": True,
    "bs1_af_dominant": 0.005,    # AD / XL / mixed genes (strong benign, -4)
    "bs1_af_recessive": 0.02,    # AR-only genes
    "pm2_af_dominant": 0.0001,   # PM2: absent or extremely rare
    "pm2_af_recessive": 0.001,
    "pm2_points": 1,             # PM2_Supporting, ClinGen SVI 2020
    "spliceai_pp3": 0.2,
    "spliceai_bp4": 0.1,
    "pm4_repeat_min_copies": 3,  # PM4 not applied when the indel unit is in a tandem repeat of >= N copies
}

# Pejaver et al. 2022 (PMID 36413997) REVEL bands, verified against
# PMC13282463 Table 1. Intervals: PP3 lower bound inclusive; BP4 upper inclusive.
REVEL_PP3 = [(0.932, 4, "PP3_Strong"), (0.773, 2, "PP3_Moderate"), (0.644, 1, "PP3")]
REVEL_BP4 = [(0.003, -8, None),          # very strong band: capped at Strong below
             (0.016, -4, "BP4_Strong"),
             (0.183, -2, "BP4_Moderate"),
             (0.290, -1, "BP4")]

NULL_CONSEQUENCES = {
    "stop_gained", "frameshift_variant",
    "splice_donor_variant", "splice_acceptor_variant",
    "transcript_ablation",
}
PM4_CONSEQUENCES = {"inframe_insertion", "inframe_deletion", "stop_lost"}
SPLICE_CONSEQUENCES = {
    "splice_region_variant", "splice_donor_5th_base_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
}

# gnomAD non-bottleneck populations (Gorzynski et al. 2022 used the same idea).
# Excluded: ASJ, FIN, AMI, MID, REMAINING (founder / small / mixed groups).
POPMAX_COLUMNS = [
    "gnomADe_AFR_AF", "gnomADe_AMR_AF", "gnomADe_EAS_AF", "gnomADe_NFE_AF", "gnomADe_SAS_AF",
    "gnomADg_AFR_AF", "gnomADg_AMR_AF", "gnomADg_EAS_AF", "gnomADg_NFE_AF", "gnomADg_SAS_AF",
]
GLOBAL_AF_COLUMNS = ["gnomADe_AF", "gnomADg_AF"]

# HLA / MHC: hypervariable, poorly represented in gnomAD, not ACMG-interpretable.
# GRCh38 MHC region (GRC definition, chr6:28,510,120-33,480,577).
MHC_CHROM = {"6", "chr6"}
MHC_START, MHC_END = 28510120, 33480577
HLA_GENE_RE = re.compile(r"^HLA-")

LABEL_ORDER = {"Pathogenic": 5, "Likely Pathogenic": 4, "VUS": 3,
               "Likely Benign": 2, "Benign": 1}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def to_float(v):
    if v is None:
        return None
    s = str(v).strip()
    if s in ("", "-", ".", "NA", "nan"):
        return None
    # VEP may write "0.01&0.02" for multiple co-located records: take the max
    try:
        return max(float(x) for x in re.split(r"[&,]", s) if x not in ("", "-", "."))
    except ValueError:
        return None


_INTRON_OFF_RE = re.compile(r"c\.[*-]?\d+([+-])(\d+)")


def intron_offset(hgvsc):
    """('+', 5) for c.123+5G>A, ('-', 19) for c.124-19A>G (first position of a range);
    None when HGVSc gives no intronic offset (exonic, UTR without offset, or empty)."""
    m = _INTRON_OFF_RE.search(hgvsc or "")
    return (m.group(1), int(m.group(2))) if m else None


def popmax_af(row):
    """Frequency used for BA1 / BS1 / PM2. If the caller supplies row['freq_override'] (e.g. the
    gnomAD GroupMax FAF95 from a direct allele-exact lookup) it is used as is.
    Otherwise: non-bottleneck popmax from gnomAD exome+genome population columns.
    Falls back to gnomAD global AF, then to VEP MAX_AF (which also includes
    1000 Genomes and bottleneck populations) when population columns are absent.
    Returns (af_or_None, source_label)."""
    if "freq_override" in row:
        return row["freq_override"], row.get("freq_override_src", "override")
    vals = [to_float(row.get(c)) for c in POPMAX_COLUMNS]
    vals = [v for v in vals if v is not None]
    if vals:
        return max(vals), "gnomAD_popmax_nonbottleneck"
    vals = [to_float(row.get(c)) for c in GLOBAL_AF_COLUMNS]
    vals = [v for v in vals if v is not None]
    if vals:
        return max(vals), "gnomAD_global"
    v = to_float(row.get("MAX_AF"))
    if v is not None:
        return v, "VEP_MAX_AF"
    return None, "absent"


def consequences(row):
    raw = row.get("Consequence", "") or ""
    return set(c for c in re.split(r"[,&]", raw) if c)


def is_hla_mhc(gene, chrom, pos):
    if gene and HLA_GENE_RE.match(gene):
        return True
    try:
        p = int(pos)
    except (TypeError, ValueError):
        return False
    return str(chrom) in MHC_CHROM and MHC_START <= p <= MHC_END


def is_last_exon(row):
    m = re.match(r"^(\d+)/(\d+)$", str(row.get("EXON", "") or ""))
    return bool(m) and m.group(1) == m.group(2)


def zygosity(gt):
    g = (gt or "").replace("|", "/")
    if g in ("", ".", "./."):
        return "unknown"
    parts = [p for p in g.split("/") if p not in ("",)]
    if len(parts) == 1:
        return "hemi" if parts[0] not in ("0", ".") else "unknown"
    if len(parts) == 2 and parts[0] == parts[1] and parts[0] not in ("0", "."):
        return "hom"
    if "0" in parts:
        return "het"
    if len(set(parts)) == 2:
        return "het_multi"          # e.g. 1/2 at a multi-allelic site
    return "unknown"


MOI_DOMINANT = {"AD", "XLD", "SD"}
MOI_RECESSIVE = {"AR", "XLR"}
_MOI_PHRASES = [("x-linked recessive", " XLR "), ("x-linked dominant", " XLD "), ("x-linked", " XL "),
                ("autosomal dominant", " AD "), ("autosomal recessive", " AR "), ("semidominant", " SD "),
                ("semi-dominant", " SD "), ("dominant", " AD "), ("recessive", " AR ")]


def moi_tokens(inh):
    """Inheritance text -> set of codes (AD, AR, XLD, XLR, XL, SD, MT, IMPRINTED ...).
    Codes are matched as whole tokens, so 'XLR' is not read as 'XL'."""
    t = " %s " % (inh or "").lower()
    for phrase, code in _MOI_PHRASES:
        t = t.replace(phrase, code)
    return {x for x in re.split(r"[^A-Za-z]+", t.upper()) if x}


def inheritance_class(inh):
    """'recessive' only when the gene is AR / XLR without any dominant mode; 'dominant' when any
    dominant mode (AD, XLD, SD) is listed; 'unknown' otherwise (plain 'XL', MT, imprinted, empty):
    unknown genes are tested against both models."""
    t = moi_tokens(inh)
    if t & MOI_DOMINANT:
        return "dominant"
    if t & MOI_RECESSIVE:
        return "recessive"
    return "unknown"


def moi_models(inh):
    """(dominant model allowed, recessive model allowed) for the inheritance fit.
    'XL' without R/D, imprinted, MT or empty -> both models (not decided by the table)."""
    t = moi_tokens(inh)
    dom, rec = bool(t & MOI_DOMINANT), bool(t & MOI_RECESSIVE)
    if not dom and not rec:
        return True, True
    return dom, rec


# ---------------------------------------------------------------------------
# Reference lists
# ---------------------------------------------------------------------------
def load_gene_list(path):
    genes = set()
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    genes.add(line.split("\t")[0].strip())
    return genes


def load_ba1_exceptions(path):
    """TSV with columns gene, hgvsc (c. notation, transcript-free), label."""
    exc = set()
    if path and os.path.isfile(path):
        with open(path) as fh:
            for i, line in enumerate(fh):
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if cols[0] == "gene":
                    continue
                if len(cols) >= 2:
                    exc.add((cols[0].strip(), cols[1].strip()))
    return exc


def hgvsc_short(hgvsc):
    """'ENST00000357654.9:c.68_69del' -> 'c.68_69del'"""
    s = hgvsc or ""
    return s.split(":", 1)[1] if ":" in s else s


# ---------------------------------------------------------------------------
# ClinVar
# ---------------------------------------------------------------------------
CLNREVSTAT_STARS = [
    ("practice_guideline", 4),
    ("reviewed_by_expert_panel", 3),
    ("criteria_provided,_multiple_submitters,_no_conflicts", 2),
    ("criteria_provided,_conflicting_classifications", 1),
    ("criteria_provided,_conflicting_interpretations", 1),
    ("criteria_provided,_single_submitter", 1),
]


def clinvar_stars(clnrevstat):
    s = (clnrevstat or "").lower().replace(" ", "_")
    for key, stars in CLNREVSTAT_STARS:
        if key in s:
            return stars
    return 0


def clinvar_germline_class(clnsig):
    """Normalise ClinVar VCF CLNSIG (germline aggregate). Secondary terms after
    '|' or ',_' (risk_factor, drug_response ...) are ignored for the class."""
    s = (clnsig or "").strip()
    if not s or s in (".", "-"):
        return "no_data"
    primary = re.split(r"[|]|,_", s)[0].lower()
    if "conflicting" in primary:
        return "conflicting"
    if primary in ("pathogenic", "pathogenic/likely_pathogenic",
                   "pathogenic/likely_pathogenic/likely_risk_allele",
                   "pathogenic/likely_risk_allele"):
        return "pathogenic"
    if primary.startswith("pathogenic") and "low_penetrance" in primary:
        return "pathogenic_low_penetrance"
    if primary.startswith("likely_pathogenic"):
        return "likely_pathogenic"
    if primary.startswith("uncertain"):
        return "VUS"
    if primary in ("benign", "benign/likely_benign"):
        return "benign"
    if primary == "likely_benign":
        return "likely_benign"
    return "other"          # drug_response, risk_factor, association, not_provided ...


def clinvar_record(clnsig, clnrevstat, source="clinvar_vcf"):
    cls = clinvar_germline_class(clnsig)
    stars = clinvar_stars(clnrevstat)
    plp = cls in ("pathogenic", "likely_pathogenic", "pathogenic_low_penetrance")
    single_no_conflict = (stars == 1 and "conflicting" not in (clnrevstat or "").lower())
    rescue = (source == "clinvar_vcf") and plp and (stars >= 2 or single_no_conflict)
    return {"clinvar": cls, "clinvar_stars": stars, "clinvar_source": source,
            "clinvar_rescue": rescue, "clinvar_raw": clnsig or ""}


def clinvar_from_vep_rsid(clin_sig):
    """VEP cache CLIN_SIG is rsID-level, unstarred, mixes alleles and somatic.
    Display-only: never rescues, shown with source 'vep_rsid_level'."""
    s = (clin_sig or "").lower()
    if not s or s in ("-", "."):
        cls = "no_data"
    elif "conflicting" in s or ("pathogenic" in s and "benign" in s):
        cls = "conflicting"
    elif "likely_pathogenic" in s:
        cls = "likely_pathogenic"
    elif "pathogenic" in s:
        cls = "pathogenic"
    elif "uncertain" in s:
        cls = "VUS"
    elif "benign" in s:
        cls = "benign"
    else:
        cls = "other"
    return {"clinvar": cls, "clinvar_stars": "", "clinvar_source": "vep_rsid_level",
            "clinvar_rescue": False, "clinvar_raw": clin_sig or ""}


# ---------------------------------------------------------------------------
# PVS1 strength (simplified Abou Tayoun 2018)
# ---------------------------------------------------------------------------
PVS1_DOWN2 = {8: 2, 4: 1, 2: 0, 1: 0}


def nmd_escape(row):
    """True when the VEP NMD plugin flags NMD escape, or (fallback) last exon."""
    if "NMD_escaping" in str(row.get("NMD", "") or ""):
        return True
    return is_last_exon(row)


def pvs1_strength(row, csq):
    if "transcript_ablation" in csq:
        return "PVS1", 8
    if csq & {"stop_gained", "frameshift_variant"}:
        if nmd_escape(row):
            return "PVS1_Moderate(NMD_escape,review)", 2
        return "PVS1", 8
    if csq & {"splice_donor_variant", "splice_acceptor_variant"}:
        return "PVS1_Strong(canonical_splice,review)", 4
    if "start_lost" in csq:
        return "PVS1_Moderate(start_lost)", 2
    return "PVS1_unverified", 0


def load_clingen_validity(path):
    """ClinGen gene-disease validity CSV (search.clinicalgenome.org/kb/gene-validity/download).
    Returns {gene: best classification} (Definitive > Strong > Moderate > Limited > other),
    or None if no file."""
    import csv
    if not path:
        return None
    if not os.path.isfile(path):
        print("  [WARN] ClinGen validity file not found: %s (PanelApp proxy used)" % path, file=sys.stderr)
        return None
    order = {"definitive": 5, "strong": 4, "moderate": 3, "limited": 2}
    genes, header, n = {}, None, 0
    with open(path, errors="replace") as fh:
        for row in csv.reader(fh):
            up = [c.strip().upper() for c in row]
            if header is None:
                if "GENE SYMBOL" in up and "CLASSIFICATION" in up:
                    header = up
                continue
            if len(row) < len(header) or row[0].startswith("+"):
                continue
            d = dict(zip(header, [c.strip() for c in row]))
            n += 1
            cls = d.get("CLASSIFICATION", "").lower()
            g = d["GENE SYMBOL"]
            if g not in genes or order.get(cls, 0) > order.get(genes[g], 0):
                genes[g] = cls
    if header is None:
        print("  [WARN] ClinGen validity file has no GENE SYMBOL/CLASSIFICATION header: %s" % path,
              file=sys.stderr)
        return None
    print("  ClinGen gene validity: %d curations, %d genes (%d Definitive/Strong)" % (
        n, len(genes), sum(1 for v in genes.values() if v in ("definitive", "strong"))), file=sys.stderr)
    return genes


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------
def points_to_label(points):
    if points >= 10:
        return "Pathogenic"
    if points >= 6:
        return "Likely Pathogenic"
    if points >= 0:
        return "VUS"
    if points >= -6:
        return "Likely Benign"
    return "Benign"


class Classifier(object):
    def __init__(self, lof_genes=None, ba1_exceptions=None, cfg=None,
                 established_genes=None, vcep_genes=None, gene_thresholds=None):
        """established_genes: genes with a Definitive/Strong gene-disease relationship
        (None = no information, treat every gene as established).
        vcep_genes: genes with a ClinGen VCEP (in-silico capped at Supporting)."""
        self.established_genes = established_genes
        self.vcep_genes = vcep_genes or set()
        self.gene_thresholds = gene_thresholds or {}     # gene -> (ba1, bs1) from VCEP specs
        self.cfg = dict(DEFAULTS)
        if cfg:
            self.cfg.update(cfg)
        self.lof_genes = lof_genes if lof_genes is not None else set()
        self.ba1_exceptions = ba1_exceptions if ba1_exceptions is not None else set()
        if not self.lof_genes:
            print("  [WARN] variant_classify: no LoF gene list loaded; PVS1 will "
                  "not be applied (recorded as PVS1_unverified(0))", file=sys.stderr)

    # -- single variant ------------------------------------------------------
    def classify(self, row, inheritance=""):
        """row: dict with VEP-style keys (SYMBOL, Consequence, EXON, HGVSc,
        REVEL, spliceai_max, gnomAD population columns, MAX_AF).
        Returns dict(acmg_class, acmg_points, acmg_evidence, popmax_af,
        popmax_source, ba1_exception)."""
        c = self.cfg
        gene = row.get("SYMBOL", "") or row.get("gene", "")
        csq = consequences(row)
        inh = inheritance_class(inheritance)
        af, af_src = popmax_af(row)
        revel = to_float(row.get("REVEL"))
        spliceai = to_float(row.get("spliceai_max"))
        ev = []
        pts = 0

        def add(code, p):
            ev.append("%s(%+d)" % (code, p))
            return p

        # BA1: stand-alone benign
        exc_key = (gene, hgvsc_short(row.get("HGVSc", "")))
        ba1_exception = exc_key in self.ba1_exceptions
        g_ba1, g_bs1 = self.gene_thresholds.get(gene, (None, None))
        ba1_thr = g_ba1 if g_ba1 is not None else c["ba1_af"]
        if af is not None and af > ba1_thr:
            if ba1_exception:
                ev.append("BA1_exception(0)")
            else:
                ev.append("BA1(stand-alone%s)" % (",VCEP" if g_ba1 is not None else ""))
                return {"acmg_class": "Benign", "acmg_points": -8,
                        "acmg_evidence": "|".join(ev),
                        "popmax_af": af, "popmax_source": af_src,
                        "ba1_exception": False}

        # BS1: more common than expected for the disorder
        if c["bs1_enabled"] and af is not None and not ba1_exception:
            thr = c["bs1_af_recessive"] if inh == "recessive" else c["bs1_af_dominant"]
            if g_bs1 is not None:
                thr = g_bs1
            if af > thr:
                pts += add("BS1(VCEP)" if g_bs1 is not None else "BS1", -4)

        # PM2: absent / extremely rare
        thr = c["pm2_af_recessive"] if inh == "recessive" else c["pm2_af_dominant"]
        if af is None or af < thr:
            pts += add("PM2", c["pm2_points"])

        # PVS1: null variant in LoF-mechanism gene (simplified Abou Tayoun tree)
        null = bool(csq & NULL_CONSEQUENCES)
        pvs1 = False
        if (null or "start_lost" in csq) and gene in self.lof_genes:
            pvs1 = True
            code, p = pvs1_strength(row, csq)
            if self.established_genes is not None and gene not in self.established_genes:
                code, p = code + "_down2(gene_validity)", PVS1_DOWN2.get(p, 0)
            if p > 0:
                pts += add(code, p)
            else:
                ev.append(code + "(0)")
        elif null:
            ev.append("PVS1_unverified(0)")

        # PM4: protein length change (not when PVS1 applied)
        if not pvs1 and (csq & PM4_CONSEQUENCES):
            # ACMG PM4 is for in-frame indels "in a nonrepeat region". Rule used here
            # (SignalScope choice): not applied when the inserted/deleted unit lies in a tandem
            # repeat of >= pm4_repeat_min_copies copies in the reference (row['tandem_repeat_copies']).
            rep = int(row.get("tandem_repeat_copies") or 0)
            if rep >= c["pm4_repeat_min_copies"]:
                ev.append("PM4_not_applied(repeat_%dx)(0)" % rep)
            else:
                pts += add("PM4", 2)

        # PS1 / PM5: other ClinVar P/LP change at the same codon (missense only)
        codon_ev = row.get("codon_evidence", "")
        if "missense_variant" in csq and codon_ev == "PS1":
            pts += add("PS1", 4)
        elif "missense_variant" in csq and codon_ev == "PM5":
            pts += add("PM5", 2)

        # In-silico, counted once
        is_missense = "missense_variant" in csq
        vcep = gene in self.vcep_genes
        if is_missense and revel is not None:
            done = False
            for thr_, p, code in REVEL_PP3:
                if revel >= thr_:
                    if vcep:
                        code, p = "PP3(VCEP_cap)", 1
                    pts += add(code, p)
                    done = True
                    break
            if not done:
                for thr_, p, code in REVEL_BP4:
                    if revel <= thr_:
                        if code is None:              # very strong band -> Strong
                            code, p = "BP4_Strong", -4
                        if vcep:
                            code, p = "BP4(VCEP_cap)", -1
                        pts += add(code, p)
                        break
            # SpliceAI can still upgrade a missense with no REVEL PP3
            if not done and spliceai is not None and spliceai >= c["spliceai_pp3"]:
                ev.append("SpliceAI_high_on_missense(review)")
        elif not null and not pvs1 and spliceai is not None:
            if spliceai >= c["spliceai_pp3"]:
                pts += add("PP3_SpliceAI", 1)
            elif spliceai <= c["spliceai_bp4"]:
                # Walker 2023: BP7 for synonymous and for intronic variants outside the
                # donor (+1..+6) / acceptor (-20..-1) regions; BP4 otherwise
                # position from the HGVSc intronic offset (donor +1..+6, acceptor -1..-20
                # are the splice regions); VEP terms only as a fallback.
                off = intron_offset(row.get("HGVSc", ""))
                if off is not None and "intron_variant" in csq:
                    in_region = (off[0] == "+" and off[1] <= 6) or (off[0] == "-" and off[1] <= 20)
                    pts += add("BP4_SpliceAI" if in_region else "BP7", -1)
                elif ("synonymous_variant" in csq or "intron_variant" in csq) and not (csq & SPLICE_CONSEQUENCES):
                    # Walker et al. 2023 (PMC10357475): "BP7 should not be applied for synonymous
                    # substitutions ... located at the first base or the last 3 bases of the exon" and
                    # "not ... for variants within a designated splice region (standard or minimal)"; the
                    # VEP splice-region terms are the conservative definition their Figure 4 allows. Such
                    # synonymous variants get BP4 (same -1 points), not BP7.
                    pts += add("BP7", -1)
                elif csq & SPLICE_CONSEQUENCES:
                    pts += add("BP4_SpliceAI", -1)
        elif "synonymous_variant" in csq and spliceai is None:
            ev.append("BP7_not_applied(no_SpliceAI)")

        return {"acmg_class": points_to_label(pts), "acmg_points": pts,
                "acmg_evidence": "|".join(ev) if ev else "none",
                "popmax_af": af, "popmax_source": af_src,
                "ba1_exception": ba1_exception}


# ---------------------------------------------------------------------------
# Tiering
# ---------------------------------------------------------------------------
def assign_tier(v, mode, disease_genes=None):
    """v: dict with gene, chrom, pos, acmg_class, clinvar_rescue, gt,
    inheritance, plus optional n_het_in_gene (for unphased compound-het hint).
    Returns (tier, flags) where tier is one of:
      'excluded_hla'      HLA/MHC: listed separately, never shortlisted
      'excluded_benign'   B/LB by points and no ClinVar rescue
      'tier1'             P/LP (or ClinVar rescue) fitting the inheritance
      'carrier'           P/LP heterozygous in an AR-only gene
      'tier2'             VUS in a disease-relevant gene
      'tier3'             everything else (collapsed in the report)"""
    flags = []
    gene = v.get("gene", "")
    if is_hla_mhc(gene, v.get("chrom"), v.get("pos")):
        return "excluded_hla", ["HLA_MHC_not_interpretable"]

    label = v.get("acmg_class", "VUS")
    rescue = bool(v.get("clinvar_rescue"))
    inh = inheritance_class(v.get("inheritance", ""))
    zyg = zygosity(v.get("gt", ""))

    if rescue and label in ("VUS", "Likely Benign", "Benign"):
        flags.append("ClinVar_P/LP_vs_ACMG_%s" % label.replace(" ", "_"))

    if label in ("Likely Benign", "Benign") and not rescue:
        return "excluded_benign", flags

    plp = label in ("Pathogenic", "Likely Pathogenic") or rescue
    if plp:
        if inh == "recessive" and zyg in ("het",):
            if int(v.get("n_het_in_gene", 1) or 1) >= 2:
                flags.append("possible_compound_het(unphased)")
                return "tier1", flags
            flags.append("carrier(het_in_AR_gene)")
            return "carrier", flags
        return "tier1", flags

    relevant = True
    if disease_genes is not None:
        relevant = gene in disease_genes
    if label == "VUS" and relevant:
        return "tier2", flags
    return "tier3", flags


TIER_RANK = {"tier1": 0, "carrier": 1, "tier2": 2, "tier3": 3,
             "excluded_benign": 8, "excluded_hla": 9}


def priority_key(v):
    """Sort key for the shortlist: tier, ACMG points, ClinVar stars."""
    try:
        stars = int(v.get("clinvar_stars") or 0)
    except ValueError:
        stars = 0
    return (TIER_RANK.get(v.get("tier"), 9),
            -int(v.get("acmg_points", 0) or 0),
            -stars)
