#!/usr/bin/env python3
"""
classify_small_variants.py - SignalScope small-variant decision tree (SNVs + indels).

Talos-style, gene-first prioritisation with long-read phasing, built for WGS scale.
Streams the VEP TSVs (never loads them whole), so memory stays at a few hundred MB.

DECISION TREE (each stage is cheap and runs before the next; counts go to funnel.tsv)
  0. Parse        VEP TSV per caller, one record per variant, best transcript chosen
                  (panel gene first, then MANE Select > canonical > first)
  1. Gene gate    variant hits a gene in the gene table (PanelApp green by default)
                  -> outside: RESEARCH track if strong evidence, else dropped
  2. HLA / MHC    -> listed separately as not interpretable
  3. Consequence  coding / splice impact (HIGH, MODERATE, or SpliceAI >= threshold)
                  or ClinVar germline P/LP at any consequence; for disease RNA genes
                  (RMRP, TERC, RNU4ATAC ...) the RNA exon and promoter also count
  4. Callers + QC depth >= 8; called by both callers (no GQ minimum then), or by one
                  caller with GQ >= 30 (GQ >= 20 when only one caller was run);
                  1-2 bp indels in a homopolymer run >= 5 (ONT weak spot) need both
                  callers and GQ >= 20
  5. Frequency    hard filter > 2% for every gene (VCEP BA1 if lower), on the gnomAD
                  v4.1 GroupMax FAF95 from an allele-exact direct lookup (--gnomad_base;
                  VEP misses gnomAD alleles without a dbSNP rsID); popmax AF of the
                  non-bottleneck populations when FAF is unavailable; ClinVar P/LP recaptured up to 5%; BA1 exceptions kept.
                  Inheritance-aware expectations (0.01% dominant, 1% recessive) only
                  set frequency_flag; a single het above 0.5% cannot fit dominant
  6. Evidence     removed when ClinVar benign/likely benign comes from an expert panel or practice
                  guideline (3-4 stars), or is well reviewed (1-2 stars) AND our own points are
                  benign/likely benign; disagreement at 1-2 stars -> kept with review_flag
                  PRIMARY: ClinVar P/LP | null variant where LoF is a disease mechanism
                           (PanelApp, ClinGen HI 3/30 vs 40, recessive gene, >= 2 reviewed
                           ClinVar null P/LP in the gene; LOEUF only as last resort)
                           | SpliceAI >= primary threshold
                           | same amino-acid change as a ClinVar P/LP (PS1)
                           | other missense at a ClinVar P/LP codon (PM5), >= 2 stars
                             AND REVEL >= 0.773
                  SUPPORTING: REVEL >= 0.644 | SpliceAI >= supporting threshold
                           | ClinVar 0-star P/LP | ClinVar VUS/conflicting
                           | PM5 without the REVEL / star condition
                           | RNA-gene exon or promoter variant
                  no label -> dropped
  7. Inheritance  dominant: one primary variant; recessive: homozygous primary, or
                  two qualifying variants with >= 1 primary (in trans if phased,
                  'phase unknown' flag otherwise; cis pairs excluded);
                  X-linked recessive in males: hemizygous counts as homozygous
                  single het primary in a recessive gene -> CARRIER list
                  supporting evidence only -> moi_fit 'supporting_evidence_only'
                  AD/AR genes: a null variant alone fits the dominant model only with
                  ClinGen HI 3 (otherwise carrier / recessive models), except the cited
                  genes in DOMINANT_NULL_OVERRIDE (COL4A3, COL4A4, Savige 2022)
                  two callers run, call by only one (not ClinVar P/LP) -> Tier 3, fit kept,
                  'single_caller_confirm'
  8. Triage       Tier1_high_priority (= clinical shortlist: primary evidence + inheritance
                  fit) | Tier2_carrier_second_hit_search (primary evidence, single het in
                  a recessive gene) | Tier3_review (supporting evidence only, or primary
                  evidence without inheritance fit). Evidence points and every criterion
                  (variant_classify.py) are shown to rank and explain, NOT as an ACMG
                  class: no Pathogenic/VUS/Benign labels are output. ClinVar separate.

INTENDED USE
  Clinical variant prioritisation to support subsequent diagnostic interpretation by
  qualified curators. Not a diagnostic report. Diagnostic use requires validation by
  the testing laboratory (NPAAC Requirements for medical testing of human genetic
  variation, 4th ed. 2025), including variant classes detected, analytical limits and
  reproducibility. Gene-specific VCEP criteria are not applied unless supplied.

OUTPUTS (in --out_dir, prefix --sample)
  .clinical_shortlist.tsv   tier1: fits gene + evidence + inheritance
  .carriers.tsv             heterozygous P/LP-type variants in recessive genes
  .research.tsv             outside the panel but strong evidence (speculative)
  .not_interpretable.tsv    HLA / MHC variants that passed QC
  .funnel.tsv               counts at every stage, per caller set
  .all_candidates.tsv       everything that reached stage 6, with the reason it stopped

TESTING OPTION
  --clinvar_mask FILE   hide these exact alleles from ClinVar (first column or
                        'variant_key' column, chrom_pos_ref/alt in VEP form). Used by the
                        spike-in test to ask "would we find it if ClinVar had never seen
                        it, while its neighbours are still known?"
"""

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import variant_classify as vc  # noqa: E402

NULL_CSQ = vc.NULL_CONSEQUENCES | {"start_lost"}
KEEP_COLS = [
    "SYMBOL", "Gene", "Feature", "Consequence", "IMPACT", "BIOTYPE", "CANONICAL", "MANE_SELECT",
    "EXON", "INTRON", "HGVSc", "HGVSp", "Existing_variation", "SIFT", "PolyPhen", "REVEL",
    "SpliceAI_pred", "MAX_AF", "MAX_AF_POPS", "gnomADe_AF", "gnomADg_AF", "CLIN_SIG",
    "STRAND", "CDS_position", "Protein_position", "Amino_acids", "Codons", "DISTANCE", "NMD",
] + vc.POPMAX_COLUMNS

# Disease genes whose product is an RNA, not a protein: their pathogenic variants
# sit in the RNA exon or the promoter (e.g. RMRP promoter insertions in
# cartilage-hair hypoplasia), so the coding filter alone would always miss them.
# Genes where a single heterozygous null variant fits the DOMINANT model although ClinGen HI
# is not 3 (HI 30 there describes the recessive form). Each entry needs a citation; the note
# is written to review_flag.
#   Savige J et al. Guidelines for genetic testing and management of Alport syndrome.
#   Clin J Am Soc Nephrol 2022;17(1):143-154 (heterozygotes: autosomal dominant Alport
#   syndrome / TBMN, monitor and treat). PanelApp Australia Haematuria_Alport: monoallelic
#   and biallelic, Green; OMIM 104200 (Alport syndrome 3, AD).
_ALPORT_NOTE = "AD_Alport_spectrum_variable_penetrance(Savige2022)"
DOMINANT_NULL_OVERRIDE = {"COL4A3": _ALPORT_NOTE, "COL4A4": _ALPORT_NOTE}
NCRNA_DISEASE_GENES = {"RMRP", "TERC", "RNU4ATAC", "RNU12", "SNORD118", "RNU7-1", "RPPH1"}
NCRNA_CSQ = {"non_coding_transcript_exon_variant", "mature_miRNA_variant"}

CODON_TABLE = dict(zip(
    [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    """The parser on its own, so other tools (for example the benchmark call-set builder)
    read the SAME QC defaults instead of copying them."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clair3_vep", required=True)
    p.add_argument("--clair3_vcf", required=True, help="normalised (split) Clair3 VCF, bgzipped")
    p.add_argument("--dv_vep", default="")
    p.add_argument("--dv_vcf", default="")
    p.add_argument("--phased_vcf", default="",
                   help="the pipeline's phasing (tree_phase: WhatsHap on the SNVs of both callers, "
                        "indels from the haplotagged reads); when given, its GT / PS replace the callers' own "
                        "phase, and variants it leaves unphased are unphased")
    p.add_argument("--gene_table", required=True, help="panelapp_fetch.py output (or same columns)")
    p.add_argument("--clinvar_vcf", default="", help="ClinVar GRCh38 VCF (bgzip + tabix)")
    p.add_argument("--clingen_dosage", default="", help="ClinGen_gene_curation_list_GRCh38.tsv (for PVS1 genes)")
    p.add_argument("--clingen_validity", default="",
                   help="ClinGen gene-disease validity CSV; genes curated below Strong get PVS1 down 2 levels")
    p.add_argument("--mim2gene", default="",
                   help="OMIM mim2gene.txt (omim.org/static/omim/data/mim2gene.txt): gene MIM number "
                        "for a direct OMIM link per variant (any gene, any panel). Absent -> no OMIM link")
    p.add_argument("--vcep_genes", default="", help="genes with a ClinGen VCEP (one per line): in-silico capped")
    p.add_argument("--gene_constraint", default="",
                   help="gnomAD v4.1 constraint_metrics.tsv: adds LOEUF, pLI, missense Z (display only)")
    p.add_argument("--ba1_exceptions", default=os.path.join(HERE, "..", "assets", "ba1_exceptions.tsv"))
    p.add_argument("--min_confidence", default="green", choices=["green", "amber"])
    p.add_argument("--sex", default="unknown", choices=["male", "female", "unknown"])
    p.add_argument("--out_dir", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--bcftools", default="bcftools")
    # thresholds
    p.add_argument("--af_filter", type=float, default=0.02,
                   help="hard popmax filter for every gene (diagnostic-lab practice); ClinVar P/LP "
                        "recaptured up to --af_clinvar")
    p.add_argument("--af_dominant", type=float, default=0.0001,
                   help="expected max popmax for a dominant disorder: above it -> frequency_flag only")
    p.add_argument("--af_recessive", type=float, default=0.01,
                   help="expected max popmax for a recessive disorder: above it -> frequency_flag only")
    p.add_argument("--af_dominant_fit", type=float, default=0.005,
                   help="a single het variant above this popmax cannot fit a dominant model "
                        "(unless ClinVar P/LP); equals BS1 for dominant genes")
    p.add_argument("--af_clinvar", type=float, default=0.05)
    p.add_argument("--min_gq", type=int, default=20)
    p.add_argument("--min_gq_both", type=int, default=0,
                   help="GQ threshold when BOTH callers make a PASS call (default 0: agreement of "
                        "the two callers replaces the GQ check; DP >= --min_dp still applies). Chosen with "
                        "a QC settings comparison on HG002")
    p.add_argument("--min_dp", type=int, default=8)
    p.add_argument("--single_caller_min_gq", type=int, default=30)
    p.add_argument("--spliceai_primary", type=float, default=0.5)
    p.add_argument("--spliceai_supporting", type=float, default=0.2)
    p.add_argument("--revel_supporting", type=float, default=0.644,
                   help="REVEL >= this is supporting missense evidence (Pejaver et al. 2022 PP3 "
                        "Supporting lower bound)")
    p.add_argument("--pm4_repeat_min_copies", type=int, default=3,
                   help="PM4 (in-frame indel) is not used when the indel unit lies in a tandem "
                        "repeat of >= this many reference copies (ACMG: 'in a nonrepeat region')")
    p.add_argument("--allow_missing_plugins", action="store_true",
                   help="testing only: do not stop when a VEP file lacks the SpliceAI_pred or REVEL column")
    p.add_argument("--pm2_points", type=int, default=1)
    p.add_argument("--pm5_primary_min_revel", type=float, default=0.773,
                   help="PM5 (other missense at a ClinVar P/LP codon) counts as PRIMARY only with "
                        "REVEL >= this and a >= 2-star ClinVar record; otherwise SUPPORTING")
    p.add_argument("--ncrna_promoter_bp", type=int, default=100,
                   help="upstream distance counted as promoter for RNA disease genes")
    p.add_argument("--clinvar_mask", default="", help="testing: alleles hidden from ClinVar")
    p.add_argument("--gnomad_base", default="",
                   help="direct, allele-exact gnomAD lookup for variants reaching the frequency stage: a "
                        "URL or folder prefix; '{ds}' and '{chrom}' are filled in, e.g. https://storage."
                        "googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/{ds}/gnomad.{ds}.v4.1.sites."
                        "chr{chrom}.vcf.bgz (ds = exomes, genomes)")
    p.add_argument("--gnomad_cache", default="",
                   help="TSV cache of gnomAD lookups (shared between samples; only public gnomAD data)")
    p.add_argument("--reference_fasta", default="", help="indexed FASTA (.fai) for homopolymer context")
    p.add_argument("--homopolymer_min", type=int, default=5,
                   help="indel touching a single-base run >= this length = low-complexity indel")
    p.add_argument("--homopolymer_min_gq", type=int, default=20,
                   help="low-complexity indels need both callers (when two were run) and GQ >= this "
                        "(chosen with a QC settings comparison on HG002)")
    p.add_argument("--homopolymer_max_indel", type=int, default=2,
                   help="guard only 1..N bp indels made of the run's own base (the ONT error mode)")
    p.add_argument("--gnomad_threads", type=int, default=8, help="parallel gnomAD queries")
    p.add_argument("--vcep_specs", default=os.path.join(HERE, "..", "assets", "vcep_specs.tsv"),
                   help="VCEP gene specifications (BA1/BS1 FAF thresholds, PP3 cap); only rows "
                        "with verified=yes are applied")
    p.add_argument("--no_dominant_null_override", action="store_true",
                   help="ignore DOMINANT_NULL_OVERRIDE (COL4A3/COL4A4); for sensitivity analyses")
    p.add_argument("--spliceai_extra", default="",
                   help="TSV from spliceai_fill.py (variant, SpliceAI_pred, status): SpliceAI computed "
                        "for variants without a precomputed score (indels); used only where VEP has none")
    p.add_argument("--clinvar_lof_min", type=int, default=2,
                   help="LoF counts as an established mechanism when ClinVar has >= N reviewed "
                        "P/LP null variants in the gene (Abou Tayoun 2018)")
    p.add_argument("--loeuf_max", type=float, default=0.6,
                   help="last-resort rule, only when ClinGen and ClinVar give no answer: a gene with "
                        "LOEUF >= this is treated as LoF-tolerant")
    return p


def parse_args(argv=None):
    return build_parser().parse_args(argv)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def strip_chr(c):
    return c[3:] if c.startswith("chr") else c


def load_gene_table(path, min_conf):
    rank = {"green": 3, "amber": 2, "red": 1}
    genes, alias = {}, {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                header = cols
                continue
            r = dict(zip(header, cols + [""] * (len(header) - len(cols))))
            if rank.get(r.get("confidence", "green"), 3) < rank[min_conf]:
                continue
            genes[r["gene"]] = r
            for a in (r.get("aliases") or "").split(","):
                if a:
                    alias.setdefault(a, r["gene"])
    return genes, alias


def load_clingen_lof(path):
    """Genes with ClinGen haploinsufficiency score 3 (dominant LoF) or 30 (AR).
    The header line is the one containing 'Haploinsufficiency Score' (it may start with '#')."""
    lof, curated = set(), {}
    if not path:
        return lof, curated
    if not os.path.isfile(path):
        print("  [WARN] ClinGen dosage file not found: %s" % path, file=sys.stderr)
        return lof, curated
    header, n_rows = None, 0
    with open(path, errors="replace") as fh:
        for line in fh:
            cols = [c.strip() for c in line.lstrip("#").rstrip("\n").split("\t")]
            if header is None:
                if any(c.lower() == "haploinsufficiency score" for c in cols):
                    header = cols
                continue
            if line.startswith("#"):
                continue
            n_rows += 1
            r = dict(zip(header, cols))
            sym = r.get("Gene Symbol") or r.get("Gene symbol") or cols[0]
            hi = r.get("Haploinsufficiency Score", "")
            if hi not in ("", "-", "Not yet evaluated"):
                curated[sym] = hi
            if hi in ("3", "30"):
                lof.add(sym)
    if header is None:
        print("  [WARN] ClinGen file has no 'Haploinsufficiency Score' column: %s" % path, file=sys.stderr)
    else:
        print("  ClinGen dosage: %d genes read, %d with HI score 3 or 30" % (n_rows, len(lof)), file=sys.stderr)
    return lof, curated


LOF_MC = ("SO:0001587", "SO:0001589", "SO:0001574", "SO:0001575",
          "nonsense", "frameshift", "splice_donor", "splice_acceptor")


def clinvar_gene_lof_evidence(clinvar_vcf, genes, alias, mask=frozenset()):
    """{gene: [n_PLP, n_PLP_null]} over rescue-level germline ClinVar records (>= 2 stars, or
    1 star without conflict). Null = nonsense / frameshift / canonical splice (INFO/MC).
    Masked alleles (testing) are not counted."""
    out = defaultdict(lambda: [0, 0])
    if not clinvar_vcf or not os.path.isfile(clinvar_vcf):
        return None
    with gzip.open(clinvar_vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#" or "athogenic" not in line:
                continue
            c = line.split("\t", 8)
            inf = c[7]
            d = {}
            for f in inf.split(";"):
                if f.startswith(("CLNSIG=", "CLNREVSTAT=", "GENEINFO=", "MC=")):
                    k, v = f.split("=", 1)
                    d[k] = v
            rec = vc.clinvar_record(d.get("CLNSIG", ""), d.get("CLNREVSTAT", ""))
            if not rec["clinvar_rescue"]:
                continue
            if mask:
                p_, r_, a_ = vep_norm(c[1], c[3], c[4].split(",")[0])
                if (strip_chr(c[0]), p_, r_, a_) in mask:
                    continue
            is_null = any(t in d.get("MC", "") for t in LOF_MC)
            for gi in d.get("GENEINFO", "").split("|"):
                g = gi.split(":")[0]
                g = g if g in genes else alias.get(g, "")
                if g:
                    out[g][0] += 1
                    out[g][1] += int(is_null)
    return dict(out)


def load_vcep_specs(path):
    """assets/vcep_specs.tsv -> {gene: {'ba1': float|None, 'bs1': float|None, 'pp3_cap': bool,
    'source': url}}; rows with verified != yes are reported and skipped."""
    out, skipped = {}, 0
    if not path or not os.path.isfile(path):
        return out
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if header is None:
                header = c
                continue
            r = dict(zip(header, c + [""] * (len(header) - len(c))))
            if r.get("verified", "").lower() != "yes":
                skipped += 1
                continue
            out[r["gene"]] = {"ba1": vc.to_float(r.get("ba1_faf")), "bs1": vc.to_float(r.get("bs1_faf")),
                              "pp3_cap": r.get("pp3_cap", "").lower() == "yes",
                              "source": r.get("source_url", "")}
    print("  VCEP specifications: %d genes applied, %d rows not yet verified (skipped)" % (len(out), skipped),
          file=sys.stderr)
    return out


def load_constraint(path):
    """gnomAD constraint TSV -> {gene: (loeuf, pli, mis_z, constraint_flags)}, MANE Select row
    first, then canonical."""
    out, rank = {}, {}
    if not path:
        return out
    if not os.path.isfile(path):
        print("  [WARN] gene constraint file not found: %s" % path, file=sys.stderr)
        return out
    with open_any(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        need = {"gene", "lof.oe_ci.upper", "lof.pLI", "mis.z_score"}
        if not need <= set(ix):
            print("  [WARN] constraint file lacks columns %s" % sorted(need - set(ix)), file=sys.stderr)
            return out
        for line in fh:
            c = line.rstrip("\n").split("\t")
            g = c[ix["gene"]]
            r_ = 2 if "mane_select" in ix and c[ix["mane_select"]].lower() == "true" else \
                1 if "canonical" in ix and c[ix["canonical"]].lower() == "true" else 0
            if g not in out or r_ > rank[g]:
                fl = c[ix["constraint_flags"]] if "constraint_flags" in ix else ""
                out[g] = (c[ix["lof.oe_ci.upper"]], c[ix["lof.pLI"]], c[ix["mis.z_score"]], fl)
                rank[g] = r_
    print("  gene constraint: %d genes" % len(out), file=sys.stderr)
    return out


class Fasta(object):
    """Minimal indexed FASTA reader (.fai), standard library only."""
    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as f:
            for line in f:
                n, ln, off, lb, lw = line.split("\t")[:5]
                self.idx[n] = (int(ln), int(off), int(lb), int(lw))

    def fetch(self, chrom, start, end):
        """1-based inclusive."""
        key = chrom if chrom in self.idx else ("chr" + chrom if "chr" + chrom in self.idx else chrom[3:])
        if key not in self.idx:
            return ""
        ln, off, lb, lw = self.idx[key]
        start, end = max(1, start), min(ln, end)
        if end < start:
            return ""
        s0, e0 = start - 1, end
        bs = off + (s0 // lb) * lw + s0 % lb
        be = off + ((e0 - 1) // lb) * lw + (e0 - 1) % lb + 1
        self.fh.seek(bs)
        return self.fh.read(be - bs).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def homopolymer_run(fasta, chrom, pos, ref, alt, max_indel=2):
    """Length of the reference single-base run in which a short indel adds or removes copies of
    that same base (the ONT error mode), e.g. AAAAAA -> AAAAA. VEP-form key: deletion removes
    ref at pos..; insertion sits between pos-1 and pos. 0 when not such an indel (SNV, complex
    or longer indel, or bases different from the run), or without FASTA."""
    if fasta is None:
        return 0
    if ref == "-" and alt != "-":
        s = alt
    elif alt == "-" and ref != "-":
        s = ref
    else:
        return 0
    if len(s) > max_indel or len(set(s)) != 1:
        return 0
    b = s[0]
    win = 30
    seq = fasta.fetch(chrom, pos - win, pos + len(s) + win)
    if not seq:
        return 0
    off = pos - win
    best = 0
    for anchor in (pos - 1, pos):
        i = anchor - off
        if i < 0 or i >= len(seq) or seq[i] != b:
            continue
        l, r = i, i
        while l > 0 and seq[l - 1] == b:
            l -= 1
        while r < len(seq) - 1 and seq[r + 1] == b:
            r += 1
        best = max(best, r - l + 1)
    return best


def tandem_repeat_copies(fasta, chrom, pos, ref, alt, win=200):
    """Number of tandem copies, in the REFERENCE, of the repeat
    unit of a pure insertion/deletion (VEP-form key; deletion = ref at pos.., insertion between
    pos-1 and pos). The unit is the shortest period of the inserted/deleted sequence. For a
    deletion the deleted copy is counted (it is in the reference); for an insertion only the
    copies already present. 0 without FASTA or for SNVs / complex changes."""
    if fasta is None:
        return 0
    if ref == "-" and alt != "-":
        s = alt
    elif alt == "-" and ref != "-":
        s = ref
    else:
        return 0
    n = len(s)
    u = s
    for k in range(1, n + 1):
        if n % k == 0 and s[:k] * (n // k) == s:
            u = s[:k]
            break
    k = len(u)
    start = max(1, pos - win)
    seq = fasta.fetch(chrom, start, pos + n + win)
    if not seq:
        return 0
    i0 = pos - start
    c, j = 0, i0
    while j >= 0 and seq[j:j + k] == u:
        c += 1
        j += k
    j = i0 - k
    while j >= 0 and seq[j:j + k] == u:
        c += 1
        j -= k
    return c


GNOMAD_POPS = ["afr", "amr", "eas", "nfe", "sas"]      # non-bottleneck, as in variant_classify
GNOMAD_FIELDS = ["AC", "AN", "AF", "nhomalt"] + ["AF_" + p for p in GNOMAD_POPS] + ["fafmax_faf95_max"]


def gnomad_lookup(keys, base, bcftools, cache_path, threads=8):
    """Allele-exact gnomAD lookup (exomes + genomes) with bcftools; works on local files or
    public URLs (only the needed index blocks are fetched). Queries run in parallel (one job per
    dataset x chromosome) and each chromosome is appended to the cache as soon as both datasets
    are done, so an interrupted run keeps its finished work.
    Returns ({key: {'popmax', 'faf', 'af', 'nhomalt', 'ac', 'datasets', 'filters'}}, absent_keys).
    Only public reference data is cached. Cache entries from an older field layout are re-queried."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    import threading
    found, absent = {}, set()
    cache = {}
    nf = len(GNOMAD_FIELDS)
    if cache_path and os.path.isfile(cache_path):
        with open(cache_path) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) != 3 or line.startswith("#"):
                    continue
                if c[1] == "found" and any(len(h.split(",")) != nf + 2 for h in c[2].split("|")):
                    continue                              # cache entry without FAF: query again
                cache[c[0]] = (c[1], c[2])

    def kstr(k):
        return "%s_%d_%s/%s" % k
    todo = [k for k in keys if kstr(k) not in cache]
    new = {}
    if todo and base:
        by_chrom = defaultdict(list)
        for k in todo:
            if k[0] in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                by_chrom[k[0]].append(k)
        tbi_dir = os.path.join(os.path.dirname(os.path.abspath(cache_path)) if cache_path else
                               tempfile.gettempdir(), "gnomad_index_cache")
        os.makedirs(tbi_dir, exist_ok=True)
        fmt = "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t" + "\t".join("%INFO/" + f for f in GNOMAD_FIELDS) + "\n"

        def job(ds, chrom, ks):
            url = base.replace("{ds}", ds).replace("{chrom}", chrom)
            tf = tempfile.NamedTemporaryFile("w", suffix=".reg", delete=False, dir=tbi_dir)
            for p_ in sorted(set(x[1] for x in ks)):
                tf.write("chr%s\t%d\t%d\n" % (chrom, max(1, p_ - 1), p_ + 1))
            tf.close()
            try:
                out = subprocess.run([bcftools, "query", "-R", tf.name, "-f", fmt, url],
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     universal_newlines=True, check=True, cwd=tbi_dir).stdout
                return ds, chrom, out, ""
            except (subprocess.CalledProcessError, OSError) as e:
                return ds, chrom, None, (getattr(e, "stderr", "") or str(e)).strip()[:200]
            finally:
                os.unlink(tf.name)

        hits = defaultdict(list)
        done = defaultdict(set)
        failed_chroms = set()
        lock = threading.Lock()
        n_fail = 0
        jobs = [(ds, ch, ks) for ch, ks in sorted(by_chrom.items()) for ds in ("exomes", "genomes")]
        print("  gnomAD: %d queries (%d chromosomes x 2 datasets), %d threads" % (
              len(jobs), len(by_chrom), threads), file=sys.stderr, flush=True)
        with ThreadPoolExecutor(max_workers=max(1, threads)) as ex:
            futs = [ex.submit(job, *j) for j in jobs]
            for n_done, f in enumerate(as_completed(futs), 1):
                ds, chrom, out, err = f.result()
                want = set(by_chrom[chrom])
                if out is None:
                    n_fail += 1
                    failed_chroms.add(chrom)
                    print("  [WARN] gnomAD %s chr%s lookup failed: %s" % (ds, chrom, err), file=sys.stderr, flush=True)
                else:
                    for line in out.splitlines():
                        c = line.split("\t")
                        if len(c) < 5 + nf:
                            continue
                        p_, r_, a_ = vep_norm(c[1], c[2], c[3])
                        k = (strip_chr(c[0]), p_, r_, a_)
                        if k in want:
                            hits[k].append([ds, c[4]] + c[5:5 + nf])
                done[chrom].add(ds)
                if len(done[chrom]) == 2:
                    if chrom not in failed_chroms:          # unknown, not absent: never cached
                        chunk = {}
                        for k in by_chrom[chrom]:
                            chunk[kstr(k)] = ("absent", "") if k not in hits else \
                                ("found", "|".join(",".join(h) for h in hits[k]))
                        new.update(chunk)
                        if cache_path:
                            with lock, open(cache_path, "a") as fh:
                                for ks_, (st, val) in chunk.items():
                                    fh.write("%s\t%s\t%s\n" % (ks_, st, val))
                    print("  gnomAD: chr%s done (%d/%d queries finished)" % (chrom, n_done, len(jobs)),
                          file=sys.stderr, flush=True)
        if n_fail:
            print("  [WARN] %d gnomAD queries failed; those variants keep the VEP frequency" % n_fail,
                  file=sys.stderr)
    cache.update(new)
    for k in keys:
        v = cache.get(kstr(k))
        if v is None:
            continue
        if v[0] == "absent":
            absent.add(k)
            continue
        pm, af, faf, nh, ac, dss, fls = None, None, None, 0, 0, [], []
        for h in v[1].split("|"):
            f = h.split(",")
            ds, filt, vals = f[0], f[1], dict(zip(GNOMAD_FIELDS, f[2:]))
            dss.append(ds)
            fls.append("%s:%s" % (ds, filt))
            for fld in ["AF_" + p for p in GNOMAD_POPS]:
                x = vc.to_float(vals.get(fld))
                if x is not None:
                    pm = x if pm is None else max(pm, x)
            x = vc.to_float(vals.get("AF"))
            if x is not None:
                af = x if af is None else max(af, x)
            x = vc.to_float(vals.get("fafmax_faf95_max"))
            if x is not None:
                faf = x if faf is None else max(faf, x)
            nh += int(vc.to_float(vals.get("nhomalt")) or 0)
            ac += int(vc.to_float(vals.get("AC")) or 0)
        found[k] = {"popmax": pm if pm is not None else af, "af": af, "faf": faf, "nhomalt": nh, "ac": ac,
                    "datasets": "+".join(dss), "filters": ";".join(fls)}
    return found, absent


def load_spliceai_extra(path):
    """spliceai_fill.py output -> {key: SpliceAI_pred} (only rows with a score)."""
    out = {}
    if not path:
        return out
    if not os.path.isfile(path):
        sys.exit("[FATAL] --spliceai_extra %s not found" % path)
    with open(path) as fh:
        for l in fh:
            if l.startswith(("#", "variant\t")):
                continue
            f = l.rstrip("\n").split("\t")
            m = UV_RE.match(f[0])
            if m and len(f) > 1 and f[1]:
                out[(strip_chr(m.group(1)), int(m.group(2)), m.group(3), m.group(4))] = f[1]
    print("  SpliceAI computed here (spliceai_fill.py): %d variants with a score" % len(out), file=sys.stderr)
    return out


def pass_call(g):
    """True for a caller record that is a real call: FILTER PASS (or '.') and a non-reference
    allele in GT."""
    if not g or g.get("filter") not in ("PASS", "."):
        return False
    return any(a not in ("0", ".", "") for a in re.split(r"[/|]", g.get("gt", "")))


def load_mim2gene(path):
    """{gene symbol: MIM number} and {Ensembl gene ID: MIM number} for OMIM entries of type
    'gene' or 'gene/phenotype' (the gene's own entry, e.g. BRCA2 -> 600185). Columns of mim2gene.txt:
    MIM Number, MIM Entry Type, Entrez Gene ID, Approved Gene Symbol (HGNC), Ensembl Gene ID."""
    by_sym, by_ens = {}, {}
    if not path:
        return by_sym, by_ens
    if not os.path.isfile(path):
        print("  [WARN] mim2gene file not found: %s (no OMIM links)" % path, file=sys.stderr)
        return by_sym, by_ens
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or c[1].strip() not in ("gene", "gene/phenotype"):
                continue
            if c[3].strip():
                by_sym.setdefault(c[3].strip(), c[0].strip())
            if len(c) > 4 and c[4].strip():
                by_ens.setdefault(c[4].strip(), c[0].strip())
    print("  OMIM mim2gene: %d gene entries" % len(by_sym), file=sys.stderr)
    return by_sym, by_ens


def load_clingen_diseases(path):
    """{gene: 'disease (MOI, classification); ...'} from the ClinGen gene-disease validity
    CSV, for the report's disease line (any gene with a ClinGen curation; display only)."""
    import csv
    out = {}
    if not path or not os.path.isfile(path):
        return out
    order = {"definitive": 0, "strong": 1, "moderate": 2, "limited": 3}
    header, rows = None, defaultdict(list)
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
            rows[d["GENE SYMBOL"]].append((order.get(d.get("CLASSIFICATION", "").lower(), 9),
                                           d.get("DISEASE LABEL", ""), d.get("MOI", ""), d.get("CLASSIFICATION", "")))
    for g, lst in rows.items():
        seen, parts = set(), []
        for _, dis, moi, cls in sorted(lst):
            if (dis, moi) in seen or not dis:
                continue
            seen.add((dis, moi))
            parts.append("%s (%s, %s)" % (dis, moi or "MOI not given", cls))
        out[g] = "; ".join(parts)
    return out


def variant_type(ref, alt):
    """SNV / MNV / insertion / deletion / delins, from the VEP-form alleles ('-' = empty)."""
    if ref == "-":
        return "insertion"
    if alt == "-":
        return "deletion"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "delins"


def vep_norm(pos, ref, alt):
    """Represent a VCF allele the way VEP prints it in Uploaded_variation:
    drop one shared leading base for indels, '-' for empty."""
    pos = int(pos)
    if len(ref) != len(alt) and ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:] or "-", alt[1:] or "-", pos + 1
    return pos, ref, alt


def read_vcf_genotypes(path):
    """{(chrom,pos,ref,alt) in VEP form: {gt, gq, dp, af, ps, filter, qual}} for a small VCF."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            fmt = c[8].split(":")
            smp = c[9].split(":")
            d = dict(zip(fmt, smp))
            for alt in c[4].split(","):
                pos, ref, a = vep_norm(c[1], c[3], alt)
                af = d.get("AF") or d.get("VAF") or ""
                out[(strip_chr(c[0]), pos, ref, a)] = {
                    "gt": d.get("GT", ""), "gq": d.get("GQ", ""), "dp": d.get("DP", ""),
                    "af": af, "ps": d.get("PS", ""), "filter": c[6], "qual": c[5],
                    # the allele as written in the VCF (anchor base kept), for gnomAD / VarSome links
                    "vcf": "%s-%s-%s-%s" % (strip_chr(c[0]), c[1], c[3], alt)}
    return out


UV_RE = re.compile(r"^(.+)_(\d+)_([^/]+)/(.+)$")


def stream_vep(path):
    """Yield (key, [rows]) per variant allele; key=(chrom,pos,ref,alt).
    VEP writes all rows of one input record (Uploaded_variation) consecutively, but at
    multi-allelic sites it alternates between alleles transcript by transcript, so rows
    are grouped by Uploaded_variation first and split by Allele afterwards.
    Only one input record is held in memory at a time."""
    idx, cur_uv, cur, header = None, None, [], None
    counts = {"records": 0, "unparsed": 0, "unmapped_allele": 0}

    def flush(uv, rows):
        counts["records"] += 1
        m = UV_RE.match(uv)
        if not m:
            counts["unparsed"] += 1
            if counts["unparsed"] <= 3:
                print("  [WARN] cannot read variant name %r in %s (VCF ID column must be '.')" % (uv, path),
                      file=sys.stderr)
            return []
        chrom, pos, ref, alts = m.group(1), int(m.group(2)), m.group(3), m.group(4)
        alt_list = alts.split("/")
        # Two naming styles exist: VEP-trimmed ('C/-', '-/G') and VCF-style with the
        # anchor base ('CCTCA/C', 'G/GGGGACGGC'; VEP >= 110). Both are mapped to the
        # same canonical key = vep_norm() of the VCF allele, which is also what the
        # genotype lookup and ClinVar lookup use. The Allele column (always trimmed
        # when all alleles share the first base) is mapped back to its alt.
        vep_style = "-" in ref or any("-" in x for x in alt_list)
        share_first = all(x[:1] == ref[:1] for x in alt_list) and ref != "-"
        # VEP prints the Allele trimmed ONLY when every allele of the record shares the
        # first base; otherwise it prints the full allele (verified: chr1:97257464 A>AG,G gives
        # 'AG' and 'G'). The printed forms are mapped first; other spellings only as a fallback
        # that never overwrites them, so an SNV is never merged with a trimmed insertion.
        amap, canons = {}, []
        for alt in alt_list:
            canon = (pos, ref, alt) if vep_style else vep_norm(pos, ref, alt)
            printed = alt if vep_style else ((alt[1:] or "-") if share_first else alt)
            canons.append((alt, printed, canon))
            amap.setdefault(printed, canon)
        for alt, printed, canon in canons:
            for form in (alt, canon[2]):
                amap.setdefault(form, canon)
        by_key = {}
        for allele, row in rows:
            canon = amap.get(allele)
            if canon is None:
                counts["unmapped_allele"] += 1      # never guessed: reported, not assigned
                continue
            by_key.setdefault(canon, []).append(row)
        return [((strip_chr(chrom),) + k, rs) for k, rs in by_key.items()]

    with open_any(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].rstrip("\n").split("\t")
                idx = {h: i for i, h in enumerate(header)}
                continue
            cols = line.rstrip("\n").split("\t")
            if header is None:
                sys.exit("[FATAL] %s has data rows before the #Uploaded_variation header" % path)
            if len(cols) != len(header):
                # VEP always writes every column, so a short or long row means a cut-off file
                # (disk full, interrupted copy): stop the run
                sys.exit("[FATAL] VEP row with %d columns, header has %d (truncated or corrupt file?) in %s: %s"
                         % (len(cols), len(header), path, line[:80]))
            uv = cols[0]
            allele = cols[idx["Allele"]] if "Allele" in idx and idx["Allele"] < len(cols) else ""
            row = {k: (cols[idx[k]] if k in idx and idx[k] < len(cols) else "") for k in KEEP_COLS}
            if uv != cur_uv:
                if cur_uv is not None:
                    for item in flush(cur_uv, cur):
                        yield item
                cur_uv, cur = uv, []
            cur.append((allele, row))
    if cur_uv is not None:
        for item in flush(cur_uv, cur):
            yield item
    n, bad = counts["records"], counts["unparsed"]
    if counts["unmapped_allele"]:
        print("  [WARN] %d VEP rows had an Allele that matches no ALT of their record in %s (not used)" % (
            counts["unmapped_allele"], path), file=sys.stderr)
    if bad:
        print("  [WARN] %d of %d VEP records had unreadable variant names in %s" % (bad, n, path), file=sys.stderr)
    if n and bad > 0.05 * n:
        sys.exit("[FATAL] %d/%d VEP records could not be parsed (Uploaded_variation must look like "
                 "chrom_pos_ref/alt; set the VCF ID column to '.' before running VEP)" % (bad, n))


def spliceai_max(rows, genes_ok=None):
    """Highest SpliceAI delta score of the variant. SpliceAI_pred entries look like
    SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|...; a variant near two genes gets one entry per gene. With genes_ok, only entries whose
    SYMBOL is in genes_ok count, so a variant never takes a neighbouring gene's score."""
    best = None
    for r in rows:
        for entry in (r.get("SpliceAI_pred") or "").split(","):
            parts = entry.split("|")
            if genes_ok is not None and parts and parts[0] not in genes_ok:
                continue
            if len(parts) >= 5:
                try:
                    m = max(float(x) for x in parts[1:5])
                    best = m if best is None else max(best, m)
                except ValueError:
                    pass
    return best


def revel_fill(r, rows):
    """REVEL scores are keyed by Ensembl transcript ID, so a missense change reported on
    a newer MANE transcript that REVEL does not list gets no score although the same change is
    scored on another transcript. When the chosen row has no REVEL, take the highest REVEL of
    the other missense rows of the same variant with the same amino-acid change (any missense
    row when Amino_acids is not given). Returns (row, source)."""
    if vc.to_float(r.get("REVEL")) is not None:
        return r, ("chosen_transcript" if "missense_variant" in (r.get("Consequence") or "") else "")
    if "missense_variant" not in (r.get("Consequence") or ""):
        return r, ""
    aa = r.get("Amino_acids") or ""
    vals = [vc.to_float(x.get("REVEL")) for x in rows
            if "missense_variant" in (x.get("Consequence") or "")
            and (not aa or aa in ("-", "") or (x.get("Amino_acids") or "") == aa)]
    vals = [v for v in vals if v is not None]
    if not vals:
        return r, "none"
    r = dict(r)
    r["REVEL"] = "%g" % max(vals)
    return r, "other_transcript_same_aa_change"


def exon_intron_label(r):
    """'exon 3/10', 'intron 2/9', or '' when the variant is in neither (upstream, downstream,
    intergenic). The full region, with the distance and the transcript, comes from
    region_annotate.py."""
    ex, it = (r.get("EXON") or "").strip(), (r.get("INTRON") or "").strip()
    if ex not in ("", "-"):
        return "exon " + ex
    if it not in ("", "-"):
        return "intron " + it
    return ""


def pick_row(rows, genes, alias):
    def gene_of(r):
        s = r.get("SYMBOL", "")
        return s if s in genes else alias.get(s, "")
    # Among panel-gene rows, use MANE Select / canonical transcripts when there are any, and
    # pick the most severe consequence among them (a variant can be upstream of one panel gene
    # and a frameshift in its neighbour, e.g. NTHL1 / TSC2 share a bidirectional promoter).
    impact_rank = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}

    def tx_rank(r):
        if r.get("MANE_SELECT", "-") not in ("-", ""):
            return 0
        return 1 if r.get("CANONICAL") == "YES" else 2
    in_panel = [r for r in rows if gene_of(r)]
    pool = in_panel or rows
    main = [r for r in pool if tx_rank(r) < 2] or pool
    best = min(main, key=lambda r: (impact_rank.get(r.get("IMPACT", ""), 4), tx_rank(r)))
    return best, gene_of(best)


# ---------------------------------------------------------------------------
# ClinVar: exact-allele lookup for a set of variants (bcftools, streaming)
# ---------------------------------------------------------------------------
def load_mask(path):
    """Variant keys (chrom,pos,ref,alt in VEP form) to hide from ClinVar (testing)."""
    mask = set()
    if not path:
        return mask
    with open(path) as fh:
        for line in fh:
            m = UV_RE.match(line.split("\t")[0].strip())
            if m:
                mask.add((strip_chr(m.group(1)), int(m.group(2)), m.group(3), m.group(4)))
    print("  ClinVar mask: %d alleles hidden from ClinVar (testing mode)" % len(mask), file=sys.stderr)
    return mask


def _clinvar_query(positions, clinvar_vcf, bcftools, pad):
    """bcftools query of ClinVar at (chrom, pos) +- pad. Returns the output lines."""
    tf = tempfile.NamedTemporaryFile("w", suffix=".regions.tsv", delete=False)
    tf.close()
    reg = sorted(set(positions), key=lambda x: (x[0], x[1]))
    with open(tf.name, "w") as fh:
        for c, p in reg:
            fh.write("%s\t%d\t%d\n" % (c, max(1, p - pad), p + pad))
    cmd = [bcftools, "query", "-R", tf.name,
           "-f", "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%ID\n", clinvar_vcf]
    try:
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, check=True).stdout
    except (subprocess.CalledProcessError, OSError) as e:
        err = getattr(e, "stderr", "") or ""
        os.unlink(tf.name)
        sys.exit("[FATAL] ClinVar lookup failed (%s). bcftools said: %s\n"
                 "Check that %s is a bgzipped, tabix-indexed ClinVar VCF with CLNSIG/CLNREVSTAT, "
                 "or run without --clinvar_vcf on purpose." % (e, err.strip()[:400], clinvar_vcf))
    if os.path.exists(tf.name):
        os.unlink(tf.name)
    return out.splitlines()


def clinvar_lookup(keys, clinvar_vcf, bcftools, mask=frozenset()):
    res = {}
    if not clinvar_vcf or not keys:
        return res
    # window of +-1 bp around each VEP position (VCF indel anchor sits at pos-1)
    lines = _clinvar_query([(k[0], k[1]) for k in keys], clinvar_vcf, bcftools, 1)
    want = set(keys) - set(mask)
    for line in lines:
        c = line.split("\t")
        if len(c) < 7:
            continue
        for alt in c[3].split(","):
            pos, ref, a = vep_norm(c[1], c[2], alt)
            k = (strip_chr(c[0]), pos, ref, a)
            if k in want:
                rec = vc.clinvar_record(c[4], c[5])
                rec["clinvar_id"] = c[6]
                res[k] = rec
    return res


# ---------------------------------------------------------------------------
# ClinVar at the same codon: PS1 (same amino-acid change, different nucleotide)
# and PM5 (different missense at a residue with a ClinVar P/LP missense)
# ---------------------------------------------------------------------------
def codon_positions(key, row):
    """For a missense SNV: (genomic positions of codon bases 0,1,2, strand, ref codon,
    ref aa, alt aa) from VEP CDS_position / Codons / Amino_acids / STRAND; None if
    not computable. Codons split by an intron are caught later by the base check."""
    if len(key[2]) != 1 or len(key[3]) != 1 or key[2] == "-" or key[3] == "-":
        return None
    try:
        cds = int(str(row.get("CDS_position", "")).split("-")[0].split("/")[0])
        strand = int(row.get("STRAND", ""))
    except ValueError:
        return None
    codons = (row.get("Codons") or "").split("/")
    aas = (row.get("Amino_acids") or "").split("/")
    if len(codons) != 2 or len(aas) != 2 or len(codons[0]) != 3:
        return None
    phase = (cds - 1) % 3
    pos = key[1]
    gpos = [pos - phase + i for i in range(3)] if strand == 1 else [pos + phase - i for i in range(3)]
    return gpos, strand, codons[0].upper(), aas[0], aas[1]


def codon_lookup(cands, clinvar_vcf, bcftools, mask=frozenset()):
    """cands: {key: row} missense SNVs. Returns {key: {'code': 'PS1'|'PM5', 'stars': n,
    'ids': 'id,id', 'aa': 'R/H'}} using rescue-level ClinVar P/LP records only."""
    res = {}
    if not clinvar_vcf or not cands:
        return res
    info, positions = {}, []
    for k, row in cands.items():
        cp = codon_positions(k, row)
        if cp:
            info[k] = cp
            positions += [(k[0], p) for p in cp[0]]
    if not info:
        return res
    by_pos = defaultdict(list)
    for line in _clinvar_query(positions, clinvar_vcf, bcftools, 0):
        c = line.split("\t")
        if len(c) < 7 or len(c[2]) != 1:
            continue
        for alt in c[3].split(","):
            if len(alt) == 1:
                by_pos[(strip_chr(c[0]), int(c[1]))].append((c[2], alt, c[4], c[5], c[6]))
    for k, (gpos, strand, rcodon, raa, aaa) in info.items():
        best = None
        for i, gp in enumerate(gpos):
            for ref, alt, clnsig, revstat, cid in by_pos.get((k[0], gp), []):
                ck = (k[0], gp, ref, alt)
                if ck == k or ck in mask:
                    continue
                rb, ab = (ref, alt) if strand == 1 else (COMP.get(ref, "N"), COMP.get(alt, "N"))
                if rcodon[i] != rb:          # codon split by an intron, or annotation mismatch
                    continue
                acodon = rcodon[:i] + ab + rcodon[i + 1:]
                aa = CODON_TABLE.get(acodon)
                if not aa or aa == "*" or aa == raa:
                    continue
                rec = vc.clinvar_record(clnsig, revstat)
                if not rec["clinvar_rescue"]:
                    continue
                code = "PS1" if aa == aaa else "PM5"
                cand = (code == "PS1", rec["clinvar_stars"], code, cid, "%s/%s" % (raa, aa))
                if best is None or cand[:2] > best[:2]:
                    best = cand
        if best:
            res[k] = {"code": best[2], "stars": best[1], "ids": best[3], "aa": best[4]}
    return res


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def check_vep_columns(paths, allow_missing):
    """Stop when a VEP file lacks the SpliceAI_pred or REVEL column: without them no variant
    could be scored. The pipeline's vep_check makes the same test; this covers direct runs."""
    for p in paths:
        if not p:
            continue
        head = ""
        with open_any(p) as fh:
            for line in fh:
                if line.startswith("#Uploaded_variation") or (line.startswith("#") and not line.startswith("##")):
                    head = line
                    break
        cols = head.lstrip("#").rstrip("\n").split("\t")
        miss = [c for c in ("SpliceAI_pred", "REVEL") if c not in cols]
        if miss:
            msg = "%s has no %s column (VEP plugin not loaded?)" % (p, " / ".join(miss))
            if allow_missing:
                print("  [WARN] " + msg + "; continuing (--allow_missing_plugins)", file=sys.stderr)
            else:
                sys.exit("[FATAL] " + msg)


PHASED = {}   # {VEP-form key: genotype} from --phased_vcf


def main():
    a = parse_args()
    if a.phased_vcf:
        PHASED.update(read_vcf_genotypes(a.phased_vcf))
        print("  phased VCF: %d alleles, %d phased" % (
            len(PHASED), sum(1 for v in PHASED.values() if "|" in v["gt"])), file=sys.stderr)
    check_vep_columns([a.clair3_vep, getattr(a, "dv_vep", "")], a.allow_missing_plugins)
    os.makedirs(a.out_dir, exist_ok=True)
    pre = os.path.join(a.out_dir, a.sample)
    genes, alias = load_gene_table(a.gene_table, a.min_confidence)
    print("[classify] sample=%s genes=%d (min %s)" % (a.sample, len(genes), a.min_confidence), file=sys.stderr)
    clingen_lof, clingen_hi = load_clingen_lof(a.clingen_dosage)
    mask = load_mask(a.clinvar_mask)
    cv_lof = clinvar_gene_lof_evidence(a.clinvar_vcf, genes, alias, mask)
    constraint = load_constraint(a.gene_constraint)

    def lof_mechanism(gene):
        """(is_mechanism, label) following Abou Tayoun 2018: is LoF an established disease
        mechanism for this gene?"""
        g_ = genes[gene]
        hi = clingen_hi.get(gene, "")
        n_plp, n_null = (cv_lof or {}).get(gene, (0, 0))
        loeuf = vc.to_float(constraint.get(gene, ("",))[0])
        if g_.get("lof_not_mechanism") == "yes":
            return False, "LoF_not_mechanism(PanelApp)"
        if hi == "40":
            return False, "ClinGen_HI_40(dosage_sensitivity_unlikely)"
        if hi in ("3", "30"):
            return True, "ClinGen_HI_%s" % hi
        if vc.inheritance_class(g_.get("moi")) == "recessive":
            return True, "recessive_gene"
        if cv_lof is not None and n_null >= a.clinvar_lof_min:
            return True, "ClinVar_null_PLP=%d" % n_null
        if cv_lof is not None and n_plp >= 5 and n_null == 0:
            return False, "no_null_among_%d_ClinVar_PLP" % n_plp
        if loeuf is not None and loeuf < a.loeuf_max:
            return True, "LOEUF=%.2f" % loeuf
        if cv_lof is not None and loeuf is not None:
            return False, "LoF_mechanism_unverified(LOEUF=%.2f,ClinVar_null=%d)" % (loeuf, n_null)
        return True, "LoF_assumed(no_ClinGen_ClinVar_or_constraint_data)"

    lof_status = {g: lof_mechanism(g) for g in genes}
    # PVS1 only where the mechanism is ESTABLISHED (curated or ClinVar evidence), not inferred
    lof_genes = set(g for g, (ok, lab_) in lof_status.items()
                    if ok and lab_.startswith(("ClinGen_HI", "recessive_gene", "ClinVar_null")))
    print("  LoF mechanism: %d genes established for PVS1 (ClinGen HI 3/30, recessive, or ClinVar "
          ">= %d null P/LP); %d more plausible (LOEUF / no data); %d not the mechanism%s" % (
              len(lof_genes), a.clinvar_lof_min,
              sum(1 for g, (ok, _) in lof_status.items() if ok and g not in lof_genes),
              sum(1 for ok, _ in lof_status.values() if not ok),
              "" if cv_lof is not None else " [no ClinVar: ClinVar rules off]"), file=sys.stderr)
    # established gene-disease relationship (for PVS1): ClinGen Definitive/Strong, or not
    # curated by ClinGen but PanelApp green. ClinGen Moderate/Limited/Disputed -> not established.
    validity = vc.load_clingen_validity(a.clingen_validity)
    diseases = load_clingen_diseases(a.clingen_validity)          # display only
    mim_sym, mim_ens = load_mim2gene(a.mim2gene)                    # display only
    green = set(g for g in genes if genes[g].get("confidence", "green") == "green")
    if validity is None:
        established = green
    else:
        established = set(g for g in genes if validity.get(g) in ("definitive", "strong")
                          or (g not in validity and g in green))
    print("  established gene-disease relationship (PVS1 full strength): %d of %d genes" % (
          len(established), len(genes)), file=sys.stderr)
    vcep = set(vc.load_gene_list(a.vcep_genes)) if a.vcep_genes else set()
    vcep_specs = load_vcep_specs(a.vcep_specs)
    vcep |= set(g for g, v in vcep_specs.items() if v["pp3_cap"])
    if vcep:
        print("  VCEP genes (in-silico capped at Supporting): %d" % len(vcep), file=sys.stderr)
    fasta = None
    if a.reference_fasta:
        if os.path.isfile(a.reference_fasta + ".fai"):
            fasta = Fasta(a.reference_fasta)
        else:
            print("  [WARN] no .fai for %s: homopolymer check off" % a.reference_fasta, file=sys.stderr)
    clf = vc.Classifier(lof_genes=lof_genes, ba1_exceptions=vc.load_ba1_exceptions(a.ba1_exceptions),
                        cfg={"pm2_points": a.pm2_points, "pm4_repeat_min_copies": a.pm4_repeat_min_copies}, established_genes=established, vcep_genes=vcep,
                        gene_thresholds=dict((g, (v["ba1"], v["bs1"])) for g, v in vcep_specs.items()))

    funnel = Counter()
    sai_extra = load_spliceai_extra(a.spliceai_extra)
    sai_other_gene = set()     # variants whose only SpliceAI >= threshold belongs to another gene
    callers = [("clair3", a.clair3_vep, a.clair3_vcf)]
    if a.dv_vep:
        callers.append(("deepvariant", a.dv_vep, a.dv_vcf))

    # ---- stages 0-3 per caller, streaming ---------------------------------------
    cand = {}           # key -> record (merged across callers)
    research = []
    research_seen = set()
    hla = []
    for cname, vep, vcf in callers:
        n_all = n_gene = 0
        for key, rows in stream_vep(vep):
            n_all += 1
            sai_src = "precomputed" if any((x.get("SpliceAI_pred") or "") not in ("", "-") for x in rows) else "none"
            if sai_src == "none" and sai_extra.get(key):
                rows = [dict(rows[0], SpliceAI_pred=sai_extra[key])] + rows[1:]
                sai_src = "computed_SpliceAI_v1.3"
            r, gene = pick_row(rows, genes, alias)
            r, revel_src = revel_fill(r, rows)
            # SpliceAI only for the gene this variant is assigned to (its symbol or a panel alias)
            sym = r.get("SYMBOL", "")
            ok = {g for g in (gene, sym) if g} | {k for k, v in alias.items() if gene and v == gene}
            sai = spliceai_max(rows, ok)
            sai_any = spliceai_max(rows)
            if sai_any is not None and sai_any >= a.spliceai_supporting and (sai is None or sai < a.spliceai_supporting):
                sai_other_gene.add(key)
            csq = vc.consequences(r)
            impact = r.get("IMPACT", "")
            coding = impact in ("HIGH", "MODERATE") or (sai is not None and sai >= a.spliceai_supporting)
            if not gene:
                # research track: outside the panel but strong evidence
                if impact == "HIGH" or (sai is not None and sai >= a.spliceai_primary):
                    if key not in research_seen:
                        research_seen.add(key)
                        research.append((key, r, sai, cname))
                continue
            n_gene += 1
            ncrna = False
            if gene in NCRNA_DISEASE_GENES:
                dist = vc.to_float(r.get("DISTANCE"))
                ncrna = bool(csq & NCRNA_CSQ) or ("upstream_gene_variant" in csq and dist is not None
                                                  and dist <= a.ncrna_promoter_bp)
                coding = coding or ncrna
            rec = cand.get(key)
            if rec is None:
                rec = {"key": key, "row": r, "gene": gene, "spliceai_max": sai, "callers": set(),
                       "coding": coding, "ncrna": ncrna, "geno": {}, "revel_src": revel_src,
                       "spliceai_src": sai_src}
                cand[key] = rec
            rec["callers"].add(cname)
        funnel["0_%s_variants_in_vep" % cname] = n_all
        funnel["1_%s_in_panel_genes" % cname] = n_gene
        print("  %s: %d variants, %d in panel genes" % (cname, n_all, n_gene), file=sys.stderr)
    funnel["1_union_in_panel_genes"] = len(cand)
    funnel["1e_spliceai_only_for_another_gene"] = len(sai_other_gene)
    funnel["1c_revel_from_other_transcript"] = sum(1 for x in cand.values() if x.get("revel_src") == "other_transcript_same_aa_change")
    funnel["1d_spliceai_scored"] = sum(1 for x in cand.values() if x.get("spliceai_max") is not None)
    funnel["1d_spliceai_computed_here"] = sum(1 for x in cand.values() if x.get("spliceai_src", "").startswith("computed"))
    funnel["1d_spliceai_none"] = sum(1 for x in cand.values() if x.get("spliceai_max") is None)
    if len(cand) > 1000 and not funnel["1d_spliceai_scored"]:
        print("  [WARN] no variant has a SpliceAI score: was the VEP SpliceAI plugin loaded? "
              "(see the VEP *_warnings.txt)", file=sys.stderr)

    # ClinVar for every in-gene variant (so non-coding ClinVar P/LP is not lost)
    cv = clinvar_lookup(list(cand.keys()), a.clinvar_vcf, a.bcftools, mask)
    funnel["1b_clinvar_allele_matches"] = len(cv)
    missense = {k: rec["row"] for k, rec in cand.items()
                if "missense_variant" in vc.consequences(rec["row"])}
    codon_cv = codon_lookup(missense, a.clinvar_vcf, a.bcftools, mask)
    funnel["1c_clinvar_same_codon_PS1_PM5"] = len(codon_cv)

    # genotypes for candidates only
    for cname, vep, vcf in callers:
        gts = read_vcf_genotypes(vcf)
        for k, rec in cand.items():
            if k in gts:
                rec["geno"][cname] = gts[k]
    # a caller supports a variant only with a PASS, non-reference call (a Clair3 LowQual
    # record does not count as support)
    for k, rec in cand.items():
        rec["vep_callers"] = set(rec["callers"])
        rec["callers"] = {c for c in rec["vep_callers"] if pass_call(rec["geno"].get(c))}
        rec["calls_not_pass"] = ",".join(sorted("%s:%s" % (c, (rec["geno"].get(c) or {}).get("filter", "no_record"))
                                                for c in rec["vep_callers"] - rec["callers"]))

    # ---- stages 2-6 -----------------------------------------------------------
    kept = []
    stop = Counter()
    all_rows = []
    for k, rec in cand.items():
        r = rec["row"]
        gene = rec["gene"]
        g = genes[gene]
        moi = g.get("moi", "unknown")
        inh = vc.inheritance_class(moi)
        if k in mask:
            rec["clin"] = vc.clinvar_record("", "", source="masked_for_testing")
        else:
            rec["clin"] = cv.get(k) or vc.clinvar_from_vep_rsid(r.get("CLIN_SIG", ""))
        rec["codon"] = codon_cv.get(k)
        clin = rec["clin"]
        reason = ""
        if vc.is_hla_mhc(gene, k[0], k[1]):
            reason = "2_hla_mhc"
        elif not rec["coding"] and not clin["clinvar_rescue"]:
            reason = "3_consequence"
        else:
            # caller support + QC (one shared function, see qc_stage)
            reason, rec["best_geno"], hp_ = qc_stage(rec["geno"], rec["callers"], len(callers),
                                                      clin["clinvar_rescue"], fasta, k, a)
            if a.phased_vcf and rec["best_geno"] is not None:
                # phase only from the pipeline's own phasing (WhatsHap on the SNVs, indels
                # placed from the haplotagged reads: indel_phase.py); quality stays from the caller.
                # A variant without a phase there loses any phase its caller wrote, so every
                # variant is compared within ONE set of phase blocks.
                ph = PHASED.get(k)
                q = dict(rec["best_geno"])
                if ph and "|" in ph["gt"] and ph["ps"] not in ("", "."):
                    q["gt"], q["ps"], q["phase_from"] = ph["gt"], ph["ps"], "signalscope_phase"
                else:
                    q["gt"], q["ps"], q["phase_from"] = (q.get("gt") or "").replace("|", "/"), "", "unphased"
                rec["best_geno"] = q
            if hp_ is not None:
                rec["homopolymer"] = hp_
        rec["pre_reason"] = reason

    # ---- direct gnomAD lookup for variants that reached the frequency stage ---------------
    freq_keys = [k for k, rec in cand.items() if not rec["pre_reason"]]
    gn_found, gn_absent = ({}, set())
    if a.gnomad_base or a.gnomad_cache:
        print("  gnomAD direct lookup: %d variants reach the frequency stage" % len(freq_keys), file=sys.stderr)
        gn_found, gn_absent = gnomad_lookup(freq_keys, a.gnomad_base, a.bcftools, a.gnomad_cache,
                                            a.gnomad_threads)
        print("  gnomAD direct lookup: %d found, %d confirmed absent, %d unknown" % (
              len(gn_found), len(gn_absent), len(freq_keys) - len(gn_found) - len(gn_absent)), file=sys.stderr)
    funnel["5a_gnomad_direct_found"] = len(gn_found)
    funnel["5a_gnomad_direct_absent"] = len(gn_absent)

    for k, rec in cand.items():
        r = rec["row"]
        gene = rec["gene"]
        g = genes[gene]
        moi = g.get("moi", "unknown")
        inh = vc.inheritance_class(moi)
        clin = rec["clin"]
        reason = rec["pre_reason"]
        if not reason:
            if True:            # stage 5: frequency
                pmax, src = vc.popmax_af(r)
                gn = gn_found.get(k)
                if gn is not None:
                    vep_pm = pmax
                    pmax, src = gn["popmax"], "gnomAD_v4.1_direct(%s)" % gn["datasets"]
                    rec["gnomad_nhomalt"] = gn["nhomalt"]
                    rec["gnomad_filters"] = gn["filters"]
                    if vep_pm is None and (pmax or 0) > 0:
                        funnel["5b_vep_frequency_missing_but_in_gnomAD"] += 1
                elif k in gn_absent:
                    if pmax is None:
                        src = "absent_gnomAD_v4.1_direct"
                    rec["gnomad_nhomalt"] = 0
                rec["popmax"], rec["popmax_src"] = pmax, src
                # frequency used for every decision: GroupMax FAF95 when the direct lookup has it
                # (gnomAD leaves it empty for very low counts: then popmax AF), else popmax AF
                if gn is not None and gn.get("faf") is not None:
                    rec["freq_used"], rec["freq_basis"] = gn["faf"], "gnomAD_GroupMax_FAF95"
                else:
                    rec["freq_used"], rec["freq_basis"] = pmax, ("popmax_AF" if pmax is not None else src)
                af = rec["freq_used"] or 0.0
                exp = a.af_recessive if inh == "recessive" else a.af_dominant
                if af > exp:
                    rec["freq_flag"] = "above_expected_for_%s(>%g)" % (
                        "recessive" if inh == "recessive" else "dominant", exp)
                thr = a.af_filter
                vb = vcep_specs.get(gene, {}).get("ba1")
                if vb is not None:
                    thr = min(thr, vb)                     # VCEP BA1 below 2%: filter there
                exc = (gene, vc.hgvsc_short(r.get("HGVSc", ""))) in clf.ba1_exceptions
                if clin["clinvar_rescue"]:
                    thr = max(thr, a.af_clinvar)
                if af > thr and not exc:
                    reason = "5_frequency"
        if reason:
            stop[reason] += 1
            if reason == "2_hla_mhc":
                hla.append(rec)
            rec["stop"] = reason
            all_rows.append(rec)
            continue

        # provisional evidence points (needed for stage 6a and for ranking)
        sai = rec["spliceai_max"]
        cd = rec.get("codon")
        crow = dict(r)
        crow["spliceai_max"] = "" if sai is None else "%.4f" % sai
        crow["codon_evidence"] = cd["code"] if cd else ""
        if vc.consequences(r) & {"inframe_insertion", "inframe_deletion"}:
            rec["repeat_copies"] = tandem_repeat_copies(fasta, k[0], k[1], k[2], k[3])
            crow["tandem_repeat_copies"] = rec["repeat_copies"]
        if rec.get("freq_used") is not None:
            crow["freq_override"], crow["freq_override_src"] = rec["freq_used"], rec["freq_basis"]
        lab = clf.classify(crow, moi)

        # ---- stage 6a: benign by BOTH well-reviewed ClinVar and our own points -> removed;
        #      ClinVar benign but our points disagree -> kept, flagged for review ----
        cv_benign = clin["clinvar_source"] == "clinvar_vcf" and clin["clinvar"] in ("benign", "likely_benign") \
            and (clin["clinvar_stars"] >= 2 or (clin["clinvar_stars"] == 1
                                                and "conflict" not in clin["clinvar_raw"].lower()))
        # an expert-panel or practice-guideline (3-4 star) Benign/Likely benign classification
        # removes the variant whatever our points
        if cv_benign and clin["clinvar_stars"] >= 3:
            stop["6_benign_expert_panel"] += 1
            rec["stop"] = "6_benign_expert_panel"
            rec.update(lab)
            all_rows.append(rec)
            continue
        if cv_benign:
            if lab["acmg_class"] in ("Benign", "Likely Benign"):
                stop["6_benign_clinvar_and_points"] += 1
                rec["stop"] = "6_benign_clinvar_and_points"
                rec.update(lab)
                all_rows.append(rec)
                continue
            rec["review_flag"] = "ClinVar_%s(%s*)_but_points_%+d" % (clin["clinvar"], clin["clinvar_stars"],
                                                                    lab["acmg_points"])

        # ---- stage 6: evidence labels ----
        csq = vc.consequences(r)
        revel = vc.to_float(r.get("REVEL"))
        prim, supp = [], []
        if clin["clinvar_rescue"]:
            prim.append("ClinVar_P/LP(%s*)" % clin["clinvar_stars"])
        if csq & NULL_CSQ:
            # null variant is PRIMARY only where loss of function is a disease mechanism
            ok, why = lof_status[gene]
            (prim if ok else supp).append("null_variant(%s)" % why)
        if sai is not None and sai >= a.spliceai_primary:
            prim.append("SpliceAI>=%.1f" % a.spliceai_primary)
        elif sai is not None and sai >= a.spliceai_supporting:
            supp.append("SpliceAI>=%.1f" % a.spliceai_supporting)
        if "missense_variant" in csq and revel is not None and revel >= a.revel_supporting:
            supp.append("REVEL>=%.3f" % a.revel_supporting)
        if clin["clinvar"] in ("pathogenic", "likely_pathogenic") and not clin["clinvar_rescue"] \
                and clin["clinvar_source"] == "clinvar_vcf":
            supp.append("ClinVar_P/LP_low_review")
        if clin["clinvar"] in ("VUS", "conflicting") and clin["clinvar_source"] == "clinvar_vcf":
            supp.append("ClinVar_%s" % clin["clinvar"])
        if cd and not clin["clinvar_rescue"]:
            lab_ = "ClinVar_%s(%s,%s*,id%s)" % ("same_AA" if cd["code"] == "PS1" else "same_codon",
                                                cd["code"], cd["stars"], cd["ids"])
            if cd["code"] == "PS1" or (cd["stars"] >= 2 and revel is not None
                                         and revel >= a.pm5_primary_min_revel):
                prim.append(lab_)
            else:
                supp.append(lab_)
        if rec.get("ncrna"):
            supp.append("RNA_gene_exon_or_promoter")
        if csq & vc.PM4_CONSEQUENCES and not (csq & NULL_CSQ):
            if rec.get("repeat_copies", 0) >= a.pm4_repeat_min_copies:
                rec["review_flag"] = (rec.get("review_flag", "") + "|" if rec.get("review_flag") else "") + \
                    "PM4_not_applied(tandem_repeat_%dx)" % rec["repeat_copies"]
            else:
                supp.append("inframe_or_stop_lost(PM4)")
        rec["primary"], rec["supporting"] = prim, supp
        if not prim and not supp:
            stop["6_no_evidence"] += 1
            rec["stop"] = "6_no_evidence"
            rec.update(lab)
            all_rows.append(rec)
            continue
        rec.update(lab)
        kept.append(rec)
    funnel.update({"stop_" + k: v for k, v in stop.items()})
    funnel["6_with_evidence"] = len(kept)

    # ---- stage 7: inheritance fit ---------------------------------------------------
    dom_override = {} if a.no_dominant_null_override else dict(DOMINANT_NULL_OVERRIDE)
    by_gene = defaultdict(list)
    for rec in kept:
        by_gene[rec["gene"]].append(rec)
    tier1, carriers, others = [], [], []
    for gene, recs in by_gene.items():
        moi = genes[gene].get("moi", "unknown")
        inh = vc.inheritance_class(moi)
        on_x = recs[0]["key"][0] in ("X", "chrX")
        for rec in recs:
            gt = rec.get("best_geno", {}).get("gt", "")
            z = vc.zygosity(gt)
            if on_x and a.sex == "male" and z == "hom":
                z = "hemi"
            rec["zygosity"] = z
            if z == "hom" and rec.get("gnomad_nhomalt") == 0:
                rec["hom_flag"] = "hom_not_seen_in_gnomAD(check_artefact_or_consanguinity)"
        prim_recs = [x for x in recs if x["primary"]]
        dominant_ok, recessive_ok = vc.moi_models(moi)     # whole-token parsing (XLR is not dominant)
        override_note = dom_override.get(gene, "")
        if override_note:
            # cited gene override: monoallelic disease is established although ClinGen HI is
            # not 3 (the HI score describes the recessive form); both models stay open
            dominant_ok = recessive_ok = True
        pair_note = ""
        if recessive_ok and len(recs) >= 2 and prim_recs:
            pair_note = _phase_relation(recs)
        # genes with both dominant and recessive disease: a null variant alone supports the
        # DOMINANT model only when haploinsufficiency is established (ClinGen HI 3); otherwise
        # the LoF evidence may come from the recessive disease (e.g. COL6A2: Bethlem AD is
        # missense / in-frame, Ullrich AR is LoF)
        mixed = dominant_ok and recessive_ok and inh != "unknown" and not override_note
        for rec in recs:
            fit = ""
            if rec["primary"]:
                dom_af_ok = (rec.get("freq_used") or 0.0) <= a.af_dominant_fit or rec["clin"]["clinvar_rescue"]
                null_only = all(p_.startswith("null_variant") for p_ in rec["primary"])
                dom_by_null_ok = not (mixed and null_only and clingen_hi.get(gene, "") != "3")
                if dominant_ok and dom_af_ok and dom_by_null_ok and \
                        rec["zygosity"] in ("het", "het_multi", "hom", "hemi"):
                    fit = "dominant"
                    if override_note and rec["zygosity"] in ("het", "het_multi"):
                        rec["review_flag"] = ((rec.get("review_flag", "") + "|") if rec.get("review_flag")
                                              else "") + override_note
                if recessive_ok and rec["zygosity"] in ("hom", "hemi"):
                    fit = fit or "recessive_homozygous"
            if not fit and recessive_ok and len(recs) >= 2 and prim_recs and pair_note != "cis":
                fit = "recessive_two_hits(%s)" % (pair_note or "phase_unknown")
            rec["moi_fit"] = fit
            # caller discordance policy: when two callers were run, a call made by only one of
            # them is never silently dropped but never ranked with concordant calls either:
            # it goes to Tier 3 with its inheritance fit kept, for confirmation (IGV / orthogonal).
            # Exception: rescue-level ClinVar P/LP alleles keep their tier, flagged.
            single = len(callers) == 2 and len(rec["callers"]) < 2
            if single and rec["clin"]["clinvar_rescue"]:
                rec["review_flag"] = ((rec.get("review_flag", "") + "|") if rec.get("review_flag") else "") + \
                    "single_caller(confirm)"
                single = False
            if single and (fit or (rec["primary"] and recessive_ok)):
                rec["moi_fit"] = (fit or "single_het_in_recessive_gene") + "|single_caller_confirm"
                rec["triage"] = "Tier3_review"
                others.append(rec)
            elif fit:
                rec["triage"] = "Tier1_high_priority"
                tier1.append(rec)
            elif rec["primary"] and recessive_ok:
                rec["moi_fit"] = "single_het_in_recessive_gene"
                rec["triage"] = "Tier2_carrier_second_hit_search"
                carriers.append(rec)
            else:
                rec["moi_fit"] = "no_inheritance_fit" if rec["primary"] else "supporting_evidence_only"
                rec["triage"] = "Tier3_review"
                others.append(rec)
    # second-hit search for Tier 2 carriers: every other rare, QC-passing variant in the same
    # gene (with or without evidence), with its phase relative to the carrier variant
    pool_by_gene = defaultdict(list)
    for x in kept + [x for x in all_rows if x.get("stop") == "6_no_evidence"]:
        pool_by_gene[x["gene"]].append(x)
    n_trans = 0
    for rec in carriers:
        cands = []
        for x in pool_by_gene.get(rec["gene"], []):
            if x is rec:
                continue
            rel = _phase_relation([rec, x]) or "phase_unknown"
            n_trans += int(rel == "trans")
            r_ = x["row"]
            cands.append("%s_%d_%s/%s(%s;%s;%s)" % (x["key"] + (
                (r_.get("HGVSp") or r_.get("Consequence") or "").split(":")[-1], rel,
                "evidence=" + ("|".join(x.get("supporting", [])) or "none"))))
        rec["second_hits"] = "; ".join(cands)
    funnel["7b_carriers_with_second_variant_in_trans"] = n_trans
    funnel["7_clinical_shortlist"] = len(tier1)
    funnel["7_carriers"] = len(carriers)
    funnel["7_evidence_but_no_fit"] = len(others)
    funnel["research_outside_panel"] = len(research)
    funnel["not_interpretable_hla"] = len(hla)

    # ---- write ----------------------------------------------------------------
    cols = ["gene", "variant", "chrom", "pos", "ref", "alt", "consequence", "impact", "hgvsc", "hgvsp",
            "exon_intron", "moi", "zygosity", "gt", "gq", "depth", "af", "phase_set", "callers",
            "popmax_af", "popmax_source", "revel", "spliceai_max",
            "clinvar", "clinvar_stars", "clinvar_source", "clinvar_raw",
            "primary_evidence", "supporting_evidence", "moi_fit",
            "triage_tier", "evidence_points", "evidence_codes", "panels", "stop_reason", "Existing_variation",
            "frequency_flag", "nmd", "loeuf", "pli", "mis_z",
            "constraint_flags", "homopolymer_run", "gnomad_nhomalt", "gnomad_filters", "hom_flag",
            "freq_used", "freq_basis", "review_flag", "lof_mechanism", "second_hit_candidates", "revel_source",
            "transcript", "calls_not_pass", "spliceai_source", "repeat_copies",
            # (display only, for curator review; not used by the tree)
            "variant_type", "vcf_variant", "clinvar_id", "sift", "polyphen", "gnomad_exome_af",
            "gnomad_genome_af", "max_af", "max_af_pops", "omim_gene", "clingen_disease"]
    def out_row(rec):
        r, k, q, c = rec["row"], rec["key"], rec.get("best_geno", {}), rec.get("clin", {})
        g = genes.get(rec["gene"], {})
        ex = r.get("EXON") or ""
        return {
            "gene": rec["gene"], "variant": "%s_%d_%s/%s" % k, "chrom": k[0], "pos": k[1],
            "ref": k[2], "alt": k[3], "consequence": r.get("Consequence"), "impact": r.get("IMPACT"),
            "hgvsc": r.get("HGVSc"), "hgvsp": r.get("HGVSp"),
            "exon_intron": exon_intron_label(r),
            "moi": g.get("moi", ""), "zygosity": rec.get("zygosity", ""), "gt": q.get("gt", ""),
            "gq": q.get("gq", ""), "depth": q.get("dp", ""), "af": q.get("af", ""), "phase_set": q.get("ps", ""),
            "callers": "+".join(sorted(rec.get("callers", []))),
            "popmax_af": "" if rec.get("popmax") is None else "%.3g" % rec["popmax"],
            "popmax_source": rec.get("popmax_src", ""), "revel": r.get("REVEL", ""),
            "spliceai_max": "" if rec.get("spliceai_max") is None else "%.2f" % rec["spliceai_max"],
            "clinvar": c.get("clinvar", ""), "clinvar_stars": c.get("clinvar_stars", ""),
            "clinvar_source": c.get("clinvar_source", ""), "clinvar_raw": c.get("clinvar_raw", ""),
            "primary_evidence": "|".join(rec.get("primary", [])),
            "supporting_evidence": "|".join(rec.get("supporting", [])),
            "moi_fit": rec.get("moi_fit", ""), "triage_tier": rec.get("triage", "not_prioritised" if rec.get("stop") else ""),
            "evidence_points": rec.get("acmg_points", ""), "evidence_codes": rec.get("acmg_evidence", ""),
            "panels": g.get("panels", ""), "stop_reason": rec.get("stop", ""),
            "Existing_variation": r.get("Existing_variation", ""),
            "frequency_flag": rec.get("freq_flag", ""), "nmd": r.get("NMD", ""),
            "loeuf": constraint.get(rec["gene"], ("", "", ""))[0],
            "pli": constraint.get(rec["gene"], ("", "", ""))[1],
            "mis_z": constraint.get(rec["gene"], ("", "", ""))[2],
            "constraint_flags": constraint.get(rec["gene"], ("", "", "", ""))[3],
            "homopolymer_run": rec.get("homopolymer", ""),
            "gnomad_nhomalt": rec.get("gnomad_nhomalt", ""),
            "gnomad_filters": rec.get("gnomad_filters", ""),
            "hom_flag": rec.get("hom_flag", ""),
            "revel_source": rec.get("revel_src", ""),
            # the VEP transcript of the row used (Feature), so the region check compares the
            # same transcript (rows without HGVSc do not name it)
            "transcript": (r.get("Feature") or "").split(".")[0],
            "calls_not_pass": rec.get("calls_not_pass", ""),
            "spliceai_source": rec.get("spliceai_src", ""),
            "repeat_copies": rec.get("repeat_copies", ""),
            "freq_used": "" if rec.get("freq_used") is None else "%.3g" % rec["freq_used"],
            "freq_basis": rec.get("freq_basis", ""),
            "review_flag": rec.get("review_flag", ""),
            "lof_mechanism": lof_status.get(rec["gene"], ("", ""))[1],
            "second_hit_candidates": rec.get("second_hits", ""),
            "variant_type": variant_type(k[2], k[3]),
            "vcf_variant": q.get("vcf", ""),
            "clinvar_id": c.get("clinvar_id", "") if c.get("clinvar_source") == "clinvar_vcf" else "",
            "sift": r.get("SIFT", ""), "polyphen": r.get("PolyPhen", ""),
            "gnomad_exome_af": r.get("gnomADe_AF", ""), "gnomad_genome_af": r.get("gnomADg_AF", ""),
            "max_af": r.get("MAX_AF", ""), "max_af_pops": r.get("MAX_AF_POPS", ""),
            "omim_gene": mim_ens.get((r.get("Gene") or "").split(".")[0]) or mim_sym.get(rec["gene"], ""),
            "clingen_disease": diseases.get(rec["gene"], ""),
        }
    def write(name, recs, sort=True):
        if sort:
            recs = sorted(recs, key=lambda x: (-int(x.get("acmg_points", 0) or 0), x["gene"]))
        path = "%s.%s.tsv" % (pre, name)
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for rec in recs:
                o = out_row(rec)
                fh.write("\t".join(str(o.get(c_, "")) for c_ in cols) + "\n")
        print("  wrote %s (%d rows)" % (path, len(recs)), file=sys.stderr)
    write("clinical_shortlist", tier1)
    write("carriers", carriers)
    write("not_interpretable", hla, sort=False)
    write("all_candidates", kept + [x for x in all_rows if x.get("stop") != "2_hla_mhc"])
    # research
    rpath = "%s.research.tsv" % pre
    with open(rpath, "w") as fh:
        fh.write("gene\tvariant\tconsequence\timpact\thgvsp\tspliceai_max\tpopmax_af\tcaller\tnote\n")
        n_research = 0
        for key, r, sai, cname in research:
            pm, _ = vc.popmax_af(r)
            if pm is not None and pm > a.af_recessive:
                continue
            if vc.is_hla_mhc(r.get("SYMBOL", ""), key[0], key[1]):
                continue   # HLA/MHC variants are not interpretable, also outside the panel
            n_research += 1
            fh.write("%s\t%s_%d_%s/%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                r.get("SYMBOL", ""), key[0], key[1], key[2], key[3], r.get("Consequence"), r.get("IMPACT"),
                r.get("HGVSp"), "" if sai is None else "%.2f" % sai, "" if pm is None else "%.3g" % pm, cname,
                "outside selected panel(s); speculative, not for clinical use"))
    funnel["research_outside_panel"] = n_research   # the rows written (rare, outside HLA/MHC)
    with open("%s.funnel.tsv" % pre, "w") as fh:
        fh.write("stage\tcount\n")
        for k_ in sorted(funnel):
            fh.write("%s\t%d\n" % (k_, funnel[k_]))
    print("  funnel: " + ", ".join("%s=%d" % (k_, funnel[k_]) for k_ in sorted(funnel)), file=sys.stderr)


def qc_stage(geno, callers_pass, n_callers_run, clinvar_rescue, fasta, key, a):
    """Stage 4 (caller support + genotype QC), one function used by the tree AND by
    the benchmark call-set builder, so the benchmarked call set follows exactly the same rules.
    geno: {caller: genotype dict from read_vcf_genotypes}; callers_pass: callers with a PASS,
    non-reference call; n_callers_run: callers given to the run. Returns
    (reason or '', best genotype dict, homopolymer run or None when not evaluated)."""
    q = None
    for cname in ("clair3", "deepvariant"):
        gg = geno.get(cname) if cname in callers_pass else None
        if gg and (q is None or _num(gg["gq"]) > _num(q["gq"])):
            q = gg
    if q is not None and "|" not in q.get("gt", ""):
        # phase from the other caller when only it is phased (panel mode: DeepVariant
        # VCFs carry WhatsHap PS tags, Clair3 VCFs do not); quality stays from q
        for cname in ("deepvariant", "clair3"):
            gg = geno.get(cname) if cname in callers_pass else None
            if gg and "|" in gg.get("gt", "") and gg.get("ps") not in ("", "."):
                q = dict(q)
                q["gt"], q["ps"], q["phase_from"] = gg["gt"], gg["ps"], cname
                break
    gq, dp = _num((q or {}).get("gq")), _num((q or {}).get("dp"))
    both = len(callers_pass) >= 2 or n_callers_run == 1
    hp = None
    if not callers_pass:
        reason = "4_no_pass_call" if geno else "4_no_genotype"
    elif q is None:
        reason = "4_no_genotype"
    elif q.get("filter") not in ("PASS", "."):
        reason = "4_filter_not_pass"
    elif gq < (a.min_gq_both if (getattr(a, "min_gq_both", None) is not None and len(callers_pass) >= 2)
               else a.min_gq) or dp < a.min_dp:
        reason = "4_low_gq_or_depth"
    elif not both and gq < a.single_caller_min_gq and not clinvar_rescue:
        reason = "4_single_caller_low_gq"
    else:
        reason = ""
        hp = homopolymer_run(fasta, key[0], key[1], key[2], key[3], a.homopolymer_max_indel)
        if hp >= a.homopolymer_min and not clinvar_rescue and (not both or gq < a.homopolymer_min_gq):
            reason = "4_indel_homopolymer_low_support"
    return reason, q or {}, hp


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return -1.0


def _phase_relation(recs):
    """For two or more variants in one gene: 'trans' / 'cis' / '' (unknown), from
    phased genotypes sharing a phase set (Clair3/WhatsHap PS tag)."""
    phased = []
    for r in recs:
        q = r.get("best_geno", {})
        gt, ps = q.get("gt", ""), q.get("ps", "")
        if "|" in gt and ps not in ("", "."):
            phased.append((ps, gt))
    by_ps = defaultdict(set)
    for ps, gt in phased:
        by_ps[ps].add(gt)
    for ps, gts in by_ps.items():
        if {"0|1", "1|0"} <= gts:
            return "trans"
        if len(gts) == 1 and len(phased) >= 2:
            return "cis"
    return ""


if __name__ == "__main__":
    main()
