# =============================================================================
# interpretation_paths.smk - file names of the decision-tree layer (Layer 3, panel mode).
# Included by the Snakefile BEFORE 'rule all' (no rules here, names only), so the
# targets can be listed; the rules are in rules/interpretation.smk (included at the end).
# =============================================================================
DT_CFG       = config.get("decision_tree", {})
TREE_ENABLED = MODE == "panel" and bool(DT_CFG.get("enabled", True))
RES          = config.get("resources", {})

TREE_DIR     = os.path.join(SAMPLE_ROOT, "interpretation")
NORM_DIR     = os.path.join(TREE_DIR, "vcf_norm")
TVEP_DIR     = os.path.join(TREE_DIR, "vep")
CLASS_DIR    = os.path.join(TREE_DIR, "classify")
SVD_DIR      = os.path.join(TREE_DIR, "sv_dosage")

NORM_VCF     = os.path.join(NORM_DIR, f"{SAMPLE}" + ".{caller}.norm.vcf.gz")
TVEP_TSV     = os.path.join(TVEP_DIR, f"{SAMPLE}" + ".{caller}.norm.vep.tsv")
GENE_TABLE   = os.path.join(TREE_DIR, f"{SAMPLE}.gene_table.tsv")
SAI_EXTRA    = os.path.join(CLASS_DIR, f"{SAMPLE}.spliceai_computed.tsv")
CANDIDATES   = os.path.join(CLASS_DIR, f"{SAMPLE}.all_candidates.tsv")
EXON_BED     = os.path.join(CLASS_DIR, f"{SAMPLE}.coding_exons_pad20.bed")
EXON_COV     = os.path.join(CLASS_DIR, f"{SAMPLE}.exon_coverage.tsv")
SV_DOSAGE    = os.path.join(SVD_DIR, f"{SAMPLE}.sv_dosage.tsv")
SV_TIERS     = os.path.join(SVD_DIR, f"{SAMPLE}.sv_tiers.tsv")
SV_PHASE     = os.path.join(SVD_DIR, f"{SAMPLE}.sv_phase.tsv")
PHASE_DIR_T  = os.path.join(TREE_DIR, "phase")
PHASED_VCF   = os.path.join(PHASE_DIR_T, f"{SAMPLE}.phased.vcf.gz")
HAPLOTAG_BAM = os.path.join(PHASE_DIR_T, f"{SAMPLE}.targeted.haplotagged.bam")
REGIONS_TSV  = os.path.join(CLASS_DIR, f"{SAMPLE}.regions.tsv")
MODEL_CHECK  = os.path.join(META_DIR, f"{SAMPLE}.clair3_model_check.txt")
TREE_TOOLS   = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.tool_versions.report.tsv")
PER_BASE     = os.path.join(COV_DIR, f"{SAMPLE}.mosdepth_report.per-base.bed.gz")
# the report is built in two steps: generate_report_html.py writes the base report, which
# report_redesign.py turns into the final report (REPORT_HTML) with the read-length histogram
REPORT_BASE  = os.path.join(REPORT_DIR, "html", f"{SAMPLE}.report.base.html")
READLEN_HIST = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.readlen_hist.tsv")

# outputs of the WGS shortlist path (snv/indel/sv_process), not built in panel mode
OLD_SHORTLIST_TARGETS = [SNV_SUMMARY, SNV_CLINVAR_SUMM, SNV_SHORTLIST_CTX, VEP_TSV, VEP_STATS,
                         DV_VEP_TSV, DV_VEP_STATS]
TREE_TARGETS = [CANDIDATES, EXON_COV, SV_DOSAGE, SV_TIERS, SV_PHASE, PHASED_VCF, REGIONS_TSV, MODEL_CHECK, REPORT_HTML, IGV_MANIFEST]
