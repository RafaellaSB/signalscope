# =============================================================================
# SignalScope — Panel-Agnostic Long-Read Sequencing Pipeline — Snakefile
#
# Oxford Nanopore R10.4.1 / GRCh38 or T2T-CHM13 / Any gene panel
#
# PIPELINE STRUCTURE
# ------------------
# All outputs are written under:
#   {samples_root}/{sample}/
#
# Input POD5 is read from (never written to):
#   {pod5_base}/{sample}/
#
# OUTPUT DIRECTORY LAYOUT PER SAMPLE
# ------------------------------------
#   raw/fastq/          — basecalled FASTQ (gzipped)
#   alignments/         — sorted BAMs and indexes
#   qc/                 — NanoPlot, flagstat, panel depth
#   coverage/           — mosdepth outputs
#   variants/
#     small_variants/
#       clair3/         — raw Clair3 output
#       deepvariant/    — DeepVariant VCF
#     structural_variants/
#       sniffles/       — Sniffles VCF.gz
#   annotation/
#     vcf_inputs/       — panel-filtered VCF ready for VEP
#     vep/              — VEP TSV + stats
#     tables/           — concordance tables
#   logs/               — per-rule stdout and stderr logs
#   tmp/                — intermediate files safe to leave
#   run_metadata/       — software versions
#
# START MODES
# -----------
# The pipeline supports three entry points via config start_from:
#   pod5  — runs Dorado then continues
#   fastq — skips Dorado; reads from config["existing_fastq"]
#   bam   — skips Dorado and alignment; reads from config["existing_bam"]
#
# LAYER STATUS
# ------------
# Layer 1 (core):     ACTIVE — basecall/align/QC/Clair3
# Layer 2 (variant analysis): SNV/indel via Clair3 + DeepVariant; SV via Sniffles2 + CuteSV; per-class concordance
# Layer 3 (clinical interpretation & reporting), panel mode: decision tree (triage tiers with evidence
#   points, no ACMG class), SV dosage, exon coverage, IGV and HTML report -> rules/interpretation.smk.
#   WGS mode uses the shortlist path (snv/indel/sv_process).
#
# USAGE
# -----
#   # Dry run first — always:
#   snakemake -n --reason --config sample=barcodeXX start_from=bam existing_bam=/abs/path/to.bam
#
#   # From existing BAM (most common — skip Dorado + alignment):
#   snakemake --config sample=barcode06 start_from=bam \
#       existing_bam=/abs/path/to/barcode06.primary.bam --cores 16
#
#   # From existing FASTQ (skip Dorado only):
#   snakemake --config sample=barcode06 start_from=fastq \
#       existing_fastq=/abs/path/to/barcode06.fastq.gz --cores 16
#
#   # From POD5 (full run including Dorado basecalling):
#   snakemake --config sample=barcode06 start_from=pod5 --cores 16
#
#   # Stop before variant calling (QC + coverage only):
#   snakemake --config sample=barcode06 start_from=bam \
#       existing_bam=/abs/path/to/barcode06.primary.bam --cores 16 --until mosdepth_coverage
#
#   # Report/preview only (assumes all outputs already exist):
#   snakemake --config sample=barcode06 start_from=bam \
#       existing_bam=/abs/path/to/barcode06.primary.bam --cores 1 \
#       --forcerun generate_report_html
#
#   # DAG visualisation:
#   snakemake --config sample=barcode06 start_from=bam \
#       existing_bam=/abs/path/to/barcode06.primary.bam --rulegraph | dot -Tpng > dag.png
# =============================================================================

configfile: "config.yaml"

import os

# stop early, with the list of settings, when a placeholder from the release config
# (a value written as <...>) has not been replaced. Settings of modules that are switched off
# (enabled: false) are not checked.
def _unfilled(cfg, prefix=""):
    out = []
    if isinstance(cfg, dict):
        if cfg.get("enabled") is False:
            return out
        for k, v in cfg.items():
            out += _unfilled(v, "%s%s." % (prefix, k))
    elif isinstance(cfg, str) and cfg.strip().startswith("<") and cfg.strip().endswith(">"):
        out.append(prefix.rstrip("."))
    return out
_PLACEHOLDERS = [k for k in _unfilled(config)
                 if not (k == "tools.dorado" and config.get("start_from", "bam") != "pod5")]
if _PLACEHOLDERS:
    raise ValueError("[ERROR] these settings still hold a <placeholder> from config.yaml; set them for "
                     "your server or in the sample config:\n  " + "\n  ".join(_PLACEHOLDERS))

# ── Auto-detect input format (SignalScope v3) ────────────────────────────
# If start_from is not explicitly set, detect from input file extension.
import os as _os
def _auto_detect_start(cfg):
    explicit = cfg.get("start_from", "")
    if explicit:
        return explicit
    ip = cfg.get("input_path", "")
    if not ip:
        if cfg.get("existing_bam", ""):
            return "bam"
        elif cfg.get("existing_fastq", ""):
            return "fastq"
        return "pod5"
    if _os.path.isdir(ip):
        return "pod5"
    ext = ip.lower()
    if ext.endswith(".pod5"):
        return "pod5"
    elif ext.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        return "fastq"
    elif ext.endswith(".bam"):
        return "bam"
    return "pod5"



# =============================================================================
# CONFIGURATION RESOLUTION
#
# All path variables are derived here. Rules reference these names only.
# Nothing is hardcoded inside rule bodies.
# =============================================================================

# SAMPLE must be provided on the command line: --config sample=barcodeXX
# The config.yaml deliberately leaves this unset to prevent accidental runs.
SAMPLE = config.get("sample", "")
if not SAMPLE:
    raise ValueError(
        "[ERROR] No sample specified.\n"
        "  Run with: snakemake --config sample=barcodeXX start_from=bam existing_bam=... --cores N"
    )

# --- Per-sample root directory -----------------------------------------------
# Every output file is under SAMPLE_ROOT. This is the only directory the
# pipeline writes to — pod5_pass/ is never touched.
SAMPLE_ROOT = os.path.join(config["samples_root"], SAMPLE)

# Sub-directories under SAMPLE_ROOT (using string constants for clarity)
RAW_DIR    = os.path.join(SAMPLE_ROOT, "raw")
ALIGN_DIR  = os.path.join(SAMPLE_ROOT, "alignments")
QC_DIR     = os.path.join(SAMPLE_ROOT, "qc")
COV_DIR    = os.path.join(SAMPLE_ROOT, "coverage")
VAR_DIR    = os.path.join(SAMPLE_ROOT, "variants")
ANNOT_DIR  = os.path.join(SAMPLE_ROOT, "annotation")
LOG_DIR    = os.path.join(SAMPLE_ROOT, "logs")
TMP_DIR    = os.path.join(SAMPLE_ROOT, "tmp")
META_DIR   = os.path.join(SAMPLE_ROOT, "run_metadata")

# Variant sub-directories
CLAIR3_DIR   = os.path.join(VAR_DIR, "small_variants", "clair3")
SNIFFLES_DIR = os.path.join(VAR_DIR, "structural_variants", "sniffles")
CUTESV_DIR   = os.path.join(VAR_DIR, "structural_variants", "cutesv")

# Annotation sub-directories
VCF_IN_DIR   = os.path.join(ANNOT_DIR, "vcf_inputs")
VEP_DIR      = os.path.join(ANNOT_DIR, "vep")
TABLES_DIR   = os.path.join(ANNOT_DIR, "tables")

# --- Tool binaries ------------------------------------------------------------
DORADO   = config["tools"]["dorado"]
MINIMAP2 = config["tools"]["minimap2"]
SAMTOOLS = config["tools"]["samtools"]
SCRIPTS = os.path.join(workflow.basedir, "scripts")   # anchored to Snakefile location
BCFTOOLS = config["tools"]["bcftools"]
SNIFFLES = config["tools"]["sniffles"]
MOSDEPTH = config["tools"]["mosdepth"]
NANOPLOT = config["tools"]["nanoplot"]
VEP_SIF = config["vep"]["sif"]
# Bind the directories that hold the data VEP needs, derived from the config
# paths (reference, cache, plugins, SIF) so the container works on any host
# without hardcoded mount points. Duplicate/empty entries are de-duplicated.
def _vep_binds():
    import os as _os
    paths = [
        _os.path.dirname(_os.path.abspath(config["ref_fasta"])),
        _os.path.abspath(config["vep"]["cache"]),
        _os.path.abspath(config["vep"].get("plugins", "") or config["vep"]["cache"]),
        _os.path.dirname(_os.path.abspath(VEP_SIF)),
        _os.path.abspath(config["samples_root"]),
    ]
    seen, binds = set(), []
    for d in paths:
        if d and d not in seen:
            seen.add(d)
            binds.append("--bind " + d)
    return " ".join(binds)
VEP = ('apptainer exec --nv --cleanenv ' + _vep_binds() + ' ' + VEP_SIF + ' vep')

# --- Clair3 -------------------------------------------------------------------
CLAIR3_RUNNER   = config["clair3"]["sif"]
CLAIR3_MODELDIR = config["clair3"]["model_dir"]
CLAIR3_MODEL    = config["clair3"]["model"]

# --- Reference: MMI index if it exists, otherwise fall back to FASTA ---------
_ref_mmi_cfg = config.get("ref_mmi", "")
ALIGN_REF = _ref_mmi_cfg if (_ref_mmi_cfg and os.path.isfile(_ref_mmi_cfg)) else config["ref_fasta"]

if ALIGN_REF == config["ref_fasta"]:
    print(f"[INFO] MMI index not found or not set — aligning against ref_fasta")
else:
    print(f"[INFO] Using MMI index: {ALIGN_REF}")

# --- Start mode --------------------------------------------------------------
START_FROM = config.get("start_from", "pod5")

if START_FROM not in ("pod5", "fastq", "bam", "vcf"):
    raise ValueError(f"start_from must be 'pod5', 'fastq', or 'bam' — got '{START_FROM}'")

# Validate required paths for each mode
if START_FROM == "fastq":
    _ef = config.get("existing_fastq", "")
    if not _ef:
        raise ValueError("start_from=fastq requires existing_fastq to be set in config")
    if not os.path.isfile(_ef):
        raise ValueError(f"existing_fastq not found: {_ef}")

if START_FROM == "bam":
    _eb = config.get("existing_bam", "")
    if not _eb:
        raise ValueError("start_from=bam requires existing_bam to be set in config")
    if not os.path.isfile(_eb):
        raise ValueError(f"existing_bam not found: {_eb}")

if START_FROM == "vcf":
    _eb = config.get("existing_bam", "")
    if not _eb:
        raise ValueError("start_from=vcf requires existing_bam")
    for _k in ["existing_clair3_vcf", "existing_dv_vcf"]:
        if not config.get(_k, ""):
            raise ValueError(f"start_from=vcf requires {_k}")

print(f"[INFO] start_from={START_FROM}  sample={SAMPLE}")

# ── Sequencing summary — optional (SignalScope v3) ───────────────────────
# ── VEP cache completeness check ─────────────────────────────────────────────
# Prevents the "cache extracted but incomplete" failure (e.g. interrupted tar):
# VEP silently drops variants on missing chromosomes, causing cryptic downstream
# errors. Fail fast here with a clear, actionable message instead.
_vep_cache = config.get("vep", {}).get("cache", "")
if _vep_cache and _os.path.isdir(_vep_cache):
    _assembly = config.get("vep", {}).get("assembly", "GRCh38")
    # find the species/version_assembly dir (e.g. homo_sapiens/116_GRCh38)
    _hs = _os.path.join(_vep_cache, "homo_sapiens")
    if _os.path.isdir(_hs):
        _subdirs = [d for d in _os.listdir(_hs) if _assembly in d]
        if _subdirs:
            _cache_chr_dir = _os.path.join(_hs, _subdirs[0])
            _expected = set(str(c) for c in range(1, 23)) | {"X", "Y"}
            _present = set(d for d in _os.listdir(_cache_chr_dir)
                           if d in _expected)
            _missing = _expected - _present
            if _missing:
                raise ValueError(
                    f"VEP cache appears INCOMPLETE at {_cache_chr_dir}\n"
                    f"  Found {len(_present)}/24 chromosome directories; "
                    f"missing: {', '.join(sorted(_missing, key=lambda x: (len(x), x)))}\n"
                    f"  This usually means the cache tarball was not fully extracted "
                    f"(e.g. interrupted). Re-extract the VEP cache completely:\n"
                    f"    tar xzf homo_sapiens_vep_*_{_assembly}.tar.gz\n"
                    f"  then re-run. See README 'External resources'."
                )
            print(f"[INFO] VEP cache complete: {len(_present)}/24 chromosomes")

_seq_summary = config.get("sequencing_summary_file", "")
HAS_SEQ_SUMMARY = bool(
    _seq_summary and _os.path.isfile(_seq_summary) and _os.path.getsize(_seq_summary) > 200
)
if not HAS_SEQ_SUMMARY:
    print("[INFO] No sequencing summary - metrics derived from BAM")


# --- Pipeline mode -----------------------------------------------------------
MODE = config.get("mode", "panel")   # "panel" (adaptive sampling) or "wgs"
if MODE not in ("panel", "wgs"):
    raise ValueError(f"[ERROR] Invalid mode '{MODE}'. Must be 'panel' or 'wgs'.")
print(f"[INFO] Pipeline mode: {MODE}")

# WGS-specific filtering overrides (fall back to panel defaults)
_cf = config.get("clinical_filtering", {})
_wgs_cf = _cf.get("wgs", {})
if MODE == "wgs":
    EFF_MAX_AF           = _wgs_cf.get("max_population_af", _cf.get("max_population_af", 0.001))
    EFF_MIN_CALLERS      = _wgs_cf.get("min_caller_support", _cf.get("min_caller_support", 2))
    EFF_MAX_SHORTLIST    = _wgs_cf.get("max_shortlist", 50)
    EFF_MAX_IGV          = _wgs_cf.get("max_igv_screenshots", 30)
    print(f"[INFO] WGS filtering: max_af={EFF_MAX_AF}, min_callers={EFF_MIN_CALLERS}, "
          f"max_shortlist={EFF_MAX_SHORTLIST}, max_igv={EFF_MAX_IGV}")
else:
    EFF_MAX_AF           = _cf.get("max_population_af", 0.01)
    EFF_MIN_CALLERS      = _cf.get("min_caller_support", 1)
    EFF_MAX_SHORTLIST    = 0    # unlimited for panel
    EFF_MAX_IGV          = 0    # unlimited for panel

# --- Active BED --------------------------------------------------------------
TEST_MODE  = config.get("test_mode", False)
ACTIVE_BED = config["test_bed"] if TEST_MODE else config.get("panel_bed", "") or []
if MODE == "wgs":
    ACTIVE_BED = []
WGS_FALLBACK_BED = _os.path.join(_os.path.dirname(workflow.snakefile), "assets", "wgs_placeholders", "whole_genome.bed") if MODE == "wgs" else ""

# --- Ensure the panel BED is contig-order + coordinate sorted (robustness) ---
# An unsorted BED can cause Clair3 to silently return zero candidates on some
# contigs. Sort ACTIVE_BED against the reference .fai order so EVERY consumer
# (Clair3, DeepVariant regions, panel filters) gets a sorted BED. Pure awk/sort,
# no bedtools dependency.
# Skipped in WGS mode (no panel BED).
import subprocess as _sp_bedsort
_raw_active_bed = ACTIVE_BED if ACTIVE_BED else []
_active_fai = config["ref_fasta"] + ".fai"
_harmonizer = _os.path.join(_os.path.dirname(workflow.snakefile), "scripts", "harmonize_and_sort_bed.py")
if ACTIVE_BED and _os.path.isfile(ACTIVE_BED) and _os.path.isfile(_active_fai):
    _sorted_active_bed = _os.path.join(
        _os.path.dirname(_raw_active_bed) or ".",
        _os.path.basename(_raw_active_bed).replace(".bed", "") + ".sorted.bed"
    )
    # Harmonise BED contig naming to the reference (chr/no-chr), then sort by
    # reference contig order. Aborts loudly if the BED can't be reconciled
    # (e.g. wrong assembly) instead of silently emitting an empty BED.
    # The result is written to a temporary file and moved into place only when the content
    # changed, so an unchanged BED keeps its timestamp and does not trigger reruns.
    import filecmp as _filecmp_bed
    _tmp_sorted_bed = _sorted_active_bed + ".tmp"
    _sp_bedsort.run(
        ["python3", _harmonizer, "--bed", _raw_active_bed,
         "--fai", _active_fai, "--out", _tmp_sorted_bed],
        check=True
    )
    if _os.path.isfile(_sorted_active_bed) and _filecmp_bed.cmp(_tmp_sorted_bed, _sorted_active_bed, shallow=False):
        _os.remove(_tmp_sorted_bed)
    else:
        _os.replace(_tmp_sorted_bed, _sorted_active_bed)
    ACTIVE_BED = _sorted_active_bed
    print(f"[INFO] Panel BED harmonised + sorted -> {ACTIVE_BED}")
elif MODE == "wgs":
    print("[INFO] WGS mode — no panel BED to sort")
elif not _os.path.isfile(_active_fai):
    print(f"[WARN] {_active_fai} not found - using panel BED unsorted (Clair3 may miss contigs)")
else:
    print(f"[WARN] Panel BED not found: {ACTIVE_BED}")

if TEST_MODE:
    print(f"[INFO] test_mode=true — region={config['test_region']}  BED={ACTIVE_BED}")

# --- Canonical file paths ----------------------------------------------------
# Defined once here. Rules refer to these constants, not f-strings inline.

FASTQ_GZ     = os.path.join(RAW_DIR, "fastq", f"{SAMPLE}.fastq.gz")
BAM_ALL      = os.path.join(ALIGN_DIR, f"{SAMPLE}.bam")
BAM_ALL_BAI  = os.path.join(ALIGN_DIR, f"{SAMPLE}.bam.bai")
BAM_PRI      = os.path.join(ALIGN_DIR, f"{SAMPLE}.primary.bam")
# Sniffles panel VCF for VEP annotation
SNIFFLES_VCF_PANEL = os.path.join(
    config["samples_root"], SAMPLE,
    "variants", "structural_variants", "sniffles", f"{SAMPLE}.sniffles.panel.vcf.gz"
)
SNIFFLES_VEP_TSV = os.path.join(
    config["samples_root"], SAMPLE, "annotation", "vep", f"{SAMPLE}.sniffles.vep.tsv"
)
BAM_PRI_BAI  = os.path.join(ALIGN_DIR, f"{SAMPLE}.primary.bam.bai")
# BAM used as input to all variant callers: the primary BAM.
_pr = config.get("pseudogene_recovery", {})
if (isinstance(_pr, dict) and _pr.get("enabled")) or config.get("pkd1_recovery", False):
    raise ValueError("[ERROR] pseudogene_recovery is not part of SignalScope; "
                     "set pseudogene_recovery.enabled: false")
BAM_FOR_CALLING     = BAM_PRI
BAM_FOR_CALLING_BAI = BAM_PRI_BAI

FLAGSTAT     = os.path.join(QC_DIR, f"{SAMPLE}.flagstat.txt")
PANEL_READS  = os.path.join(QC_DIR, f"{SAMPLE}.panel_read_count.txt")
DEPTH_SUMM   = os.path.join(QC_DIR, f"{SAMPLE}.panel_depth_summary.txt")
NANOPLOT_HTML = os.path.join(QC_DIR, "nanoplot", f"{SAMPLE}_NanoPlot-report.html")
NANOPLOT_TXT  = os.path.join(QC_DIR, "nanoplot", f"{SAMPLE}_NanoStats.txt")

MOSDEPTH_SUMMARY = os.path.join(COV_DIR, f"{SAMPLE}.mosdepth.summary.txt")
MOSDEPTH_GLOBAL  = os.path.join(COV_DIR, f"{SAMPLE}.mosdepth.global.dist.txt")
MOSDEPTH_REGION  = os.path.join(COV_DIR, f"{SAMPLE}.mosdepth.region.dist.txt")
MOSDEPTH_BED_GZ  = os.path.join(COV_DIR, f"{SAMPLE}.regions.bed.gz")
# REPORT mosdepth output — uses REPORT BED (COL4A3 and COL4A4 separate)
MOSDEPTH_REPORT_BED_GZ = os.path.join(COV_DIR, f"{SAMPLE}.mosdepth_report.regions.bed.gz")
MOSDEPTH_BED_CSI = os.path.join(COV_DIR, f"{SAMPLE}.regions.bed.gz.csi")

CLAIR3_VCF_RAW   = os.path.join(CLAIR3_DIR, "merge_output.vcf.gz")
CLAIR3_VCF_PANEL = os.path.join(VCF_IN_DIR, f"{SAMPLE}.clair3.panel.vcf.gz")
CLAIR3_VCF_TBI   = os.path.join(VCF_IN_DIR, f"{SAMPLE}.clair3.panel.vcf.gz.tbi")


SNIFFLES_VCF       = os.path.join(SNIFFLES_DIR, f"{SAMPLE}.sniffles.vcf.gz")
SNIFFLES_TBI       = os.path.join(SNIFFLES_DIR, f"{SAMPLE}.sniffles.vcf.gz.tbi")
SNIFFLES_VCF_PANEL = os.path.join(SNIFFLES_DIR, f"{SAMPLE}.sniffles.panel.vcf.gz")
SNIFFLES_TBI_PANEL = os.path.join(SNIFFLES_DIR, f"{SAMPLE}.sniffles.panel.vcf.gz.tbi")

# DeepVariant small variant calling (indel second caller)
DEEPVARIANT_DIR  = os.path.join(VAR_DIR, "small_variants", "deepvariant")
DV_VCF_GZ        = os.path.join(DEEPVARIANT_DIR, f"{SAMPLE}.deepvariant.vcf.gz")
DV_VCF_TBI       = os.path.join(DEEPVARIANT_DIR, f"{SAMPLE}.deepvariant.vcf.gz.tbi")
DV_VCF_PANEL     = os.path.join(VCF_IN_DIR,      f"{SAMPLE}.deepvariant.panel.vcf.gz")
DV_VCF_PANEL_TBI = os.path.join(VCF_IN_DIR,      f"{SAMPLE}.deepvariant.panel.vcf.gz.tbi")
DV_VEP_TSV       = os.path.join(VEP_DIR,         f"{SAMPLE}.deepvariant.vep.tsv")
DV_VEP_STATS     = os.path.join(VEP_DIR,         f"{SAMPLE}.deepvariant.vep_stats.html")


CONCORDANCE_TSV  = os.path.join(TABLES_DIR, f"{SAMPLE}.caller_concordance.tsv")
CONCORDANCE_SUMM = os.path.join(TABLES_DIR, f"{SAMPLE}.concordance_summary.txt")

VEP_TSV          = os.path.join(VEP_DIR, f"{SAMPLE}.clair3.vep.tsv")
VEP_STATS        = os.path.join(VEP_DIR, f"{SAMPLE}.clair3.vep_stats.html")

VERSIONS_TXT     = os.path.join(META_DIR, "software_versions.txt")

# Report sub-directories
REPORT_DIR       = os.path.join(SAMPLE_ROOT, "report")
REPORT_DATA_DIR  = os.path.join(REPORT_DIR, "data")
REPORT_PLOT_DIR  = os.path.join(REPORT_DIR, "plots")

# Coverage report outputs
COV_BY_GENE      = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.coverage_by_gene.tsv")
COV_FINAL        = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.coverage.final.tsv")
COV_SUMMARY      = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.coverage_summary.tsv")
COV_BELOW10      = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.genes_below_10x.tsv")
COV_BELOW20      = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.genes_below_20x.tsv")
COV_PLOT         = os.path.join(REPORT_PLOT_DIR,  f"{SAMPLE}.coverage.png")

# QC / run summary report outputs
QC_THRESHOLDS_TSV  = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.qc_thresholds.tsv")
RUN_SUMMARY_TSV    = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.run_summary.tsv")
RUN_SUMMARY_DISP   = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.run_summary_display.tsv")
SHARED_TARGETS_TSV   = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.shared_targets.tsv")
QC_THRESHOLDS_PNG  = os.path.join(REPORT_PLOT_DIR,  f"{SAMPLE}.qc_thresholds.png")
REPORT_HTML        = os.path.join(REPORT_DIR, "html", f"{SAMPLE}.report.html")

# Logo — embedded in report header if file exists at this path
LOGO_PATH = os.path.join(REPORT_DIR, "logo",
                config.get("report", {}).get("logo_filename", "SignalScope_logo.png"))

# SNV report sub-directory and output files
REPORT_SNV_DIR      = os.path.join(REPORT_DIR, "snv")
SNV_SUMMARY         = os.path.join(REPORT_SNV_DIR, f"{SAMPLE}.snv.summary.tsv")
SNV_CLINVAR_SUMM    = os.path.join(REPORT_SNV_DIR, f"{SAMPLE}.snv.clinvar_summary.normalized.tsv")
SNV_SHORTLIST_CTX   = os.path.join(REPORT_SNV_DIR,
                          f"{SAMPLE}.snv.shortlist.with_gene_context.header.tsv")

# IGV screenshots — manifest TSV (base64 PNGs), embedded in the HTML report
IGV_DIR          = os.path.join(REPORT_DIR, "igv")
IGV_MANIFEST     = os.path.join(IGV_DIR, f"{SAMPLE}.igv_manifest.tsv")

# Sequencing performance (whole vs targeted) and tool versions
TARGETED_BAM       = os.path.join(QC_DIR, "targeted", f"{SAMPLE}.targeted.bam")
TARGETED_BAI       = os.path.join(QC_DIR, "targeted", f"{SAMPLE}.targeted.bam.bai")
TARGETED_NANOSTATS = os.path.join(QC_DIR, "targeted", "nanoplot", f"{SAMPLE}_NanoStats.txt")
SEQ_PERF_TSV       = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.sequencing_performance.tsv")
TOOL_VERSIONS_TSV  = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.tool_versions.tsv")

# AVI annotation sentinel (written by annotate_avi; consumed by report rules)
AVI_ENABLED = config.get("avi", {}).get("enabled", False) and MODE == "wgs"
AVI_SENTINEL = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.avi_annotated.done") if AVI_ENABLED else []

# Methylation module (modkit pileup + clinical screening)
_meth_cfg = config.get("methylation", {})
METH_ENABLED = _meth_cfg.get("enabled", False) and MODE == "wgs"
METH_DIR = os.path.join(SAMPLE_ROOT, "methylation")
METH_BEDMETHYL = os.path.join(METH_DIR, f"{SAMPLE}.cpg.combined.bed.gz")
METH_BEDMETHYL_TBI = METH_BEDMETHYL + ".tbi"
METH_CLINICAL_TSV = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.clinical_methylation.tsv")
METH_SENTINEL = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.meth_annotated.done") if METH_ENABLED else []

# WhatsHap phasing module
_phase_cfg = config.get("phasing", {})
PHASE_ENABLED = _phase_cfg.get("enabled", False) and MODE == "wgs"
PHASE_DIR = os.path.join(SAMPLE_ROOT, "phasing")
PHASED_VCF = os.path.join(PHASE_DIR, f"{SAMPLE}.phased.vcf.gz")
PHASE_STATS = os.path.join(PHASE_DIR, f"{SAMPLE}.phase_stats.tsv")
COMPOUND_HETS_TSV = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.compound_hets.tsv")
PHASE_SENTINEL = os.path.join(REPORT_DATA_DIR, f"{SAMPLE}.phase_annotated.done") if PHASE_ENABLED else []


# =============================================================================
# HELPER: create_sample_dirs
#
# Ensures the full directory skeleton exists before any rule writes output.
# Called via a shell checkpoint at the start of the pipeline.
# All rules also use mkdir -p in their own shell blocks for safety.
# =============================================================================

def get_all_sample_dirs(sample_root):
    """Return the full list of directories that should exist for a sample."""
    return [
        os.path.join(sample_root, "raw", "pod5"),
        os.path.join(sample_root, "raw", "fastq"),
        os.path.join(sample_root, "alignments"),
        os.path.join(sample_root, "qc", "nanoplot"),
        os.path.join(sample_root, "coverage"),
        os.path.join(sample_root, "variants", "small_variants", "clair3"),
        os.path.join(sample_root, "variants", "small_variants", "deepvariant"),
        os.path.join(sample_root, "variants", "structural_variants", "sniffles"),
        os.path.join(sample_root, "annotation", "vcf_inputs"),
        os.path.join(sample_root, "annotation", "vep"),
        os.path.join(sample_root, "annotation", "clinvar"),
        os.path.join(sample_root, "annotation", "tables"),
        os.path.join(sample_root, "logs"),
        os.path.join(sample_root, "tmp"),
        os.path.join(sample_root, "failed_attempts"),
        os.path.join(sample_root, "run_metadata"),
        os.path.join(sample_root, "report", "data"),
        os.path.join(sample_root, "report", "plots"),
        os.path.join(sample_root, "report", "scripts"),
        os.path.join(sample_root, "report", "html"),
        os.path.join(sample_root, "report", "logs"),
        os.path.join(sample_root, "report", "notes"),
        os.path.join(sample_root, "report", "snv"),
        os.path.join(sample_root, "methylation"),
        os.path.join(sample_root, "phasing"),
        os.path.join(sample_root, "qc", "targeted"),
        os.path.join(sample_root, "qc", "targeted", "nanoplot"),
    ]


# =============================================================================
# RULE ALL — Active pipeline targets
#
# All layers are active: core processing, variant analysis, and clinical interpretation/reporting.
#
# Start-mode logic:
#   pod5  → FASTQ_GZ is a pipeline output (from dorado_basecall)
#   fastq → FASTQ_GZ is a symlink/copy created by the handoff rule
#   bam   → BAM_PRI is handed off directly; FASTQ_GZ is not required
#
# Note on ancient() inputs in generate_report_html:
#   SNV outputs (snv_summary, snv_clinvar, snv_shortlist),
#   seq_perf_tsv, and tool_versions_tsv are declared with ancient().
#   This means if those files already exist, running --forcerun generate_report_html
#   will NOT re-trigger snv_process, build_sequencing_performance, or
#   record_tool_versions. This allows safe report previewing without
#   rerunning expensive variant callers or targeted BAM analysis.
# =============================================================================

# Pre-compute start-mode-dependent targets so rule all contains only plain names
# (Snakemake's parser rejects inline ternary expressions inside input: blocks)
# --- Build target list based on mode ------------------------------------------
_ALL_TARGETS = [
    BAM_PRI,
    BAM_PRI_BAI,
    NANOPLOT_HTML,
    FLAGSTAT,
    PANEL_READS,
    DEPTH_SUMM,
    MOSDEPTH_SUMMARY,
    MOSDEPTH_GLOBAL,
    VERSIONS_TXT,
    CLAIR3_VCF_RAW,
    CLAIR3_VCF_TBI,
    DV_VCF_GZ,
    DV_VCF_TBI,
    DV_VEP_TSV,
    DV_VEP_STATS,
    SNIFFLES_VCF,
    SNIFFLES_TBI,
    CONCORDANCE_TSV,
    CONCORDANCE_SUMM,
    VEP_TSV,
    VEP_STATS,
    COV_PLOT,
    RUN_SUMMARY_TSV,
    RUN_SUMMARY_DISP,
    REPORT_HTML,
    IGV_MANIFEST,
    SNV_SUMMARY,
    SNV_CLINVAR_SUMM,
    SNV_SHORTLIST_CTX,
    SEQ_PERF_TSV,
    TOOL_VERSIONS_TSV,
]

# Targets only needed when starting from pod5 or fastq
if START_FROM not in ("bam", "vcf"):
    _ALL_TARGETS.extend([FASTQ_GZ, BAM_ALL, BAM_ALL_BAI])

# Panel-only targets (skipped in WGS mode)
if MODE == "panel":
    _ALL_TARGETS += [
        MOSDEPTH_REGION,
        MOSDEPTH_BED_GZ,
        MOSDEPTH_BED_CSI,
        CLAIR3_VCF_PANEL,
        DV_VCF_PANEL,
        DV_VCF_PANEL_TBI,
        SNIFFLES_VCF_PANEL,
        SNIFFLES_TBI_PANEL,
        COV_BY_GENE,
        COV_FINAL,
        QC_THRESHOLDS_TSV,
        QC_THRESHOLDS_PNG,
        SHARED_TARGETS_TSV,
    ]

# WGS-only optional modules
if AVI_ENABLED:
    _ALL_TARGETS.append(AVI_SENTINEL)
if METH_ENABLED:
    _ALL_TARGETS.extend([METH_BEDMETHYL, METH_CLINICAL_TSV, METH_SENTINEL])
if PHASE_ENABLED:
    _ALL_TARGETS.extend([PHASED_VCF, PHASE_STATS, COMPOUND_HETS_TSV, PHASE_SENTINEL])

# Layer 3, panel mode: decision tree, SV prioritisation and report
include: "rules/interpretation_paths.smk"
_ALL_TARGETS.append(MODEL_CHECK)
if TREE_ENABLED:
    _ALL_TARGETS = [t for t in _ALL_TARGETS if t not in OLD_SHORTLIST_TARGETS] + \
                   [t for t in TREE_TARGETS if t not in _ALL_TARGETS]

rule all:
    input: _ALL_TARGETS



# =============================================================================
# LAYER 1 — CORE PIPELINE
# =============================================================================

# -----------------------------------------------------------------------------
# Rule: create_sample_structure
#
# Creates the full directory skeleton for this sample before anything else runs.
# This is a lightweight rule — it touches a sentinel file so Snakemake can
# track it as a dependency.
# -----------------------------------------------------------------------------
rule create_sample_structure:
    output:
        sentinel = os.path.join(META_DIR, ".structure_created")
    run:
        import shutil as _shutil
        for d in get_all_sample_dirs(SAMPLE_ROOT):
            os.makedirs(d, exist_ok=True)
        # Copy the SignalScope logo into the report so the HTML report shows branding.
        _logo_name = config.get("report", {}).get("logo_filename", "SignalScope_logo.png")
        _logo_src  = os.path.join(os.path.dirname(workflow.snakefile), "assets", _logo_name)
        _logo_dst  = os.path.join(SAMPLE_ROOT, "report", "logo", _logo_name)
        if os.path.isfile(_logo_src):
            os.makedirs(os.path.dirname(_logo_dst), exist_ok=True)
            _shutil.copyfile(_logo_src, _logo_dst)
        # Write sentinel
        with open(output.sentinel, "w") as f:
            f.write(f"Sample directory structure created for: {SAMPLE}\n")


# =============================================================================
# INPUT SOURCE RULES — mutually exclusive by START_FROM
#
# Each output path has exactly ONE rule that can produce it.
# Only the block matching START_FROM is parsed into the DAG.
# The other two blocks are never evaluated by Snakemake.
#
#   START_FROM=pod5   → dorado_basecall  produces raw/fastq/{sample}.fastq.gz
#   START_FROM=fastq  → fastq_handoff    produces raw/fastq/{sample}.fastq.gz
#   START_FROM=bam    → bam_handoff      produces alignments/{sample}.primary.bam
#
# filter_primary (pod5/fastq only) and align (pod5/fastq only) are also
# gated here so no rule tries to produce primary.bam when start_from=bam.
# =============================================================================

if START_FROM == "pod5":

    # -------------------------------------------------------------------------
    # Rule: dorado_basecall  [active only when start_from=pod5]
    #
    # Runs Dorado SUP on all POD5 files for this barcode.
    # Pipes stdout through gzip — no intermediate uncompressed file written.
    # Chemistry confirmed: dna_r10.4.1_e8.2_400bps_sup@v4.2.0
    # -------------------------------------------------------------------------
    rule dorado_basecall:
        input:
            sentinel = os.path.join(META_DIR, ".structure_created"),
            pod5_dir = os.path.join(config["pod5_base"], "{sample}")
        output:
            fastq_gz = os.path.join(
                config["samples_root"], "{sample}", "raw", "fastq", "{sample}.fastq.gz"
            )
        params:
            dorado      = DORADO,
            model       = config["dorado"]["model"],
            extra_flags = config["dorado"]["extra_flags"]
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "dorado_basecall.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "dorado_basecall.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] dorado_basecall started" > {log.out}
            echo "  sample   : {wildcards.sample}" >> {log.out}
            echo "  pod5_dir : {input.pod5_dir}"   >> {log.out}
            echo "  model    : {params.model}"     >> {log.out}

            {params.dorado} basecaller \
                {params.model} \
                {input.pod5_dir} \
                {params.extra_flags} \
                2>> {log.err} \
            | gzip -c > {output.fastq_gz}

            READ_COUNT=$(zcat {output.fastq_gz} | grep -c '^@' || true)
            echo "  reads written : $READ_COUNT"                           >> {log.out}
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] dorado_basecall complete" >> {log.out}
            """

    # -------------------------------------------------------------------------
    # Rule: align  [pod5 variant — input is dorado output]
    # -------------------------------------------------------------------------
    rule align:
        input:
            fastq_gz = os.path.join(
                config["samples_root"], "{sample}", "raw", "fastq", "{sample}.fastq.gz"
            ),
            ref = ALIGN_REF
        output:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.bam"
            )
        params:
            minimap2    = MINIMAP2,
            samtools    = SAMTOOLS,
            test_mode   = TEST_MODE,
            test_region = config.get("test_region", "")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "align.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "align.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bam})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] align started"  > {log.out}
            echo "  ref       : {input.ref}"       >> {log.out}
            echo "  test_mode : {params.test_mode}" >> {log.out}

            if [ "{params.test_mode}" = "True" ]; then
                {params.minimap2} -t {threads} -ax map-ont {input.ref} {input.fastq_gz} \
                    2>> {log.err} \
                | {params.samtools} sort -@ {threads} \
                    2>> {log.err} \
                | {params.samtools} view -@ {threads} -b "{params.test_region}" \
                    -o {output.bam} 2>> {log.err}
            else
                {params.minimap2} -t {threads} -ax map-ont {input.ref} {input.fastq_gz} \
                    2>> {log.err} \
                | {params.samtools} sort -@ {threads} -o {output.bam} 2>> {log.err}
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] align complete" >> {log.out}
            """

    # -------------------------------------------------------------------------
    # Rule: filter_primary  [pod5 variant]
    # -------------------------------------------------------------------------
    rule filter_primary:
        input:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.bam"
            )
        output:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.primary.bam"
            )
        params:
            samtools = SAMTOOLS
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] filter_primary started" > {log.out}
            echo "  no flag filter: secondary and supplementary alignments kept (split reads are SV evidence)"  >> {log.out}
            {params.samtools} view -@ {threads} -b {input.bam} \
                2>> {log.err} \
            | {params.samtools} sort -@ {threads} -o {output.bam} \
                2>> {log.err}
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] filter_primary complete" >> {log.out}
            """


elif START_FROM == "fastq":

    # -------------------------------------------------------------------------
    # Rule: fastq_handoff  [active only when start_from=fastq]
    #
    # Produces the canonical raw/fastq/{sample}.fastq.gz from an existing
    # FASTQ file without rerunning Dorado.
    #
    # Handles two cases automatically, detected from the source file extension:
    #
    #   source ends with .fastq.gz
    #     -> create a symlink (no data copied, instant)
    #
    #   source ends with .fastq  (e.g. barcode04.sup.Q10.fastq)
    #     -> gzip it into the canonical location
    #     -> original file is NOT modified or moved
    #
    # After this rule, every downstream rule sees only the canonical path.
    # -------------------------------------------------------------------------
    rule fastq_handoff:
        input:
            sentinel    = os.path.join(META_DIR, ".structure_created"),
            existing_fq = config["existing_fastq"]
        output:
            fastq_gz = os.path.join(
                config["samples_root"], "{sample}", "raw", "fastq", "{sample}.fastq.gz"
            )
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "fastq_handoff.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "fastq_handoff.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.fastq_gz})

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastq_handoff started" > {log.out}
            echo "  source : {input.existing_fq}" >> {log.out}
            echo "  target : {output.fastq_gz}"   >> {log.out}

            SOURCE="{input.existing_fq}"

            if [[ "$SOURCE" == *.fastq.gz ]]; then
                # Already gzipped — symlink directly, no data movement
                echo "  action : symlink (.fastq.gz source)" >> {log.out}
                ln -sf "$SOURCE" {output.fastq_gz}

            elif [[ "$SOURCE" == *.fastq ]]; then
                # Plain FASTQ — compress into canonical location
                # Original file is left untouched
                echo "  action : gzip (.fastq source)" >> {log.out}
                gzip -c "$SOURCE" > {output.fastq_gz} 2>> {log.err}
                READ_COUNT=$(zcat {output.fastq_gz} | grep -c "^@" || true)
                echo "  reads  : $READ_COUNT"         >> {log.out}

            else
                echo "[ERROR] existing_fastq must end in .fastq or .fastq.gz" >&2
                echo "  received: $SOURCE" >&2
                exit 1
            fi

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] fastq_handoff complete" >> {log.out}
            """

    # -------------------------------------------------------------------------
    # Rule: align  [fastq variant — input is the handed-off FASTQ symlink]
    # -------------------------------------------------------------------------
    rule align:
        input:
            fastq_gz = os.path.join(
                config["samples_root"], "{sample}", "raw", "fastq", "{sample}.fastq.gz"
            ),
            ref = ALIGN_REF
        output:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.bam"
            )
        params:
            minimap2    = MINIMAP2,
            samtools    = SAMTOOLS,
            test_mode   = TEST_MODE,
            test_region = config.get("test_region", "")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "align.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "align.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bam})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] align started"  > {log.out}
            echo "  ref       : {input.ref}"        >> {log.out}
            echo "  test_mode : {params.test_mode}" >> {log.out}

            if [ "{params.test_mode}" = "True" ]; then
                {params.minimap2} -t {threads} -ax map-ont {input.ref} {input.fastq_gz} \
                    2>> {log.err} \
                | {params.samtools} sort -@ {threads} \
                    2>> {log.err} \
                | {params.samtools} view -@ {threads} -b "{params.test_region}" \
                    -o {output.bam} 2>> {log.err}
            else
                {params.minimap2} -t {threads} -ax map-ont {input.ref} {input.fastq_gz} \
                    2>> {log.err} \
                | {params.samtools} sort -@ {threads} -o {output.bam} 2>> {log.err}
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] align complete" >> {log.out}
            """

    # -------------------------------------------------------------------------
    # Rule: filter_primary  [fastq variant]
    # -------------------------------------------------------------------------
    rule filter_primary:
        input:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.bam"
            )
        output:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.primary.bam"
            )
        params:
            samtools = SAMTOOLS
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] filter_primary started" > {log.out}
            echo "  no flag filter: secondary and supplementary alignments kept (split reads are SV evidence)"  >> {log.out}
            {params.samtools} view -@ {threads} -b {input.bam} \
                2>> {log.err} \
            | {params.samtools} sort -@ {threads} -o {output.bam} \
                2>> {log.err}
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] filter_primary complete" >> {log.out}
            """


elif START_FROM == "bam":

    # -------------------------------------------------------------------------
    # Rule: bam_handoff  [active only when start_from=bam]
    #
    # Symlinks an existing primary-filtered BAM into the canonical location.
    # The BAM must already be filtered (primary only) and indexed.
    # align and filter_primary are NOT defined in this mode.
    # -------------------------------------------------------------------------
    rule bam_handoff:
        input:
            sentinel     = os.path.join(META_DIR, ".structure_created"),
            existing_bam = config["existing_bam"]
        output:
            bam = os.path.join(
                config["samples_root"], "{sample}", "alignments", "{sample}.bam"
            )
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "bam_handoff.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "bam_handoff.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bam})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] bam_handoff started" > {log.out}
            echo "  source (all-reads): {input.existing_bam}" >> {log.out}
            ln -sf {input.existing_bam} {output.bam}
            if [ -f "{input.existing_bam}.bai" ]; then
                ln -sf {input.existing_bam}.bai {output.bam}.bai
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] bam_handoff complete" >> {log.out}
            """

    rule filter_primary:
        input:
            bam = os.path.join(config["samples_root"], "{sample}", "alignments", "{sample}.bam")
        output:
            bam = os.path.join(config["samples_root"], "{sample}", "alignments", "{sample}.primary.bam")
        params:
            samtools = SAMTOOLS
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            echo "filter_primary started" > {log.out}
            {params.samtools} view -@ {threads} -b {input.bam} 2>> {log.err} | {params.samtools} sort -@ {threads} -o {output.bam} 2>> {log.err}
            echo "filter_primary complete" >> {log.out}
            """





elif START_FROM == "vcf":
    rule bam_handoff:
        input:
            sentinel = os.path.join(META_DIR, ".structure_created"),
            existing_bam = config["existing_bam"]
        output:
            bam = os.path.join(config["samples_root"], "{sample}", "alignments", "{sample}.bam")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "bam_handoff.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "bam_handoff.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bam})
            ln -sf {input.existing_bam} {output.bam}
            [ -f "{input.existing_bam}.bai" ] && ln -sf {input.existing_bam}.bai {output.bam}.bai || true
            """

    rule filter_primary:
        input:
            bam = os.path.join(config["samples_root"], "{sample}", "alignments", "{sample}.bam")
        output:
            bam = os.path.join(config["samples_root"], "{sample}", "alignments", "{sample}.primary.bam")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "filter_primary.stderr.log")
        shell:
            """
            ln -sf $(realpath {input.bam}) {output.bam}
            """

# =============================================================================
# SHARED RULES — active for all start modes that produce a BAM
# =============================================================================

# -----------------------------------------------------------------------------
# Rule: index_all_bam
#
# Indexes the sorted all-reads BAM.
# Only relevant for pod5 and fastq modes (where an all-reads BAM exists).
# In bam mode this rule is never triggered (no all-reads BAM is produced).
# -----------------------------------------------------------------------------
rule index_all_bam:
    input:
        bam = os.path.join(
            config["samples_root"], "{sample}", "alignments", "{sample}.bam"
        )
    output:
        bai = os.path.join(
            config["samples_root"], "{sample}", "alignments", "{sample}.bam.bai"
        )
    params:
        samtools = SAMTOOLS
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "index_all_bam.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "index_all_bam.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing all-reads BAM" > {log.out}
        {params.samtools} index -@ {threads} {input.bam} 2>> {log.err}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: index_primary_bam
#
# Indexes the primary-only BAM.
# Required by Clair3, mosdepth, and all downstream tools.
# -----------------------------------------------------------------------------
rule index_primary_bam:
    input:
        bam = os.path.join(
            config["samples_root"], "{sample}", "alignments", "{sample}.primary.bam"
        )
    output:
        bai = os.path.join(
            config["samples_root"], "{sample}", "alignments", "{sample}.primary.bam.bai"
        )
    params:
        samtools = SAMTOOLS
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "index_primary_bam.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "index_primary_bam.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Indexing primary BAM" > {log.out}
        {params.samtools} index -@ {threads} {input.bam} 2>> {log.err}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: nanoplot_qc
#
# Alignment QC visualisations using NanoPlot run on the primary BAM.
# Using --bam produces alignment-aware stats (identity, coverage etc.)
# and works correctly regardless of start_from mode.
# -----------------------------------------------------------------------------
rule nanoplot_qc:
    input:
        bam = BAM_PRI,
        bai = BAM_PRI_BAI
    output:
        html    = os.path.join(
            config["samples_root"], "{sample}", "qc", "nanoplot", "{sample}_NanoPlot-report.html"
        ),
        summary = os.path.join(
            config["samples_root"], "{sample}", "qc", "nanoplot", "{sample}_NanoStats.txt"
        )
    params:
        nanoplot = NANOPLOT,
        outdir   = os.path.join(config["samples_root"], "{sample}", "qc", "nanoplot"),
        prefix   = "{sample}_"
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "nanoplot.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "nanoplot.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] nanoplot_qc started" > {log.out}
        echo "  BAM    : {input.bam}"     >> {log.out}
        echo "  outdir : {params.outdir}" >> {log.out}

        {params.nanoplot} \
            --bam {input.bam} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --threads {threads} \
            > {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] nanoplot_qc complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: basic_qc
#
# Three alignment QC checks from one rule (all share the same BAM):
#   1. flagstat    — overall alignment statistics
#   2. panel_reads — reads overlapping the CKD panel BED
#   3. depth_summ  — mean and max depth across panel positions
# -----------------------------------------------------------------------------
if MODE == "wgs":
    rule basic_qc:
        """WGS mode: flagstat only (no panel BED for read count / depth)."""
        input:
            bam = BAM_PRI,
            bai = BAM_PRI_BAI
        output:
            flagstat    = os.path.join(config["samples_root"], "{sample}", "qc", "{sample}.flagstat.txt"),
            panel_reads = os.path.join(config["samples_root"], "{sample}", "qc", "{sample}.panel_read_count.txt"),
            depth_summ  = os.path.join(config["samples_root"], "{sample}", "qc", "{sample}.panel_depth_summary.txt")
        params:
            samtools = SAMTOOLS
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "basic_qc.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "basic_qc.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.flagstat})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] basic_qc (WGS mode) started" > {log.out}
            echo "  BAM : {input.bam}" >> {log.out}

            echo "--- flagstat ---" >> {log.out}
            {params.samtools} flagstat -@ {threads} {input.bam} \
                > {output.flagstat} 2>> {log.err}
            cat {output.flagstat} >> {log.out}

            echo "--- panel read count (WGS: n/a) ---" >> {log.out}
            echo "n/a (WGS mode)" > {output.panel_reads}

            echo "--- panel depth summary (WGS: n/a) ---" >> {log.out}
            echo "n/a (WGS mode)" > {output.depth_summ}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] basic_qc complete" >> {log.out}
            """
else:
    rule basic_qc:
        """Panel mode: flagstat + panel read count + panel depth summary."""
        input:
            bam = BAM_PRI,
            bai = BAM_PRI_BAI,
            bed = ACTIVE_BED if ACTIVE_BED else []
        output:
            flagstat    = os.path.join(config["samples_root"], "{sample}", "qc", "{sample}.flagstat.txt"),
            panel_reads = os.path.join(config["samples_root"], "{sample}", "qc", "{sample}.panel_read_count.txt"),
            depth_summ  = os.path.join(config["samples_root"], "{sample}", "qc", "{sample}.panel_depth_summary.txt")
        params:
            samtools = SAMTOOLS
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "basic_qc.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "basic_qc.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail

            mkdir -p $(dirname {output.flagstat})

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] basic_qc started" > {log.out}
            echo "  BAM : {input.bam}" >> {log.out}
            echo "  BED : {input.bed}" >> {log.out}

            echo "--- flagstat ---" >> {log.out}
            {params.samtools} flagstat -@ {threads} {input.bam} \
                > {output.flagstat} 2>> {log.err}
            cat {output.flagstat} >> {log.out}

            echo "--- panel read count ---" >> {log.out}
            {params.samtools} view -@ {threads} -L {input.bed} -c {input.bam} \
                > {output.panel_reads} 2>> {log.err}
            echo "  reads on panel: $(cat {output.panel_reads})" >> {log.out}

            echo "--- panel depth summary ---" >> {log.out}
            {params.samtools} depth -b {input.bed} {input.bam} \
                2>> {log.err} \
            | awk '{{sum+=$3; n++; if($3>max) max=$3}} \
                   END{{print "mean=" sum/n, "max=" max, "positions=" n}}' \
                > {output.depth_summ}
            cat {output.depth_summ} >> {log.out}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] basic_qc complete" >> {log.out}
            """


# -----------------------------------------------------------------------------
# Rule: mosdepth_coverage
# Panel mode: per-region coverage with --by BED.
# WGS mode: genome-wide per-chromosome coverage (no --by).
# -----------------------------------------------------------------------------
if MODE == "wgs":
    rule mosdepth_coverage:
        input:
            bam = BAM_PRI,
            bai = BAM_PRI_BAI
        output:
            summary = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth.summary.txt"),
            global_ = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth.global.dist.txt"),
        params:
            mosdepth = MOSDEPTH,
            prefix   = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "mosdepth.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "mosdepth.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {params.prefix})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth_coverage (WGS) started" > {log.out}
            echo "  BAM    : {input.bam}"     >> {log.out}

            {params.mosdepth} \
                --no-per-base \
                --threads {threads} \
                {params.prefix} \
                {input.bam} \
                >> {log.out} 2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth_coverage (WGS) complete" >> {log.out}
            """
else:
    rule mosdepth_coverage:
        input:
            bam = BAM_PRI,
            bai = BAM_PRI_BAI,
            bed = config.get("panel_bed", "") or [] if MODE == "panel" else []
        output:
            summary = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth.summary.txt"),
            global_ = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth.global.dist.txt"),
            region  = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth.region.dist.txt"),
            bed_gz  = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.regions.bed.gz"),
            bed_csi = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.regions.bed.gz.csi")
        params:
            mosdepth = MOSDEPTH,
            prefix   = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "mosdepth.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "mosdepth.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {params.prefix})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth_coverage started" > {log.out}
            echo "  BAM    : {input.bam}"     >> {log.out}

            {params.mosdepth} \
                --by {input.bed} \
                --threads {threads} \
                {params.prefix} \
                {input.bam} \
                > {log.out} 2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth_coverage complete" >> {log.out}
            """


# -----------------------------------------------------------------------------
# Rule: mosdepth_coverage_report
# Panel mode: second mosdepth run with REPORT BED (separate COL4A3/COL4A4).
# WGS mode: skipped (outputs are placeholders pointing to the main mosdepth).
# -----------------------------------------------------------------------------
if MODE == "wgs":
    rule mosdepth_coverage_report:
        input:
            summary_main = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth.summary.txt")
        output:
            summary = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.mosdepth.summary.txt"),
            global_ = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.mosdepth.global.dist.txt"),
            region  = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.mosdepth.region.dist.txt"),
            bed_gz  = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.regions.bed.gz"),
            bed_csi = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.regions.bed.gz.csi")
        shell:
            """
            cp {input.summary_main} {output.summary}
            touch {output.global_} {output.region} {output.bed_gz} {output.bed_csi}
            """
else:
    rule mosdepth_coverage_report:
        input:
            bam = BAM_PRI,
            bai = BAM_PRI_BAI,
            bed = config.get("report_bed","") if config.get("report_bed","") else []
        output:
            summary = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.mosdepth.summary.txt"),
            global_ = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.mosdepth.global.dist.txt"),
            region  = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.mosdepth.region.dist.txt"),
            bed_gz  = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.regions.bed.gz"),
            bed_csi = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.regions.bed.gz.csi"),
            per_base = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report.per-base.bed.gz")
        params:
            mosdepth = MOSDEPTH,
            prefix   = os.path.join(config["samples_root"], "{sample}", "coverage", "{sample}.mosdepth_report")
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs", "mosdepth_report.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs", "mosdepth_report.stderr.log")
        threads: config["threads"]
        shell:
            """
        set -euo pipefail

        mkdir -p $(dirname {params.prefix})

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth_coverage_report started" > {log.out}
        echo "  BAM    : {input.bam}"     >> {log.out}
        echo "  prefix : {params.prefix}" >> {log.out}
        echo "  BED    : {input.bed}"     >> {log.out}

        {params.mosdepth} \
            --by {input.bed} \
            --threads {threads} \
            {params.prefix} \
            {input.bam} \
            > {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth_coverage_report complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: record_versions
#
# Captures tool version strings for reproducibility auditing.
# Runs automatically with every pipeline execution.
# -----------------------------------------------------------------------------
rule record_versions:
    input:
        sentinel = os.path.join(META_DIR, ".structure_created")
    output:
        txt = os.path.join(
            config["samples_root"], "{sample}", "run_metadata", "software_versions.txt"
        )
    params:
        dorado   = DORADO,
        minimap2 = MINIMAP2,
        samtools = SAMTOOLS,
        bcftools = BCFTOOLS,
        clair3   = CLAIR3_RUNNER,
        mosdepth = MOSDEPTH,
        nanoplot = NANOPLOT,
        model    = CLAIR3_MODEL
    log:
        err = os.path.join(config["samples_root"], "{sample}", "logs", "record_versions.stderr.log")
    shell:
        """
        set -euo pipefail
        {{
            echo "# CKD pipeline software versions"
            echo "# Sample       : {wildcards.sample}"
            echo "# Date         : $(date '+%Y-%m-%d %H:%M:%S')"
            echo "# Clair3 model : {params.model}"
            echo ""
            echo "=== dorado ==="
            {params.dorado} --version 2>&1 || echo "unavailable"
            echo ""
            echo "=== minimap2 ==="
            {params.minimap2} --version 2>&1 || echo "unavailable"
            echo ""
            echo "=== samtools ==="
            {params.samtools} --version 2>&1 | head -2
            echo ""
            echo "=== bcftools ==="
            {params.bcftools} --version 2>&1 | head -2
            echo ""
            echo "=== clair3 ==="
            bash {params.clair3} --version 2>&1 || echo "unavailable"
            echo ""
            echo "=== mosdepth ==="
            {params.mosdepth} --version 2>&1 || echo "unavailable"
            echo ""
            echo "=== nanoplot ==="
            {params.nanoplot} --version 2>&1 || echo "unavailable"
            echo ""
            echo "=== snakemake ==="
            snakemake --version 2>&1 || echo "unavailable"
        }} > {output.txt} 2>> {log.err}
        """


# -----------------------------------------------------------------------------
# Rule: clair3_call
#
# Small variant calling (SNV + indel) on the primary BAM.
# All parameters confirmed from the validated barcode06 manual run.
# Model: r1041_e82_400bps_sup_v420 (chemistry-matched to R10.4.1)
# -----------------------------------------------------------------------------


rule clair3_call:
    input:
        bam = BAM_FOR_CALLING,
        bai = BAM_FOR_CALLING_BAI,
        model_check = MODEL_CHECK,
        ref = config["ref_fasta"],
        bed = ACTIVE_BED if ACTIVE_BED else []
    output:
        vcf = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "small_variants", "clair3", "merge_output.vcf.gz"
        )
    params:
        sif        = config["clair3"]["sif"],
        model_dir  = CLAIR3_MODELDIR,
        outdir     = os.path.join(
            config["samples_root"], "{sample}", "variants", "small_variants", "clair3"
        ),
        existing_vcf = config.get("existing_clair3_vcf", "")
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "clair3.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "clair3.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail

        mkdir -p {params.outdir}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] clair3_call started" > {log.out}
        echo "  sample : {wildcards.sample}"  >> {log.out}
        echo "  BAM    : {input.bam}"         >> {log.out}
        echo "  BED    : {input.bed}"         >> {log.out}
        # model chosen or checked by clair3_model_check (clair3.model: auto or a name)
        MODEL_NAME=$(awk '/^model: /{{print $2; exit}}' {input.model_check})
        [ -n "$MODEL_NAME" ] || {{ echo "[ERROR] no model in {input.model_check}" >> {log.err}; exit 1; }}
        echo "  model  : $MODEL_NAME" >> {log.out}

        # VCF mode: symlink existing VCF
        if [ -n "{params.existing_vcf}" ]; then
            mkdir -p $(dirname {output.vcf})
            ln -sf $(realpath {params.existing_vcf}) {output.vcf}
            [ -f "{params.existing_vcf}.tbi" ] && ln -sf $(realpath {params.existing_vcf}.tbi) {output.vcf}.tbi || tabix -p vcf {output.vcf}
            exit 0
        fi

        ABS_REF=$(realpath {input.ref})
        ABS_BAM=$(realpath {input.bam})
        ABS_OUT=$(realpath {params.outdir})
        ABS_MODEL=$(realpath {params.model_dir}/$MODEL_NAME)
        export APPTAINER_TMPDIR=$ABS_OUT/apptainer_tmp
        export APPTAINER_CACHEDIR=$ABS_OUT/apptainer_tmp
        mkdir -p $ABS_OUT/apptainer_tmp
        BED_BIND=""
        BED_FLAG=""
        if [ -n "{input.bed}" ]; then
            ABS_BED=$(realpath {input.bed})
            BED_BIND="--bind $(dirname $ABS_BED)"
            BED_FLAG="--bed_fn=$ABS_BED"
        fi
        apptainer exec \
            --nv \
            --no-home \
            --cleanenv \
            --home /tmp \
            --bind $(dirname $ABS_REF) \
            $BED_BIND \
            --bind $(dirname $ABS_BAM) \
            --bind $ABS_OUT \
            --bind $ABS_MODEL \
            {params.sif} \
            /opt/bin/run_clair3.sh \
                --bam_fn=$ABS_BAM \
                --ref_fn=$ABS_REF \
                --threads={threads} \
                --platform=ont \
                --model_path=$ABS_MODEL \
                --output=$ABS_OUT \
                $BED_FLAG \
                --use_gpu \
                --qual=2 \
                --snp_min_af=0.08 \
                --indel_min_af=0.15 \
                --min_mq=5 \
                --min_coverage=2 \
            >> {log.out} 2>> {log.err}

        VARIANT_COUNT=$(zcat {output.vcf} | grep -vc '^#' || echo "unknown")
        echo "  raw variants called : $VARIANT_COUNT" >> {log.out}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] clair3_call complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: panel_filter_clair3
#
# Restricts Clair3 VCF to panel BED regions.
# Output goes to annotation/vcf_inputs/ — the authoritative input for VEP.
# Always uses full panel_bed (not ACTIVE_BED) so test-mode runs are filtered
# to the real panel, not the test region.
# -----------------------------------------------------------------------------
# In WGS mode, panel_filter rules are not in the DAG (their outputs are
# excluded from rule all). If somehow triggered, they pass through.
rule panel_filter_clair3:
    input:
        vcf = CLAIR3_VCF_RAW,
        bed = config.get("panel_bed", "") or [] if MODE == "panel" else []
    output:
        vcf = os.path.join(
            config["samples_root"], "{sample}",
            "annotation", "vcf_inputs", "{sample}.clair3.panel.vcf.gz"
        ),
        tbi = os.path.join(
            config["samples_root"], "{sample}",
            "annotation", "vcf_inputs", "{sample}.clair3.panel.vcf.gz.tbi"
        )
    params:
        bcftools = BCFTOOLS
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "panel_filter_clair3.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "panel_filter_clair3.stderr.log")
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.vcf})

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] panel_filter_clair3 started" > {log.out}
        echo "  input : {input.vcf}" >> {log.out}
        echo "  BED   : {input.bed}" >> {log.out}

        if [ -n "{input.bed}" ]; then
            {params.bcftools} view -R {input.bed} {input.vcf} 2>> {log.err} | {params.bcftools} sort 2>> {log.err} | bgzip -c > {output.vcf} 2>> {log.err}
        else
            {params.bcftools} sort {input.vcf} 2>> {log.err} | bgzip -c > {output.vcf} 2>> {log.err}
        fi
        tabix -p vcf {output.vcf} 2>> {log.err}

        VARIANT_COUNT=$(zcat {output.vcf} | grep -vc '^#' || echo "unknown")
        echo "  variants after panel filter : $VARIANT_COUNT" >> {log.out}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] panel_filter_clair3 complete" >> {log.out}
        """



# =============================================================================
# Rule: cutesv_call
#
# Calls structural variants using CuteSV 2.1.3 (settings: config 'cutesv' block).
# Run genome-wide on the primary BAM, then filtered to panel BED.
# FORMAT: GT:DR:DV:PL:GQ — proper genotype and read support fields.
# Self-filters low-confidence calls with q5 flag.
# Replaces SVIM which lacked FORMAT columns (GT=N/A, depth=0).
# =============================================================================
# Rule: sniffles_call
#
# Calls structural variants using Sniffles2.
# Run genome-wide on the primary BAM, then panel-filtered.
# =============================================================================
rule sniffles_call:
    input:
        bam = BAM_FOR_CALLING,
        bai = BAM_FOR_CALLING_BAI,
        bed = config.get("panel_bed","") if config.get("panel_bed","") else []
    output:
        vcf_gz = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "structural_variants", "sniffles", "{sample}.sniffles.vcf.gz"
        ),
        tbi = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "structural_variants", "sniffles", "{sample}.sniffles.vcf.gz.tbi"
        ),
        panel = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "structural_variants", "sniffles", "{sample}.sniffles.panel.vcf.gz"
        ),
        panel_tbi = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "structural_variants", "sniffles", "{sample}.sniffles.panel.vcf.gz.tbi"
        )
    params:
        sniffles   = SNIFFLES,
        bcftools   = BCFTOOLS,
        existing_vcf = config.get("existing_sniffles_vcf", ""),
        minsvlen   = config["sniffles"]["minsvlen"],
        minsupport = config["sniffles"]["minsupport"],
        workdir    = os.path.join(
            config["samples_root"], "{sample}", "variants", "structural_variants", "sniffles"
        ),
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "sniffles.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "sniffles.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p {params.workdir}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] sniffles_call started" > {log.out}

        if [ -n "{params.existing_vcf}" ]; then
            mkdir -p $(dirname {output.vcf_gz})
            ln -sf $(realpath {params.existing_vcf}) {output.vcf_gz}
            [ -f "{params.existing_vcf}.tbi" ] && ln -sf $(realpath {params.existing_vcf}.tbi) {output.tbi} || tabix -p vcf {output.vcf_gz}
            cp {output.vcf_gz} {output.panel} && tabix -p vcf {output.panel}
            echo "  symlinked from {params.existing_vcf}" >> {log.out}
            exit 0
        fi
        echo "  sample : {wildcards.sample}" >> {log.out}
        echo "  BAM    : {input.bam}"        >> {log.out}

        {params.sniffles} \
            --input {input.bam} \
            --vcf {params.workdir}/{wildcards.sample}.sniffles.raw.vcf \
            --threads {threads} \
            --minsvlen {params.minsvlen} \
            --minsupport {params.minsupport} \
            --output-rnames \
            2>> {log.err}

        bgzip -c {params.workdir}/{wildcards.sample}.sniffles.raw.vcf > {output.vcf_gz} 2>> {log.err}
        tabix -p vcf {output.vcf_gz} 2>> {log.err}
        rm -f {params.workdir}/{wildcards.sample}.sniffles.raw.vcf

        if [ -n "{input.bed}" ]; then
            {params.bcftools} view --regions-file {input.bed} {output.vcf_gz} -O z -o {output.panel} 2>> {log.err}
        else
            cp {output.vcf_gz} {output.panel}
        fi
        tabix -p vcf {output.panel} 2>> {log.err}

        RAW=$(zcat {output.vcf_gz} | grep -vc '^#' || true)
        PANEL=$(zcat {output.panel} | grep -vc '^#' || true)
        echo "  SVs genome-wide : $RAW"   >> {log.out}
        echo "  SVs panel       : $PANEL" >> {log.out}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] sniffles_call complete" >> {log.out}
        """
# =============================================================================
rule cutesv_call:
    input:
        bam = BAM_FOR_CALLING,
        bai = BAM_FOR_CALLING_BAI,
        ref = config["ref_fasta"],
        bed = config.get("panel_bed","") if config.get("panel_bed","") else []
    output:
        vcf_gz    = os.path.join(config["samples_root"], "{sample}", "variants", "structural_variants", "cutesv", "{sample}.cutesv.vcf.gz"),
        tbi       = os.path.join(config["samples_root"], "{sample}", "variants", "structural_variants", "cutesv", "{sample}.cutesv.vcf.gz.tbi"),
        panel     = os.path.join(config["samples_root"], "{sample}", "variants", "structural_variants", "cutesv", "{sample}.cutesv.panel.vcf.gz"),
        panel_tbi = os.path.join(config["samples_root"], "{sample}", "variants", "structural_variants", "cutesv", "{sample}.cutesv.panel.vcf.gz.tbi")
    params:
        cutesv  = config["tools"]["cutesv"],
        bcftools = BCFTOOLS,
        existing_vcf = config.get("existing_cutesv_vcf", ""),
        workdir  = os.path.join(CUTESV_DIR, "tmp"),
        # the settings of the SV benchmark (cuteSV default minimum support 10; ONT
        # clustering options) with min size 50 bp and SVs up to 5 Mb (config cutesv: block)
        min_support = config.get("cutesv", {}).get("min_support", 10),
        min_size    = config.get("cutesv", {}).get("min_size", 50),
        max_size    = config.get("cutesv", {}).get("max_size", 5000000),
        bias_ins    = config.get("cutesv", {}).get("max_cluster_bias_INS", 100),
        ratio_ins   = config.get("cutesv", {}).get("diff_ratio_merging_INS", 0.3),
        bias_del    = config.get("cutesv", {}).get("max_cluster_bias_DEL", 100),
        ratio_del   = config.get("cutesv", {}).get("diff_ratio_merging_DEL", 0.3)
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "cutesv.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "cutesv.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p {params.workdir}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] cutesv_call started" > {log.out}
        echo "  sample : {wildcards.sample}" >> {log.out}
        echo "  BAM    : {input.bam}"        >> {log.out}

        {params.cutesv} \
            {input.bam} \
            {input.ref} \
            {params.workdir}/{wildcards.sample}.cutesv.raw.vcf \
            {params.workdir} \
            --threads {threads} \
            --sample  {wildcards.sample} \
            --genotype \
            --report_readid \
            --min_support {params.min_support} \
            --min_size    {params.min_size} \
            -L {params.max_size} \
            --max_cluster_bias_INS {params.bias_ins} --diff_ratio_merging_INS {params.ratio_ins} \
            --max_cluster_bias_DEL {params.bias_del} --diff_ratio_merging_DEL {params.ratio_del} \
            >> {log.out} 2>> {log.err}

        bgzip -c {params.workdir}/{wildcards.sample}.cutesv.raw.vcf > {output.vcf_gz} 2>> {log.err}
        tabix -p vcf {output.vcf_gz} 2>> {log.err}

        if [ -n "{input.bed}" ]; then
            {params.bcftools} view --regions-file {input.bed} --apply-filters PASS {output.vcf_gz} -O z -o {output.panel} 2>> {log.err}
        else
            {params.bcftools} view --apply-filters PASS {output.vcf_gz} -O z -o {output.panel} 2>> {log.err}
        fi
        tabix -p vcf {output.panel} 2>> {log.err}

        RAW=$(zcat {output.vcf_gz} | grep -vc '^#' || true)
        PANEL=$(zcat {output.panel} | grep -vc '^#' || true)
        echo "  SVs genome-wide : $RAW"   >> {log.out}
        echo "  SVs panel PASS  : $PANEL" >> {log.out}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] cutesv_call complete" >> {log.out}
        """


rule caller_concordance:
    input:
        clair3   = CLAIR3_VCF_PANEL,
        deepvariant = DV_VCF_PANEL   # panel-filtered — matches Clair3 scope
    output:
        table   = os.path.join(
            config["samples_root"], "{sample}", "annotation", "tables",
            "{sample}.caller_concordance.tsv"
        ),
        summary = os.path.join(
            config["samples_root"], "{sample}", "annotation", "tables",
            "{sample}.concordance_summary.txt"
        )
    params:
        bcftools = BCFTOOLS
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "concordance.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "concordance.stderr.log")
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.table})

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] caller_concordance started" > {log.out}

        TMPDIR=$(mktemp -d)

        {params.bcftools} isec \
            -p $TMPDIR \
            {input.clair3} \
            {input.deepvariant} \
            2>> {log.err}

        # bcftools isec layout:
        #   0000.vcf — Clair3 private
        #   0002.vcf — Clair3 records at shared sites

        CLAIR3_TOTAL=$(zcat {input.clair3}   | grep -vc '^#' || echo 0)
        DV_TOTAL=$(zcat {input.deepvariant} | grep -vc '^#' || echo 0)
        SHARED=$(grep -vc '^#' $TMPDIR/0002.vcf   || echo 0)
        CLAIR3_ONLY=$(grep -vc '^#' $TMPDIR/0000.vcf  || echo 0)
        DV_ONLY=$(grep -vc '^#' $TMPDIR/0001.vcf || echo 0)

        {{
            echo "# Caller concordance — {wildcards.sample}"
            echo "Clair3_total    $CLAIR3_TOTAL"
            echo "DeepVariant_total  $DV_TOTAL"
            echo "Shared          $SHARED"
            echo "Clair3_only     $CLAIR3_ONLY"
            echo "DeepVariant_only   $DV_ONLY"
        }} | tee {output.summary} >> {log.out}

        echo -e "CHROM\tPOS\tREF\tALT\tIN_CLAIR3\tIN_DEEPVARIANT\tCALLER_SUPPORT" \
            > {output.table}

        grep -v '^#' $TMPDIR/0000.vcf \
            | awk '{{print $1"\t"$2"\t"$4"\t"$5"\t1\t0\t1"}}' \
            >> {output.table} || true

        grep -v '^#' $TMPDIR/0002.vcf \
            | awk '{{print $1"\t"$2"\t"$4"\t"$5"\t1\t1\t2"}}' \
            >> {output.table} || true

        grep -v '^#' $TMPDIR/0001.vcf \
            | awk '{{print $1"\t"$2"\t"$4"\t"$5"\t0\t1\t1"}}' \
            >> {output.table} || true

        rm -rf $TMPDIR
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] caller_concordance complete" >> {log.out}
        """


# =============================================================================
# LAYER 3 — CLINICAL INTERPRETATION & REPORTING
# VEP annotation, ACMG/AMP classification, variant shortlisting, and HTML report.
# =============================================================================

# -----------------------------------------------------------------------------
# Rule: vep_annotate
#
# Annotates the panel-filtered Clair3 VCF using Ensembl VEP offline.
# Output is tab-separated, one row per transcript per variant.
# Transcript selection and consequence prioritisation happen downstream in snv_process.
# -----------------------------------------------------------------------------
rule vep_annotate:
    input:
        vcf = CLAIR3_VCF_PANEL,
        ref = config["ref_fasta"]
    output:
        tsv   = os.path.join(
            config["samples_root"], "{sample}", "annotation", "vep",
            "{sample}.clair3.vep.tsv"
        ),
        stats = os.path.join(
            config["samples_root"], "{sample}", "annotation", "vep",
            "{sample}.clair3.vep_stats.html"
        )
    params:
        vep      = VEP,
        cache    = config["vep"]["cache"],
        plugins  = config["vep"]["plugins"],
        assembly = config["vep"]["assembly"],
        flags    = config["vep"]["extra_flags"]
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "vep.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "vep.stderr.log")
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {output.tsv})

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] vep_annotate started" > {log.out}
        echo "  input  : {input.vcf}"    >> {log.out}
        echo "  cache  : {params.cache}" >> {log.out}

        {params.vep} \
            --input_file {input.vcf} \
            --output_file {output.tsv} \
            --stats_file {output.stats} \
            --fasta {input.ref} \
            --cache \
            --dir_cache {params.cache} \
            --dir_plugins {params.plugins} \
            --assembly {params.assembly} \
            {params.flags} \
            > {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] vep_annotate complete" >> {log.out}
        """



# =============================================================================
# Rule: vep_annotate_sniffles
# Runs VEP on Sniffles panel VCF to annotate 10-49bp indels for the
# indel table (two-caller concordance: Clair3 + Sniffles).
# =============================================================================
rule vep_annotate_sniffles:
    input:
        vcf = SNIFFLES_VCF_PANEL,
        ref = config["ref_fasta"]
    output:
        tsv   = os.path.join(
            config["samples_root"], "{sample}", "annotation", "vep",
            "{sample}.sniffles.vep.tsv"
        ),
        stats = os.path.join(
            config["samples_root"], "{sample}", "annotation", "vep",
            "{sample}.sniffles.vep_stats.html"
        )
    params:
        vep      = VEP,
        cache    = config["vep"]["cache"],
        plugins  = config["vep"]["plugins"],
        assembly = config["vep"]["assembly"],
        flags    = config["vep"]["extra_flags"]
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs",
                           "vep_sniffles.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs",
                           "vep_sniffles.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] vep_annotate_sniffles started" > {log.out}
        {params.vep} \
            --input_file {input.vcf} \
            --output_file {output.tsv} \
            --stats_file {output.stats} \
            --fasta {input.ref} \
            --cache \
            --dir_cache {params.cache} \
            --dir_plugins {params.plugins} \
            --assembly {params.assembly} \
            {params.flags} \
            > {log.out} 2>> {log.err}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] vep_annotate_sniffles complete" >> {log.out}
        """


# =============================================================================
# Rule: deepvariant_call
#
# Runs DeepVariant (ONT_R104 model) via Apptainer on the primary BAM.
# Restricted to the panel BED — reduces runtime significantly.
# Apptainer must be available on the execution host (deepvariant_sif in config.yaml).
#
# Output VCF is bgzipped + tabixed in a post-processing step because
# DeepVariant's --output_vcf writes a plain .vcf.gz that is already bgzipped
# but NOT tabixed; the tbi is created here to allow downstream bcftools view.
#
# Mount strategy: bind-mount the run directory and the reference directory
# separately so all absolute paths resolve correctly inside the container.
# =============================================================================
rule deepvariant_call:
    input:
        bam = BAM_FOR_CALLING,
        bai = BAM_FOR_CALLING_BAI,
        ref = config["ref_fasta"],
        bed = config.get("deepvariant_bed", "") or []
    output:
        vcf = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "small_variants", "deepvariant", "{sample}.deepvariant.vcf.gz"
        ),
        tbi = os.path.join(
            config["samples_root"], "{sample}",
            "variants", "small_variants", "deepvariant", "{sample}.deepvariant.vcf.gz.tbi"
        )
    params:
        outdir     = os.path.join(
            config["samples_root"], "{sample}", "variants", "small_variants", "deepvariant"
        ),
        dv_version = "1.10.0",
        sif        = config["deepvariant_sif"],
        samtools   = SAMTOOLS,
        existing_vcf = config.get("existing_dv_vcf", ""),
        # DeepVariant small model, off by default (deepvariant.small_model: false).
        # Benchmark (GIAB v5.0q, HG002): at 20x-40x it changes F1 by <= 0.003; at 10x it lowers
        # SNV recall to 0.49, so the panel default keeps it off.
        small_model = "" if (config.get("deepvariant") or {}).get("small_model", False) else "--disable_small_model"
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs", "deepvariant.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs", "deepvariant.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] deepvariant_call started" > {log.out}
        echo "  sample  : {wildcards.sample}" >> {log.out}
        echo "  BAM     : {input.bam}"        >> {log.out}
        echo "  BED     : {input.bed}"        >> {log.out}
        echo "  model   : ONT_R104 (DeepVariant {params.dv_version}) {params.small_model}" >> {log.out}
        echo "  runner  : Apptainer {params.sif}" >> {log.out}

        mkdir -p {params.outdir}/tmp
        mkdir -p {params.outdir}/dv_logs

        if [ -n "{params.existing_vcf}" ]; then
            mkdir -p $(dirname {output.vcf})
            ln -sf $(realpath {params.existing_vcf}) {output.vcf}
            [ -f "{params.existing_vcf}.tbi" ] && ln -sf $(realpath {params.existing_vcf}.tbi) {output.tbi} || tabix -p vcf {output.vcf}
            exit 0
        fi

        ABS_BAM=$(realpath {input.bam})
        ABS_OUT=$(realpath {params.outdir})
        ABS_REF=$(realpath {input.ref})

        BED_BIND=""
        REGIONS_FLAG=""
        if [ -n "{input.bed}" ]; then
            ABS_BED=$(realpath {input.bed})
            # the BED file itself is passed to --regions (run_deepvariant accepts
            # "paths to BED/BEDPE files"), so every interval is called.
            echo "  regions : $ABS_BED ($(grep -vc '^#' $ABS_BED) intervals)" >> {log.out}
            BED_BIND="--bind $(dirname $ABS_BED)"
            REGIONS_FLAG="--regions=$ABS_BED"
        else
            echo "  regions : whole genome (WGS mode)" >> {log.out}
        fi

        export APPTAINER_TMPDIR=$ABS_OUT/apptainer_tmp
        export APPTAINER_CACHEDIR=$ABS_OUT/apptainer_tmp
        mkdir -p $ABS_OUT/apptainer_tmp

        apptainer exec \
            --nv \
            --no-home \
            --cleanenv \
            --home /tmp \
            --bind $(dirname $ABS_REF) \
            $BED_BIND \
            --bind $(dirname $ABS_BAM) \
            --bind $ABS_OUT \
            {params.sif} \
            run_deepvariant \
                --model_type=ONT_R104 \
                --ref=$ABS_REF \
                --reads=$ABS_BAM \
                $REGIONS_FLAG \
                --output_vcf=$ABS_OUT/{wildcards.sample}.deepvariant.vcf.gz \
                --num_shards={threads} \
                --intermediate_results_dir=$ABS_OUT/tmp \
                --logging_dir=$ABS_OUT/dv_logs \
                {params.small_model} \

            >> {log.out} 2>> {log.err}

        # Ensure tabix index exists
        # DV 1.10.0 creates index automatically — just verify files exist
        ls $ABS_OUT/{wildcards.sample}.deepvariant.vcf.gz $ABS_OUT/{wildcards.sample}.deepvariant.vcf.gz.tbi >> {log.out}

        VARIANT_COUNT=$(zcat {output.vcf} | grep -vc '^#' || echo "unknown")
        echo "  raw variants called : $VARIANT_COUNT" >> {log.out}

        # safety check: with a regions BED, the calls must not all come from one
        # chromosome when the BED spans several (a sign that only part of the BED was used)
        if [ -n "{input.bed}" ]; then
            BED_CHR=$(grep -v '^#' $ABS_BED | cut -f1 | sort -u | wc -l)
            VCF_CHR=$(zcat {output.vcf} | grep -v '^#' | cut -f1 | sort -u | wc -l)
            echo "  chromosomes: $BED_CHR in the BED, $VCF_CHR with DeepVariant calls" >> {log.out}
            if [ "$BED_CHR" -ge 3 ] && [ "$VCF_CHR" -le 1 ]; then
                echo "[FATAL] DeepVariant calls cover $VCF_CHR chromosome(s) but the BED spans $BED_CHR: the regions were not all processed" | tee -a {log.out} >> {log.err}
                exit 1
            fi
        fi
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] deepvariant_call complete" >> {log.out}
        """


# =============================================================================
# Rule: panel_filter_deepvariant
#
# Restricts DeepVariant VCF to panel BED regions and keeps only PASS variants.
# Mirrors panel_filter_clair3 — output goes to annotation/vcf_inputs/.
# Always uses full panel_bed (not ACTIVE_BED) so test-mode runs are correctly
# filtered to the real panel region.
# =============================================================================
rule panel_filter_deepvariant:
    input:
        vcf = DV_VCF_GZ,
        bed = config.get("panel_bed","") if config.get("panel_bed","") else []
    output:
        vcf = os.path.join(
            config["samples_root"], "{sample}",
            "annotation", "vcf_inputs", "{sample}.deepvariant.panel.vcf.gz"
        ),
        tbi = os.path.join(
            config["samples_root"], "{sample}",
            "annotation", "vcf_inputs", "{sample}.deepvariant.panel.vcf.gz.tbi"
        )
    params:
        bcftools = BCFTOOLS
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs",
                           "panel_filter_deepvariant.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs",
                           "panel_filter_deepvariant.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.vcf})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] panel_filter_deepvariant started" > {log.out}

        if [ -n "{input.bed}" ]; then
            {params.bcftools} view -R {input.bed} -f PASS -O z -o {output.vcf} {input.vcf} 2>> {log.err}
        else
            {params.bcftools} view -f PASS -O z -o {output.vcf} {input.vcf} 2>> {log.err}
        fi
        {params.bcftools} index --tbi {output.vcf} 2>> {log.err}

        VARIANT_COUNT=$(zcat {output.vcf} | grep -vc '^#' || echo "unknown")
        echo "  PASS panel variants : $VARIANT_COUNT" >> {log.out}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] panel_filter_deepvariant complete" >> {log.out}
        """


# =============================================================================
# Rule: vep_annotate_deepvariant
#
# Annotates the panel-filtered DeepVariant VCF using Ensembl VEP offline.
# Same flags as vep_annotate (Clair3) — consistent annotation for concordance.
# Output VEP TSV is consumed by indel_process (--deepvariant_vep).
# =============================================================================
rule vep_annotate_deepvariant:
    input:
        vcf = DV_VCF_PANEL,
        ref = config["ref_fasta"]
    output:
        tsv   = os.path.join(
            config["samples_root"], "{sample}", "annotation", "vep",
            "{sample}.deepvariant.vep.tsv"
        ),
        stats = os.path.join(
            config["samples_root"], "{sample}", "annotation", "vep",
            "{sample}.deepvariant.vep_stats.html"
        )
    params:
        vep      = VEP,
        cache    = config["vep"]["cache"],
        plugins  = config["vep"]["plugins"],
        assembly = config["vep"]["assembly"],
        flags    = config["vep"]["extra_flags"]
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs",
                           "vep_deepvariant.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs",
                           "vep_deepvariant.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.tsv})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] vep_annotate_deepvariant started" > {log.out}
        echo "  input  : {input.vcf}"    >> {log.out}
        echo "  cache  : {params.cache}" >> {log.out}

        {params.vep} \
            --input_file {input.vcf} \
            --output_file {output.tsv} \
            --stats_file {output.stats} \
            --fasta {input.ref} \
            --cache \
            --dir_cache {params.cache} \
            --dir_plugins {params.plugins} \
            --assembly {params.assembly} \
            {params.flags} \
            > {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] vep_annotate_deepvariant complete" >> {log.out}
        """


# =============================================================================

# -----------------------------------------------------------------------------
# Rule: coverage_build_tables
#
# Steps 1 + 2 of the coverage module:
#   1. Parse mosdepth regions BED → coverage_by_gene.tsv (adds target_size_bp)
#   2. Expand semicolon-joined gene labels, join against panel gene reference,
#      fill missing genes with depth=0 → coverage_complete.final.tsv
#
# Validates: gene count matches panel BED in the final TSV.
# -----------------------------------------------------------------------------
rule coverage_build_tables:
    input:
        regions_bed = MOSDEPTH_REPORT_BED_GZ,
        # the panel gene list; when not given, the report BED (gene name in column 4)
        gene_ref    = config.get("gene_intervals_ref", "") or config.get("report_bed", "") or config.get("panel_bed", "")
    output:
        by_gene = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.coverage_by_gene.tsv"
        ),
        final   = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.coverage.final.tsv"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "coverage_build_tables.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "coverage_build_tables.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.by_gene})
        mkdir -p $(dirname {log.out})

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_build_tables started" > {log.out}
        echo "  sample     : {wildcards.sample}"    >> {log.out}
        echo "  regions_bed: {input.regions_bed}"   >> {log.out}
        echo "  gene_ref   : {input.gene_ref}"      >> {log.out}

        python3 {SCRIPTS}/coverage_build_tables.py \
            --regions_bed {input.regions_bed} \
            --gene_ref    {input.gene_ref} \
            --out_by_gene {output.by_gene} \
            --out_final   {output.final} \
            --sample      {wildcards.sample} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_build_tables complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: coverage_summarise
#
# Steps 3 + 4 of the coverage module:
#   3. Compute summary metrics (mean depth, median, genes below 10x / 20x)
#      → coverage_summary.tsv
#   4. Write low-coverage gene tables
#      → genes_below_10x.tsv, genes_below_20x.tsv
#
# All outputs are derived from the master table only (no re-reading BED).
# -----------------------------------------------------------------------------
rule coverage_summarise:
    input:
        final = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.coverage.final.tsv"
        )
    output:
        summary  = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.coverage_summary.tsv"
        ),
        below10  = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.genes_below_10x.tsv"
        ),
        below20  = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.genes_below_20x.tsv"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "coverage_summarise.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "coverage_summarise.stderr.log")
    shell:
        """
        set -euo pipefail
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_summarise started" > {log.out}

        python3 {SCRIPTS}/coverage_summarise.py \
            --final_tsv   {input.final} \
            --out_summary {output.summary} \
            --out_below10 {output.below10} \
            --out_below20 {output.below20} \
            --sample      {wildcards.sample} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_summarise complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: coverage_plot
#
# Step 5 of the coverage module:
#   Horizontal bar plot of per-gene mean depth, colour-coded by tier.
#   Genes below 10x are labelled directly on the plot.
#   Validated design: ggplot2, red/amber/green tiers, dashed thresholds.
#
# Requires R with: ggplot2, readr, dplyr
# -----------------------------------------------------------------------------
# Per-gene coverage plot: rule coverage_plot_vertical (below)


# =============================================================================
# REPORT MODULE — QC + RUN SUMMARY + HTML
#
# Rule execution order:
#   nanoplot_qc
#     → parse_nanostats_qc       (NanoStats.txt → qc_quality_thresholds.tsv)
#     → plot_qc_quality_thresholds (TSV → PNG)
#
#   [sequencing_summary_file in config]
#   + nanoplot_qc
#     → build_run_summary        (seq_summary + NanoStats → run_summary TSVs)
#
#   mosdepth_coverage
#     → generate_shared_target_flags (regions BED → shared_target_flags.tsv)
#     → coverage_plot_vertical       (coverage TSV + flags → vertical PNG)
#
#   all of the above
#     → generate_report_html     (HTML report, self-contained with embedded plots)
# =============================================================================


# -----------------------------------------------------------------------------
# Rule: parse_nanostats_qc
#
# Extracts Q10–Q30 percentage values from NanoStats.txt.
# NanoStats.txt is the ONLY authoritative source for quality thresholds.
# Do NOT derive these from BAM or any other source.
# -----------------------------------------------------------------------------
rule parse_nanostats_qc:
    input:
        nanostats = os.path.join(
            config["samples_root"], "{sample}", "qc", "nanoplot", "{sample}_NanoStats.txt"
        )
    output:
        tsv = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.qc_thresholds.tsv"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "parse_nanostats_qc.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "parse_nanostats_qc.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] parse_nanostats_qc started" > {log.out}

        python3 {SCRIPTS}/parse_nanostats.py \
            --nanostats {input.nanostats} \
            --out_tsv   {output.tsv} \
            --sample    {wildcards.sample} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] parse_nanostats_qc complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: plot_qc_quality_thresholds
#
# Bar plot of Q10–Q30 quality thresholds. Locked colour scheme and axis limits.
# Requires R with: ggplot2, readr
# -----------------------------------------------------------------------------
rule plot_qc_quality_thresholds:
    input:
        tsv = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.qc_thresholds.tsv"
        )
    output:
        png = os.path.join(
            config["samples_root"], "{sample}", "report", "plots",
            "{sample}.qc_thresholds.png"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "plot_qc_quality_thresholds.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "plot_qc_quality_thresholds.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.png})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] plot_qc_quality_thresholds started" > {log.out}

        python3 {SCRIPTS}/plot_qc_thresholds.py \
            --input  {input.tsv} \
            --output {output.png} \
            >> {log.out} 2>> {log.err}

        if [ ! -s "{output.png}" ]; then
            echo "[ERROR] output PNG missing or empty" >&2; exit 1
        fi
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] plot_qc_quality_thresholds complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: build_run_summary
#
# Combines per-barcode metrics from the run sequencing_summary file with
# NanoPlot read-quality metrics into machine-readable and display TSVs.
#
# Sources (authoritative, do NOT mix):
#   sequencing_summary_*.txt → Number_of_reads, Total_yield_Gb, Run_duration_hours
#   NanoStats.txt            → Mean_read_length_bp, Read_N50_bp, Mean_read_quality
# -----------------------------------------------------------------------------
rule build_run_summary:
    input:
        nanostats   = os.path.join(
            config["samples_root"], "{sample}", "qc", "nanoplot", "{sample}_NanoStats.txt"
        )
    params:
        seq_summary = config["sequencing_summary_file"]
    output:
        summary = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.run_summary.tsv"
        ),
        display = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.run_summary_display.tsv"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "build_run_summary.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "build_run_summary.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.summary}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] build_run_summary started" > {log.out}

        python3 {SCRIPTS}/build_run_summary.py \
            --seq_summary  {params.seq_summary} \
            --nanostats    {input.nanostats} \
            --out_summary  {output.summary} \
            --out_display  {output.display} \
            --sample       {wildcards.sample} \
            --barcode_col  25 \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] build_run_summary complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: generate_shared_target_flags
#
# Detects semicolon-joined gene labels in the mosdepth BED (e.g. COL4A3;COL4A4)
# and writes a flags TSV used by the coverage plot and HTML report.
# Rule: do NOT divide depth by 2; annotate with * only.
# -----------------------------------------------------------------------------
rule generate_shared_target_flags:
    input:
        regions_bed = MOSDEPTH_BED_GZ if MODE == "panel" else []
    output:
        flags = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.shared_targets.tsv"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "shared_target_flags.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "shared_target_flags.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.flags}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] generate_shared_target_flags started" > {log.out}

        python3 {SCRIPTS}/generate_shared_target_flags.py \
            --regions_bed {input.regions_bed} \
            --out_flags   {output.flags} \
            --sample      {wildcards.sample} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] generate_shared_target_flags complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: coverage_plot_vertical
#
# Vertical bar plot of per-gene mean depth.
# Genes sharing a merged target interval are annotated with *.
# Depth labels shown for genes below 10x.
# Softened palette: #D97B66 / #E7BE73 / #7FB6A4
# Requires R with: ggplot2, readr, dplyr
# -----------------------------------------------------------------------------
if MODE == "wgs":
    rule coverage_plot_vertical:
        """WGS mode: per-chromosome barplot from mosdepth summary."""
        input:
            mosdepth_summary = os.path.join(
                config["samples_root"], "{sample}", "coverage",
                "{sample}.mosdepth.summary.txt"
            )
        output:
            png = os.path.join(
                config["samples_root"], "{sample}", "report", "plots",
                "{sample}.coverage.png"
            )
        log:
            out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "coverage_plot_vertical.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "coverage_plot_vertical.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.png})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_plot_wgs started" > {log.out}

            python3 {SCRIPTS}/plot_coverage_wgs.py \
                --input  {input.mosdepth_summary} \
                --output {output.png} \
                --sample {wildcards.sample} \
                >> {log.out} 2>> {log.err}

            if [ ! -s "{output.png}" ]; then
                echo "[ERROR] output PNG missing or empty" >&2; exit 1
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_plot_wgs complete" >> {log.out}
            """
else:
    rule coverage_plot_vertical:
        """Panel mode: per-gene barplot from coverage TSV."""
        input:
            final = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.coverage.final.tsv"
            ),
            flags = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.shared_targets.tsv"
            )
        output:
            png = os.path.join(
                config["samples_root"], "{sample}", "report", "plots",
                "{sample}.coverage.png"
            )
        log:
            out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "coverage_plot_vertical.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "coverage_plot_vertical.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.png})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_plot_vertical started" > {log.out}

            python3 {SCRIPTS}/plot_coverage.py \
                --input  {input.final} \
                --flags  {input.flags} \
                --output {output.png} \
                --sample {wildcards.sample} \
                >> {log.out} 2>> {log.err}

            if [ ! -s "{output.png}" ]; then
                echo "[ERROR] output PNG missing or empty" >&2; exit 1
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] coverage_plot_vertical complete" >> {log.out}
            """


# -----------------------------------------------------------------------------
# Rule: generate_report_html
#
# Generates the SignalScope per-sample HTML report.
# The report is fully self-contained — plots are embedded as base64 data URIs.
#
# Layout:
#   Header: SignalScope™ + metadata bar (sample / run ID / date / version)
#   Top row: [Sequencing Summary table] [QC quality plot]
#   Full width: Coverage vertical bar plot
#   Footer: coverage note + shared-target note
# -----------------------------------------------------------------------------
rule generate_igv_screenshots:
    # WGS mode (panel mode uses tree_igv in rules/interpretation.smk): screenshots from the
    # BAM the variants were called on (the primary BAM).
    input:
        bam             = BAM_FOR_CALLING,
        bai             = BAM_FOR_CALLING_BAI,
        snv_shortlist   = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.shortlist.with_gene_context.header.tsv"
        )),
        indel_shortlist = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "indel",
            "{sample}.indel.shortlist.with_gene_context.header.tsv"
        )),
        sv_shortlist    = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "sv",
            "{sample}.sv.shortlist.with_gene_context.header.tsv"
        ))
    output:
        manifest = os.path.join(
            config["samples_root"], "{sample}", "report", "igv",
            "{sample}.igv_manifest.tsv"
        )
    params:
        ref_fasta = config["ref_fasta"],
        igv_genome = config.get("igv", {}).get("genome", "hg38"),
        igv_sh    = config["tools"]["igv"],
        flank     = config.get("igv", {}).get("flank", 200),
        igv_mem   = config.get("igv", {}).get("memory", "4g")
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "generate_igv_screenshots.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "generate_igv_screenshots.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.manifest}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] generate_igv_screenshots started" > {log.out}

        python3 {SCRIPTS}/generate_igv_screenshots.py \
            --snv_shortlist   {input.snv_shortlist} \
            --indel_shortlist {input.indel_shortlist} \
            --sv_shortlist    {input.sv_shortlist} \
            --bam             {input.bam} \
            --ref_fasta       {params.ref_fasta} \
            --igv_genome      {params.igv_genome} \
            --out_dir         $(dirname {output.manifest}) \
            --sample          {wildcards.sample} \
            --igv_sh          {params.igv_sh} \
            --flank           {params.flank} \
            --igv_memory      {params.igv_mem} \
            --preflight       clean \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] generate_igv_screenshots complete" >> {log.out}
        """


# =============================================================================
# METHYLATION MODULE  [WGS only, opt-in via config methylation.enabled]
#
# Requires BAM with MM/ML tags (Dorado --modified-bases 5mCG_5hmCG).
# 1. modkit_pileup: genome-wide CpG methylation BED
# 2. check_clinical_methylation: screen curated clinical loci
# 3. extract_variant_methylation: per-variant methylation context (in-place)
# =============================================================================

if METH_ENABLED:

    rule modkit_pileup:
        """Run modkit pileup to extract CpG methylation fractions from BAM."""
        input:
            bam = BAM_PRI,
            bai = BAM_PRI_BAI,
        output:
            bed_gz  = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed.gz"
            ),
            bed_tbi = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed.gz.tbi"
            ),
        params:
            ref_fasta = config["ref_fasta"],
            modkit    = _meth_cfg.get("modkit_bin", "modkit"),
            raw_bed   = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed"
            ),
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs",
                               "modkit_pileup.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs",
                               "modkit_pileup.stderr.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {params.raw_bed}) $(dirname {log.out})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] modkit_pileup started" > {log.out}

            {params.modkit} pileup {input.bam} {params.raw_bed} \
                --cpg \
                --combine-strands \
                --combine-mods \
                --ref {params.ref_fasta} \
                --threads {threads} \
                --log-filepath {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] modkit pileup done, compressing..." >> {log.out}

            bgzip -c {params.raw_bed} > {output.bed_gz}
            tabix -p bed {output.bed_gz}
            rm -f {params.raw_bed}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] modkit_pileup complete" >> {log.out}
            """


    rule check_clinical_methylation:
        """Screen curated clinical methylation loci for abnormal patterns."""
        input:
            bed_gz = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed.gz"
            ),
            bed_tbi = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed.gz.tbi"
            ),
        output:
            tsv = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.clinical_methylation.tsv"
            ),
        params:
            clinical_bed = _meth_cfg.get("clinical_bed",
                os.path.join(os.path.dirname(workflow.snakefile), "assets",
                             "clinical_methylation_loci.bed")),
            min_coverage = _meth_cfg.get("min_coverage", 5),
        log:
            out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "check_clinical_methylation.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "check_clinical_methylation.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log.out})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] check_clinical_methylation started" > {log.out}

            python3 {SCRIPTS}/check_clinical_methylation.py \
                --bedmethyl    {input.bed_gz} \
                --clinical_bed {params.clinical_bed} \
                --output       {output.tsv} \
                --min_coverage {params.min_coverage} \
                >> {log.out} 2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] check_clinical_methylation complete" >> {log.out}
            """


    rule extract_variant_methylation:
        """Add per-variant methylation context to shortlisted SNVs and indels."""
        input:
            snv_shortlist = os.path.join(
                config["samples_root"], "{sample}", "report", "snv",
                "{sample}.snv.shortlist.with_gene_context.header.tsv"
            ),
            indel_shortlist = os.path.join(
                config["samples_root"], "{sample}", "report", "indel",
                "{sample}.indel.shortlist.with_gene_context.header.tsv"
            ),
            bed_gz = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed.gz"
            ),
            bed_tbi = os.path.join(
                config["samples_root"], "{sample}", "methylation",
                "{sample}.cpg.combined.bed.gz.tbi"
            ),
        output:
            sentinel = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.meth_annotated.done"
            ),
        params:
            window       = _meth_cfg.get("variant_window", 500),
            min_coverage = _meth_cfg.get("min_coverage", 5),
        log:
            out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "extract_variant_methylation.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "extract_variant_methylation.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.sentinel}) $(dirname {log.out})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] extract_variant_methylation started" > {log.out}

            # Annotate SNV shortlist in-place
            python3 {SCRIPTS}/extract_variant_methylation.py \
                --shortlist  {input.snv_shortlist} \
                --bedmethyl  {input.bed_gz} \
                --output     {input.snv_shortlist} \
                --window     {params.window} \
                --min_coverage {params.min_coverage} \
                >> {log.out} 2>> {log.err}

            # Annotate indel shortlist in-place
            python3 {SCRIPTS}/extract_variant_methylation.py \
                --shortlist  {input.indel_shortlist} \
                --bedmethyl  {input.bed_gz} \
                --output     {input.indel_shortlist} \
                --window     {params.window} \
                --min_coverage {params.min_coverage} \
                >> {log.out} 2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] extract_variant_methylation complete" >> {log.out}
            touch {output.sentinel}
            """


# =============================================================================
# WHATSHAP PHASING MODULE  [WGS only, opt-in via config phasing.enabled]
#
# 1. whatshap_phase: phase Clair3 VCF using long reads
# 2. whatshap_stats: phase block statistics
# 3. annotate_phasing: add phase info to shortlists + flag compound hets
# =============================================================================

if PHASE_ENABLED:

    rule whatshap_phase:
        """Phase variants using WhatsHap and long-read BAM."""
        input:
            vcf = CLAIR3_VCF_RAW,
            bam = BAM_PRI,
            bai = BAM_PRI_BAI,
        output:
            vcf = os.path.join(
                config["samples_root"], "{sample}", "phasing",
                "{sample}.phased.vcf.gz"
            ),
        params:
            ref_fasta = config["ref_fasta"],
            whatshap   = _phase_cfg.get("whatshap_bin", "whatshap"),
            raw_vcf    = os.path.join(
                config["samples_root"], "{sample}", "phasing",
                "{sample}.phased.vcf"
            ),
        log:
            out = os.path.join(config["samples_root"], "{sample}", "logs",
                               "whatshap_phase.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "logs",
                               "whatshap_phase.stderr.log")
        threads: 4
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {params.raw_vcf}) $(dirname {log.out})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] whatshap_phase started" > {log.out}

            {params.whatshap} phase \
                --indels \
                --ignore-read-groups \
                --reference={params.ref_fasta} \
                -o {params.raw_vcf} \
                {input.vcf} \
                {input.bam} \
                2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Phasing done, compressing..." >> {log.out}
            bgzip -c {params.raw_vcf} > {output.vcf}
            tabix -p vcf {output.vcf}
            rm -f {params.raw_vcf}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] whatshap_phase complete" >> {log.out}
            """


    rule whatshap_stats:
        """Compute phase block statistics from phased VCF."""
        input:
            vcf = os.path.join(
                config["samples_root"], "{sample}", "phasing",
                "{sample}.phased.vcf.gz"
            ),
        output:
            tsv = os.path.join(
                config["samples_root"], "{sample}", "phasing",
                "{sample}.phase_stats.tsv"
            ),
        params:
            whatshap = _phase_cfg.get("whatshap_bin", "whatshap"),
        log:
            err = os.path.join(config["samples_root"], "{sample}", "logs",
                               "whatshap_stats.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
            {params.whatshap} stats \
                --tsv={output.tsv} \
                {input.vcf} \
                2> {log.err}
            """


    rule annotate_phasing:
        """Add phase info to shortlists and flag compound het candidates."""
        input:
            phased_vcf      = os.path.join(
                config["samples_root"], "{sample}", "phasing",
                "{sample}.phased.vcf.gz"
            ),
            snv_shortlist   = os.path.join(
                config["samples_root"], "{sample}", "report", "snv",
                "{sample}.snv.shortlist.with_gene_context.header.tsv"
            ),
            indel_shortlist = os.path.join(
                config["samples_root"], "{sample}", "report", "indel",
                "{sample}.indel.shortlist.with_gene_context.header.tsv"
            ),
        output:
            compound_hets = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.compound_hets.tsv"
            ),
            sentinel = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.phase_annotated.done"
            ),
        log:
            out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "annotate_phasing.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "annotate_phasing.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.compound_hets}) $(dirname {log.out})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotate_phasing started" > {log.out}

            python3 {SCRIPTS}/annotate_phasing.py \
                --phased_vcf          {input.phased_vcf} \
                --shortlist_snv       {input.snv_shortlist} \
                --shortlist_indel     {input.indel_shortlist} \
                --output_snv          {input.snv_shortlist} \
                --output_indel        {input.indel_shortlist} \
                --output_compound_hets {output.compound_hets} \
                >> {log.out} 2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotate_phasing complete" >> {log.out}
            touch {output.sentinel}
            """


# -----------------------------------------------------------------------------
# Rule: annotate_avi  [WGS only, opt-in via config avi.enabled]
#
# Queries AlphaGenome Atlas for precomputed AVI scores on shortlisted
# SNVs and indels. Adds avi_phred, avi_raw, and avi_shap_* columns
# IN-PLACE to the shortlist TSVs. Writes a sentinel when done.
#
# Runs inside the "alphagenome" conda env (Python 3.10 + alphagenome SDK).
# Requires ALPHAGENOME_API_KEY env var (free, non-commercial).
# ~50 variants per sample, takes < 1 min on the Atlas precomputed endpoint.
# -----------------------------------------------------------------------------
if AVI_ENABLED:
    rule annotate_avi:
        input:
            snv_shortlist   = os.path.join(
                config["samples_root"], "{sample}", "report", "snv",
                "{sample}.snv.shortlist.with_gene_context.header.tsv"
            ),
            indel_shortlist = os.path.join(
                config["samples_root"], "{sample}", "report", "indel",
                "{sample}.indel.shortlist.with_gene_context.header.tsv"
            ),
        output:
            sentinel = os.path.join(
                config["samples_root"], "{sample}", "report", "data",
                "{sample}.avi_annotated.done"
            ),
        params:
            api_key_env = config.get("avi", {}).get("api_key_env", "ALPHAGENOME_API_KEY"),
            max_workers = config.get("avi", {}).get("max_workers", 4),
        log:
            out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "annotate_avi.stdout.log"),
            err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                               "annotate_avi.stderr.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.sentinel}) $(dirname {log.out})
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotate_avi started" > {log.out}

            python3 {SCRIPTS}/annotate_avi.py \
                --shortlist_snv   {input.snv_shortlist} \
                --shortlist_indel {input.indel_shortlist} \
                --output_snv      {input.snv_shortlist} \
                --output_indel    {input.indel_shortlist} \
                --api_key_env     {params.api_key_env} \
                --max_workers     {params.max_workers} \
                >> {log.out} 2>> {log.err}

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] annotate_avi complete" >> {log.out}
            touch {output.sentinel}
            """


rule generate_report_html:
    input:
        run_display  = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.run_summary_display.tsv"
        ),
        qc_plot      = os.path.join(
            config["samples_root"], "{sample}", "report", "plots",
            "{sample}.qc_thresholds.png"
        ) if MODE == "panel" else [],
        cov_plot     = os.path.join(
            config["samples_root"], "{sample}", "report", "plots",
            "{sample}.coverage.png"
        ),
        shared_flags = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.shared_targets.tsv"
        ) if MODE == "panel" else [],
        snv_summary  = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.summary.tsv"
        )),
        snv_clinvar  = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.clinvar_summary.normalized.tsv"
        )),
        snv_shortlist = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.shortlist.with_gene_context.header.tsv"
        )),
        indel_summary  = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "indel",
            "{sample}.indel.summary.tsv"
        )),
        indel_shortlist = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "indel",
            "{sample}.indel.shortlist.with_gene_context.header.tsv"
        )),
        sv_summary     = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "sv",
            "{sample}.sv.summary.tsv"
        )),
        sv_shortlist   = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "sv",
            "{sample}.sv.shortlist.with_gene_context.header.tsv"
        )),
        mosdepth_summary = ancient(os.path.join(
            config["samples_root"], "{sample}", "coverage",
            "{sample}.mosdepth.summary.txt"
        )),
        seq_perf_tsv     = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.sequencing_performance.tsv"
        )),
        tool_versions_tsv = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.tool_versions.tsv"
        )),
        igv_manifest = os.path.join(
            config["samples_root"], "{sample}", "report", "igv",
            "{sample}.igv_manifest.tsv"
        ),
        avi_sentinel = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.avi_annotated.done"
        )) if AVI_ENABLED else [],
        meth_sentinel = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.meth_annotated.done"
        )) if METH_ENABLED else [],
        meth_clinical = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.clinical_methylation.tsv"
        )) if METH_ENABLED else [],
        phase_sentinel = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.phase_annotated.done"
        )) if PHASE_ENABLED else [],
        compound_hets = ancient(os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.compound_hets.tsv"
        )) if PHASE_ENABLED else [],
    output:
        html = os.path.join(
            config["samples_root"], "{sample}", "report", "html",
            "{sample}.report.html"
        )
    params:
        mode           = MODE,
        run_id         = config["run_id"],
        report_date    = config.get("report_date", __import__("datetime").date.today().strftime("%d %B %Y")),
        report_version = config["report"]["version"],
        platform       = config["report"]["platform"],
        library        = config["report"]["library"],
        targeting      = config["report"]["targeting"],
        report_bed     = config.get("report_bed", "") or [],
        basecalling    = config["report"]["basecalling"],
        alignment      = config["report"]["alignment"],
        variants       = config["report"]["variants"],
        logo_path      = lambda wildcards: os.path.join(
            config["samples_root"], wildcards.sample, "report", "logo",
            config.get("report", {}).get("logo_filename", "SignalScope_logo.png")
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "generate_report_html.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "generate_report_html.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] generate_report_html started" > {log.out}

        python3 {SCRIPTS}/generate_report_html.py \
            --mode           {params.mode} \
            --sample         {wildcards.sample} \
            --run_id         "{params.run_id}" \
            --report_date    "{params.report_date}" \
            --run_display    {input.run_display} \
            $([ -n "{input.qc_plot}" ] && echo "--qc_plot {input.qc_plot}" || echo "") \
            --cov_plot       {input.cov_plot} \
            $([ -n "{input.shared_flags}" ] && echo "--shared_flags {input.shared_flags}" || echo "") \
            --report_version "{params.report_version}" \
            --platform       "{params.platform}" \
            --library        "{params.library}" \
            --targeting      "{params.targeting}" \
            --panel_bed      "{params.report_bed}" \
            --basecalling    "{params.basecalling}" \
            --alignment      "{params.alignment}" \
            --variants       "{params.variants}" \
            --snv_summary    {input.snv_summary} \
            --snv_clinvar    {input.snv_clinvar} \
            --snv_shortlist      {input.snv_shortlist} \
            --indel_summary      {input.indel_summary} \
            --indel_shortlist    {input.indel_shortlist} \
            --sv_summary         {input.sv_summary} \
            --sv_shortlist       {input.sv_shortlist} \
            --mosdepth_summary   {input.mosdepth_summary} \
            --seq_perf_tsv       {input.seq_perf_tsv} \
            --tool_versions_tsv  {input.tool_versions_tsv} \
            --logo_path          "{params.logo_path}" \
            --igv_dir            $(dirname {input.snv_shortlist})/../igv \
            $([ -n "{input.meth_clinical}" ] && echo "--meth_clinical {input.meth_clinical}" || echo "") \
            $([ -n "{input.compound_hets}" ] && echo "--compound_hets {input.compound_hets}" || echo "") \
            --out_html           {output.html} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] generate_report_html complete" >> {log.out}
        """


# =============================================================================
# SNV REPORT MODULE
#
# Implements the validated barcode06 SNV workflow per sample.
#
# Rule: snv_process
#            gene context reference
#   Outputs: snv.summary.tsv, snv.clinvar_summary.normalized.tsv,
#            snv.shortlist.with_gene_context.header.tsv
#            (plus intermediate snv.merged.tsv and snv.annotated.tsv)
#
# The generate_report_html rule already accepts these as optional inputs
# and renders the SNV section when they are present.
# =============================================================================

# -----------------------------------------------------------------------------
# Rule: snv_process
#
# Steps performed (see scripts/snv_process.py for full detail):
#   1.  Extract SNPs from Clair3 VCF
#   4.  Join VEP annotation (best transcript: MANE_SELECT > CANONICAL > first)
#   5.  Normalise ClinVar categories
#   6.  Apply shortlist filter: IMPACT HIGH/MODERATE + MAX_AF < threshold
#   7.  Join gene context (disease_summary, renal_relevance, inheritance,
#       phenotype_group) from shared gene context reference
#   8.  Write all output TSVs
# -----------------------------------------------------------------------------
rule snv_process:
    input:
        clair3_vcf       = CLAIR3_VCF_RAW,
        deepvariant_vcf  = DV_VCF_PANEL,
        panel_bed        = config.get("panel_bed", "") or [],
        vep_tsv          = VEP_TSV,
        dv_vep_tsv       = DV_VEP_TSV,
        gene_context_ref = config.get("gene_context_ref", "") or [],
        gene_features    = config.get("gene_features_ref", "") or []
    output:
        summary    = os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.summary.tsv"
        ),
        clinvar    = os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.clinvar_summary.normalized.tsv"
        ),
        shortlist  = os.path.join(
            config["samples_root"], "{sample}", "report", "snv",
            "{sample}.snv.shortlist.with_gene_context.header.tsv"
        )
    params:
        out_dir            = os.path.join(config["samples_root"], "{sample}", "report", "snv"),
        max_af             = EFF_MAX_AF,
        mode               = MODE,
        min_caller_support = EFF_MIN_CALLERS,
        max_shortlist      = EFF_MAX_SHORTLIST,
        fallback_bed       = WGS_FALLBACK_BED,
        fallback_gene_context = os.path.join(os.path.dirname(workflow.snakefile), "assets", "wgs_placeholders", "empty_gene_context.tsv") if MODE == "wgs" else ""
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "snv_process.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "snv_process.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p {params.out_dir} $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] snv_process started" > {log.out}
        echo "  sample         : {wildcards.sample}"          >> {log.out}
        echo "  mode           : {params.mode}"               >> {log.out}
        echo "  clair3_vcf     : {input.clair3_vcf}"          >> {log.out}
        echo "  deepvariant_vcf: {input.deepvariant_vcf}"     >> {log.out}
        echo "  vep_tsv        : {input.vep_tsv}"             >> {log.out}
        echo "  max_af         : {params.max_af}"             >> {log.out}

        PB=""; if [ -n "{input.panel_bed}" ]; then PB="--panel_bed {input.panel_bed}"; fi; if [ -z "$PB" ] && [ -n "{params.fallback_bed}" ]; then PB="--panel_bed {params.fallback_bed}"; fi
        GCR=""; if [ -n "{input.gene_context_ref}" ]; then GCR="--gene_context_ref {input.gene_context_ref}"; fi; if [ -z "$GCR" ] && [ -n "{params.fallback_gene_context}" ]; then GCR="--gene_context_ref {params.fallback_gene_context}"; fi
        GF=""; if [ -n "{input.gene_features}" ]; then GF="--gene_features {input.gene_features}"; fi

        python3 {SCRIPTS}/snv_process.py \
            --clair3_vcf       {input.clair3_vcf} \
            --deepvariant_vcf  {input.deepvariant_vcf} \
            $PB \
            --vep_tsv          {input.vep_tsv} \
            --deepvariant_vep  {input.dv_vep_tsv} \
            $GCR $GF \
            --out_dir          {params.out_dir} \
            --sample           {wildcards.sample} \
            --max_af           {params.max_af} \
            --mode             {params.mode} \
            --min_caller_support {params.min_caller_support} \
            --max_shortlist    {params.max_shortlist} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] snv_process complete" >> {log.out}
        """


# =============================================================================
# SEQUENCING PERFORMANCE + TOOL VERSIONS MODULE
# =============================================================================

# -----------------------------------------------------------------------------
# Rule: create_targeted_bam
#
# Filters the primary BAM to reads overlapping the panel BED.
# This isolates adaptive-sampling-accepted reads from rejected short reads,
# giving a clean view of actual sequencing performance on the target.
# -----------------------------------------------------------------------------

# =============================================================================
# Rule: indel_process
# =============================================================================
rule indel_process:
    input:
        clair3_vcf    = os.path.join(VCF_IN_DIR, f"{SAMPLE}.clair3.panel.vcf.gz"),
        vep_tsv       = os.path.join(VEP_DIR, f"{SAMPLE}.clair3.vep.tsv"),
        dv_vcf        = DV_VCF_PANEL,
        dv_vep_tsv    = DV_VEP_TSV,
        gene_context  = config.get("gene_context_ref", "") or [],
        gene_features = config.get("gene_features_ref", "") or [],
        panel_bed     = ACTIVE_BED
    output:
        summary   = os.path.join(
            config["samples_root"], "{sample}", "report", "indel",
            "{sample}.indel.summary.tsv"),
        shortlist = os.path.join(
            config["samples_root"], "{sample}", "report", "indel",
            "{sample}.indel.shortlist.with_gene_context.header.tsv")
    params:
        mode               = MODE,
        max_af             = EFF_MAX_AF,
        min_caller_support = EFF_MIN_CALLERS,
        max_shortlist      = EFF_MAX_SHORTLIST,
        fallback_bed       = WGS_FALLBACK_BED,
        fallback_gene_context = os.path.join(os.path.dirname(workflow.snakefile), "assets", "wgs_placeholders", "empty_gene_context.tsv") if MODE == "wgs" else ""
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "indel_process.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "indel_process.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.summary}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] indel_process started" > {log.out}
        PB=""; if [ -n "{input.panel_bed}" ]; then PB="--panel_bed {input.panel_bed}"; fi; if [ -z "$PB" ] && [ -n "{params.fallback_bed}" ]; then PB="--panel_bed {params.fallback_bed}"; fi
        GC=""; if [ -n "{input.gene_context}" ]; then GC="--gene_context {input.gene_context}"; fi; if [ -z "$GC" ] && [ -n "{params.fallback_gene_context}" ]; then GC="--gene_context {params.fallback_gene_context}"; fi
        GF=""; if [ -n "{input.gene_features}" ]; then GF="--gene_features {input.gene_features}"; fi

        python3 {SCRIPTS}/indel_process.py \
            --clair3_vcf       {input.clair3_vcf} \
            --vep_tsv          {input.vep_tsv} \
            --deepvariant_vcf  {input.dv_vcf} \
            --deepvariant_vep  {input.dv_vep_tsv} \
            $GC $GF $PB \
            --out_dir          $(dirname {output.summary}) \
            --sample           {wildcards.sample} \
            --max_af           {params.max_af} \
            --mode             {params.mode} \
            --min_caller_support {params.min_caller_support} \
            --max_shortlist    {params.max_shortlist} \
            >> {log.out} 2>> {log.err}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] indel_process complete" >> {log.out}
        """

# =============================================================================
# Rule: sv_process
# =============================================================================
rule sv_process:
    input:
        sniffles_vcf  = os.path.join(SNIFFLES_DIR, f"{SAMPLE}.sniffles.panel.vcf.gz"),
        cutesv_vcf    = os.path.join(CUTESV_DIR,   f"{SAMPLE}.cutesv.panel.vcf.gz"),
        gene_features = config.get("gene_features_ref", "") or [],
        gene_context  = config.get("gene_context_ref", "") or [],
        exon_bounds   = config.get("exon_boundaries_ref", "") or []
    output:
        summary   = os.path.join(
            config["samples_root"], "{sample}", "report", "sv",
            "{sample}.sv.summary.tsv"),
        shortlist = os.path.join(
            config["samples_root"], "{sample}", "report", "sv",
            "{sample}.sv.shortlist.with_gene_context.header.tsv")
    params:
        mode               = MODE,
        min_caller_support = EFF_MIN_CALLERS,
        max_shortlist      = EFF_MAX_SHORTLIST
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "sv_process.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "sv_process.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.summary}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] sv_process started" > {log.out}
        GF_FLAG=""; if [ -n "{input.gene_features}" ]; then GF_FLAG="--gene_features {input.gene_features}"; fi
        GC_FLAG=""; if [ -n "{input.gene_context}" ]; then GC_FLAG="--gene_context {input.gene_context}"; fi
        EB_FLAG=""; if [ -n "{input.exon_bounds}" ]; then EB_FLAG="--exon_boundaries {input.exon_bounds}"; fi

        python3 {SCRIPTS}/sv_process.py \
            --sniffles_vcf    {input.sniffles_vcf} \
            --cutesv_vcf      {input.cutesv_vcf} \
            $GF_FLAG $GC_FLAG $EB_FLAG \
            --out_dir         $(dirname {output.summary}) \
            --sample          {wildcards.sample} \
            --mode            {params.mode} \
            --min_caller_support {params.min_caller_support} \
            --max_shortlist   {params.max_shortlist} \
            >> {log.out} 2>> {log.err}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] sv_process complete" >> {log.out}
        """

# In WGS mode this rule is not in the DAG.
rule create_targeted_bam:
    input:
        bam = BAM_PRI,
        bai = BAM_PRI_BAI,
        bed = config.get("panel_bed","") if config.get("panel_bed","") else []
    output:
        bam = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted",
            "{sample}.targeted.bam"
        ),
        bai = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted",
            "{sample}.targeted.bam.bai"
        )
    params:
        samtools = SAMTOOLS
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs",
                           "create_targeted_bam.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs",
                           "create_targeted_bam.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] create_targeted_bam started" > {log.out}
        echo "  BED : {input.bed}" >> {log.out}

        if [ -n "{input.bed}" ]; then
            {params.samtools} view -@ {threads} -b -L {input.bed} {input.bam} > {output.bam} 2>> {log.err}
        else
            ln -sf $(realpath {input.bam}) {output.bam}
        fi

        {params.samtools} index -@ {threads} {output.bam} 2>> {log.err}

        READ_COUNT=$({params.samtools} flagstat {output.bam} | awk 'NR==1{{print $1}}')
        echo "  Targeted reads: $READ_COUNT" >> {log.out}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] create_targeted_bam complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: targeted_nanoplot
#
# Runs NanoPlot on the BED-filtered BAM to get targeted read statistics.
# Uses --bam mode, consistent with the whole-sample QC run.
# -----------------------------------------------------------------------------
rule targeted_nanoplot:
    input:
        bam = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted",
            "{sample}.targeted.bam"
        ),
        bai = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted",
            "{sample}.targeted.bam.bai"
        )
    output:
        nanostats = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted", "nanoplot",
            "{sample}_NanoStats.txt"
        )
    params:
        nanoplot = NANOPLOT,
        outdir   = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted", "nanoplot"
        ),
        prefix   = "{sample}_"
    log:
        out = os.path.join(config["samples_root"], "{sample}", "logs",
                           "targeted_nanoplot.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "logs",
                           "targeted_nanoplot.stderr.log")
    threads: config["threads"]
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] targeted_nanoplot started" > {log.out}

        {params.nanoplot} \
            --bam    {input.bam} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --threads {threads} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] targeted_nanoplot complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: build_sequencing_performance
#
# Combines whole-dataset metrics (from run_summary.tsv) with targeted-region
# metrics (from targeted BAM + NanoStats) into a single comparison TSV.
#
# Output: {sample}.sequencing_performance.tsv
# Columns: metric | whole_dataset | targeted_regions
# -----------------------------------------------------------------------------
rule build_sequencing_performance:
    input:
        run_summary      = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.run_summary.tsv"
        ),
        targeted_bam     = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted",
            "{sample}.targeted.bam"
        ),
        targeted_nanostats = os.path.join(
            config["samples_root"], "{sample}", "qc", "targeted", "nanoplot",
            "{sample}_NanoStats.txt"
        )
    params:
        seq_summary  = config.get("sequencing_summary_file", ""),
    output:
        tsv = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.sequencing_performance.tsv"
        )
    log:
        out = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "build_sequencing_performance.stdout.log"),
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "build_sequencing_performance.stderr.log")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log.out})
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] build_sequencing_performance started" > {log.out}

        python3 {SCRIPTS}/build_sequencing_performance.py \
            --run_summary        {input.run_summary} \
            --targeted_bam       {input.targeted_bam} \
            --targeted_nanostats {input.targeted_nanostats} \
            --seq_summary        {params.seq_summary} \
            --out_tsv            {output.tsv} \
            --sample             {wildcards.sample} \
            >> {log.out} 2>> {log.err}

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] build_sequencing_performance complete" >> {log.out}
        """


# -----------------------------------------------------------------------------
# Rule: record_tool_versions
#
# Captures version strings and model info for all key pipeline tools.
# Output: {sample}.tool_versions.tsv
# Columns: tool | version | model | path
#
# Errors in version capture are non-fatal (tool may emit to stderr).
# -----------------------------------------------------------------------------
rule record_tool_versions:
    input:
        model_check = MODEL_CHECK
    output:
        tsv = os.path.join(
            config["samples_root"], "{sample}", "report", "data",
            "{sample}.tool_versions.tsv"
        )
    params:
        dorado       = DORADO,
        start_from   = START_FROM,
        bam_path     = BAM_PRI,
        minimap2     = MINIMAP2,
        samtools     = SAMTOOLS,

        sniffles     = SNIFFLES,
        cutesv       = config["tools"]["cutesv"],
        clair3_sif   = config["clair3"]["sif"],
        vep_sif      = config["vep"]["sif"],
        dv_vcf_glob  = os.path.join(config["samples_root"], "{sample}", "variants", "small_variants", "deepvariant", "{sample}.deepvariant.vcf.gz"),
        dorado_model = config["dorado"]["model"],
        clair3_model = CLAIR3_MODEL,
        dv_sif_path  = config["deepvariant_sif"],
        dv_small     = "small model on" if (config.get("deepvariant") or {}).get("small_model", False) else "--disable_small_model",
    log:
        err = os.path.join(config["samples_root"], "{sample}", "report", "logs",
                           "record_tool_versions.stderr.log")
    shell:
        """
        # Explicitly disable strict mode for this rule — version commands
        # may exit non-zero and we want to capture "unavailable" gracefully.
        set +euo pipefail
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})

        echo -e "tool\tversion\tmodel\tpath" > {output.tsv}

        # Dorado + Minimap2: extract from BAM header when start_from=bam
        if [ "{params.start_from}" = "bam" ] && [ -f "{params.bam_path}" ]; then
            HDR=$(samtools view -H {params.bam_path} 2>/dev/null)
            # Basecaller: look for a dorado @PG or basecall_model in @RG; else report as not recorded
            BC_MODEL=$(echo "$HDR" | grep "^@RG" | grep -o "basecall_model=[^\t ]*" | head -1 | cut -d= -f2)
            DORADO_PG=$(echo "$HDR" | grep "^@PG" | grep -i "dorado\\|basecaller\\|guppy" | head -1)
            if [ -n "$DORADO_PG" ]; then
                BC_PN=$(echo "$DORADO_PG" | grep -o "PN:[^\t ]*" | cut -d: -f2)
                BC_VER=$(echo "$DORADO_PG" | grep -o "VN:[^\t ]*" | cut -d: -f2)
                echo -e "basecaller\t$BC_PN $BC_VER\t${{BC_MODEL:-see config}}\tExtracted from BAM header" >> {output.tsv}
            else
                echo -e "basecaller\tnot recorded in BAM\t${{BC_MODEL:-{params.dorado_model}}}\tBAM has no basecaller provenance; model from config" >> {output.tsv}
            fi
            # Aligner: the minimap2 (or other) @PG line
            AL_PG=$(echo "$HDR" | grep "^@PG" | grep -i "minimap2\\|bwa\\|winnowmap" | head -1)
            if [ -n "$AL_PG" ]; then
                AL_PN=$(echo "$AL_PG" | grep -o "PN:[^\t ]*" | cut -d: -f2)
                AL_VER=$(echo "$AL_PG" | grep -o "VN:[^\t ]*" | cut -d: -f2)
                echo -e "aligner\t$AL_PN $AL_VER\t-\tExtracted from BAM header" >> {output.tsv}
            else
                echo -e "aligner\tnot recorded in BAM\t-\tBAM has no aligner provenance" >> {output.tsv}
            fi
        else
            DORADO_VER=$({{ {params.dorado} --version 2>&1 || echo "unavailable"; }} | head -1)
            echo -e "dorado\t$DORADO_VER\t$(basename {params.dorado_model})\t{params.dorado}" >> {output.tsv}
            MM2_VER=$({{ {params.minimap2} --version 2>&1 || echo "unavailable"; }} | head -1)
            echo -e "minimap2\t$MM2_VER\t-\t{params.minimap2}" >> {output.tsv}
        fi

        # Samtools: emits version to stdout
        SAM_VER=$({{ {params.samtools} --version 2>&1 || echo "unavailable"; }} | head -1)
        echo -e "samtools\t$SAM_VER\t-\t{params.samtools}" >> {output.tsv}

        # Clair3 v2: query version from the Apptainer container
        CLAIR3_VER=$({{ apptainer exec {params.clair3_sif} /opt/bin/run_clair3.sh --version 2>&1 || echo "unavailable"; }} | head -1)
        CLAIR3_USED=$(awk '/^model: /{{print $2; exit}}' {input.model_check})
        echo -e "clair3\t$CLAIR3_VER\t${{CLAIR3_USED:-{params.clair3_model}}}\t{params.clair3_sif}" >> {output.tsv}

        # Sniffles: emits version to stdout or stderr depending on build
        SNIFFLES_VER=$({{ {params.sniffles} --version 2>&1 || echo "unavailable"; }} | head -1)
        echo -e "sniffles\t$SNIFFLES_VER\t-\t{params.sniffles}" >> {output.tsv}

        # CuteSV: emits version to stdout
        CUTESV_VER=$({{ {params.cutesv} --version 2>&1 || echo "unavailable"; }} | head -1)
        echo -e "cutesv\t$CUTESV_VER\t-\t{params.cutesv}" >> {output.tsv}

        # DeepVariant: runs in Apptainer; extract version from output VCF header
        DV_VCF={params.dv_vcf_glob}
        if [ -f "$DV_VCF" ]; then
            DV_VER=$(zcat "$DV_VCF" 2>/dev/null | grep "^##DeepVariant_version=" | head -1 | cut -d= -f2 || echo "unavailable")
        else
            DV_VER="unavailable"
        fi
        echo -e "deepvariant\t$DV_VER\tONT_R104 ({params.dv_small})\t{params.dv_sif_path}" >> {output.tsv}

        # VEP: query version from the Apptainer container
        VEP_VER=$({{ apptainer exec {params.vep_sif} vep 2>&1 | grep -i "ensembl-vep" | head -1 | grep -o "[0-9]*\\.[0-9]*" || echo "unavailable"; }})
        echo -e "vep\t$VEP_VER\tGRCh38 cache\t{params.vep_sif}" >> {output.tsv}
        echo "  Written: {output.tsv}" 2>> {log.err}
        """


# =============================================================================
# LAYER 3 (panel mode): decision tree, SV dosage, exon coverage, IGV, report
# =============================================================================
include: "rules/interpretation.smk"
