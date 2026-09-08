# SignalScope

Panel-agnostic long-read clinical genomics pipeline for Oxford Nanopore adaptive sampling.

SignalScope takes Oxford Nanopore (ONT) sequencing data through a complete, automated workflow, from raw signal to an interactive clinical HTML report, for any targeted gene panel. It was developed for hereditary-disease gene panels and validated on both a chronic kidney disease (CKD) panel and a hereditary cancer panel.

## Overview

SignalScope is organised in three layers:

1. **Signal acquisition and alignment.** Dorado basecalling (SUP or HAC), alignment with minimap2 to GRCh38, and per-target coverage with mosdepth.
2. **Variant analysis.** Three independent variant-class tracks, each with two complementary callers:
   * SNV / indel: Clair3 and DeepVariant
   * Structural variants: Sniffles2 and CuteSV
   * Per-class concordance is computed at the panel level.
3. **Clinical interpretation and reporting.** VEP annotation, ACMG/AMP classification, variant shortlisting, and a self-contained HTML report with per-variant evidence and embedded IGV screenshots.

The pipeline is panel-agnostic: the panel BED is the single source of truth. Point it at any panel and it adapts, including automatic detection of pseudogene-recovery targets.

## Repository layout

```
signalscope/
├── run_signalscope.sh          # runner (wraps Snakemake with robust defaults)
├── setup_igv.sh                # one-time download of IGV + hg38 genome bundle
├── environment.yml             # conda environment (general tools)
├── config/
│   ├── config.template.yaml    # copy this and fill in your paths
│   ├── config.ckd.yaml         # worked example (CKD panel)
│   └── config.cancer.yaml      # worked example (hereditary cancer panel)
├── pipeline/
│   ├── Snakefile               # workflow definition
│   └── scripts/                # analysis and reporting scripts
└── resources/
    ├── panels/                 # panel BED files
    ├── gene_context/           # per-gene disease-context tables
    └── reporting/              # gene intervals, exon boundaries, gene features
```

## Requirements

* Linux with a CUDA-capable GPU (for Dorado basecalling and GPU variant callers)
* [conda](https://docs.conda.io/) / [mamba](https://mamba.readthedocs.io/)
* [Apptainer](https://apptainer.org/) (or Singularity) for the containerised callers
* `wget` and `unzip` (for the setup script)

## Installation

### 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate signalscope
```

This installs the general-purpose tools (minimap2, samtools, bcftools, Sniffles2, CuteSV, mosdepth, NanoPlot, Snakemake).

### 2. Pull the variant-caller containers

The deep-learning callers and the annotation engine run from pinned Apptainer images:

```bash
apptainer pull docker://hkubal/clair3:v2.0.1_gpu
apptainer pull docker://google/deepvariant:1.10.0-gpu
apptainer pull docker://ensemblorg/ensembl-vep:release_116.0
```

### 3. Download IGV and the hg38 genome bundle

IGV read-pileup screenshots use the official IGV desktop distribution (run headless). Download it once:

```bash
bash setup_igv.sh
```

### 4. Provide the reference data (not bundled)

The following are user-provided or downloaded separately, and their paths are set in the config file:

| Resource | Notes |
|---|---|
| Reference genome (GRCh38 `.fa` + `.fai`) | Contig naming is auto-harmonised to your panel BED. |
| Clair3 PyTorch model | Matched to the basecaller (for example, `clair3_pytorch_sup_v500` for Dorado `sup@v5.0.0`). |
| VEP cache (GRCh38, release 116) | About 26 GB. Extract the tarball completely; an incomplete extraction is detected and reported at startup. |
| VEP plugin data (optional) | REVEL and SpliceAI, if used. |
| Dorado and a basecalling model | Needed only when starting from POD5. |

## Quick start

Copy the template config, set your paths, then run:

```bash
cp config/config.template.yaml config/my_panel.yaml
# edit config/my_panel.yaml: reference, container SIF paths, VEP cache, panel BED, Clair3 model

# from raw signal (POD5):
./run_signalscope.sh --config config/my_panel.yaml --sample SAMPLE --start-from pod5

# from an aligned BAM:
./run_signalscope.sh --config config/my_panel.yaml --sample SAMPLE \
    --start-from bam --bam /path/to/sample.bam
```

The runner applies robust defaults (greedy scheduler, safe resumption of interrupted runs) so that a single command works reliably on any host.

## Adding a new panel

SignalScope adapts to any panel from three inputs:

1. A **panel BED** listing the target regions (one entry per gene in column 4).
2. A **gene-context table** (tab-separated) with each gene's disease, inheritance, relevance, and phenotype group. Generate one for any panel with `pipeline/scripts/build_gene_context.py` from the [GenCC](https://thegencc.org/) gene-disease assertions.
3. A **YAML config** defining references, tool paths, and thresholds.

No disease-specific logic is embedded in the core pipeline. Genes with high-identity pseudogenes (for example, *PKD1* or *PMS2*) are handled by a competitive-realignment recovery module; new pseudogene-affected genes are added by listing their coordinates in the recovery database, with no code changes.

## Output

For each sample, SignalScope produces a single self-contained HTML report containing:

* Run and sample metadata, sequencing performance, and adaptive-sampling enrichment
* Per-gene coverage across the panel
* Pseudogene-recovery summaries where applicable
* Shortlisted SNVs, indels, and structural variants, each with dual-caller concordance, pathogenicity annotation, gene-disease context, and an embedded IGV read-pileup screenshot

All plots and screenshots are embedded as base64 resources, so the report is portable as a single file with no external dependencies.

## Citation

If you use SignalScope, please cite:

> Barichello RS, Mallett AJ, Schmitz U. SignalScope: an end-to-end framework from nanopore signal to clinician-ready genomic reports. (Manuscript in preparation.)

## License

Released under the MIT License. See [LICENSE](LICENSE).
