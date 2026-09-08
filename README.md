# SignalScope

**Panel-agnostic long-read clinical genomics pipeline for Oxford Nanopore adaptive sampling.**

SignalScope takes Oxford Nanopore (ONT) sequencing data through a complete, automated
workflow — from raw signal to an interactive clinical HTML report — for any targeted
gene panel. It was developed for hereditary-disease gene panels and validated on both a
chronic kidney disease (CKD) panel and a hereditary cancer panel.

---

## Overview

SignalScope is organised in three layers:

1. **Signal acquisition & alignment** — Dorado basecalling (SUP or HAC), alignment with
   minimap2 to GRCh38 (or T2T-CHM13), and per-target coverage with mosdepth.
2. **Variant analysis** — three independent variant-class tracks, each with two
   complementary callers:
   - SNV / indel: **Clair3** + **DeepVariant**
   - Structural variants: **Sniffles2** + **CuteSV**
   - Per-class concordance is computed at the panel level.
3. **Clinical interpretation & reporting** — VEP annotation, ACMG/AMP classification,
   variant shortlisting, and a self-contained HTML report with per-variant evidence and
   embedded IGV screenshots.

The pipeline is **panel-agnostic**: the panel BED is the single source of truth. Point it
at any panel and it adapts, including automatic detection of pseudogene-recovery targets.

---

## Repository layout

```
signalscope/
├── environment.yml            # conda environment (tool versions)
├── environment.lock.yml       # fully pinned lock file
├── pipeline/
│   ├── Snakefile              # the workflow
│   ├── scripts/               # processing & reporting scripts
│   └── assets/                # report logo
├── config/
│   ├── config.ckd.yaml        # worked example: CKD panel (SUP basecalling)
│   └── config.cancer.yaml     # worked example: hereditary cancer panel (HAC basecalling)
├── resources/
│   ├── panels/                # panel BEDs
│   ├── gene_context/          # per-gene clinical context tables
│   └── reporting/             # gene intervals, features, exon boundaries
└── models/                    # Clair3 + Dorado models (chemistry-matched)
```

---

## Installation

### 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate signalscope
```

This installs the general-purpose tools (minimap2, samtools, bcftools, Sniffles2, CuteSV,
mosdepth, NanoPlot, Snakemake, IGV).

### 2. External resources (not bundled — obtain separately)

The general-purpose tools install via the conda environment (step 1). The three
deep-learning / large-data components run from pinned Apptainer (Singularity)
containers, and the reference data are user-provided. Pull the containers and set
all paths in your config file.

**Containers** (pull once with Apptainer):

```bash
apptainer pull docker://hkubal/clair3:v2.0.1_gpu           # Clair3 v2.0.1 (GPU)
apptainer pull docker://google/deepvariant:1.10.0-gpu      # DeepVariant 1.10.0 (GPU)
apptainer pull docker://ensemblorg/ensembl-vep:release_116.0   # Ensembl VEP 116
```

| Resource | Notes |
|---|---|
| **Reference genome** (GRCh38 `.fa` + `.fai`) | User-provided. Naming (`chr` vs no-`chr`) is auto-harmonised to your panel BED. |
| **Clair3 v2.0.1** container + PyTorch model | Model must match the basecaller: Dorado `sup@v4.2.0` → `clair3_pytorch_sup_v420`; `hac@v4.3.0` → `clair3_pytorch_hac_v430`. Matched PyTorch models are shipped in `models/`. |
| **DeepVariant 1.10.0** container | Run via Apptainer with `--disable_small_model` (set by the pipeline). ONT_R104 model. |
| **Ensembl VEP 116** container + cache (GRCh38) | Cache is ~26 GB. **Extract the tarball completely** — an interrupted extraction silently drops variants on missing chromosomes (the pipeline now checks for this and fails fast). Optional plugin data: REVEL, SpliceAI. |
| **Dorado** + basecalling model | From ONT; needed only for `start_from=pod5`. GPU required for basecalling. |

---

## Usage

Run from the repository root (so `resources/` paths resolve). Provide the sample name and
entry point on the command line.

**From an aligned BAM** (most common):

```bash
conda activate signalscope
snakemake \
  --configfile config/config.ckd.yaml \
  --config sample=SAMPLE start_from=bam existing_bam=/path/to/aligned.bam \
  --directory . --cores 16
```

Entry points (`start_from`):
- `pod5` — full run from raw POD5 (includes Dorado basecalling; GPU)
- `fastq` — from basecalled FASTQ
- `bam` — from an aligned BAM (with secondary alignments retained, for pseudogene recovery)

The final report is written to `samples/{sample}/report/html/{sample}.report.html`.

---

## Applying SignalScope to a new panel

1. Provide your panel as a BED file (`resources/panels/`).
2. Generate a per-gene context table with the included script:
   ```bash
   python pipeline/scripts/build_gene_context_gencc.py \
       --gencc gencc-submissions.tsv --bed your_panel.bed \
       --relevance_col relevance --out your_gene_context
   ```
   (GenCC submissions export: https://search.thegencc.org/download)
3. Point a copy of a config file at your panel, reference, and matched models.
4. Run as above.

---

## Basecaller / model matching

The Clair3 model **must** match the Dorado basecalling model, which depends on the ONT
chemistry. SignalScope ships the two validated pairs:

| Basecaller | Clair3 model |
|---|---|
| Dorado `sup@v4.2.0` | `r1041_e82_400bps_sup_v420` |
| Dorado `hac@v4.3.0` | `r1041_e82_400bps_hac_v430` |

Set `dorado.model` and `clair3.model` in the config to match your data.

---

## Citation

_(Manuscript in preparation.)_

Key methods and resources: GenCC (DiStefano et al., Genet Med 2022); Clair3 (Zheng et al.,
Nat Comput Sci 2022); DeepVariant (Poplin et al., Nat Biotechnol 2018); Sniffles2 (Smolka
et al., Nat Biotechnol 2024); CuteSV (Jiang et al., Genome Biol 2020); minimap2 (Li,
Bioinformatics 2018); Ensembl VEP (McLaren et al., Genome Biol 2016).

---

## License

_(To be added.)_
