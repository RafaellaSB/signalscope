#!/usr/bin/env python3
# Update config.template.yaml for the 3-container architecture (Clair3 v2, DeepVariant, VEP 116).
# Run from signalscope/config/
import sys

p = "config.template.yaml"
s = open(p).read()

# 1. tools block: vep is now a container (remove from tools), keep dorado external
old_tools = '''  # External (not in conda env) — set your paths:
  dorado:   ""               # [SET THIS] ONT Dorado binary
  vep:      ""               # [SET THIS] Ensembl VEP'''
new_tools = '''  # External (not in conda env) — needed only for start_from=pod5 (basecalling):
  dorado:   ""               # [SET THIS if start_from=pod5] ONT Dorado binary'''
s = s.replace(old_tools, new_tools)

# 2. DeepVariant container comment
old_dv = '''# --- DeepVariant container (EXTERNAL — pull the pinned version) -------------
deepvariant_sif: ""          # [SET THIS] e.g. deepvariant_1.10.0-gpu.sif'''
new_dv = '''# --- Variant-caller containers (pull the pinned Apptainer/Singularity images) ---
# SNV/indel callers run from GPU containers (see README "External resources").
deepvariant_sif: ""          # [SET THIS] deepvariant_1.10.0-gpu.sif (apptainer pull docker://google/deepvariant:1.10.0-gpu)'''
s = s.replace(old_dv, new_dv)

# 3. Clair3: now a container (sif), not a runner script
old_clair3 = '''clair3:
  runner:    ""              # [SET THIS] run_clair3.sh
  model_dir: "models"
  model:     ""              # [SET THIS] matched to your basecaller'''
new_clair3 = '''clair3:
  sif:       ""              # [SET THIS] clair3_v2.0.1_gpu.sif (apptainer pull docker://hkubal/clair3:v2.0.1_gpu)
  model_dir: "models"        # shipped PyTorch models live here (relative to repo root)
  model:     ""              # [SET THIS] PyTorch model matching your basecaller:
                             #   Dorado sup@v4.2.0 -> clair3_pytorch_sup_v420
                             #   Dorado hac@v4.3.0 -> clair3_pytorch_hac_v430'''
s = s.replace(old_clair3, new_clair3)

# 4. VEP: now a container + cache
old_vep = '''vep:
  cache:    ""               # [SET THIS] Ensembl VEP cache (~20 GB, version-matched)
  plugins:  ""               # [SET THIS] VEP Plugins dir
  assembly: "GRCh38"
  extra_flags: "--everything --offline --format vcf --tab"
  # Add REVEL / SpliceAI plugin paths to extra_flags if available.'''
new_vep = '''vep:
  sif:      ""               # [SET THIS] ensembl-vep_release_116.0.sif (apptainer pull docker://ensemblorg/ensembl-vep:release_116.0)
  cache:    ""               # [SET THIS] VEP 116 GRCh38 cache dir (~26 GB; extract FULLY - see README)
  plugins:  ""               # [SET THIS] VEP Plugins dir (for REVEL/SpliceAI, optional)
  assembly: "GRCh38"
  extra_flags: "--everything --offline --format vcf --tab"
  # To add REVEL / SpliceAI, append to extra_flags:
  #   --plugin REVEL,file=/path/to/revel.tsv.gz
  #   --plugin SpliceAI,snv=/path/to/spliceai.vcf.gz,indel=/path/to/spliceai.vcf.gz'''
s = s.replace(old_vep, new_vep)

open(p, "w").write(s)
print("OK: template updated for 3-container architecture (Clair3 v2, DeepVariant, VEP 116)")
