#!/usr/bin/env bash
#
# SignalScope pipeline runner
# ---------------------------
# Wraps Snakemake with the settings the pipeline needs so every run is
# reproducible and robust. Always uses the greedy scheduler (the ILP/pulp
# scheduler in Snakemake >=9.10 can crash with an IndexError) and enables
# safe resumption of interrupted runs. Runs Snakemake inside the project's
# conda environment.
#
# Usage:
#   ./run_signalscope.sh --config config/config.ckd.yaml --sample barcode06 \
#                        --start-from pod5
#   ./run_signalscope.sh --config config/config.cancer.yaml --sample barcode01 \
#                        --start-from bam --bam /path/to/sample.bam
#
# Options:
#   --config FILE       Config YAML (required)
#   --sample NAME       Sample name (required)
#   --start-from STAGE  pod5 | fastq | bam   (default: bam)
#   --bam FILE          Existing BAM (required if --start-from bam)
#   --fastq FILE        Existing FASTQ (required if --start-from fastq)
#   --cores N           CPU cores (default: 16)
#   --dir PATH          Working directory / repo root (default: script's dir)
#   --env NAME          Conda env with Snakemake >=9 (default: signalscope)
#   --dry-run           Show the plan without running
#   --                  Everything after -- is passed straight to snakemake
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$SCRIPT_DIR"

CONFIG=""
SAMPLE=""
START_FROM="bam"
BAM=""
FASTQ=""
CORES=16
CONDA_ENV="signalscope"
DRYRUN=""
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)     CONFIG="$2"; shift 2;;
    --sample)     SAMPLE="$2"; shift 2;;
    --start-from) START_FROM="$2"; shift 2;;
    --bam)        BAM="$2"; shift 2;;
    --fastq)      FASTQ="$2"; shift 2;;
    --cores)      CORES="$2"; shift 2;;
    --dir)        WORKDIR="$2"; shift 2;;
    --env)        CONDA_ENV="$2"; shift 2;;
    --dry-run)    DRYRUN="-n"; shift;;
    --)           shift; EXTRA=("$@"); break;;
    *) echo "Unknown option: $1" >&2; exit 2;;
  esac
done

[[ -z "$CONFIG" ]] && { echo "ERROR: --config is required" >&2; exit 2; }
[[ -z "$SAMPLE" ]] && { echo "ERROR: --sample is required" >&2; exit 2; }
[[ ! -f "$CONFIG" ]] && { echo "ERROR: config not found: $CONFIG" >&2; exit 2; }

CONFIG_ARGS=( "sample=$SAMPLE" "start_from=$START_FROM" )
if [[ "$START_FROM" == "bam" ]]; then
  [[ -z "$BAM" ]] && { echo "ERROR: --start-from bam requires --bam" >&2; exit 2; }
  [[ ! -f "$BAM" ]] && { echo "ERROR: BAM not found: $BAM" >&2; exit 2; }
  CONFIG_ARGS+=( "existing_bam=$BAM" )
elif [[ "$START_FROM" == "fastq" ]]; then
  [[ -z "$FASTQ" ]] && { echo "ERROR: --start-from fastq requires --fastq" >&2; exit 2; }
  CONFIG_ARGS+=( "existing_fastq=$FASTQ" )
fi

echo "[run_signalscope] env=$CONDA_ENV config=$CONFIG sample=$SAMPLE start_from=$START_FROM cores=$CORES"
echo "[run_signalscope] scheduler=greedy (ILP scheduler disabled for robustness)"

# Run Snakemake INSIDE the project conda env (which has Snakemake >=9).
# Robustness flags always applied:
#   --scheduler greedy   avoid the Snakemake >=9.10 ILP/pulp scheduler crash
#   --rerun-incomplete   safely resume an interrupted run
#   --keep-going         finish independent jobs even if one fails
exec conda run --no-capture-output -n "$CONDA_ENV" snakemake \
  --snakefile "$WORKDIR/pipeline/Snakefile" \
  --scheduler greedy \
  --rerun-incomplete \
  --keep-going \
  --configfile "$CONFIG" \
  --config "${CONFIG_ARGS[@]}" \
  --directory "$WORKDIR" \
  --cores "$CORES" \
  $DRYRUN \
  "${EXTRA[@]}"
