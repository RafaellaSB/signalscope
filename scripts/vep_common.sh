# =============================================================================
# vep_common.sh - shared VEP plugin settings and checks, sourced by every VEP step so all
# runs use the same plugins.
#
# SpliceAI: the VEP 116 SpliceAI plugin needs BOTH snv= and indel=; with only snv= it does not
# load, VEP only warns and the output has no SpliceAI column. When no indel score file is
# installed, the SNV file is also given as indel=, so the plugin loads and indel lookups find
# nothing; indels are then scored by spliceai_fill.py.
# =============================================================================
: "${PDATA:?PDATA is not set: the pipeline sets it from config.yaml (resources: vep_plugin_data)}"
SAI_SNV=$PDATA/spliceai/spliceai_scores.masked.snv.hg38.vcf.gz
SAI_INDEL=$PDATA/spliceai/spliceai_scores.masked.indel.hg38.vcf.gz
if [ -s "$SAI_INDEL" ]; then
  SAI_NOTE="SpliceAI: SNV + indel scores"
else
  SAI_INDEL=$SAI_SNV
  SAI_NOTE="SpliceAI: SNV scores only (no indel file; SNV file passed as indel= so the plugin loads)"
fi
SPLICEAI_ARG="--plugin SpliceAI,snv=$SAI_SNV,indel=$SAI_INDEL"
REVEL_ARG="--plugin REVEL,file=$PDATA/revel/new_tabbed_revel_grch38.tsv.gz"

# true when a VEP tab file has the plugin columns (a file without them is rebuilt)
vep_has_plugins() {
  local h
  h=$(grep -m1 "^#Uploaded_variation" "$1" 2>/dev/null) || return 1
  echo "$h" | grep -q "SpliceAI_pred" && echo "$h" | grep -q "REVEL"
}

# remove the previous run's warnings file before VEP starts: VEP writes it only when there
# are warnings, so an old one would otherwise survive a clean run and be read by vep_check
vep_prepare() { rm -f "${1}_warnings.txt" "$1.done"; }

# stop the run when any plugin failed to load or a plugin column is missing
vep_check() {
  local out=$1 w=${1}_warnings.txt
  if [ -s "$w" ] && grep -q -i "Failed to instantiate plugin" "$w"; then
    echo "[FATAL] a VEP plugin did not load ($w):"; grep -i "Failed to instantiate plugin" "$w" | sort -u
    rm -f "$out.done"; exit 1
  fi
  if ! vep_has_plugins "$out"; then
    echo "[FATAL] $out has no SpliceAI_pred or REVEL column"; rm -f "$out.done"; exit 1
  fi
  local n s
  n=$(grep -vc "^#" "$out")
  s=$(awk -F'\t' '/^#Uploaded/{for(i=1;i<=NF;i++) if($i=="SpliceAI_pred") c=i; next}
                  !/^#/ && c && $c!="-" && $c!="" {k++} END{print k+0}' "$out")
  echo "  VEP check OK: plugins loaded; $s of $n rows carry a SpliceAI score ($SAI_NOTE)"
}
