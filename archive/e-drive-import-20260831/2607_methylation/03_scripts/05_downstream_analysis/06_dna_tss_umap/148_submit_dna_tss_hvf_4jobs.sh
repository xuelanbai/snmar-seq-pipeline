#!/usr/bin/env bash
# Step 148: preview or submit four fixed high-variable promoter sensitivity jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/06_dna_tss_umap"
TEMPLATE="$ROOT/147_submit_dna_tss_hvf_template_32core.csub"
SOURCE_SUMMARY="$PROJECT/08_qc_summary/38_dna_tss_promoters/dna_tss_complete_summary.tsv"
JOB_ROOT="$PROJECT/08_qc_summary/39_dna_tss_umap/hvf_job_scripts"
CONFIRM="${CONFIRM_SUBMIT_DNA_TSS_HVF:-0}"
FORCE="${FORCE:-0}"

[[ "$CONFIRM" =~ ^[01]$ && "$FORCE" =~ ^[01]$ ]] || exit 2
for path in "$TEMPLATE" "$SOURCE_SUMMARY"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
[[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$SOURCE_SUMMARY")" == "validated" ]] || {
  echo "[ERROR] run and validate step 128 first" >&2
  exit 2
}

specifications=(
  "20000 1000 HVF1000_F20_C100"
  "20000 2000 HVF2000_F20_C100"
  "50000 1000 HVF1000_F20_C100"
  "50000 2000 HVF2000_F20_C100"
)

mkdir -p "$JOB_ROOT" "$PROJECT/09_logs/05_downstream_analysis/dna_tss_umap/hvf"
submitted=0
for specification in "${specifications[@]}"; do
  read -r minimum_total maximum_features run_tag <<< "$specification"
  depth_tag="D$((minimum_total / 1000))k"
  run_id="ortholog_single_gene_${depth_tag}_TSSctx100_${run_tag}"
  validated_summary="$PROJECT/08_qc_summary/39_dna_tss_umap/$run_id/dna_tss_umap_summary.tsv"
  if [[ "$FORCE" == "0" && -s "$validated_summary" ]] && \
    [[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$validated_summary")" == "validated" ]] && \
    [[ "$(awk -F '\t' '$1=="run_id"{print $2; exit}' "$validated_summary")" == "$run_id" ]]; then
    echo "[SKIP] validated $run_id"
    continue
  fi
  job="$JOB_ROOT/${depth_tag}_${run_tag}.csub"
  sed \
    -e "s/__MIN_TOTAL_CALLS__/$minimum_total/g" \
    -e "s/__MAX_FEATURES__/$maximum_features/g" \
    -e "s/__RUN_TAG__/$run_tag/g" \
    -e "s/__FORCE__/$FORCE/g" \
    "$TEMPLATE" > "$job"
  grep -Eq '__MIN_TOTAL_CALLS__|__MAX_FEATURES__|__RUN_TAG__|__FORCE__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  if [[ "$CONFIRM" == "1" ]]; then
    csub < "$job"
    submitted=$((submitted + 1))
  else
    echo "[DRY RUN] csub < $job"
  fi
done

if [[ "$CONFIRM" == "1" ]]; then
  echo "[OK] submitted=$submitted"
else
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_TSS_HVF=1"
fi
