#!/usr/bin/env bash
# Step 152: preview or submit D20/D50 reference-aware promoter UMAPs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/06_dna_tss_umap"
TEMPLATE="$ROOT/151_submit_dna_tss_reference_aware_template_32core.csub"
SOURCE_SUMMARY="$PROJECT/08_qc_summary/38_dna_tss_promoters/dna_tss_complete_summary.tsv"
JOB_ROOT="$PROJECT/08_qc_summary/39_dna_tss_umap/reference_aware_job_scripts"
CONFIRM="${CONFIRM_SUBMIT_DNA_TSS_REFERENCE_AWARE:-0}"
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
  "20000 D20k"
  "50000 D50k"
)
mkdir -p "$JOB_ROOT" "$PROJECT/09_logs/05_downstream_analysis/dna_tss_umap/reference_aware"
submitted=0
for specification in "${specifications[@]}"; do
  read -r minimum_total depth_tag <<< "$specification"
  run_id="promoter_cluster_${depth_tag}_TSSctx25_reference_aware"
  validated_summary="$PROJECT/08_qc_summary/39_dna_tss_umap/$run_id/dna_tss_umap_summary.tsv"
  if [[ "$FORCE" == "0" && -s "$validated_summary" ]] && \
    [[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$validated_summary")" == "validated" ]] && \
    [[ "$(awk -F '\t' '$1=="run_id"{print $2; exit}' "$validated_summary")" == "$run_id" ]]; then
    echo "[SKIP] validated $run_id"
    continue
  fi
  job="$JOB_ROOT/${depth_tag}.csub"
  sed \
    -e "s/__MIN_TOTAL_CALLS__/$minimum_total/g" \
    -e "s/__DEPTH_TAG__/$depth_tag/g" \
    -e "s/__FORCE__/$FORCE/g" \
    "$TEMPLATE" > "$job"
  grep -Eq '__MIN_TOTAL_CALLS__|__DEPTH_TAG__|__FORCE__' "$job" && {
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
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_TSS_REFERENCE_AWARE=1"
fi
