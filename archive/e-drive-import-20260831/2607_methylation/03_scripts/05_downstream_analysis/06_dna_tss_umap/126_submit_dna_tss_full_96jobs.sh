#!/usr/bin/env bash
# Step 126: preview or submit 96 resumable physical-well TSS jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/06_dna_tss_umap"
TEMPLATE="$ROOT/125_submit_dna_tss_well_template_8core.csub"
INDEX_SUMMARY="$PROJECT/08_qc_summary/38_dna_tss_promoters/dna_tss_promoter_index_summary.tsv"
TASKS="$PROJECT/08_qc_summary/20_single_cytosine_inputs/dna_single_cytosine_tasks.tsv"
JOB_ROOT="$PROJECT/08_qc_summary/38_dna_tss_promoters/job_scripts"
START_WELL="${START_WELL:-1}"
END_WELL="${END_WELL:-96}"
CONFIRM="${CONFIRM_SUBMIT_DNA_TSS:-0}"
FORCE="${FORCE:-0}"

for name in START_WELL END_WELL CONFIRM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_WELL >= 1 && END_WELL <= 96 && START_WELL <= END_WELL && \
  CONFIRM <= 1 && FORCE <= 1 )) || exit 2
for path in "$TEMPLATE" "$INDEX_SUMMARY" "$TASKS"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
[[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$INDEX_SUMMARY")" == "validated" ]] || {
  echo "[ERROR] run and validate step 124 first" >&2
  exit 2
}
rows="$(awk -F '\t' 'NR > 1 {count++} END {print count+0}' "$TASKS")"
(( rows == 768 )) || {
  echo "[ERROR] expected 768 single-cytosine tasks, found $rows" >&2
  exit 2
}

mkdir -p "$JOB_ROOT" "$PROJECT/09_logs/05_downstream_analysis/dna_tss_umap"
submitted=0
for ((index=START_WELL; index<=END_WELL; index++)); do
  job="$JOB_ROOT/well_$(printf '%02d' "$index").csub"
  sed \
    -e "s/__WELL_INDEX__/$index/g" \
    -e "s/__FORCE__/$FORCE/g" \
    "$TEMPLATE" > "$job"
  grep -Eq '__WELL_INDEX__|__FORCE__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  else
    echo "[DRY RUN] csub < $job"
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted=$submitted"
else
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_TSS=1"
fi
