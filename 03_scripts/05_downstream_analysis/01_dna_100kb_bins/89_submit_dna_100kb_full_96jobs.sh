#!/usr/bin/env bash
# Step 89: preview/submit 96 resumable sample-level 100 kb jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
BASE="$SCRIPT_ROOT/05_downstream_analysis/01_dna_100kb_bins"
TEMPLATE="$BASE/88_submit_dna_100kb_sample_template_8core.csub"
TASKS="$PROJECT/08_qc_summary/20_single_cytosine_inputs/dna_single_cytosine_tasks.tsv"
JOB_ROOT="$PROJECT/08_qc_summary/27_dna_100kb/job_scripts"
START_SAMPLE="${START_SAMPLE:-1}"
END_SAMPLE="${END_SAMPLE:-96}"
CONFIRM="${CONFIRM_SUBMIT_DNA_100KB:-0}"
FORCE="${FORCE:-0}"

for name in START_SAMPLE END_SAMPLE CONFIRM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_SAMPLE >= 1 && END_SAMPLE <= 96 && START_SAMPLE <= END_SAMPLE && \
  CONFIRM <= 1 && FORCE <= 1 )) || exit 2
for path in "$TEMPLATE" "$TASKS"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
rows="$(awk -F '\t' 'NR > 1 {count++} END {print count+0}' "$TASKS")"
(( rows == 768 )) || {
  echo "[ERROR] expected 768 source tasks, found $rows" >&2
  exit 2
}

mkdir -p "$JOB_ROOT" "$PROJECT/09_logs/05_downstream_analysis/dna_100kb"
submitted=0
for ((index=START_SAMPLE; index<=END_SAMPLE; index++)); do
  job="$JOB_ROOT/sample_$(printf '%02d' "$index").csub"
  sed \
    -e "s/__SAMPLE_INDEX__/$index/g" \
    -e "s/__FORCE__/$FORCE/g" \
    "$TEMPLATE" > "$job"
  grep -Eq '__SAMPLE_INDEX__|__FORCE__' "$job" && {
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
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_100KB=1"
fi
