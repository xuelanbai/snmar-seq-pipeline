#!/usr/bin/env bash
# Step 70: preview/submit 96 context-aware cytosine sample validators.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
TASKS="$PROJECT/08_qc_summary/20_single_cytosine_inputs/dna_single_cytosine_tasks.tsv"
TEMPLATE="$POST_DIR/69_submit_dna_single_cytosine_validation_template_8core.csub"
JOB_ROOT="$PROJECT/08_qc_summary/20_single_cytosine_inputs/dna_single_cytosine_validation_jobs"
START_INDEX="${START_INDEX:-1}"
END_INDEX="${END_INDEX:-96}"
CONFIRM="${CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_VALIDATION:-0}"
FORCE="${FORCE:-0}"
VALIDATION_SIGNATURE="snmar_single_cytosine_validation_20260801"

for name in START_INDEX END_INDEX CONFIRM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_INDEX >= 1 && END_INDEX <= 96 && START_INDEX <= END_INDEX && \
  CONFIRM <= 1 && FORCE <= 1 )) || exit 2
for path in "$TASKS" "$TEMPLATE"; do
  [[ -s "$path" ]] || exit 2
done

mkdir -p "$JOB_ROOT"
submitted=0
skipped=0
for ((index=START_INDEX; index<=END_INDEX; index++)); do
  sample="$(awk -F '\t' -v wanted="$index" \
    'NR > 1 && $2 == wanted {print $3; exit}' "$TASKS")"
  [[ -n "$sample" ]] || exit 2
  summary="$PROJECT/08_qc_summary/21_dna_single_cytosine_validation/$sample/$sample.single_cytosine.validation_summary.tsv"
  if (( FORCE == 0 )) && [[ -s "$summary" ]] && \
    awk -F '\t' -v signature="$VALIDATION_SIGNATURE" \
      '$1 == "status" && $2 == "validated" {status=1} \
       $1 == "pipeline_signature" && $2 == signature {version=1} \
       END {exit !(status && version)}' \
      "$summary"; then
    echo "[SKIP] $index $sample"
    skipped=$((skipped + 1))
    continue
  fi
  job="$JOB_ROOT/sample_$(printf '%03d' "$index").csub"
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
    echo "[DRY RUN] csub < $job  # $sample"
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted=$submitted skipped_validated=$skipped"
else
  echo "[DRY RUN] no jobs submitted; set "
  echo "CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_VALIDATION=1"
fi
