#!/usr/bin/env bash
# Step 57: preview/submit one 32-CPU full TAPS job per DNA sample.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
INPUT_ROOT="$PROJECT/08_qc_summary/17_downstream_full_inputs"
MANIFEST="$INPUT_ROOT/dna_taps_cpg_full_manifest.tsv"
TEMPLATE="$POST_DIR/56_submit_dna_taps_cpg_full_template_32core.csub"
JOB_ROOT="$INPUT_ROOT/dna_taps_cpg_full_job_scripts"
OUTPUT_ROOT="$PROJECT/07_postmapping/dna_taps_cpg_full"
START_INDEX="${START_INDEX:-1}"
END_INDEX="${END_INDEX:-96}"
CONFIRM="${CONFIRM_SUBMIT_DNA_TAPS_FULL:-0}"
FORCE="${FORCE:-0}"

for name in START_INDEX END_INDEX CONFIRM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_INDEX >= 1 && END_INDEX <= 96 && START_INDEX <= END_INDEX && \
  CONFIRM <= 1 && FORCE <= 1 )) || exit 2
for path in "$MANIFEST" "$TEMPLATE"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
rows="$(awk -F '\t' 'NR > 1 {count++} END {print count+0}' "$MANIFEST")"
(( rows == 96 )) || {
  echo "[ERROR] expected 96 DNA manifest rows, found $rows" >&2
  exit 2
}

mkdir -p "$JOB_ROOT"
submitted=0
skipped=0
for ((index=START_INDEX; index<=END_INDEX; index++)); do
  row="$(awk -F '\t' -v wanted="$index" \
    'NR > 1 && $1 == wanted {print; exit}' "$MANIFEST")"
  [[ -n "$row" ]] || exit 2
  IFS=$'\t' read -r observed run_id sample group retained input_bam <<< "$row"
  [[ "$observed" == "$index" ]] || exit 2
  summary="$OUTPUT_ROOT/$sample/$sample.taps_cpg_full.summary.tsv"
  if (( FORCE == 0 )) && [[ -s "$summary" ]] && \
    awk -F '\t' -v bam="$input_bam" \
      '$1 == "status" && $2 == "validated" {ok=1}
       $1 == "input_bam" && $2 == bam {same=1}
       END {exit !(ok && same)}' "$summary"; then
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
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_TAPS_FULL=1"
fi
