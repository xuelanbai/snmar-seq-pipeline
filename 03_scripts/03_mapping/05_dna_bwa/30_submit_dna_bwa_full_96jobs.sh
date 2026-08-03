#!/usr/bin/env bash
# Step 30: render and submit one 64-CPU BWA job per DNA sample.
#
# This launcher is a dry run unless CONFIRM_SUBMIT_DNA_BWA=1.
# START_INDEX and END_INDEX allow staged submission, for example 1-24.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
MAPPING_DIR="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}/03_mapping"
INPUT_QC="${DNA_BWA_INPUT_QC_ROOT:-$PROJECT/08_qc_summary/10_dna_bwa_inputs}"
SAMPLE_TABLE="$INPUT_QC/dna_bwa_samples.tsv"
PROBE_QC="${DNA_BWA_PROBE_QC_ROOT:-$PROJECT/08_qc_summary/11_dna_bwa_probe}"
PROBE_SUMMARY="$PROBE_QC/dna_bwa_probe_summary.tsv"
TEMPLATE="$MAPPING_DIR/05_dna_bwa/29_submit_dna_bwa_sample_template_64core.csub"
RENDERED_ROOT="$PROJECT/09_logs/03_mapping/dna_bwa_submit_scripts"
JOB_LOG_ROOT="$PROJECT/09_logs/03_mapping/dna_bwa_full_jobs"
INDEX_PREFIX="${DNA_BWA_INDEX_PREFIX:-$PROJECT/02_reference/dna_human_pig_bwa/human_pig.prefixed}"
INDEX_DONE="$(dirname "$INDEX_PREFIX")/.snmar_bwa_index_done.tsv"
START_INDEX="${START_INDEX:-1}"
END_INDEX="${END_INDEX:-96}"
CONFIRM="${CONFIRM_SUBMIT_DNA_BWA:-0}"

for name in START_INDEX END_INDEX CONFIRM; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] $name must be a non-negative integer: $value" >&2
    exit 2
  }
done
(( START_INDEX >= 1 && END_INDEX <= 96 && START_INDEX <= END_INDEX )) || {
  echo "[ERROR] require 1 <= START_INDEX <= END_INDEX <= 96" >&2
  exit 2
}
(( CONFIRM == 0 || CONFIRM == 1 )) || {
  echo "[ERROR] CONFIRM_SUBMIT_DNA_BWA must be 0 or 1" >&2
  exit 2
}
for path in "$SAMPLE_TABLE" "$PROBE_SUMMARY" "$TEMPLATE" "$INDEX_DONE"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
probe_status="$(awk -F '\t' '$1=="status"{print $2; exit}' "$PROBE_SUMMARY")"
probe_runs="$(awk -F '\t' '$1=="completed_runs"{print $2; exit}' "$PROBE_SUMMARY")"
probe_pairs="$(awk -F '\t' '$1=="expected_read_pairs"{print $2; exit}' "$PROBE_SUMMARY")"
probe_markdup="$(awk -F '\t' '$1=="global_duplicate_marking"{print $2; exit}' "$PROBE_SUMMARY")"
[[ "$probe_status" == "validated" && "$probe_runs" == "8" && \
   "$probe_pairs" == "8000000" && \
   "$probe_markdup" == "disabled_multi_cell_CB_aware_required" ]] || {
  echo "[ERROR] DNA BWA probe summary has not passed the required code-level checks" >&2
  exit 2
}
sample_rows=$(( $(wc -l < "$SAMPLE_TABLE") - 1 ))
(( sample_rows == 96 )) || {
  echo "[ERROR] expected 96 samples in manifest, found $sample_rows" >&2
  exit 2
}
command -v csub >/dev/null 2>&1 || {
  echo "[ERROR] csub is unavailable" >&2
  exit 127
}

mkdir -p "$RENDERED_ROOT" "$JOB_LOG_ROOT"
echo "[INFO] requested sample-index range: $START_INDEX-$END_INDEX"
if (( CONFIRM == 0 )); then
  echo "[DRY RUN] no jobs will be submitted"
  echo "[DRY RUN] set CONFIRM_SUBMIT_DNA_BWA=1 to submit this range"
fi

submitted=0
for ((index = START_INDEX; index <= END_INDEX; index++)); do
  sample="$(awk -F '\t' -v wanted="$index" 'NR == wanted + 1 {print $2; exit}' "$SAMPLE_TABLE")"
  [[ -n "$sample" ]] || {
    echo "[ERROR] no sample for index $index" >&2
    exit 2
  }
  rendered="$RENDERED_ROOT/sample_$(printf '%03d' "$index").csub"
  sed "s/__SAMPLE_INDEX__/$index/g" "$TEMPLATE" > "$rendered"
  if grep -q '__SAMPLE_INDEX__' "$rendered"; then
    echo "[ERROR] unresolved template placeholder in $rendered" >&2
    exit 2
  fi
  echo "[PLAN] index=$index sample=$sample job_file=$rendered"
  if (( CONFIRM == 1 )); then
    csub < "$rendered"
    submitted=$((submitted + 1))
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted $submitted independent 64-CPU DNA BWA jobs"
else
  echo "[DRY RUN] rendered $((END_INDEX - START_INDEX + 1)) job files; submitted 0"
fi
