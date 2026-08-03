#!/usr/bin/env bash
# Step 44: render and optionally submit full DNA CB-aware dedup jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}/04_postmapping"
MANIFEST="$PROJECT/08_qc_summary/14_postmapping_full_inputs/dna_cb_markdup_full_manifest.tsv"
PROBE_SUMMARY="$PROJECT/08_qc_summary/13_postmapping_probe/dna_cb_markdup/dna_cb_markdup_probe_summary.tsv"
TEMPLATE="$POST_DIR/02_dna_deduplication/43_submit_dna_cb_markdup_full_template_64core.csub"
RENDERED="$PROJECT/09_logs/04_postmapping/dna_cb_markdup_full/jobs"
JOB_LOG_ROOT="$PROJECT/09_logs/04_postmapping/dna_cb_markdup_full"
START_INDEX="${START_INDEX:-1}"
END_INDEX="${END_INDEX:-96}"
CONFIRM="${CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL:-0}"
DIRECT_DNA_FULL="${DIRECT_DNA_FULL:-1}"

for name in START_INDEX END_INDEX CONFIRM DIRECT_DNA_FULL; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_INDEX >= 1 && END_INDEX <= 96 && START_INDEX <= END_INDEX )) || {
  echo "[ERROR] require 1 <= START_INDEX <= END_INDEX <= 96" >&2
  exit 2
}
(( CONFIRM == 0 || CONFIRM == 1 )) || exit 2
(( DIRECT_DNA_FULL == 0 || DIRECT_DNA_FULL == 1 )) || exit 2
for path in "$MANIFEST" "$TEMPLATE"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
rows=$(( $(wc -l < "$MANIFEST") - 1 ))
(( rows == 96 )) || {
  echo "[ERROR] expected 96 DNA full rows, found $rows" >&2
  exit 2
}
if (( DIRECT_DNA_FULL == 0 )); then
  [[ -s "$PROBE_SUMMARY" ]] || {
    echo "[ERROR] missing/empty Probe summary: $PROBE_SUMMARY" >&2
    exit 2
  }
  awk -F '\t' '
    $1 == "status" && $2 == "validated" {status=1}
    $1 == "pipeline_signature" &&
      $2 == "snmar_dna_cb_markdup_corrected_20260731" {signature=1}
    END {exit !(status && signature)}
  ' "$PROBE_SUMMARY" || {
    echo "[ERROR] corrected DNA CB-aware Probe is not validated" >&2
    exit 2
  }
else
  echo "[INFO] DIRECT_DNA_FULL=1: user-authorized full run without Probe"
fi
command -v csub >/dev/null 2>&1 || {
  echo "[ERROR] csub is unavailable" >&2
  exit 127
}

mkdir -p "$RENDERED" "$JOB_LOG_ROOT"
submitted=0
for ((index = START_INDEX; index <= END_INDEX; index++)); do
  sample="$(awk -F '\t' -v wanted="$index" \
    'NR > 1 && $1 == wanted {print $3; exit}' "$MANIFEST")"
  [[ -n "$sample" ]] || exit 2
  job="$RENDERED/sample_$(printf '%03d' "$index").csub"
  sed "s/__SAMPLE_INDEX__/$index/g" "$TEMPLATE" > "$job"
  grep -q '__SAMPLE_INDEX__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  echo "[PLAN] index=$index sample=$sample job=$job"
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted $submitted DNA CB-aware dedup jobs"
else
  echo "[DRY RUN] submitted 0; set CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=1"
fi
