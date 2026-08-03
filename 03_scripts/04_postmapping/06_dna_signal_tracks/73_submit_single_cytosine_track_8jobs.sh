#!/usr/bin/env bash
# Step 73: preview/submit eight resumable global WCG/GCH track jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}/04_postmapping"
TEMPLATE="$POST_DIR/06_dna_signal_tracks/72_submit_single_cytosine_track_template_8core.csub"
JOB_ROOT="$PROJECT/08_qc_summary/20_single_cytosine_inputs/dna_single_cytosine_track_jobs"
START_SHARD="${START_SHARD:-1}"
END_SHARD="${END_SHARD:-8}"
CONFIRM="${CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_TRACKS:-0}"
FORCE="${FORCE:-0}"

for name in START_SHARD END_SHARD CONFIRM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_SHARD >= 1 && END_SHARD <= 8 && START_SHARD <= END_SHARD && \
  CONFIRM <= 1 && FORCE <= 1 )) || exit 2
[[ -s "$TEMPLATE" ]] || exit 2
mkdir -p "$JOB_ROOT"
submitted=0
skipped=0
for ((shard=START_SHARD; shard<=END_SHARD; shard++)); do
  padded="$(printf '%02d' "$shard")"
  summary="$PROJECT/08_qc_summary/22_dna_single_cytosine_tracks/dna_single_cytosine.track_shard_${padded}.summary.tsv"
  if (( FORCE == 0 )) && [[ -s "$summary" ]] && \
    awk -F '\t' '$1 == "status" && $2 == "validated" {ok=1} END {exit !ok}' \
      "$summary"; then
    echo "[SKIP] shard=$shard"
    skipped=$((skipped + 1))
    continue
  fi
  job="$JOB_ROOT/shard_${padded}.csub"
  sed \
    -e "s/__SHARD__/$shard/g" \
    -e "s/__FORCE__/$FORCE/g" \
    "$TEMPLATE" > "$job"
  grep -Eq '__SHARD__|__FORCE__' "$job" && {
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
  echo "[OK] submitted=$submitted skipped_validated=$skipped"
else
  echo "[DRY RUN] no jobs submitted; set "
  echo "CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_TRACKS=1"
fi
