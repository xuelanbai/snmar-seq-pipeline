#!/usr/bin/env bash
# Step 35 worker: mark/remove DNA PCR duplicates within CB, never across cells.
#
# DNA has no molecular UMI.  A duplicate family is therefore defined by:
#   CB + reference/mate references + unclipped paired-template endpoints
#   + pair orientation
# as implemented by samtools fixmate/markdup in template mode.  Read length
# alone is not a duplicate key.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 RUN_ID INPUT_BAM" >&2
  exit 2
fi
RUN_ID="$1"
INPUT_BAM="$2"

: "${PROJECT:?PROJECT is required}"
: "${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
: "${QC_ROOT:?QC_ROOT is required}"
: "${TMP_ROOT:?TMP_ROOT is required}"
: "${SAMTOOLS_BIN:?SAMTOOLS_BIN is required}"

PYTHON_BIN="${PYTHON_BIN:-$(command -v python || true)}"
POST_DIR="${POST_DIR:-${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}/04_postmapping}"
FORCE="${FORCE:-0}"
NAME_SORT_THREADS="${NAME_SORT_THREADS:-20}"
FIXMATE_THREADS="${FIXMATE_THREADS:-4}"
COORD_SORT_THREADS="${COORD_SORT_THREADS:-32}"
MARKDUP_THREADS="${MARKDUP_THREADS:-12}"
FILTER_THREADS="${FILTER_THREADS:-8}"
NAME_SORT_MEMORY="${NAME_SORT_MEMORY:-1G}"
COORD_SORT_MEMORY="${COORD_SORT_MEMORY:-1G}"
KEEP_MARKED_BAM="${KEEP_MARKED_BAM:-0}"
PIPELINE_SIGNATURE="snmar_dna_cb_markdup_corrected_20260731"

[[ "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]] || {
  echo "[ERROR] unsafe RUN_ID: $RUN_ID" >&2
  exit 2
}
[[ -s "$INPUT_BAM" && -s "$INPUT_BAM.bai" ]] || {
  echo "[ERROR] input BAM or BAI missing/empty: $INPUT_BAM" >&2
  exit 2
}
[[ -x "$SAMTOOLS_BIN" && -x "$PYTHON_BIN" ]] || exit 127
[[ "$FORCE" =~ ^[01]$ ]] || exit 2
for name in NAME_SORT_THREADS FIXMATE_THREADS COORD_SORT_THREADS \
  MARKDUP_THREADS FILTER_THREADS; do
  value="${!name}"
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || {
    echo "[ERROR] $name must be a positive integer: $value" >&2
    exit 2
  }
done
[[ "$KEEP_MARKED_BAM" =~ ^[01]$ ]] || {
  echo "[ERROR] KEEP_MARKED_BAM must be 0 or 1" >&2
  exit 2
}
(( NAME_SORT_THREADS + FIXMATE_THREADS + COORD_SORT_THREADS <= 64 )) || {
  echo "[ERROR] streaming sort/fixmate thread total exceeds 64" >&2
  exit 2
}
(( MARKDUP_THREADS + FILTER_THREADS <= 64 )) || {
  echo "[ERROR] markdup/filter thread total exceeds 64" >&2
  exit 2
}
[[ -s "$POST_DIR/90_optional_validation/38_summarize_dna_cb_markdup_probe.py" ]] || exit 2

run_dir="$OUTPUT_ROOT/$RUN_ID"
tmp_dir="$TMP_ROOT/${RUN_ID}.$$"
final_bam="$run_dir/${RUN_ID}.cb_dedup.bam"
final_bai="$final_bam.bai"
marked_bam="$run_dir/${RUN_ID}.cb_marked.bam"
marked_bai="$marked_bam.bai"
stats="$QC_ROOT/by_run/${RUN_ID}.markdup.stats.txt"
done_file="$QC_ROOT/by_run/${RUN_ID}.done.tsv"
tmp_coord="$tmp_dir/${RUN_ID}.fixmate.coord.bam"
tmp_bam="$run_dir/.${RUN_ID}.cb_dedup.bam.tmp.$$"
tmp_marked="$run_dir/.${RUN_ID}.cb_marked.bam.tmp.$$"
tmp_stats="$stats.tmp.$$"

mkdir -p "$run_dir" "$QC_ROOT/by_run" "$tmp_dir"

cleanup() {
  rc=$?
  rm -f "$tmp_bam" "$tmp_bam.bai" "$tmp_marked" \
    "$tmp_marked.bai" "$tmp_stats" "$tmp_coord" "$tmp_coord.bai"
  if [[ -d "$tmp_dir" ]]; then
    find "$tmp_dir" -mindepth 1 -maxdepth 1 -type f -delete
    rmdir "$tmp_dir" 2>/dev/null || true
  fi
  exit "$rc"
}
trap cleanup EXIT INT TERM

if (( FORCE == 0 )) && [[ -s "$final_bam" && -s "$final_bai" && \
  -s "$stats" && -s "$done_file" ]] && \
  awk -F '\t' -v signature="$PIPELINE_SIGNATURE" '
    $1 == "pipeline_signature" && $2 == signature {signed=1}
    $1 == "status" && $2 == "completed" {complete=1}
    END {exit !(signed && complete)}
  ' "$done_file"; then
  "$SAMTOOLS_BIN" quickcheck -v "$final_bam"
  if (( KEEP_MARKED_BAM == 1 )); then
    [[ -s "$marked_bam" && -s "$marked_bai" ]] || {
      echo "[ERROR] KEEP_MARKED_BAM=1 but marked BAM is missing" >&2
      exit 2
    }
    "$SAMTOOLS_BIN" quickcheck -v "$marked_bam"
  fi
  echo "[INFO] corrected CB-aware dedup output exists; skipping: $RUN_ID"
  exit 0
fi
if (( FORCE == 1 )); then
  rm -f "$final_bam" "$final_bai" "$marked_bam" "$marked_bai" \
    "$stats" "$done_file"
fi

"$SAMTOOLS_BIN" quickcheck -v "$INPUT_BAM"
echo "[INFO] CB-aware duplicate removal: $RUN_ID"
echo "[INFO] qname-sort=$NAME_SORT_THREADS fixmate=$FIXMATE_THREADS coordinate-sort=$COORD_SORT_THREADS"

set +e
set -o pipefail
# In samtools sort, -l is BAM compression level.
"$SAMTOOLS_BIN" sort \
  -n -@ "$NAME_SORT_THREADS" -m "$NAME_SORT_MEMORY" -u \
  -T "$tmp_dir/qname" "$INPUT_BAM" |
  "$SAMTOOLS_BIN" fixmate \
    -m -u -@ "$FIXMATE_THREADS" - - |
  "$SAMTOOLS_BIN" sort \
    -@ "$COORD_SORT_THREADS" -m "$COORD_SORT_MEMORY" -l 1 \
    -T "$tmp_dir/coordinate" -o "$tmp_coord" -
pipeline_rc=$?
set -e
(( pipeline_rc == 0 )) || {
  echo "[ERROR] name-sort/fixmate/coordinate-sort pipeline failed" >&2
  exit "$pipeline_rc"
}
[[ -s "$tmp_coord" ]] || {
  echo "[ERROR] fixmate coordinate BAM is empty" >&2
  exit 2
}

if (( KEEP_MARKED_BAM == 1 )); then
  "$SAMTOOLS_BIN" markdup \
    -@ "$MARKDUP_THREADS" -c -s -m t \
    --barcode-tag CB \
    -f "$tmp_stats" \
    "$tmp_coord" "$tmp_marked"
  "$SAMTOOLS_BIN" quickcheck -v "$tmp_marked"
  "$SAMTOOLS_BIN" view \
    -@ "$FILTER_THREADS" -b -F 1024 \
    -o "$tmp_bam" "$tmp_marked"
else
  set +e
  set -o pipefail
  "$SAMTOOLS_BIN" markdup \
    -@ "$MARKDUP_THREADS" -c -s -m t \
    --barcode-tag CB \
    -f "$tmp_stats" \
    "$tmp_coord" - |
    "$SAMTOOLS_BIN" view \
      -@ "$FILTER_THREADS" -b -F 1024 \
      -o "$tmp_bam" -
  markdup_rc=$?
  set -e
  (( markdup_rc == 0 )) || {
    echo "[ERROR] corrected markdup/filter pipeline failed" >&2
    exit "$markdup_rc"
  }
fi
"$SAMTOOLS_BIN" quickcheck -v "$tmp_bam"
mv -f "$tmp_bam" "$final_bam"
mv -f "$tmp_stats" "$stats"
"$SAMTOOLS_BIN" index -@ "$MARKDUP_THREADS" "$final_bam"
if (( KEEP_MARKED_BAM == 1 )); then
  mv -f "$tmp_marked" "$marked_bam"
  "$SAMTOOLS_BIN" index -@ "$MARKDUP_THREADS" "$marked_bam"
else
  rm -f "$marked_bam" "$marked_bai"
fi

tmp_done="$done_file.tmp.$$"
{
  printf 'metric\tvalue\n'
  printf 'status\tcompleted\n'
  printf 'pipeline_signature\t%s\n' "$PIPELINE_SIGNATURE"
  printf 'run_id\t%s\n' "$RUN_ID"
  printf 'input_bam\t%s\n' "$INPUT_BAM"
  printf 'output_bam\t%s\n' "$final_bam"
  printf 'markdup_stats\t%s\n' "$stats"
  printf 'duplicate_scope\twithin_CB_and_reference_species\n'
  printf 'duplicate_key\tCB+reference+unclipped_template_endpoints+orientation\n'
  printf 'markdup_mode\ttemplate\n'
  printf 'markdup_expected_max_read_length\tdefault_300\n'
  printf 'duplicate_flag_filter\t0x400\n'
  printf 'duplicates_removed\t1\n'
  printf 'marked_bam_retained\t%s\n' "$KEEP_MARKED_BAM"
  if (( KEEP_MARKED_BAM == 1 )); then
    printf 'marked_bam\t%s\n' "$marked_bam"
  fi
  printf 'completed_at\t%s\n' "$(date -Is)"
} > "$tmp_done"
mv -f "$tmp_done" "$done_file"

echo "[OK] CB-aware deduplicated BAM: $final_bam"
