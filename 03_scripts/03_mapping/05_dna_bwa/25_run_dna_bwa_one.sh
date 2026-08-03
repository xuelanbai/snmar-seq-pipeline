#!/usr/bin/env bash
# Step 25 helper: map one DNA probe lane or one complete 8-lane sample.
#
# Required environment:
#   RUN_ID SAMPLE R1_FILES R2_FILES EXPECTED_PAIRS OUTPUT_ROOT QC_ROOT
#   BWA_BIN SAMTOOLS_BIN PYTHON_BIN INDEX_PREFIX MAPPING_DIR TMP_ROOT
#
# R1_FILES and R2_FILES are comma-separated, ordered path lists. In full mode
# all lane FASTQs are concatenated in the same order and mapped as one sample,
# while retaining CB tags for later cell-aware processing.

set -euo pipefail

required_vars=(
  RUN_ID SAMPLE R1_FILES R2_FILES EXPECTED_PAIRS OUTPUT_ROOT QC_ROOT
  BWA_BIN SAMTOOLS_BIN PYTHON_BIN INDEX_PREFIX MAPPING_DIR TMP_ROOT
)
for name in "${required_vars[@]}"; do
  [[ -n "${!name:-}" ]] || {
    echo "[ERROR] missing required environment variable: $name" >&2
    exit 2
  }
done

MODE="${MODE:-full}"
MAX_PAIRS="${MAX_PAIRS:-0}"
FORCE="${FORCE:-0}"
REVALIDATE_EXISTING="${REVALIDATE_EXISTING:-0}"
BWA_THREADS="${BWA_THREADS:-44}"
COORD_SORT_THREADS="${COORD_SORT_THREADS:-12}"
DECOMPRESS_THREADS="${DECOMPRESS_THREADS:-2}"
VALIDATE_THREADS="${VALIDATE_THREADS:-4}"
SORT_MEMORY="${SORT_MEMORY:-1G}"
LIBRARY="${LIBRARY:-snMAR_DNA}"
PIGZ_BIN="${PIGZ_BIN:-$(command -v pigz || true)}"

for name in MAX_PAIRS FORCE REVALIDATE_EXISTING BWA_THREADS COORD_SORT_THREADS \
  DECOMPRESS_THREADS VALIDATE_THREADS; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] $name must be a non-negative integer: $value" >&2
    exit 2
  }
done
[[ "$MODE" == "full" || "$MODE" == "probe" ]] || {
  echo "[ERROR] MODE must be full or probe: $MODE" >&2
  exit 2
}
(( EXPECTED_PAIRS >= 1 && FORCE <= 1 && REVALIDATE_EXISTING <= 1 && \
   BWA_THREADS >= 1 && \
   COORD_SORT_THREADS >= 1 && \
   DECOMPRESS_THREADS >= 1 && VALIDATE_THREADS >= 1 )) || exit 2
if [[ "$MODE" == "full" ]]; then
  (( MAX_PAIRS == 0 )) || {
    echo "[ERROR] full mode requires MAX_PAIRS=0" >&2
    exit 2
  }
  [[ -n "$PIGZ_BIN" && -x "$PIGZ_BIN" ]] || {
    echo "[ERROR] pigz is required for full-mode streaming" >&2
    exit 127
  }
else
  (( MAX_PAIRS == EXPECTED_PAIRS && MAX_PAIRS >= 1 )) || {
    echo "[ERROR] probe mode requires MAX_PAIRS=EXPECTED_PAIRS>=1" >&2
    exit 2
  }
fi

for binary_name in BWA_BIN SAMTOOLS_BIN PYTHON_BIN; do
  binary="${!binary_name}"
  [[ -x "$binary" ]] || {
    echo "[ERROR] unavailable executable $binary_name=$binary" >&2
    exit 127
  }
done
for suffix in amb ann bwt pac sa; do
  [[ -s "$INDEX_PREFIX.$suffix" ]] || {
    echo "[ERROR] missing/empty BWA index component: $INDEX_PREFIX.$suffix" >&2
    exit 2
  }
done
[[ -s "$MAPPING_DIR/05_dna_bwa/24_validate_dna_bwa_bam.py" ]] || exit 2
[[ "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]] || {
  echo "[ERROR] unsafe RUN_ID: $RUN_ID" >&2
  exit 2
}
[[ "$SAMPLE" =~ ^[A-Za-z0-9_.-]+$ ]] || {
  echo "[ERROR] unsafe SAMPLE: $SAMPLE" >&2
  exit 2
}

IFS=',' read -r -a r1_files <<< "$R1_FILES"
IFS=',' read -r -a r2_files <<< "$R2_FILES"
(( ${#r1_files[@]} == ${#r2_files[@]} && ${#r1_files[@]} >= 1 )) || {
  echo "[ERROR] R1/R2 file-list lengths differ or are empty" >&2
  exit 2
}
if [[ "$MODE" == "full" ]]; then
  (( ${#r1_files[@]} == 8 )) || {
    echo "[ERROR] full sample must contain exactly 8 lane pairs; found ${#r1_files[@]}" >&2
    exit 2
  }
else
  (( ${#r1_files[@]} == 1 )) || {
    echo "[ERROR] probe run must contain exactly one temporary FASTQ pair" >&2
    exit 2
  }
fi
for path in "${r1_files[@]}" "${r2_files[@]}"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty input FASTQ: $path" >&2
    exit 2
  }
done

run_dir="$OUTPUT_ROOT/$RUN_ID"
qc_dir="$QC_ROOT/by_run"
log_dir="$run_dir/logs"
tmp_dir="$TMP_ROOT/${RUN_ID}.$$"
final_bam="$run_dir/${RUN_ID}.sorted.bam"
final_bai="$final_bam.bai"
validation_json="$qc_dir/${RUN_ID}.validation.json"
done_file="$run_dir/.snmar_dna_bwa_done.tsv"
tmp_bam="$run_dir/.${RUN_ID}.sorted.bam.tmp.$$"
bwa_log="$log_dir/${RUN_ID}.bwa.stderr.log"

mkdir -p "$run_dir" "$qc_dir" "$log_dir" "$tmp_dir"

if (( FORCE == 0 )) && [[ -s "$final_bam" && -s "$final_bai" && -s "$validation_json" && -s "$done_file" ]]; then
  "$SAMTOOLS_BIN" quickcheck -v "$final_bam"
  if (( REVALIDATE_EXISTING == 1 )); then
    echo "[INFO] validated output exists; running requested full revalidation: $RUN_ID"
    "$PYTHON_BIN" "$MAPPING_DIR/05_dna_bwa/24_validate_dna_bwa_bam.py" \
      --bam "$final_bam" \
      --sample "$SAMPLE" \
      --run-id "$RUN_ID" \
      --expected-pairs "$EXPECTED_PAIRS" \
      --samtools "$SAMTOOLS_BIN" \
      --threads "$VALIDATE_THREADS" \
      --out "$validation_json"
  else
    echo "[INFO] validated output and done marker exist; quickcheck passed, skipping: $RUN_ID"
  fi
  exit 0
fi

stream_pids=()
cleanup() {
  rc=$?
  for pid in "${stream_pids[@]}"; do
    kill "$pid" 2>/dev/null || true
  done
  rm -f "$tmp_bam" "$tmp_bam.bai"
  if [[ -d "$tmp_dir" ]]; then
    rm -f "$tmp_dir/R1.fifo" "$tmp_dir/R2.fifo"
    find "$tmp_dir" -mindepth 1 -maxdepth 1 -type f -delete
    rmdir "$tmp_dir" 2>/dev/null || true
  fi
  exit "$rc"
}
trap cleanup EXIT INT TERM

if (( FORCE == 1 )); then
  rm -f "$final_bam" "$final_bai" "$validation_json" "$done_file"
fi

input_r1="${r1_files[0]}"
input_r2="${r2_files[0]}"
if [[ "$MODE" == "full" ]]; then
  input_r1="$tmp_dir/R1.fifo"
  input_r2="$tmp_dir/R2.fifo"
  mkfifo "$input_r1" "$input_r2"

  # Only the source mismatch-vector tag is renamed. sed reads the complete
  # decompressed stream and therefore does not truncate either mate.
  (
    "$PIGZ_BIN" -dc -p "$DECOMPRESS_THREADS" "${r1_files[@]}" |
      LC_ALL=C sed $'s/\tXE:Z:/\tZE:Z:/g' > "$input_r1"
  ) &
  stream_pids+=("$!")
  (
    "$PIGZ_BIN" -dc -p "$DECOMPRESS_THREADS" "${r2_files[@]}" |
      LC_ALL=C sed $'s/\tXE:Z:/\tZE:Z:/g' > "$input_r2"
  ) &
  stream_pids+=("$!")
fi

printf -v rg_line '@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPL:ILLUMINA' \
  "$RUN_ID" "$SAMPLE" "$LIBRARY"
echo "[INFO] run=$RUN_ID mode=$MODE pairs=$EXPECTED_PAIRS lanes=${#r1_files[@]}"
echo "[INFO] BWA=$BWA_THREADS coordinate-sort=$COORD_SORT_THREADS"
echo "[INFO] global duplicate marking is disabled because this is a multi-cell CB-tagged BAM"

set +e
set -o pipefail
"$BWA_BIN" mem \
  -t "$BWA_THREADS" \
  -M \
  -C \
  -R "$rg_line" \
  "$INDEX_PREFIX" \
  "$input_r1" "$input_r2" \
  2> "$bwa_log" |
  "$SAMTOOLS_BIN" sort \
    -@ "$COORD_SORT_THREADS" -m "$SORT_MEMORY" -O bam \
    -T "$tmp_dir/coordinate_sort" -o "$tmp_bam" -
pipeline_rc=$?

stream_rc=0
for pid in "${stream_pids[@]}"; do
  wait "$pid" || stream_rc=$?
done
stream_pids=()
set -e

if (( pipeline_rc != 0 || stream_rc != 0 )); then
  echo "[ERROR] mapping pipeline failed: pipeline_rc=$pipeline_rc stream_rc=$stream_rc" >&2
  echo "[ERROR] inspect BWA log: $bwa_log" >&2
  exit 1
fi

[[ -s "$tmp_bam" ]] || {
  echo "[ERROR] mapping pipeline produced no BAM: $tmp_bam" >&2
  exit 2
}
mv -f "$tmp_bam" "$final_bam"
"$SAMTOOLS_BIN" index -@ "$VALIDATE_THREADS" "$final_bam"

"$PYTHON_BIN" "$MAPPING_DIR/05_dna_bwa/24_validate_dna_bwa_bam.py" \
  --bam "$final_bam" \
  --sample "$SAMPLE" \
  --run-id "$RUN_ID" \
  --expected-pairs "$EXPECTED_PAIRS" \
  --samtools "$SAMTOOLS_BIN" \
  --threads "$VALIDATE_THREADS" \
  --out "$validation_json"

tmp_done="$done_file.tmp.$$"
{
  printf 'metric\tvalue\n'
  printf 'status\tvalidated\n'
  printf 'mode\t%s\n' "$MODE"
  printf 'run_id\t%s\n' "$RUN_ID"
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'lane_pairs\t%s\n' "${#r1_files[@]}"
  printf 'expected_read_pairs\t%s\n' "$EXPECTED_PAIRS"
  printf 'bam\t%s\n' "$final_bam"
  printf 'validation_json\t%s\n' "$validation_json"
  printf 'source_tag_rename\tXE_to_ZE\n'
  printf 'global_duplicate_marking\tdisabled_multi_cell_CB_aware_required\n'
  printf 'completed_at\t%s\n' "$(date -Is)"
} > "$tmp_done"
mv -f "$tmp_done" "$done_file"

echo "[OK] validated DNA BWA BAM: $final_bam"
