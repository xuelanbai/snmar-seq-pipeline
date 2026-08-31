#!/usr/bin/env bash
# Step 06: orchestrate one demultiplexing mode.
set -euo pipefail

MODE="${1:-probe}"

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_ROOT="$PIPELINE_DIR"
source "$PIPELINE_DIR/01_config.sh"

INDEX_DIR="$PIPELINE_DIR/index"
PAIR_TABLE="$QC_ROOT/pairs.tsv"
UNDETERMINED_TABLE="$QC_ROOT/pairs_undetermined.tsv"

case "$MODE" in
  probe)
    MAX_READS="$PROBE_READS"
    MAX_SEGMENT_MISMATCH=1
    MAX_TOTAL_MISMATCH=2
    WRITE=0
    RUN_NAME="probe"
    ;;
  full-stats)
    MAX_READS=0
    MAX_SEGMENT_MISMATCH=1
    MAX_TOTAL_MISMATCH=2
    WRITE=0
    RUN_NAME="full_stats"
    ;;
  write-exact)
    MAX_READS=0
    MAX_SEGMENT_MISMATCH=0
    MAX_TOTAL_MISMATCH=0
    WRITE=1
    RUN_NAME="write_exact"
    ;;
  write-rescue)
    MAX_READS=0
    MAX_SEGMENT_MISMATCH=1
    MAX_TOTAL_MISMATCH=2
    WRITE=1
    RUN_NAME="write_rescue"
    ;;
  *)
    echo "Usage: bash 06_run_demultiplex.sh {probe|full-stats|write-exact|write-rescue}" >&2
    exit 2
    ;;
esac

STAT_DIR="$QC_ROOT/${RUN_NAME}_by_pair"
BC_DIR="$QC_ROOT/${RUN_NAME}_barcode_counts"
UNKNOWN_DIR="$QC_ROOT/${RUN_NAME}_unknown"
FASTQ_DIR="$SPLIT_ROOT/${RUN_NAME}"
LOG_DIR="$LOG_ROOT/${RUN_NAME}"

mkdir -p "$QC_ROOT" "$STAT_DIR" "$BC_DIR" "$UNKNOWN_DIR" "$FASTQ_DIR" "$LOG_DIR"

"$PYTHON" "$PIPELINE_DIR/03_discover_fastq_pairs.py" \
  --raw-root "$RAW_ROOT" \
  --out "$PAIR_TABLE" \
  --undetermined-out "$UNDETERMINED_TABLE"

run_one() {
  local pair_id="$1"
  local sample="$2"
  local r1="$3"
  local r2="$4"
  local done_file="$LOG_DIR/${pair_id}.done"
  local stat_file="$STAT_DIR/${pair_id}.stat.tsv"
  local bc_file="$BC_DIR/${pair_id}.barcode_count.tsv"
  local unknown_file="$UNKNOWN_DIR/${pair_id}.unknown.json"
  local out1="$FASTQ_DIR/$sample/${pair_id}.assigned_R1.fq.gz"
  local out2="$FASTQ_DIR/$sample/${pair_id}.assigned_R2.fq.gz"

  if [[ "$FORCE" -eq 0 && -s "$done_file" && -s "$stat_file" && -s "$bc_file" && -s "$unknown_file" ]]; then
    if [[ "$WRITE" -eq 0 || ( -s "$out1" && -s "$out2" ) ]]; then
      echo "[INFO] skip completed $pair_id"
      return 0
    fi
  fi

  local cmd=(
    "$PYTHON" "$PIPELINE_DIR/04_barcode_qc_demultiplex.py"
    --pair-id "$pair_id"
    --sample "$sample"
    --r1 "$r1"
    --r2 "$r2"
    --t5 "$T5_TABLE"
    --t7 "$T7_TABLE"
    --i5 "$I5_TABLE"
    --i7 "$I7_TABLE"
    --stat "$stat_file"
    --bc-count "$bc_file"
    --unknown "$unknown_file"
    --max-reads "$MAX_READS"
    --max-segment-mismatch "$MAX_SEGMENT_MISMATCH"
    --max-total-mismatch "$MAX_TOTAL_MISMATCH"
  )
  if [[ "$WRITE" -eq 1 ]]; then
    mkdir -p "$FASTQ_DIR/$sample"
    cmd+=(
      --write
      --out-r1 "$out1"
      --out-r2 "$out2"
    )
  fi
  rm -f "$done_file"
  echo "[INFO] start $pair_id" > "$LOG_DIR/${pair_id}.log"
  "${cmd[@]}" >> "$LOG_DIR/${pair_id}.log" 2>&1
  echo "[INFO] done $pair_id" >> "$LOG_DIR/${pair_id}.log"
  date -Is > "$done_file"
}

export -f run_one
export PYTHON PIPELINE_DIR INDEX_DIR STAT_DIR BC_DIR UNKNOWN_DIR FASTQ_DIR LOG_DIR
export T5_TABLE T7_TABLE I5_TABLE I7_TABLE
export MAX_READS MAX_SEGMENT_MISMATCH MAX_TOTAL_MISMATCH WRITE FORCE

tail -n +2 "$PAIR_TABLE" |
  while IFS=$'\t' read -r pair_id sample r1 r2; do
    printf '%s\0%s\0%s\0%s\0' "$pair_id" "$sample" "$r1" "$r2"
  done |
  xargs -0 -P "$JOBS" -n 4 bash -c 'run_one "$@"' _

"$PYTHON" "$PIPELINE_DIR/05_merge_demultiplex_stats.py" \
  --stat-dir "$STAT_DIR" \
  --unknown-dir "$UNKNOWN_DIR" \
  --summary "$QC_ROOT/${RUN_NAME}_summary.tsv" \
  --per-pair "$QC_ROOT/${RUN_NAME}_per_pair.tsv" \
  --unknown-summary "$QC_ROOT/${RUN_NAME}_top_unknown_barcodes.tsv" \
  --pair-table "$PAIR_TABLE"

echo "[INFO] completed mode=$MODE"
