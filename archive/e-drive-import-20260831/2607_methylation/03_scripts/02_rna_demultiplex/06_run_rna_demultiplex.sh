#!/usr/bin/env bash
# Step 06: orchestrate one RNA strict-demultiplexing mode.
set -euo pipefail

MODE="${1:-probe}"
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$PIPELINE_DIR/01_config.sh"

PAIR_TABLE="$QC_ROOT/pairs.tsv"

case "$MODE" in
  probe)
    MAX_READS="$PROBE_READS"
    WRITE=0
    RUN_NAME="probe_strict10"
    ;;
  full-strict10)
    MAX_READS=0
    WRITE=0
    RUN_NAME="full_strict10"
    ;;
  write-strict10)
    MAX_READS=0
    WRITE=1
    RUN_NAME="write_strict10"
    ;;
  *)
    echo "Usage: bash 06_run_rna_demultiplex.sh {probe|full-strict10|write-strict10}" >&2
    exit 2
    ;;
esac

STAT_DIR="$QC_ROOT/${RUN_NAME}_by_pair"
BC_DIR="$QC_ROOT/${RUN_NAME}_barcode_counts"
UNKNOWN_DIR="$QC_ROOT/${RUN_NAME}_unknown"
LOG_DIR="$LOG_ROOT/${RUN_NAME}"

mkdir -p "$QC_ROOT" "$STAT_DIR" "$BC_DIR" "$UNKNOWN_DIR" "$LOG_DIR"
if [[ "$WRITE" -eq 1 ]]; then
  mkdir -p "$SPLIT_ROOT"
fi

"$PYTHON" "$PIPELINE_DIR/03_discover_rna_fastq_pairs.py" \
  --raw-root "$RAW_ROOT" \
  --out "$PAIR_TABLE" \
  --expected-samples "$EXPECTED_SAMPLES"

run_one() {
  local pair_id="$1"
  local sample="$2"
  local r1="$3"
  local r2="$4"
  local done_file="$LOG_DIR/${pair_id}.done"
  local stat_file="$STAT_DIR/${pair_id}.stat.tsv"
  local bc_file="$BC_DIR/${pair_id}.barcode_count.tsv"
  local unknown_file="$UNKNOWN_DIR/${pair_id}.unknown.json"
  local out1="$SPLIT_ROOT/$sample/${pair_id}.assigned_R1.fq.gz"
  local out2="$SPLIT_ROOT/$sample/${pair_id}.assigned_R2.fq.gz"

  if [[ "$FORCE" -eq 0 && -s "$done_file" && -s "$stat_file" && \
        -s "$bc_file" && -s "$unknown_file" ]]; then
    if [[ "$WRITE" -eq 0 || ( -s "$out1" && -s "$out2" ) ]]; then
      echo "[INFO] skip completed $pair_id"
      return 0
    fi
  fi

  local cmd=(
    "$PYTHON" "$PIPELINE_DIR/04_rna_strict_qc_demultiplex.py"
    --pair-id "$pair_id"
    --sample "$sample"
    --r1 "$r1"
    --r2 "$r2"
    --t7-linker "$T7_LINKER_TABLE"
    --t5 "$T5_TABLE"
    --i7 "$I7_TABLE"
    --i5 "$I5_TABLE"
    --stat "$stat_file"
    --bc-count "$bc_file"
    --unknown "$unknown_file"
    --max-reads "$MAX_READS"
    --gzip-level "$GZIP_LEVEL"
  )
  if [[ "$WRITE" -eq 1 ]]; then
    mkdir -p "$SPLIT_ROOT/$sample"
    cmd+=(--write --out-r1 "$out1" --out-r2 "$out2")
  fi

  rm -f "$done_file"
  echo "[INFO] start $pair_id" > "$LOG_DIR/${pair_id}.log"
  "${cmd[@]}" >> "$LOG_DIR/${pair_id}.log" 2>&1
  echo "[INFO] done $pair_id" >> "$LOG_DIR/${pair_id}.log"
  date -Is > "$done_file"
}

export -f run_one
export PYTHON PIPELINE_DIR STAT_DIR BC_DIR UNKNOWN_DIR LOG_DIR SPLIT_ROOT
export T7_LINKER_TABLE T5_TABLE I7_TABLE I5_TABLE
export MAX_READS WRITE GZIP_LEVEL FORCE

tail -n +2 "$PAIR_TABLE" |
  while IFS=$'\t' read -r pair_id sample r1 r2; do
    printf '%s\0%s\0%s\0%s\0' "$pair_id" "$sample" "$r1" "$r2"
  done |
  xargs -0 -P "$JOBS" -n 4 bash -c 'run_one "$@"' _

"$PYTHON" "$PIPELINE_DIR/05_merge_rna_stats.py" \
  --stat-dir "$STAT_DIR" \
  --unknown-dir "$UNKNOWN_DIR" \
  --summary "$QC_ROOT/${RUN_NAME}_summary.tsv" \
  --per-pair "$QC_ROOT/${RUN_NAME}_per_pair.tsv" \
  --unknown-summary "$QC_ROOT/${RUN_NAME}_top_unknown_barcodes.tsv" \
  --pair-table "$PAIR_TABLE"

echo "[INFO] completed RNA mode=$MODE"
