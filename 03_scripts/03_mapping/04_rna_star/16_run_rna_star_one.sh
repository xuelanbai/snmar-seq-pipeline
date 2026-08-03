#!/usr/bin/env bash
# Step 16 worker: STAR-map one encoded RNA R1 and restore BAM barcode/UMI tags.

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 MODE PAIR_ID SAMPLE CLEAN_R1 EXPECTED_READS" >&2
  exit 2
fi
MODE="$1"
PAIR_ID="$2"
SAMPLE="$3"
CLEAN_R1="$4"
EXPECTED_READS="$5"

for name in MODE PAIR_ID SAMPLE CLEAN_R1 EXPECTED_READS; do
  value="${!name}"
  [[ -n "$value" ]] || {
    echo "[ERROR] empty argument: $name" >&2
    exit 2
  }
done
[[ "$MODE" == "probe" || "$MODE" == "full" ]] || {
  echo "[ERROR] MODE must be probe or full: $MODE" >&2
  exit 2
}
[[ "$PAIR_ID" != *[[:space:]]* && "$SAMPLE" != *[[:space:]]* && \
   "$CLEAN_R1" != *[[:space:]]* ]] || {
  echo "[ERROR] whitespace is not supported in worker arguments" >&2
  exit 2
}
[[ "$EXPECTED_READS" =~ ^[1-9][0-9]*$ ]] || {
  echo "[ERROR] invalid EXPECTED_READS: $EXPECTED_READS" >&2
  exit 2
}

: "${PROJECT:?PROJECT is required}"
: "${MAPPING_DIR:?MAPPING_DIR is required}"
: "${INDEX_DIR:?INDEX_DIR is required}"
: "${OUTPUT_ROOT:?OUTPUT_ROOT is required}"
: "${LOG_ROOT:?LOG_ROOT is required}"
: "${TMP_ROOT:?TMP_ROOT is required}"
: "${STAR_BIN:?STAR_BIN is required}"
: "${SAMTOOLS_BIN:?SAMTOOLS_BIN is required}"
: "${PYTHON_BIN:?PYTHON_BIN is required}"
: "${STAR_VERSION:?STAR_VERSION is required}"
: "${SAMTOOLS_VERSION:?SAMTOOLS_VERSION is required}"

TAG_CODEC="$MAPPING_DIR/04_rna_star/15_rna_fastq_sam_tags.py"
[[ -s "$TAG_CODEC" ]] || {
  echo "[ERROR] missing tag codec: $TAG_CODEC" >&2
  exit 2
}
[[ -s "$CLEAN_R1" ]] || {
  echo "[ERROR] missing/empty clean R1: $CLEAN_R1" >&2
  exit 2
}

STAR_THREADS="${STAR_THREADS:-8}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-4}"
BAM_SORT_THREADS="${BAM_SORT_THREADS:-2}"
BAM_SORT_RAM="${BAM_SORT_RAM:-8000000000}"
GENOME_LOAD="${GENOME_LOAD:-LoadAndKeep}"
FORCE="${FORCE:-0}"

for name in STAR_THREADS SAMTOOLS_THREADS BAM_SORT_THREADS BAM_SORT_RAM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || {
    echo "[ERROR] $name must be a non-negative integer: $value" >&2
    exit 2
  }
done
(( STAR_THREADS >= 1 && SAMTOOLS_THREADS >= 1 && \
   BAM_SORT_THREADS >= 1 && BAM_SORT_RAM >= 1000000000 )) || exit 2
(( FORCE == 0 || FORCE == 1 )) || exit 2

OUT_DIR="$OUTPUT_ROOT/$SAMPLE/$PAIR_ID"
FINAL_BAM="$OUT_DIR/$PAIR_ID.rna_joint.sorted.bam"
FINAL_BAI="$FINAL_BAM.bai"
FINAL_QC="$OUT_DIR/$PAIR_ID.bam_qc.tsv"
FINAL_FLAGSTAT="$OUT_DIR/$PAIR_ID.flagstat.txt"
FINAL_LOG="$LOG_ROOT/$PAIR_ID.log"
DONE_FILE="$LOG_ROOT/$PAIR_ID.done"
FAILED_FILE="$LOG_ROOT/$PAIR_ID.failed"
LOG_TMP="$LOG_ROOT/$PAIR_ID.log.tmp.$$"

mkdir -p "$OUTPUT_ROOT" "$LOG_ROOT" "$TMP_ROOT"

manifest_value() {
  local key="$1"
  awk -F $'\t' -v key="$key" '$1 == key {print $2; exit}' "$DONE_FILE"
}
qc_value() {
  local key="$1"
  awk -F $'\t' -v key="$key" '$1 == key {print $2; exit}' "$FINAL_QC"
}

completed_ok() {
  [[ -s "$DONE_FILE" && -s "$FINAL_BAM" && -s "$FINAL_BAI" && \
     -s "$FINAL_QC" && -s "$OUT_DIR/star.Log.final.out" && \
     -s "$OUT_DIR/star.ReadsPerGene.out.tab" ]] || return 1
  [[ "$(manifest_value status)" == "complete" ]] || return 1
  [[ "$(manifest_value pipeline_signature)" == \
     "snmar_rna_star_joint_v10" ]] || return 1
  [[ "$(manifest_value mode)" == "$MODE" ]] || return 1
  [[ "$(manifest_value pair_id)" == "$PAIR_ID" ]] || return 1
  [[ "$(manifest_value sample)" == "$SAMPLE" ]] || return 1
  [[ "$(manifest_value clean_R1_path)" == "$CLEAN_R1" ]] || return 1
  [[ "$(manifest_value clean_R1_bytes)" == \
     "$(stat -c '%s' "$CLEAN_R1")" ]] || return 1
  [[ "$(manifest_value clean_R1_mtime)" == \
     "$(stat -c '%Y' "$CLEAN_R1")" ]] || return 1
  [[ "$(manifest_value expected_R1_reads)" == "$EXPECTED_READS" ]] || return 1
  [[ "$(manifest_value star_version)" == "$STAR_VERSION" ]] || return 1
  [[ "$(manifest_value samtools_version)" == "$SAMTOOLS_VERSION" ]] || return 1
  [[ "$(manifest_value bam_bytes)" == \
     "$(stat -c '%s' "$FINAL_BAM")" ]] || return 1
  [[ "$(qc_value pipeline_signature)" == \
     "snmar_rna_bam_validate_v10" ]] || return 1
  [[ "$(qc_value primary_records)" == "$EXPECTED_READS" ]] || return 1
  "$SAMTOOLS_BIN" quickcheck -v "$FINAL_BAM" >/dev/null 2>&1 || return 1
}

if (( FORCE == 0 )) && completed_ok; then
  echo "[SKIP] validated $MODE RNA STAR output: $PAIR_ID"
  exit 0
fi

rm -f -- "$FAILED_FILE"
exec >"$LOG_TMP" 2>&1

TMP_JOB=""
SUCCESS=0
finish() {
  rc=$?
  if [[ -n "$TMP_JOB" && -d "$TMP_JOB" ]]; then
    case "$TMP_JOB" in
      "$TMP_ROOT"/*) rm -rf -- "$TMP_JOB" ;;
      *) echo "[WARN] not removing unsafe temp path: $TMP_JOB" >&2 ;;
    esac
  fi
  if (( SUCCESS == 1 && rc == 0 )); then
    rm -f -- "$FAILED_FILE"
  else
    {
      printf 'status\tfailed\n'
      printf 'pipeline_signature\tsnmar_rna_star_joint_v10\n'
      printf 'mode\t%s\n' "$MODE"
      printf 'pair_id\t%s\n' "$PAIR_ID"
      printf 'sample\t%s\n' "$SAMPLE"
      printf 'exit_code\t%s\n' "$rc"
      printf 'failed_utc\t%s\n' "$(date -u +%FT%TZ)"
    } >"$FAILED_FILE.tmp.$$"
    mv -f -- "$FAILED_FILE.tmp.$$" "$FAILED_FILE"
  fi
  if [[ -f "$LOG_TMP" ]]; then
    mv -f -- "$LOG_TMP" "$FINAL_LOG"
  fi
  exit "$rc"
}
trap finish EXIT

echo "[INFO] mode=$MODE pair_id=$PAIR_ID sample=$SAMPLE"
echo "[INFO] clean_R1=$CLEAN_R1"
echo "[INFO] expected_R1_reads=$EXPECTED_READS"
echo "[INFO] STAR=$STAR_VERSION"
echo "[INFO] samtools=$SAMTOOLS_VERSION"
echo "[INFO] genome_load=$GENOME_LOAD"

TMP_JOB="$(mktemp -d "$TMP_ROOT/${PAIR_ID}.XXXXXX")"
ENCODED_FASTQ="$TMP_JOB/$PAIR_ID.tags_in_qname.fastq"
ENCODE_STATS="$TMP_JOB/$PAIR_ID.encode_stats.tsv"
RESTORE_STATS="$TMP_JOB/$PAIR_ID.restore_stats.tsv"
TMP_PREFIX="$TMP_JOB/star."
ENCODED_BAM="$TMP_JOB/star.Aligned.sortedByCoord.out.bam"
RESTORED_BAM="$TMP_JOB/$PAIR_ID.rna_joint.sorted.bam"
RESTORED_BAI="$RESTORED_BAM.bai"
BAM_QC="$TMP_JOB/$PAIR_ID.bam_qc.tsv"
FLAGSTAT="$TMP_JOB/$PAIR_ID.flagstat.txt"

"$PYTHON_BIN" "$TAG_CODEC" encode-fastq \
  --input "$CLEAN_R1" \
  --output "$ENCODED_FASTQ" \
  --stats "$ENCODE_STATS"
encoded_reads=$(
  awk -F $'\t' '$1 == "encoded_reads" {print $2; exit}' "$ENCODE_STATS"
)
[[ "$encoded_reads" == "$EXPECTED_READS" ]] || {
  echo "[ERROR] encoded reads=$encoded_reads, expected=$EXPECTED_READS" >&2
  exit 65
}

"$STAR_BIN" \
  --runMode alignReads \
  --runThreadN "$STAR_THREADS" \
  --genomeDir "$INDEX_DIR" \
  --genomeLoad "$GENOME_LOAD" \
  --readFilesIn "$ENCODED_FASTQ" \
  --outFileNamePrefix "$TMP_PREFIX" \
  --outTmpDir "$TMP_JOB/STARtmp" \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortingThreadN "$BAM_SORT_THREADS" \
  --limitBAMsortRAM "$BAM_SORT_RAM" \
  --outSAMunmapped Within \
  --outSAMattributes NH HI AS nM NM MD jM jI \
  --outSAMattrRGline \
    "ID:$PAIR_ID" "SM:$SAMPLE" "LB:snMAR_RNA" "PL:ILLUMINA" \
  --outSAMmapqUnique 60 \
  --outSAMprimaryFlag OneBestScore \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --quantMode GeneCounts

[[ -s "$ENCODED_BAM" ]] || {
  echo "[ERROR] STAR did not produce sorted BAM" >&2
  exit 65
}
for path in \
  "$TMP_JOB/star.Log.final.out" \
  "$TMP_JOB/star.ReadsPerGene.out.tab"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty STAR result: $path" >&2
    exit 65
  }
done
[[ -e "$TMP_JOB/star.SJ.out.tab" ]] || {
  echo "[ERROR] STAR did not create SJ.out.tab" >&2
  exit 65
}

"$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" -h "$ENCODED_BAM" |
  "$PYTHON_BIN" "$TAG_CODEC" restore-sam \
    --input - \
    --output - \
    --stats "$RESTORE_STATS" |
  "$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" -b \
    -o "$RESTORED_BAM" -

"$SAMTOOLS_BIN" quickcheck -v "$RESTORED_BAM"
"$SAMTOOLS_BIN" index -@ "$SAMTOOLS_THREADS" "$RESTORED_BAM"
"$SAMTOOLS_BIN" flagstat -@ "$SAMTOOLS_THREADS" \
  "$RESTORED_BAM" >"$FLAGSTAT"
"$SAMTOOLS_BIN" view -@ "$SAMTOOLS_THREADS" -h "$RESTORED_BAM" |
  "$PYTHON_BIN" "$TAG_CODEC" validate-sam \
    --input - \
    --expected-primary-reads "$EXPECTED_READS" \
    --check-primary-tags 1000 \
    --stats "$BAM_QC"

[[ -s "$RESTORED_BAM" && -s "$RESTORED_BAI" && \
   -s "$BAM_QC" && -s "$FLAGSTAT" ]] || {
  echo "[ERROR] missing restored BAM validation output" >&2
  exit 65
}

mkdir -p "$OUT_DIR"
mv -f -- "$RESTORED_BAM" "$FINAL_BAM"
mv -f -- "$RESTORED_BAI" "$FINAL_BAI"
mv -f -- "$BAM_QC" "$FINAL_QC"
mv -f -- "$FLAGSTAT" "$FINAL_FLAGSTAT"
for suffix in Log.final.out Log.out Log.progress.out ReadsPerGene.out.tab \
  SJ.out.tab; do
  source_path="$TMP_JOB/star.$suffix"
  if [[ -e "$source_path" ]]; then
    mv -f -- "$source_path" "$OUT_DIR/star.$suffix"
  fi
done

{
  printf 'status\tcomplete\n'
  printf 'pipeline_signature\tsnmar_rna_star_joint_v10\n'
  printf 'mode\t%s\n' "$MODE"
  printf 'pair_id\t%s\n' "$PAIR_ID"
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'clean_R1_path\t%s\n' "$CLEAN_R1"
  printf 'clean_R1_bytes\t%s\n' "$(stat -c '%s' "$CLEAN_R1")"
  printf 'clean_R1_mtime\t%s\n' "$(stat -c '%Y' "$CLEAN_R1")"
  printf 'expected_R1_reads\t%s\n' "$EXPECTED_READS"
  printf 'star_version\t%s\n' "$STAR_VERSION"
  printf 'samtools_version\t%s\n' "$SAMTOOLS_VERSION"
  printf 'star_threads\t%s\n' "$STAR_THREADS"
  printf 'genome_load\t%s\n' "$GENOME_LOAD"
  printf 'bam_path\t%s\n' "$FINAL_BAM"
  printf 'bam_bytes\t%s\n' "$(stat -c '%s' "$FINAL_BAM")"
  printf 'bam_index_bytes\t%s\n' "$(stat -c '%s' "$FINAL_BAI")"
  printf 'primary_records\t%s\n' "$(qc_value primary_records)"
  printf 'primary_mapped_records\t%s\n' \
    "$(qc_value primary_mapped_records)"
  printf 'primary_human_records\t%s\n' \
    "$(qc_value primary_human_records)"
  printf 'primary_pig_records\t%s\n' \
    "$(qc_value primary_pig_records)"
  printf 'checked_primary_tag_records\t%s\n' \
    "$(qc_value checked_primary_tag_records)"
  printf 'completed_utc\t%s\n' "$(date -u +%FT%TZ)"
} >"$DONE_FILE.tmp.$$"
mv -f -- "$DONE_FILE.tmp.$$" "$DONE_FILE"

SUCCESS=1
echo "[OK] completed $MODE RNA STAR mapping: $PAIR_ID"
