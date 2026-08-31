#!/usr/bin/env bash
# Step 67: preview/submit 768 resumable DNA context-aware cytosine jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
INPUT_ROOT="$PROJECT/08_qc_summary/20_single_cytosine_inputs"
TASKS="$INPUT_ROOT/dna_single_cytosine_tasks.tsv"
TEMPLATE="$POST_DIR/66_submit_dna_single_cytosine_task_template_8core.csub"
JOB_ROOT="$INPUT_ROOT/dna_single_cytosine_job_scripts"
START_TASK="${START_TASK:-1}"
END_TASK="${END_TASK:-768}"
CONFIRM="${CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE:-0}"
FORCE="${FORCE:-0}"

for name in START_TASK END_TASK CONFIRM FORCE; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( START_TASK >= 1 && END_TASK <= 768 && START_TASK <= END_TASK && \
  CONFIRM <= 1 && FORCE <= 1 )) || exit 2
for path in "$TASKS" "$TEMPLATE"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
rows="$(awk -F '\t' 'NR > 1 {count++} END {print count+0}' "$TASKS")"
(( rows == 768 )) || {
  echo "[ERROR] expected 768 task rows, found $rows" >&2
  exit 2
}

mkdir -p "$JOB_ROOT"
submitted=0
skipped=0
for ((index=START_TASK; index<=END_TASK; index++)); do
  row="$(awk -F '\t' -v wanted="$index" \
    'NR > 1 && $1 == wanted {print; exit}' "$TASKS")"
  [[ -n "$row" ]] || exit 2
  IFS=$'\t' read -r observed sample_index sample group shard input_bam \
    reference_fasta contig_file cell_index contig_index output_dir min_mapq \
    min_baseq exclude_read_ends <<< "$row"
  [[ "$observed" == "$index" ]] || exit 2
  prefix="$sample.shard_$(printf '%02d' "$shard")"
  summary="$output_dir/$prefix.summary.tsv"
  if (( FORCE == 0 )) && [[ -s "$summary" ]] && \
    awk -F '\t' -v bam="$input_bam" -v mq="$min_mapq" \
      -v bq="$min_baseq" -v ends="$exclude_read_ends" \
      '$1 == "status" && $2 == "validated" {ok=1}
       $1 == "input_bam" && $2 == bam {same=1}
       $1 == "minimum_mapq" && $2 == mq {mapq=1}
       $1 == "minimum_baseq" && $2 == bq {baseq=1}
       $1 == "excluded_bases_each_read_end" && $2 == ends {trim=1}
       END {exit !(ok && same && mapq && baseq && trim)}' "$summary"; then
    echo "[SKIP] $index $sample shard=$shard"
    skipped=$((skipped + 1))
    continue
  fi
  job="$JOB_ROOT/task_$(printf '%03d' "$index").csub"
  sed \
    -e "s/__TASK_INDEX__/$index/g" \
    -e "s/__FORCE__/$FORCE/g" \
    "$TEMPLATE" > "$job"
  grep -Eq '__TASK_INDEX__|__FORCE__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  else
    echo "[DRY RUN] csub < $job  # $sample shard=$shard"
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted=$submitted skipped_validated=$skipped"
else
  echo "[DRY RUN] no jobs submitted; set "
  echo "CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE=1"
fi
