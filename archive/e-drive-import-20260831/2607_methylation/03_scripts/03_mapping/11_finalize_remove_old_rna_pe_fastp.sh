#!/usr/bin/env bash
# Step 11: remove obsolete paired-end RNA fastp code/results after v9 passes.
#
# Default: dry-run only.
# Execute:
#   CONFIRM_REMOVE_OLD_RNA_PE=1 \
#     bash 03_scripts/03_mapping/11_finalize_remove_old_rna_pe_fastp.sh
#
# The strict10 split input, including assigned_R2 FASTQ files, is NOT removed.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
CONFIRM_REMOVE_OLD_RNA_PE="${CONFIRM_REMOVE_OLD_RNA_PE:-0}"
NEW_QC="$PROJECT/08_qc_summary/07_rna_r1_fastp_full"
NEW_SUMMARY="$NEW_QC/rna_r1_fastp_full_summary.tsv"
NEW_CLEAN="$PROJECT/05_fastp/rna_r1_strict"
NEW_LOG="$PROJECT/09_logs/03_mapping/rna_r1_fastp_full_by_file"
DEMUX="$PROJECT/08_qc_summary/02_rna_demultiplex/write_strict10_per_pair.tsv"

[[ -s "$NEW_SUMMARY" ]] || {
  echo "[ERROR] missing new R1-only summary: $NEW_SUMMARY" >&2
  exit 2
}
[[ -s "$DEMUX" ]] || {
  echo "[ERROR] missing RNA demultiplex table: $DEMUX" >&2
  exit 2
}

completed=$(
  awk -F $'\t' '$1 == "completed_R1_FASTQ_files" {print $2; exit}' \
    "$NEW_SUMMARY"
)
summary_input=$(
  awk -F $'\t' \
    '$1 == "demultiplex_assigned_XD0_R1_reads" {print $2; exit}' \
    "$NEW_SUMMARY"
)
demux_input=$(
  awk -F $'\t' \
    'NR > 1 && $3 == "assigned_XD0" {sum += $4} END {print sum + 0}' \
    "$DEMUX"
)
clean_count=$(
  find "$NEW_CLEAN" -type f -name '*.clean_R1.fq.gz' | wc -l
)
failed_count=$(
  find "$NEW_LOG" -maxdepth 1 -type f -name '*.failed' | wc -l
)

[[ "$completed" == "96" ]] || {
  echo "[ERROR] completed_R1_FASTQ_files=$completed, expected 96" >&2
  exit 65
}
[[ "$summary_input" == "$demux_input" ]] || {
  echo "[ERROR] new summary input=$summary_input, demux input=$demux_input" >&2
  exit 65
}
[[ "$clean_count" == "96" ]] || {
  echo "[ERROR] clean R1 files=$clean_count, expected 96" >&2
  exit 65
}
[[ "$failed_count" == "0" ]] || {
  echo "[ERROR] new R1-only workflow still has $failed_count failed files" >&2
  exit 65
}

targets=(
  "$PROJECT/03_scripts/03_mapping/08_discover_rna_assigned_pairs.py"
  "$PROJECT/03_scripts/03_mapping/09_prepare_rna_dynamic_front_trim.py"
  "$PROJECT/03_scripts/03_mapping/10_summarize_rna_fastp_probe.py"
  "$PROJECT/03_scripts/03_mapping/11_submit_rna_fastp_probe_64core.csub"
  "$PROJECT/03_scripts/03_mapping/12_summarize_rna_fastp_full.py"
  "$PROJECT/03_scripts/03_mapping/13_submit_rna_fastp_full_64core.csub"
  "$PROJECT/03_scripts/03_mapping/14_finalize_remove_legacy.sh"
  "$PROJECT/05_fastp/rna_probe"
  "$PROJECT/05_fastp/rna_strict"
  "$PROJECT/05_fastp/rna_strict_tmp"
  "$PROJECT/08_qc_summary/05_rna_fastp_probe"
  "$PROJECT/08_qc_summary/06_rna_fastp_full"
  "$PROJECT/09_logs/03_mapping/rna_fastp_probe_by_pair"
  "$PROJECT/09_logs/03_mapping/rna_fastp_full_by_pair"
)

echo "[INFO] v9 R1-only outputs validated."
echo "[INFO] obsolete paired-end RNA fastp targets:"
for target in "${targets[@]}"; do
  case "$target" in
    "$PROJECT"/*) ;;
    *)
      echo "[ERROR] unsafe target outside project: $target" >&2
      exit 65
      ;;
  esac
  if [[ -e "$target" ]]; then
    printf 'REMOVE\t%s\n' "$target"
  else
    printf 'ABSENT\t%s\n' "$target"
  fi
done

if [[ "$CONFIRM_REMOVE_OLD_RNA_PE" != "1" ]]; then
  echo "[DRY-RUN] nothing removed."
  echo "[DRY-RUN] rerun with CONFIRM_REMOVE_OLD_RNA_PE=1 after review."
  exit 0
fi

for target in "${targets[@]}"; do
  [[ -e "$target" ]] || continue
  rm -rf -- "$target"
  echo "[REMOVED] $target"
done

echo "[INFO] obsolete paired-end RNA fastp code/results removed."
echo "[INFO] strict10 split inputs, including R2, were preserved."
