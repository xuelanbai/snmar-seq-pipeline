#!/usr/bin/env bash
# Step 20: safely remove RNA STAR probe outputs after full mapping validates.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
FULL_SUMMARY="${RNA_STAR_FULL_SUMMARY:-$PROJECT/08_qc_summary/09_rna_star_full/rna_star_full_summary.tsv}"
FULL_ROOT="${RNA_STAR_FULL_ROOT:-$PROJECT/06_mapping/rna_star_joint}"
FULL_LOG="$PROJECT/09_logs/03_mapping/rna_star_full_by_file"
CONFIRM_REMOVE_RNA_STAR_PROBE="${CONFIRM_REMOVE_RNA_STAR_PROBE:-0}"

[[ "$CONFIRM_REMOVE_RNA_STAR_PROBE" =~ ^[01]$ ]] || {
  echo "[ERROR] CONFIRM_REMOVE_RNA_STAR_PROBE must be 0 or 1" >&2
  exit 2
}
[[ -s "$FULL_SUMMARY" ]] || {
  echo "[ERROR] missing full RNA STAR summary: $FULL_SUMMARY" >&2
  exit 65
}

summary_value() {
  local key="$1"
  awk -F $'\t' -v key="$key" '$1 == key {print $2; exit}' "$FULL_SUMMARY"
}
completed="$(summary_value completed_R1_files)"
input_reads="$(summary_value STAR_input_R1_reads)"
bam_count=$(
  find "$FULL_ROOT" -type f -name '*.rna_joint.sorted.bam' | wc -l
)
done_count=$(
  find "$FULL_LOG" -maxdepth 1 -type f -name '*.done' | wc -l
)
failed_count=$(
  find "$FULL_LOG" -maxdepth 1 -type f -name '*.failed' | wc -l
)

[[ "$completed" == "96" && "$bam_count" == "96" && \
   "$done_count" == "96" && "$failed_count" == "0" ]] || {
  echo "[ERROR] full validation failed: completed=$completed BAM=$bam_count " \
       "done=$done_count failed=$failed_count" >&2
  exit 65
}
[[ "$input_reads" == "24737553" ]] || {
  echo "[ERROR] full STAR input=$input_reads, expected=24737553" >&2
  exit 65
}

targets=(
  "$PROJECT/06_mapping/rna_star_joint_probe"
  "$PROJECT/09_logs/03_mapping/rna_star_probe_by_file"
)

echo "[INFO] full RNA STAR outputs validated."
echo "[INFO] removable probe targets:"
for target in "${targets[@]}"; do
  case "$target" in
    "$PROJECT"/06_mapping/*|"$PROJECT"/09_logs/03_mapping/*) ;;
    *)
      echo "[ERROR] unsafe probe target: $target" >&2
      exit 65
      ;;
  esac
  if [[ -e "$target" ]]; then
    printf 'REMOVE\t%s\n' "$target"
  else
    printf 'ABSENT\t%s\n' "$target"
  fi
done

if (( CONFIRM_REMOVE_RNA_STAR_PROBE == 0 )); then
  echo "[DRY-RUN] nothing removed."
  echo "[DRY-RUN] set CONFIRM_REMOVE_RNA_STAR_PROBE=1 after review."
  exit 0
fi

for target in "${targets[@]}"; do
  [[ -e "$target" ]] || continue
  rm -rf -- "$target"
  echo "[REMOVED] $target"
done
echo "[OK] RNA STAR probe BAM/log outputs removed; probe QC summary retained."
