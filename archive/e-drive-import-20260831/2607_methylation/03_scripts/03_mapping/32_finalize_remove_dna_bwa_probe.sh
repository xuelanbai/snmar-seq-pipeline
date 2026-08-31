#!/usr/bin/env bash
# Step 32: remove DNA BWA probe BAM/log outputs only after full validation.
#
# Dry run:
#   bash 03_scripts/03_mapping/32_finalize_remove_dna_bwa_probe.sh
#
# Confirmed:
#   CONFIRM_REMOVE_DNA_BWA_PROBE=1 \
#     bash 03_scripts/03_mapping/32_finalize_remove_dna_bwa_probe.sh

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
PROBE_ROOT="${DNA_BWA_PROBE_ROOT:-$PROJECT/06_mapping/dna_bwa_joint_probe}"
PROBE_TMP_ROOT="${DNA_BWA_PROBE_TMP_ROOT:-$PROJECT/06_mapping/dna_bwa_tmp_probe}"
FULL_QC="${DNA_BWA_FULL_QC_ROOT:-$PROJECT/08_qc_summary/12_dna_bwa_full}"
FULL_SUMMARY="$FULL_QC/dna_bwa_full_summary.tsv"
FULL_PER_SAMPLE="$FULL_QC/dna_bwa_full_per_sample.tsv"
CONFIRM="${CONFIRM_REMOVE_DNA_BWA_PROBE:-0}"

[[ "$CONFIRM" =~ ^[01]$ ]] || {
  echo "[ERROR] CONFIRM_REMOVE_DNA_BWA_PROBE must be 0 or 1" >&2
  exit 2
}
[[ -s "$FULL_SUMMARY" && -s "$FULL_PER_SAMPLE" ]] || {
  echo "[ERROR] full DNA BWA summary is missing; probe cleanup is blocked" >&2
  exit 2
}
status="$(awk -F '\t' '$1=="status"{print $2; exit}' "$FULL_SUMMARY")"
runs="$(awk -F '\t' '$1=="completed_runs"{print $2; exit}' "$FULL_SUMMARY")"
[[ "$status" == "validated" && "$runs" == "96" ]] || {
  echo "[ERROR] full summary is not validated for all 96 samples" >&2
  exit 2
}

project_resolved="$(cd "$PROJECT" && pwd -P)"
for target in "$PROBE_ROOT" "$PROBE_TMP_ROOT"; do
  if [[ -e "$target" ]]; then
    parent="$(cd "$(dirname "$target")" && pwd -P)"
    resolved="$parent/$(basename "$target")"
    case "$resolved" in
      "$project_resolved"/06_mapping/dna_bwa_joint_probe|\
      "$project_resolved"/06_mapping/dna_bwa_tmp_probe)
        ;;
      *)
        echo "[ERROR] cleanup target is outside the two allowed project paths: $resolved" >&2
        exit 2
        ;;
    esac
    if (( CONFIRM == 0 )); then
      echo "[DRY RUN] would remove: $resolved"
    else
      rm -rf -- "$resolved"
      echo "[REMOVED] $resolved"
    fi
  else
    echo "[INFO] target does not exist: $target"
  fi
done

if (( CONFIRM == 0 )); then
  echo "[DRY RUN] no files were removed"
  echo "[DRY RUN] set CONFIRM_REMOVE_DNA_BWA_PROBE=1 after reviewing the targets"
else
  echo "[OK] DNA probe BAM/log outputs removed; probe QC tables were retained"
fi
