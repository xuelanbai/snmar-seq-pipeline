#!/usr/bin/env bash
# Step 13: remove only validated Probe/full-stats intermediate directories.
set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
PROJECT="$(cd "$PROJECT" && pwd)"
SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED_SCRIPT_ROOT="$PROJECT/03_scripts/01_demultiplex"
QC_ROOT="$PROJECT/08_qc_summary/01_demultiplex"
LOG_ROOT="$PROJECT/09_logs/01_demultiplex"
SPLIT_ROOT="$PROJECT/04_fastq_split"

if [[ "$SCRIPT_ROOT" != "$EXPECTED_SCRIPT_ROOT" ]]; then
  echo "[ERROR] Refusing cleanup outside $EXPECTED_SCRIPT_ROOT" >&2
  echo "[ERROR] Current script directory: $SCRIPT_ROOT" >&2
  exit 2
fi

required_summaries=(
  "$QC_ROOT/probe_summary.tsv"
  "$QC_ROOT/probe_per_pair.tsv"
  "$QC_ROOT/probe_top_unknown_barcodes.tsv"
  "$QC_ROOT/full_stats_summary.tsv"
  "$QC_ROOT/full_stats_per_pair.tsv"
  "$QC_ROOT/full_stats_top_unknown_barcodes.tsv"
)

for required_file in "${required_summaries[@]}"; do
  if [[ ! -s "$required_file" ]]; then
    echo "[ERROR] Required summary is missing or empty: $required_file" >&2
    exit 3
  fi
done

targets=(
  "$QC_ROOT/probe_by_pair"
  "$QC_ROOT/probe_barcode_counts"
  "$QC_ROOT/probe_unknown"
  "$QC_ROOT/full_stats_by_pair"
  "$QC_ROOT/full_stats_barcode_counts"
  "$QC_ROOT/full_stats_unknown"
  "$LOG_ROOT/probe"
  "$LOG_ROOT/full_stats"
  "$SPLIT_ROOT/probe"
  "$SPLIT_ROOT/full_stats"
)

archives=(
  "$PROJECT/snmar_2607_project_code_v3.zip"
  "$PROJECT/snmar_2607_project_code_v4.zip"
  "$PROJECT/snmar_2607_project_code_v5.zip"
)

echo "[INFO] Required aggregate summaries are present."
echo "[INFO] Cleanup targets:"
for target in "${targets[@]}"; do
  case "$target" in
    "$QC_ROOT"/probe_*|"$QC_ROOT"/full_stats_*|"$LOG_ROOT"/probe|"$LOG_ROOT"/full_stats|"$SPLIT_ROOT"/probe|"$SPLIT_ROOT"/full_stats)
      ;;
    *)
      echo "[ERROR] Unsafe cleanup target: $target" >&2
      exit 4
      ;;
  esac
  if [[ -e "$target" ]]; then
    du -sh -- "$target" 2>/dev/null || true
  else
    echo "[ABSENT] $target"
  fi
done

for archive in "${archives[@]}"; do
  if [[ -e "$archive" ]]; then
    du -h -- "$archive" 2>/dev/null || true
  fi
done

if [[ "${CONFIRM_CLEANUP:-NO}" != "YES" ]]; then
  echo "[DRY_RUN] Nothing was deleted."
  echo "[DRY_RUN] After copying the six summary TSV files locally, run:"
  echo "CONFIRM_CLEANUP=YES bash $SCRIPT_ROOT/13_cleanup_test_intermediates.sh"
  exit 0
fi

for target in "${targets[@]}"; do
  if [[ -e "$target" ]]; then
    rm -rf -- "$target"
    echo "[REMOVED] $target"
  fi
done

for archive in "${archives[@]}"; do
  if [[ -f "$archive" ]]; then
    rm -f -- "$archive"
    echo "[REMOVED] $archive"
  fi
done

for required_file in "${required_summaries[@]}"; do
  if [[ ! -s "$required_file" ]]; then
    echo "[ERROR] Required summary disappeared during cleanup: $required_file" >&2
    exit 5
  fi
done

echo "[CLEANUP_OK] Probe/full-stats aggregate summaries were retained."
echo "[CLEANUP_OK] write_rescue outputs and logs were not touched."
