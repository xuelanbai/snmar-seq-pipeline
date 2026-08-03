#!/usr/bin/env bash
# Step 02: create the project output directories.
set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/01_config.sh"

mkdir -p \
  "${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}/01_demultiplex" \
  "$SPLIT_ROOT" \
  "$QC_ROOT" \
  "$LOG_ROOT"

chmod +x \
  "$SCRIPT_ROOT/02_setup_project_dirs.sh" \
  "$SCRIPT_ROOT/03_discover_fastq_pairs.py" \
  "$SCRIPT_ROOT/04_barcode_qc_demultiplex.py" \
  "$SCRIPT_ROOT/05_merge_demultiplex_stats.py" \
  "$SCRIPT_ROOT/06_run_demultiplex.sh" \
  "$SCRIPT_ROOT/07_verify_install.py" \
  "$SCRIPT_ROOT/08_finalize_remove_legacy.sh" \
  "$SCRIPT_ROOT/12_submit_split_rescue_32core.csub" \
  "$SCRIPT_ROOT/13_cleanup_test_intermediates.sh"

echo "[INFO] PROJECT=$PROJECT"
echo "[INFO] RAW_ROOT=$RAW_ROOT"
echo "[INFO] SCRIPT_ROOT=$SCRIPT_ROOT"
echo "[INFO] SPLIT_ROOT=$SPLIT_ROOT"
echo "[INFO] QC_ROOT=$QC_ROOT"
echo "[INFO] LOG_ROOT=$LOG_ROOT"
echo "[INFO] project directories are ready"
