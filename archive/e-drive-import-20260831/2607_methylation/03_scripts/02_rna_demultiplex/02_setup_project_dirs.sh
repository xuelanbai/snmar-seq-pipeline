#!/usr/bin/env bash
# Step 02: create RNA output directories and set executable bits.
set -euo pipefail

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_ROOT/01_config.sh"

mkdir -p \
  "$SCRIPT_ROOT" \
  "$SPLIT_ROOT" \
  "$QC_ROOT" \
  "$LOG_ROOT"

chmod +x \
  "$SCRIPT_ROOT/02_setup_project_dirs.sh" \
  "$SCRIPT_ROOT/03_discover_rna_fastq_pairs.py" \
  "$SCRIPT_ROOT/04_rna_strict_qc_demultiplex.py" \
  "$SCRIPT_ROOT/05_merge_rna_stats.py" \
  "$SCRIPT_ROOT/06_run_rna_demultiplex.sh" \
  "$SCRIPT_ROOT/07_verify_install.py" \
  "$SCRIPT_ROOT/08_submit_probe_32core.csub" \
  "$SCRIPT_ROOT/09_submit_full_stats_64core.csub" \
  "$SCRIPT_ROOT/10_submit_split_exact_32core.csub" \
  "$SCRIPT_ROOT/11_probe_rna_structure_offsets.py" \
  "$SCRIPT_ROOT/12_submit_structure_offsets_32core.csub"

echo "[INFO] PROJECT=$PROJECT"
echo "[INFO] RAW_ROOT=$RAW_ROOT"
echo "[INFO] SCRIPT_ROOT=$SCRIPT_ROOT"
echo "[INFO] SPLIT_ROOT=$SPLIT_ROOT"
echo "[INFO] QC_ROOT=$QC_ROOT"
echo "[INFO] LOG_ROOT=$LOG_ROOT"
