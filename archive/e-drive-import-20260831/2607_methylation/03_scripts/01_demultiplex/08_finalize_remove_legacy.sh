#!/usr/bin/env bash
# Step 08: verify the new installation, then remove only known legacy filenames.
set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
PROJECT="$(cd "$PROJECT" && pwd)"
SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED_SCRIPT_ROOT="$PROJECT/03_scripts/01_demultiplex"

if [[ "$SCRIPT_ROOT" != "$EXPECTED_SCRIPT_ROOT" ]]; then
  echo "[ERROR] Refusing cleanup outside $EXPECTED_SCRIPT_ROOT" >&2
  echo "[ERROR] Current script directory: $SCRIPT_ROOT" >&2
  exit 2
fi

python3 "$SCRIPT_ROOT/07_verify_install.py"
python3 "$PROJECT/03_scripts/02_rna_demultiplex/07_verify_install.py"

for required_mapping_file in \
  "$PROJECT/03_scripts/03_mapping/00_README.md" \
  "$PROJECT/03_scripts/03_mapping/01_stream_strict_fastq.py" \
  "$PROJECT/03_scripts/03_mapping/02_discover_dna_assigned_pairs.py" \
  "$PROJECT/03_scripts/03_mapping/03_prepare_dna_xd0_front_trim.py" \
  "$PROJECT/03_scripts/03_mapping/04_summarize_dna_fastp_probe.py" \
  "$PROJECT/03_scripts/03_mapping/05_submit_dna_fastp_probe_64core.csub" \
  "$PROJECT/03_scripts/03_mapping/06_summarize_dna_fastp_full.py" \
  "$PROJECT/03_scripts/03_mapping/07_submit_dna_fastp_full_64core.csub"; do
  if [[ ! -s "$required_mapping_file" ]]; then
    echo "[ERROR] Missing new mapping file: $required_mapping_file" >&2
    exit 2
  fi
done

legacy_files=(
  "$SCRIPT_ROOT/README.md"
  "$SCRIPT_ROOT/T7_VALIDATION.txt"
  "$SCRIPT_ROOT/config.sh"
  "$SCRIPT_ROOT/setup_project_dirs.sh"
  "$SCRIPT_ROOT/discover_pairs.py"
  "$SCRIPT_ROOT/snmar_barcode_qc.py"
  "$SCRIPT_ROOT/merge_stats.py"
  "$SCRIPT_ROOT/run_step1.sh"
  "$SCRIPT_ROOT/submit_01_probe.32core.csub"
  "$SCRIPT_ROOT/submit_02_full_stats.64core.csub"
  "$SCRIPT_ROOT/submit_03_split_exact.32core.csub"
  "$PROJECT/03_scripts/02_mapping/MAPPING_RESOURCE_PLAN.md"
  "$PROJECT/03_scripts/02_mapping/00_MAPPING_RESOURCE_PLAN.md"
  "$PROJECT/03_scripts/02_rna_demultiplex/13_probe_rna_t5_polyt_offsets.py"
  "$PROJECT/03_scripts/02_rna_demultiplex/14_submit_t5_polyt_offsets_32core.csub"
  "$PROJECT/snmar_2607_project_code_v2.zip"
  "$PROJECT/snmar_2607_project_code_v3.zip"
  "$PROJECT/snmar_2607_project_code_v4.zip"
  "$PROJECT/snmar_dna_step1_bundle.zip"
)

removed=0
for legacy_file in "${legacy_files[@]}"; do
  if [[ -e "$legacy_file" ]]; then
    rm -f -- "$legacy_file"
    echo "[REMOVED] $legacy_file"
    removed=$((removed + 1))
  fi
done

legacy_mapping_dir="$PROJECT/03_scripts/02_mapping"
if [[ -d "$legacy_mapping_dir" ]]; then
  if find "$legacy_mapping_dir" -mindepth 1 -print -quit | grep -q .; then
    echo "[ERROR] Legacy mapping directory contains an unknown file; not deleting:" >&2
    find "$legacy_mapping_dir" -mindepth 1 -maxdepth 2 -print >&2
    exit 3
  fi
  rmdir -- "$legacy_mapping_dir"
  echo "[REMOVED] $legacy_mapping_dir"
  removed=$((removed + 1))
fi

legacy_cache="$SCRIPT_ROOT/__pycache__"
if [[ -d "$legacy_cache" ]]; then
  rm -rf -- "$legacy_cache"
  echo "[REMOVED] $legacy_cache"
  removed=$((removed + 1))
fi

for legacy_file in "${legacy_files[@]}"; do
  if [[ -e "$legacy_file" ]]; then
    echo "[ERROR] Legacy file still exists: $legacy_file" >&2
    exit 3
  fi
done

echo "[CLEAN_OK] Removed $removed legacy item(s); DNA, RNA and mapping code use the current numbering."
