#!/usr/bin/env bash
# Install-time cleanup: remove only explicitly obsolete matrix/UMAP code.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
CONFIRM="${CONFIRM_REMOVE_OBSOLETE_MATRIX_CODE:-0}"
[[ "$CONFIRM" =~ ^[01]$ ]] || exit 2

required=(
  63_prepare_single_cytosine_full.py
  64_submit_single_cytosine_prepare_8core.csub
  65_call_dna_single_cytosine_task.py
  66_submit_dna_single_cytosine_task_template_8core.csub
  67_submit_dna_single_cytosine_768jobs.sh
  68_validate_dna_single_cytosine_sample.py
  69_submit_dna_single_cytosine_validation_template_8core.csub
  70_submit_dna_single_cytosine_validation_96jobs.sh
  71_aggregate_single_cytosine_track_shard.py
  72_submit_single_cytosine_track_template_8core.csub
  73_submit_single_cytosine_track_8jobs.sh
  74_summarize_dna_single_cytosine_full.py
  75_submit_dna_single_cytosine_summary_8core.csub
  76_validate_rna_complete_matrix.py
  77_submit_rna_complete_matrix_8core.csub
  78_summarize_complete_modalities.py
  79_submit_complete_modalities_summary_4core.csub
  80_check_complete_matrix_progress.py
)
for name in "${required[@]}"; do
  [[ -s "$POST_DIR/$name" ]] || {
    echo "[ERROR] v19 prerequisite missing/empty: $POST_DIR/$name" >&2
    exit 2
  }
done

obsolete=(
  63_prepare_complete_matrices.py
  64_build_dna_100kb_complete_matrix.py
  65_validate_rna_complete_matrix.py
  66_submit_prepare_complete_matrices_8core.csub
  67_submit_dna_100kb_complete_matrix_32core.csub
  68_submit_rna_complete_matrix_validation_8core.csub
  69_submit_complete_matrices.sh
  70_summarize_complete_matrices.py
  71_submit_complete_matrix_summary_4core.csub
  72_prepare_dna_100kb_umap_inputs.py
  73_build_dna_100kb_residual_matrix.py
  74_run_dna_100kb_umap.py
  75_submit_dna_100kb_umap_prepare_8core.csub
  76_submit_dna_100kb_umap_template_32core.csub
  77_submit_dna_100kb_umap_3jobs.sh
  78_summarize_dna_100kb_umap.py
  79_submit_dna_100kb_umap_summary_4core.csub
  63_prepare_single_cpg_full.py
  64_submit_single_cpg_prepare_8core.csub
  65_call_dna_single_cpg_task.py
  66_submit_dna_single_cpg_task_template_8core.csub
  67_submit_dna_single_cpg_768jobs.sh
  68_validate_dna_single_cpg_sample.py
  69_submit_dna_single_cpg_validation_template_8core.csub
  70_submit_dna_single_cpg_validation_96jobs.sh
  71_aggregate_single_cpg_track_shard.py
  72_submit_single_cpg_track_template_8core.csub
  73_submit_single_cpg_track_8jobs.sh
  74_summarize_dna_single_cpg_full.py
  75_submit_dna_single_cpg_summary_8core.csub
)

found=0
for name in "${obsolete[@]}"; do
  path="$POST_DIR/$name"
  if [[ -e "$path" ]]; then
    found=$((found + 1))
    if (( CONFIRM == 1 )); then
      rm -f -- "$path"
      echo "[REMOVED] $path"
    else
      echo "[DRY RUN] remove $path"
    fi
  fi
done
for archive in \
  "$PROJECT/snmar_2607_project_code_v17.zip" \
  "$PROJECT/snmar_2607_project_code_v18.zip"; do
  if [[ -e "$archive" ]]; then
    found=$((found + 1))
    if (( CONFIRM == 1 )); then
      rm -f -- "$archive"
      echo "[REMOVED] $archive"
    else
      echo "[DRY RUN] remove $archive"
    fi
  fi
done

for output in \
  "$PROJECT/07_postmapping/complete_matrices/dna_100kb" \
  "$PROJECT/08_qc_summary/20_complete_matrices/dna_100kb" \
  "$PROJECT/07_postmapping/complete_matrices/dna_single_cpg" \
  "$PROJECT/08_qc_summary/20_single_cpg_inputs" \
  "$PROJECT/08_qc_summary/21_dna_single_cpg_validation" \
  "$PROJECT/08_qc_summary/22_dna_single_cpg_tracks" \
  "$PROJECT/08_qc_summary/23_dna_single_cpg_complete"; do
  if [[ -e "$output" ]]; then
    echo "[NOTICE] obsolete result directory retained: $output"
  fi
done

if (( CONFIRM == 0 )); then
  echo "[DRY RUN] candidates=$found; set "
  echo "CONFIRM_REMOVE_OBSOLETE_MATRIX_CODE=1"
else
  echo "[OK] removed obsolete code/archive candidates=$found"
fi
