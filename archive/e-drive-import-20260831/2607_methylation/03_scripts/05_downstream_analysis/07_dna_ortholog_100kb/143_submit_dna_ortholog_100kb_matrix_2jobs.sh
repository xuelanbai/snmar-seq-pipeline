#!/usr/bin/env bash
# Step 143: preview or submit D20 and D50 matrix-export jobs.

set -euo pipefail
PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/07_dna_ortholog_100kb"
CONFIRM="${CONFIRM_SUBMIT_DNA_ORTHO100_MATRIX:-0}"
jobs=(
  "$ROOT/141_submit_dna_ortholog_100kb_matrix_D20_32core.csub"
  "$ROOT/142_submit_dna_ortholog_100kb_matrix_D50_32core.csub"
)
for job in "${jobs[@]}"; do [[ -s "$job" ]] || exit 2; done
if [[ "$CONFIRM" != "1" ]]; then
  for job in "${jobs[@]}"; do echo "[DRY RUN] csub < $job"; done
  echo "[DRY RUN] set CONFIRM_SUBMIT_DNA_ORTHO100_MATRIX=1"
  exit 0
fi
for job in "${jobs[@]}"; do csub < "$job"; done
echo "[OK] submitted D20 and D50 shared-orthobin matrix jobs"

