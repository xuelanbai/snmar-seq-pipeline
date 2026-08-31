#!/usr/bin/env bash
# Step 132: preview or submit the D20 and D50 TSS UMAP jobs together.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/06_dna_tss_umap"
CONFIRM="${CONFIRM_SUBMIT_DNA_TSS_UMAP:-0}"

jobs=(
  "$ROOT/130_submit_dna_tss_umap_D20_32core.csub"
  "$ROOT/131_submit_dna_tss_umap_D50_32core.csub"
)
for job in "${jobs[@]}"; do
  [[ -s "$job" ]] || {
    echo "[ERROR] missing/empty job: $job" >&2
    exit 2
  }
done

if [[ "$CONFIRM" != "1" ]]; then
  for job in "${jobs[@]}"; do
    echo "[DRY RUN] csub < $job"
  done
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_TSS_UMAP=1"
  exit 0
fi

for job in "${jobs[@]}"; do
  csub < "$job"
done
echo "[OK] submitted D20 and D50 shared-promoter TSS UMAP jobs"
