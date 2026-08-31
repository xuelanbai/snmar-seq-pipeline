#!/usr/bin/env bash
# Step 106: preview or submit the D20 and D50 DNA UMAP jobs together.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/04_dna_umap"
CONFIRM="${CONFIRM_SUBMIT_DNA_UMAP:-0}"

jobs=(
  "$ROOT/104_submit_dna_100kb_umap_D20_32core.csub"
  "$ROOT/105_submit_dna_100kb_umap_D50_32core.csub"
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
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_UMAP=1"
  exit 0
fi

for job in "${jobs[@]}"; do
  csub < "$job"
done
echo "[OK] submitted D20 and D50 DNA UMAP jobs"
