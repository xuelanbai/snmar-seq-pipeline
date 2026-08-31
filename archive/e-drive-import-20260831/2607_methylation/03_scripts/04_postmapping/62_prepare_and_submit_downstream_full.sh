#!/usr/bin/env bash
# Step 62: validate inputs and submit direct full DNA/RNA downstream jobs.
#
# Dry-run by default. Confirmation submits all 96 DNA jobs plus RNA filtering.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
CONFIRM="${CONFIRM_SUBMIT_DOWNSTREAM_FULL:-0}"
MIN_EXACT_UMIS="${MIN_EXACT_UMIS:-10}"
MIN_GENES="${MIN_GENES:-5}"

for name in CONFIRM MIN_EXACT_UMIS MIN_GENES; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( CONFIRM <= 1 && MIN_EXACT_UMIS >= 1 && MIN_GENES >= 1 )) || exit 2

source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate "${ALIGN_ENV:-snTAPS_align}"
PYTHON_BIN="$(command -v "${PYTHON:-python}" || true)"
[[ -x "$PYTHON_BIN" ]] || exit 127

cd "$PROJECT"
"$PYTHON_BIN" "$POST_DIR/53_prepare_downstream_full.py" \
  --project "$PROJECT" \
  --rna-min-umis "$MIN_EXACT_UMIS" \
  --rna-min-genes "$MIN_GENES"

echo "[PLAN] DNA: 96 jobs x 32 CPUs; each sample uses 8 x 4-thread contig shards"
echo "[PLAN] RNA: one 8-CPU full filtered-matrix job"

export START_INDEX=1
export END_INDEX=96
if (( CONFIRM == 1 )); then
  export CONFIRM_SUBMIT_DNA_TAPS_FULL=1
  bash "$POST_DIR/57_submit_dna_taps_cpg_full_96jobs.sh"
  MIN_EXACT_UMIS="$MIN_EXACT_UMIS" MIN_GENES="$MIN_GENES" \
    csub < "$POST_DIR/61_submit_rna_filtered_matrix_full_8core.csub"
  echo "[OK] submitted all full DNA TAPS jobs and the full RNA matrix job"
else
  export CONFIRM_SUBMIT_DNA_TAPS_FULL=0
  bash "$POST_DIR/57_submit_dna_taps_cpg_full_96jobs.sh"
  echo "[DRY RUN] csub < 03_scripts/04_postmapping/61_submit_rna_filtered_matrix_full_8core.csub"
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DOWNSTREAM_FULL=1"
fi
