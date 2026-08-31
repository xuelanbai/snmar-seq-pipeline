#!/usr/bin/env bash
# Step 52: prepare manifests and submit DNA/RNA full dedup jobs together.
#
# Dry-run by default. RNA always previews/submits all four balanced shards.
# DNA can be staged with DNA_START_INDEX and DNA_END_INDEX.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
CONFIRM="${CONFIRM_SUBMIT_FULL_DEDUP:-0}"
DNA_START_INDEX="${DNA_START_INDEX:-1}"
DNA_END_INDEX="${DNA_END_INDEX:-24}"

for name in CONFIRM DNA_START_INDEX DNA_END_INDEX; do
  value="${!name}"
  [[ "$value" =~ ^[0-9]+$ ]] || exit 2
done
(( CONFIRM == 0 || CONFIRM == 1 )) || exit 2
(( DNA_START_INDEX >= 1 && DNA_END_INDEX <= 96 && \
  DNA_START_INDEX <= DNA_END_INDEX )) || exit 2

source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate "${ALIGN_ENV:-snTAPS_align}"
PYTHON_BIN="$(command -v "${PYTHON:-python}" || true)"
[[ -x "$PYTHON_BIN" ]] || exit 127

cd "$PROJECT"
"$PYTHON_BIN" "$POST_DIR/42_prepare_postmapping_full_manifests.py" \
  --project "$PROJECT"

if (( CONFIRM == 1 )); then
  export CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=1
  export CONFIRM_SUBMIT_RNA_REVERSE_UMI_FULL=1
else
  export CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=0
  export CONFIRM_SUBMIT_RNA_REVERSE_UMI_FULL=0
fi
export START_INDEX="$DNA_START_INDEX"
export END_INDEX="$DNA_END_INDEX"

bash "$POST_DIR/44_submit_dna_cb_markdup_full_96jobs.sh"
bash "$POST_DIR/48_submit_rna_reverse_genefull_full_4jobs.sh"

if (( CONFIRM == 1 )); then
  echo "[OK] submitted DNA range $DNA_START_INDEX-$DNA_END_INDEX and all RNA shards"
else
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_FULL_DEDUP=1"
fi

