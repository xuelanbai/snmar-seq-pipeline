#!/usr/bin/env bash
# Step 146: preview or submit 2 depth sets x 4 matrix representations.

set -euo pipefail
PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
ROOT="$SCRIPT_ROOT/05_downstream_analysis/07_dna_ortholog_100kb"
TEMPLATE="$ROOT/145_submit_dna_ortholog_100kb_seurat_template_16core.csub"
JOB_ROOT="$PROJECT/08_qc_summary/42_dna_ortholog_100kb_seurat/job_scripts"
CONFIRM="${CONFIRM_SUBMIT_DNA_ORTHO100_SEURAT:-0}"
FORCE="${FORCE:-0}"
REGRESS_DEPTH="${REGRESS_DEPTH:-0}"

[[ "$CONFIRM" =~ ^[01]$ && "$FORCE" =~ ^[01]$ && "$REGRESS_DEPTH" =~ ^[01]$ ]] || exit 2
[[ -s "$TEMPLATE" ]] || exit 2
runs=(
  "Dtotal20k_orthoWCG100_call5_top2000"
  "Dtotal50k_orthoWCG100_call5_top2000"
)
methods=(ratio_paper beta_posterior pearson_residual coverage_control)
mkdir -p "$JOB_ROOT" "$PROJECT/09_logs/05_downstream_analysis/dna_ortholog_100kb"
submitted=0
for run_id in "${runs[@]}"; do
  summary="$PROJECT/08_qc_summary/41_dna_ortholog_100kb_seurat_input/$run_id/dna_ortholog_100kb_seurat_input_summary.tsv"
  [[ -s "$summary" ]] || {
    echo "[ERROR] missing matrix-export summary: $summary" >&2
    exit 2
  }
  depth="${run_id#Dtotal}"
  depth="${depth%%_*}"
  for method in "${methods[@]}"; do
    tag="${depth}_${method}"
    job="$JOB_ROOT/${tag}.csub"
    sed \
      -e "s/__RUN_ID__/$run_id/g" \
      -e "s/__METHOD__/$method/g" \
      -e "s/__JOB_TAG__/$tag/g" \
      -e "s/__FORCE__/$FORCE/g" \
      -e "s/__REGRESS_DEPTH__/$REGRESS_DEPTH/g" \
      "$TEMPLATE" > "$job"
    grep -Eq '__RUN_ID__|__METHOD__|__JOB_TAG__|__FORCE__|__REGRESS_DEPTH__' "$job" && {
      echo "[ERROR] unresolved template token: $job" >&2
      exit 2
    }
    if [[ "$CONFIRM" == "1" ]]; then
      csub < "$job"
      submitted=$((submitted + 1))
    else
      echo "[DRY RUN] csub < $job"
    fi
  done
done
if [[ "$CONFIRM" == "1" ]]; then
  echo "[OK] submitted=$submitted"
else
  echo "[DRY RUN] no jobs submitted; set CONFIRM_SUBMIT_DNA_ORTHO100_SEURAT=1"
fi

