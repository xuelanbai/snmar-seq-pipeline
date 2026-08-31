#!/usr/bin/env bash
# Step 115: render and optionally submit four RNA evidence-recount jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
MANIFEST="$PROJECT/08_qc_summary/14_postmapping_full_inputs/rna_reverse_genefull_full_manifest.tsv"
AUDIT_SUMMARY="$PROJECT/08_qc_summary/34_rna_cell_line_audit/rna_read_umi_cellline_audit_summary.tsv"
TEMPLATE="$SCRIPT_ROOT/05_downstream_analysis/05_rna_umap/114_submit_rna_expression_evidence_template_32core.csub"
RENDERED="$PROJECT/09_logs/05_downstream_analysis/rna_expression_evidence/jobs"
CONFIRM="${CONFIRM_SUBMIT_RNA_EVIDENCE:-0}"

[[ "$CONFIRM" =~ ^[01]$ ]] || exit 2
for path in "$MANIFEST" "$AUDIT_SUMMARY" "$TEMPLATE"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done
rows=$(( $(wc -l < "$MANIFEST") - 1 ))
(( rows == 96 )) || {
  echo "[ERROR] expected 96 RNA rows, found $rows" >&2
  exit 2
}
[[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$AUDIT_SUMMARY")" == "validated" ]] || {
  echo "[ERROR] RNA read/UMI audit is not validated" >&2
  exit 2
}
command -v csub >/dev/null 2>&1 || {
  echo "[ERROR] csub is unavailable" >&2
  exit 127
}

mkdir -p "$RENDERED" \
  "$PROJECT/09_logs/05_downstream_analysis/rna_expression_evidence"
submitted=0
for shard in 1 2 3 4; do
  job="$RENDERED/shard_${shard}.csub"
  sed "s/__SHARD_INDEX__/$shard/g" "$TEMPLATE" > "$job"
  grep -q '__SHARD_INDEX__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  files="$(awk -F '\t' -v wanted="$shard" \
    'NR > 1 && $2 == wanted {count++} END {print count+0}' "$MANIFEST")"
  reads="$(awk -F '\t' -v wanted="$shard" \
    'NR > 1 && $2 == wanted {sum += $6} END {print sum+0}' "$MANIFEST")"
  echo "[PLAN] shard=$shard files=$files STAR_input_reads=$reads job=$job"
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted $submitted RNA expression-evidence jobs"
else
  echo "[DRY RUN] submitted 0; set CONFIRM_SUBMIT_RNA_EVIDENCE=1"
fi
