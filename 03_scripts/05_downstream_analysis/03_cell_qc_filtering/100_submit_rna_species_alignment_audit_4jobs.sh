#!/usr/bin/env bash
# Step 100: render and optionally submit four RNA species audit jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
SCRIPT_ROOT="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}"
MANIFEST="$PROJECT/08_qc_summary/14_postmapping_full_inputs/rna_reverse_genefull_full_manifest.tsv"
THRESHOLD_SUMMARY="$PROJECT/08_qc_summary/30_cell_qc_threshold_audit/cell_qc_threshold_audit_summary.tsv"
TEMPLATE="$SCRIPT_ROOT/05_downstream_analysis/03_cell_qc_filtering/99_submit_rna_species_alignment_audit_template_32core.csub"
RENDERED="$PROJECT/09_logs/05_downstream_analysis/rna_species_residual_audit/jobs"
CONFIRM="${CONFIRM_SUBMIT_RNA_SPECIES_AUDIT:-0}"

[[ "$CONFIRM" =~ ^[01]$ ]] || exit 2
for path in "$MANIFEST" "$THRESHOLD_SUMMARY" "$TEMPLATE"; do
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
[[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$THRESHOLD_SUMMARY")" \
  == "validated" ]] || {
  echo "[ERROR] cell-QC threshold audit is not validated" >&2
  exit 2
}
command -v csub >/dev/null 2>&1 || {
  echo "[ERROR] csub is unavailable" >&2
  exit 127
}

mkdir -p "$RENDERED" \
  "$PROJECT/09_logs/05_downstream_analysis/rna_species_residual_audit"
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
  primary="$(awk -F '\t' -v wanted="$shard" \
    'NR > 1 && $2 == wanted {sum += $7} END {print sum+0}' "$MANIFEST")"
  echo "[PLAN] shard=$shard files=$files primary_mapped=$primary job=$job"
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted $submitted RNA species alignment audit jobs"
else
  echo "[DRY RUN] submitted 0; set CONFIRM_SUBMIT_RNA_SPECIES_AUDIT=1"
fi
