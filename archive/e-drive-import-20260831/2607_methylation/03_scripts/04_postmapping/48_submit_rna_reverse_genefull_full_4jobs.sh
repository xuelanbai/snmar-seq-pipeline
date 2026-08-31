#!/usr/bin/env bash
# Step 48: render and optionally submit four RNA exact-UMI shard jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
MANIFEST="$PROJECT/08_qc_summary/14_postmapping_full_inputs/rna_reverse_genefull_full_manifest.tsv"
PROBE_SUMMARY="$PROJECT/08_qc_summary/13_postmapping_probe/rna_genefull/rna_genefull_probe_summary.tsv"
TEMPLATE="$POST_DIR/47_submit_rna_reverse_genefull_full_template_32core.csub"
RENDERED="$PROJECT/09_logs/04_postmapping/rna_reverse_genefull_full/jobs"
JOB_LOG_ROOT="$PROJECT/09_logs/04_postmapping/rna_reverse_genefull_full"
CONFIRM="${CONFIRM_SUBMIT_RNA_REVERSE_UMI_FULL:-0}"

[[ "$CONFIRM" =~ ^[01]$ ]] || exit 2
for path in "$MANIFEST" "$PROBE_SUMMARY" "$TEMPLATE"; do
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
[[ "$(awk -F '\t' '$1=="status"{print $2; exit}' "$PROBE_SUMMARY")" \
  == "validated" ]] || {
  echo "[ERROR] RNA GeneFull Probe is not validated" >&2
  exit 2
}
command -v csub >/dev/null 2>&1 || {
  echo "[ERROR] csub is unavailable" >&2
  exit 127
}

mkdir -p "$RENDERED" "$JOB_LOG_ROOT"
submitted=0
for shard in 1 2 3 4; do
  job="$RENDERED/shard_${shard}.csub"
  sed "s/__SHARD_INDEX__/$shard/g" "$TEMPLATE" > "$job"
  grep -q '__SHARD_INDEX__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  expected_reads="$(awk -F '\t' -v wanted="$shard" \
    'NR > 1 && $2 == wanted {sum += $6} END {print sum+0}' "$MANIFEST")"
  echo "[PLAN] shard=$shard files=24 input_reads=$expected_reads job=$job"
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted $submitted RNA reverse GeneFull exact-UMI jobs"
else
  echo "[DRY RUN] submitted 0; set CONFIRM_SUBMIT_RNA_REVERSE_UMI_FULL=1"
fi

