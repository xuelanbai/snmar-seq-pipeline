#!/usr/bin/env bash
# Step 37: render and optionally submit four independent DNA Probe jobs.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="$PROJECT/03_scripts/04_postmapping"
TEMPLATE="$POST_DIR/36_submit_dna_cb_markdup_probe_template_64core.csub"
MANIFEST="$PROJECT/08_qc_summary/13_postmapping_probe/dna_cb_markdup_probe_manifest.tsv"
RENDERED="$PROJECT/09_logs/04_postmapping/dna_cb_markdup_probe/jobs"
CONFIRM="${CONFIRM_SUBMIT_DNA_CB_MARKDUP:-0}"

[[ "$CONFIRM" =~ ^[01]$ ]] || exit 2
[[ -s "$TEMPLATE" && -s "$MANIFEST" ]] || exit 2
rows=$(( $(wc -l < "$MANIFEST") - 1 ))
(( rows == 4 )) || {
  echo "[ERROR] expected four DNA Probe rows, found $rows" >&2
  exit 2
}
command -v csub >/dev/null 2>&1 || {
  echo "[ERROR] csub is unavailable" >&2
  exit 127
}

mkdir -p "$RENDERED" \
  "$PROJECT/09_logs/04_postmapping/dna_cb_markdup_probe"
submitted=0
for index in 1 2 3 4; do
  sample="$(awk -F '\t' -v wanted="$index" \
    'NR > 1 && $1 == wanted {print $3; exit}' "$MANIFEST")"
  [[ -n "$sample" ]] || exit 2
  job="$RENDERED/probe_${index}.csub"
  sed "s/__PROBE_INDEX__/$index/g" "$TEMPLATE" > "$job"
  grep -q '__PROBE_INDEX__' "$job" && {
    echo "[ERROR] unresolved template token: $job" >&2
    exit 2
  }
  echo "[PLAN] index=$index sample=$sample job=$job"
  if (( CONFIRM == 1 )); then
    csub < "$job"
    submitted=$((submitted + 1))
  fi
done

if (( CONFIRM == 1 )); then
  echo "[OK] submitted $submitted independent DNA CB-aware markdup jobs"
else
  echo "[DRY RUN] submitted 0; set CONFIRM_SUBMIT_DNA_CB_MARKDUP=1"
fi

