#!/usr/bin/env bash
# Step 34: verify software and input capabilities before post-mapping Probe.

set -euo pipefail

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
POST_DIR="${SNMAR_SCRIPT_ROOT:-$PROJECT/03_scripts}/04_postmapping"
QC_ROOT="$PROJECT/08_qc_summary/13_postmapping_probe"
JOINT_GTF="${JOINT_GTF:-$PROJECT/02_reference/rna_human_pig_joint/human_pig.prefixed.gtf}"

SAMTOOLS_BIN="$(command -v "${SAMTOOLS:-samtools}" || true)"
PYTHON_BIN="$(command -v "${PYTHON:-python}" || true)"
for name in SAMTOOLS_BIN PYTHON_BIN; do
  value="${!name}"
  [[ -n "$value" && -x "$value" ]] || {
    echo "[ERROR] unavailable executable: $name" >&2
    exit 127
  }
done

"$PYTHON_BIN" -c 'import csv, gzip, json; print("[OK] Python standard library")'
"$SAMTOOLS_BIN" --version
markdup_help="$("$SAMTOOLS_BIN" markdup --help 2>&1 || true)"
if ! grep -q -- '--barcode-tag' <<< "$markdup_help"; then
  echo "[ERROR] this samtools markdup lacks --barcode-tag" >&2
  exit 2
fi

for path in \
  "$POST_DIR/90_optional_validation/33_prepare_postmapping_probe.py" \
  "$POST_DIR/02_dna_deduplication/35_run_dna_cb_markdup_one.sh" \
  "$POST_DIR/90_optional_validation/38_summarize_dna_cb_markdup_probe.py" \
  "$POST_DIR/90_optional_validation/39_count_rna_genefull_umi_probe.py" \
  "$POST_DIR/90_optional_validation/40_summarize_rna_genefull_umi_probe.py" \
  "$QC_ROOT/dna_cb_markdup_probe_manifest.tsv" \
  "$QC_ROOT/rna_genefull_probe_manifest.tsv" \
  "$JOINT_GTF"; do
  [[ -s "$path" ]] || {
    echo "[ERROR] missing/empty prerequisite: $path" >&2
    exit 2
  }
done

echo "[OK] post-mapping Probe requirements passed"
