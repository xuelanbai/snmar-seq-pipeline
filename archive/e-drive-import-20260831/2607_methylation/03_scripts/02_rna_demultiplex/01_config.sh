#!/usr/bin/env bash
# Step 01: shared paths and runtime settings for RNA demultiplexing.

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
RAW_ROOT="${RAW_ROOT:-${PROJECT}/00raw_data/rna/CellLine-snMAR-RNA}"

SCRIPT_ROOT="${SCRIPT_ROOT:-${PROJECT}/03_scripts/02_rna_demultiplex}"
SPLIT_ROOT="${SPLIT_ROOT:-${PROJECT}/04_fastq_split/rna_strict10}"
QC_ROOT="${QC_ROOT:-${PROJECT}/08_qc_summary/02_rna_demultiplex}"
LOG_ROOT="${LOG_ROOT:-${PROJECT}/09_logs/02_rna_demultiplex}"

PYTHON="${PYTHON:-python3}"
JOBS="${JOBS:-32}"
PROBE_READS="${PROBE_READS:-200000}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-96}"
GZIP_LEVEL="${GZIP_LEVEL:-4}"
FORCE="${FORCE:-0}"

T7_LINKER_TABLE="${T7_LINKER_TABLE:-${SCRIPT_ROOT}/index/T7_linker.tsv}"
T5_TABLE="${T5_TABLE:-${SCRIPT_ROOT}/index/T5.tsv}"
I7_TABLE="${I7_TABLE:-${SCRIPT_ROOT}/index/I7.tsv}"
I5_TABLE="${I5_TABLE:-${SCRIPT_ROOT}/index/I5.tsv}"
