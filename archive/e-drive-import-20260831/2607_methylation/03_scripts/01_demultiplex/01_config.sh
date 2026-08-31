#!/usr/bin/env bash
# Step 01: shared paths and runtime settings.

# snMAR 2607 project-wide paths.
# Environment variables supplied at submission time take precedence.

PROJECT="${PROJECT:-/share/org/YZWL/yzwl_wujx/Projects/2607_methylation}"
RAW_ROOT="${RAW_ROOT:-${PROJECT}/00raw_data/dna/zhaoqiulin0720095/X101SC26014480-Z01-J026_multipath}"

SCRIPT_ROOT="${SCRIPT_ROOT:-${PROJECT}/03_scripts/01_demultiplex}"
SPLIT_ROOT="${SPLIT_ROOT:-${PROJECT}/04_fastq_split}"
QC_ROOT="${QC_ROOT:-${PROJECT}/08_qc_summary/01_demultiplex}"
LOG_ROOT="${LOG_ROOT:-${PROJECT}/09_logs/01_demultiplex}"

PYTHON="${PYTHON:-python3}"
JOBS="${JOBS:-32}"
PROBE_READS="${PROBE_READS:-200000}"
FORCE="${FORCE:-0}"

# Current batch: only T7-1..T7-96 are valid.
T5_TABLE="${T5_TABLE:-${SCRIPT_ROOT}/index/T5.tsv}"
T7_TABLE="${T7_TABLE:-${SCRIPT_ROOT}/index/T7.current96.tsv}"
I5_TABLE="${I5_TABLE:-${SCRIPT_ROOT}/index/I5.tsv}"
I7_TABLE="${I7_TABLE:-${SCRIPT_ROOT}/index/I7.tsv}"
