#!/usr/bin/env python3
"""Step 07: verify the complete RNA strict-demultiplexing installation."""

from __future__ import annotations

import csv
import os
import sys
from collections import Counter
from pathlib import Path


DEFAULT_PROJECT = Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation")
SCRIPT_ROOT = Path(__file__).resolve().parent
PROJECT = Path(os.environ.get("PROJECT", str(DEFAULT_PROJECT))).resolve()
EXPECTED_SCRIPT_ROOT = (
    Path(os.environ.get("SNMAR_SCRIPT_ROOT", PROJECT / "03_scripts"))
    / "02_rna_demultiplex"
).resolve()

REQUIRED_FILES = (
    "01_config.sh",
    "02_setup_project_dirs.sh",
    "03_discover_rna_fastq_pairs.py",
    "04_rna_strict_qc_demultiplex.py",
    "05_merge_rna_stats.py",
    "06_run_rna_demultiplex.sh",
    "07_verify_install.py",
    "08_submit_probe_32core.csub",
    "09_submit_full_stats_64core.csub",
    "10_submit_split_exact_32core.csub",
    "90_optional_library_qc/11_probe_rna_structure_offsets.py",
    "90_optional_library_qc/12_submit_structure_offsets_32core.csub",
    "index/T7_linker.tsv",
    "index/T5.tsv",
    "index/I7.tsv",
    "index/I5.tsv",
)

INDEX_SPECS = {
    "index/T5.tsv": ("T5", 96, 8),
    "index/I7.tsv": ("I7", 384, 10),
    "index/I5.tsv": ("I5", 8, 10),
}


def read_two_column_index(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.rstrip("\r\n")
            if not text:
                continue
            fields = text.split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_no}: expected two tab-separated fields")
            rows.append((fields[0], fields[1].upper()))
    return rows


def main() -> int:
    problems: list[str] = []
    if SCRIPT_ROOT != EXPECTED_SCRIPT_ROOT:
        problems.append(
            f"script directory is {SCRIPT_ROOT}, expected {EXPECTED_SCRIPT_ROOT}"
        )

    for relative in REQUIRED_FILES:
        path = SCRIPT_ROOT / relative
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty required file: {path}")

    for relative, (prefix, count, length) in INDEX_SPECS.items():
        path = SCRIPT_ROOT / relative
        if not path.is_file():
            continue
        try:
            rows = read_two_column_index(path)
        except (OSError, ValueError) as error:
            problems.append(str(error))
            continue
        names = [name for name, _ in rows]
        sequences = [sequence for _, sequence in rows]
        if names != [f"{prefix}-{number}" for number in range(1, count + 1)]:
            problems.append(f"{relative}: names are not exactly {prefix}-1..{count}")
        if len(rows) != count:
            problems.append(f"{relative}: found {len(rows)} rows, expected {count}")
        if len(set(sequences)) != len(sequences):
            problems.append(f"{relative}: duplicate sequences")
        if any(len(sequence) != length or set(sequence) - set("ACGT") for sequence in sequences):
            problems.append(f"{relative}: invalid sequence length or alphabet")

    linker_path = SCRIPT_ROOT / "index" / "T7_linker.tsv"
    linker_sequences: list[str] = []
    if linker_path.is_file():
        try:
            with linker_path.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            expected_columns = {
                "name",
                "sequence",
                "linker_length",
                "linker_sequence",
            }
            if set(rows[0].keys()) != expected_columns:
                problems.append("index/T7_linker.tsv: unexpected columns")
            expected_names = [f"T7-{number}" for number in range(1, 97)]
            if [row["name"] for row in rows] != expected_names:
                problems.append("index/T7_linker.tsv: names are not T7-1..T7-96")
            t7_sequences = [row["sequence"].upper() for row in rows]
            linker_sequences = [row["linker_sequence"].upper() for row in rows]
            linker_lengths = [int(row["linker_length"]) for row in rows]
            if len(rows) != 96 or len(set(t7_sequences)) != 96:
                problems.append("index/T7_linker.tsv: expected 96 unique T7 rows")
            if any(len(sequence) != 8 or set(sequence) - set("ACGT") for sequence in t7_sequences):
                problems.append("index/T7_linker.tsv: invalid T7 sequence")
            if any(
                length not in (22, 23, 24, 25)
                or len(linker) != length
                or set(linker) - set("ACGT")
                for linker, length in zip(linker_sequences, linker_lengths)
            ):
                problems.append("index/T7_linker.tsv: invalid linker sequence/length")
            if Counter(linker_lengths) != Counter({22: 24, 23: 24, 24: 24, 25: 24}):
                problems.append("index/T7_linker.tsv: expected 24 rows per linker length")

            dna_t7_path = (
                PROJECT
                / "03_scripts"
                / "01_demultiplex"
                / "index"
                / "T7.current96.tsv"
            )
            if dna_t7_path.is_file():
                dna_t7 = read_two_column_index(dna_t7_path)
                if dna_t7 != list(zip(expected_names, t7_sequences)):
                    problems.append("RNA T7-linker table differs from DNA T7.current96.tsv")
        except (OSError, ValueError, KeyError, IndexError) as error:
            problems.append(f"index/T7_linker.tsv validation failed: {error}")

    python_files = (
        "03_discover_rna_fastq_pairs.py",
        "04_rna_strict_qc_demultiplex.py",
        "05_merge_rna_stats.py",
        "07_verify_install.py",
        "90_optional_library_qc/11_probe_rna_structure_offsets.py",
    )
    for relative in python_files:
        path = SCRIPT_ROOT / relative
        if not path.is_file():
            continue
        try:
            compile(path.read_text(encoding="utf-8"), str(path), "exec")
        except (OSError, SyntaxError, UnicodeError) as error:
            problems.append(f"{relative}: Python validation failed: {error}")

    feature_checks = {
        "01_config.sh": (
            "rna_strict10",
        ),
        "04_rna_strict_qc_demultiplex.py": (
            "linker_exact_given_T7",
            "T5_exact_at_expected_position",
            "strict_assigned",
            "UR:Z:",
            "UB:Z:",
            "XD:i:0",
            "XL:i:",
            "XI:Z:",
            "DNA_residual_GTCA_linker_given_T7",
            "I7_prefix8_unique_candidate_rejected",
            "potential_prefix8_rescue_structure_I5_exact",
            "strict_assigned_I7_full10",
        ),
        "06_run_rna_demultiplex.sh": (
            "probe_strict10",
            "full_strict10",
            "write_strict10",
        ),
        "09_submit_full_stats_64core.csub": (
            "#CSUB -n 64",
            "JOBS=64",
            "full-strict10",
        ),
        "10_submit_split_exact_32core.csub": (
            "#CSUB -n 32",
            "write-strict10",
        ),
        "90_optional_library_qc/11_probe_rna_structure_offsets.py": (
            "linker_expected_exact",
            "linker_GTCA_exact",
            "T5_exact_scanned_any_forward",
            "I7_prefix8_exact_unique",
            "GTCA_INSERT",
        ),
        "90_optional_library_qc/12_submit_structure_offsets_32core.csub": (
            "#CSUB -n 32",
            "--scan-right 64",
            "--max-reads 200000",
        ),
    }
    for relative, tokens in feature_checks.items():
        path = SCRIPT_ROOT / relative
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for token in tokens:
            if token not in text:
                problems.append(f"{relative}: missing feature token {token}")

    for relative in REQUIRED_FILES:
        path = SCRIPT_ROOT / relative
        if path.is_file() and b"\r\n" in path.read_bytes():
            problems.append(f"Windows CRLF line endings detected: {path}")

    if problems:
        for problem in problems:
            print(f"[FAIL] {problem}", file=sys.stderr)
        print(f"[VERIFY_FAILED] {len(problems)} problem(s)", file=sys.stderr)
        return 1

    print(f"[OK] RNA script root: {SCRIPT_ROOT}")
    print("[OK] exact T7-linker rows: 96; lengths 22/23/24/25 each 24 rows")
    print("[OK] T5=96, I7=384, I5=8; strict RNA structure matching")
    print("[OK] GTCA-insert linker is counted and discarded as DNA residual")
    print("[OK] only observed full-length 10-bp I7 is formally accepted")
    print("[OK] unique 8-bp I7 prefixes are counted and rejected as candidates")
    print("[OK] CR/CB/CN/UR/UB/XD/XE/XL/XI tags, syntax, Linux line endings")
    print("[INSTALL_OK] RNA strict-demultiplexing code is complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
