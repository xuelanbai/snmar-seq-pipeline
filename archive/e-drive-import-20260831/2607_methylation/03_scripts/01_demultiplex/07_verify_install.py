#!/usr/bin/env python3
"""Step 07: verify the numbered installation before legacy cleanup."""

from __future__ import annotations

import os
import sys
from pathlib import Path


DEFAULT_PROJECT = Path(
    "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
)
SCRIPT_ROOT = Path(__file__).resolve().parent
PROJECT = Path(os.environ.get("PROJECT", str(DEFAULT_PROJECT))).resolve()
EXPECTED_SCRIPT_ROOT = (PROJECT / "03_scripts" / "01_demultiplex").resolve()

REQUIRED_FILES = (
    "00_README.md",
    "01_config.sh",
    "02_setup_project_dirs.sh",
    "03_discover_fastq_pairs.py",
    "04_barcode_qc_demultiplex.py",
    "05_merge_demultiplex_stats.py",
    "06_run_demultiplex.sh",
    "07_verify_install.py",
    "08_finalize_remove_legacy.sh",
    "09_submit_probe_32core.csub",
    "10_submit_full_stats_64core.csub",
    "11_submit_split_exact_32core.csub",
    "12_submit_split_rescue_32core.csub",
    "13_cleanup_test_intermediates.sh",
    "index/I5.tsv",
    "index/I7.tsv",
    "index/T5.tsv",
    "index/T7.current96.tsv",
    "reference/00_T7_VALIDATION.txt",
)

INDEX_SPECS = {
    "index/I5.tsv": ("I5", 8, 10),
    "index/I7.tsv": ("I7", 384, 10),
    "index/T5.tsv": ("T5", 96, 8),
    "index/T7.current96.tsv": ("T7", 96, 8),
}


def fail(message: str) -> None:
    print(f"[FAIL] {message}", file=sys.stderr)


def read_index(path: Path) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}: line {line_number} does not have two tab-separated fields"
                )
            rows.append((fields[0], fields[1].upper()))
    return rows


def main() -> int:
    problems: list[str] = []

    if SCRIPT_ROOT != EXPECTED_SCRIPT_ROOT:
        problems.append(
            f"script directory is {SCRIPT_ROOT}, expected {EXPECTED_SCRIPT_ROOT}"
        )

    for required_dir in (PROJECT / "00raw_data", PROJECT / "03_scripts"):
        if not required_dir.is_dir():
            problems.append(f"missing project directory: {required_dir}")

    for relative in REQUIRED_FILES:
        path = SCRIPT_ROOT / relative
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty required file: {path}")

    for relative, (prefix, expected_count, sequence_length) in INDEX_SPECS.items():
        path = SCRIPT_ROOT / relative
        if not path.is_file():
            continue
        try:
            rows = read_index(path)
        except (OSError, ValueError) as error:
            problems.append(str(error))
            continue

        expected_names = [f"{prefix}-{number}" for number in range(1, expected_count + 1)]
        names = [name for name, _ in rows]
        sequences = [sequence for _, sequence in rows]

        if len(rows) != expected_count:
            problems.append(
                f"{relative}: found {len(rows)} rows, expected {expected_count}"
            )
        if names != expected_names:
            problems.append(f"{relative}: barcode names are not exactly {prefix}-1..{expected_count}")
        if len(set(sequences)) != len(sequences):
            problems.append(f"{relative}: duplicate barcode sequences detected")
        invalid = [
            sequence
            for sequence in sequences
            if len(sequence) != sequence_length or set(sequence) - set("ACGT")
        ]
        if invalid:
            problems.append(
                f"{relative}: {len(invalid)} sequence(s) have an invalid length or alphabet"
            )

    python_files = (
        "03_discover_fastq_pairs.py",
        "04_barcode_qc_demultiplex.py",
        "05_merge_demultiplex_stats.py",
        "07_verify_install.py",
    )
    for relative in python_files:
        path = SCRIPT_ROOT / relative
        if not path.is_file():
            continue
        try:
            source = path.read_text(encoding="utf-8")
            compile(source, str(path), "exec")
        except (OSError, SyntaxError, UnicodeError) as error:
            problems.append(f"{relative}: Python validation failed: {error}")

    feature_checks = {
        "04_barcode_qc_demultiplex.py": (
            "CR:Z:",
            "CB:Z:",
            "CN:Z:",
            "XD:i:",
            "XE:Z:",
            "assigned_XD0",
            "assigned_XD1",
            "assigned_XD2",
        ),
        "12_submit_split_rescue_32core.csub": (
            "#CSUB -n 32",
            "JOBS=32",
            "write-rescue",
        ),
        "13_cleanup_test_intermediates.sh": (
            "CONFIRM_CLEANUP",
            "full_stats_summary.tsv",
            "write_rescue outputs and logs were not touched",
        ),
    }
    for relative, required_tokens in feature_checks.items():
        path = SCRIPT_ROOT / relative
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for token in required_tokens:
            if token not in text:
                problems.append(f"{relative}: required feature token is missing: {token}")

    text_files = [
        SCRIPT_ROOT / relative
        for relative in REQUIRED_FILES
        if (SCRIPT_ROOT / relative).is_file()
    ]
    for path in text_files:
        if b"\r\n" in path.read_bytes():
            problems.append(f"Windows CRLF line endings detected: {path}")

    if problems:
        for problem in problems:
            fail(problem)
        print(f"[VERIFY_FAILED] {len(problems)} problem(s)", file=sys.stderr)
        return 1

    print(f"[OK] project root: {PROJECT}")
    print(f"[OK] numbered script root: {SCRIPT_ROOT}")
    print("[OK] index rows: T5=96, T7=96, I7=384, I5=8")
    print("[OK] CR/CB/CN/XD/XE audit tags and rescue submission")
    print("[OK] Python syntax and Linux line endings")
    print("[INSTALL_OK] New numbered code is complete; legacy cleanup is now allowed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
