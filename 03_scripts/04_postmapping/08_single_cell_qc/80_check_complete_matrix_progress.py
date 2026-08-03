#!/usr/bin/env python3
"""Step 80: report resumable DNA/RNA complete-matrix progress."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


VALIDATION_SIGNATURE = "snmar_single_cytosine_validation_20260801"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def is_validated(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    rows = read_tsv(path)
    return any(
        row.get("metric") == "status" and row.get("value") == "validated"
        for row in rows
    )


def is_current_validation(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    rows = read_tsv(path)
    metrics = {row.get("metric"): row.get("value") for row in rows}
    return (
        metrics.get("status") == "validated"
        and metrics.get("pipeline_signature") == VALIDATION_SIGNATURE
    )


def atomic_tsv(
    path: Path, fields: list[str], rows: list[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
    prepare = input_root / "single_cytosine_prepare_summary.tsv"
    tasks = read_tsv(tasks_path) if tasks_path.is_file() else []
    caller_completed = 0
    samples: dict[int, str] = {}
    for task in tasks:
        sample_index = int(task["sample_index"])
        sample = task["sample"]
        samples[sample_index] = sample
        shard = int(task["shard"])
        prefix = f"{sample}.shard_{shard:02d}"
        summary = Path(task["output_dir"]) / f"{prefix}.summary.tsv"
        caller_completed += int(is_validated(summary))
    validation_completed = 0
    for sample_index, sample in samples.items():
        _ = sample_index
        summary = (
            project
            / "08_qc_summary/21_dna_single_cytosine_validation"
            / sample
            / f"{sample}.single_cytosine.validation_summary.tsv"
        )
        validation_completed += int(is_current_validation(summary))
    track_completed = 0
    for shard in range(1, 9):
        summary = (
            project
            / "08_qc_summary/22_dna_single_cytosine_tracks"
            / f"dna_single_cytosine.track_shard_{shard:02d}.summary.tsv"
        )
        track_completed += int(is_validated(summary))
    dna_final = (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete"
        / "dna_single_cytosine_complete_summary.tsv"
    )
    rna_final = (
        project
        / "08_qc_summary/20_complete_matrices/rna_expression"
        / "rna_complete_expression_matrix_summary.tsv"
    )
    modalities_final = (
        project
        / "08_qc_summary/24_complete_modalities"
        / "complete_modalities_summary.tsv"
    )
    bigwig_final = (
        project
        / "08_qc_summary/24_dna_single_cytosine_bigwig"
        / "dna_single_base_bigwig_summary.tsv"
    )
    cell_qc_final = (
        project
        / "08_qc_summary/25_dna_single_cell_qc"
        / "dna_single_cell_coverage_qc_summary.tsv"
    )
    expected_tasks = 768 if tasks else 0
    expected_samples = 96 if tasks else 0
    rows: list[dict[str, object]] = [
        {
            "stage": "63-64_prepare",
            "completed": int(is_validated(prepare)),
            "expected": 1,
            "ready": int(is_validated(prepare)),
        },
        {
            "stage": "65-67_DNA_single_cytosine_callers",
            "completed": caller_completed,
            "expected": expected_tasks,
            "ready": int(
                expected_tasks == 768 and caller_completed == expected_tasks
            ),
        },
        {
            "stage": "68-70_DNA_sample_validation",
            "completed": validation_completed,
            "expected": expected_samples,
            "ready": int(
                expected_samples == 96
                and validation_completed == expected_samples
            ),
        },
        {
            "stage": "71-73_DNA_WCG_GCH_pseudobulk_tracks",
            "completed": track_completed,
            "expected": 8,
            "ready": int(track_completed == 8),
        },
        {
            "stage": "74-75_DNA_complete_VDS",
            "completed": int(is_validated(dna_final)),
            "expected": 1,
            "ready": int(is_validated(dna_final)),
        },
        {
            "stage": "76-77_RNA_complete_matrices",
            "completed": int(is_validated(rna_final)),
            "expected": 1,
            "ready": int(is_validated(rna_final)),
        },
        {
            "stage": "78-79_complete_modalities",
            "completed": int(is_validated(modalities_final)),
            "expected": 1,
            "ready": int(is_validated(modalities_final)),
        },
        {
            "stage": "81-82_DNA_single_base_BigWigs",
            "completed": int(is_validated(bigwig_final)),
            "expected": 1,
            "ready": int(is_validated(bigwig_final)),
        },
        {
            "stage": "83-84_DNA_single_cell_coverage_QC",
            "completed": int(is_validated(cell_qc_final)),
            "expected": 1,
            "ready": int(is_validated(cell_qc_final)),
        },
    ]
    output = args.output or input_root / "complete_matrix_progress.tsv"
    atomic_tsv(output, list(rows[0]), rows)
    for row in rows:
        print(
            f"{row['stage']}: {row['completed']}/{row['expected']} "
            f"ready={row['ready']}"
        )
    print(f"[OK] progress table: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
