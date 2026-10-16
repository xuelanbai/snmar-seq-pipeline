#!/usr/bin/env python3
"""Step 80: report resumable DNA/RNA complete-matrix progress."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


VALIDATION_SIGNATURE = "snmar_single_cytosine_validation_20260801"
RNA_SPECIES_ALIGNMENT_SIGNATURE = (
    "snmar_rna_species_alignment_audit_20260806"
)


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


def is_validated_single_row(path: Path, signature: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    rows = read_tsv(path)
    return (
        len(rows) == 1
        and rows[0].get("status") == "validated"
        and rows[0].get("pipeline_signature") == signature
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
    wells: dict[int, str] = {}
    for task in tasks:
        sample_index = int(task["sample_index"])
        sample = task["sample"]
        wells[sample_index] = sample
        shard = int(task["shard"])
        prefix = f"{sample}.shard_{shard:02d}"
        summary = Path(task["output_dir"]) / f"{prefix}.summary.tsv"
        caller_completed += int(is_validated(summary))
    validation_completed = 0
    for sample_index, sample in wells.items():
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
    igv_bigwig_final = (
        project
        / "08_qc_summary/26_dna_igv_species_bigwig"
        / "dna_igv_species_bigwig_summary.tsv"
    )
    dna_100kb_sample_root = (
        project / "08_qc_summary/27_dna_100kb/sample_matrices"
    )
    dna_100kb_samples = sum(
        is_validated(path)
        for path in dna_100kb_sample_root.glob(
            "*.dna_100kb.summary.tsv"
        )
    )
    dna_100kb_final = (
        project
        / "08_qc_summary/27_dna_100kb"
        / "dna_100kb_complete_summary.tsv"
    )
    cross_well_final = (
        project
        / "08_qc_summary/28_cross_well_overview"
        / "snmar_96well_batch_overall_summary.tsv"
    )
    cb_overlap_final = (
        project
        / "08_qc_summary/29_cb_overlap_96wells"
        / "cb_overlap_audit_summary.tsv"
    )
    threshold_audit_final = (
        project
        / "08_qc_summary/30_cell_qc_threshold_audit"
        / "cell_qc_threshold_audit_summary.tsv"
    )
    rna_species_alignment_root = (
        project / "08_qc_summary/31_rna_species_residual_audit/by_file"
    )
    rna_species_alignment_completed = sum(
        is_validated_single_row(path, RNA_SPECIES_ALIGNMENT_SIGNATURE)
        for path in rna_species_alignment_root.glob(
            "*.rna_species_alignment_audit.tsv"
        )
    )
    rna_species_audit_final = (
        project
        / "08_qc_summary/31_rna_species_residual_audit"
        / "rna_species_residual_audit_summary.tsv"
    )
    dna_umap_root = project / "08_qc_summary/32_dna_100kb_umap"
    dna_umap_completed = sum(
        is_validated(
            dna_umap_root / threshold / "dna_100kb_umap_summary.tsv"
        )
        for threshold in ("D_total_20k", "D_total_50k")
    )
    dna_umap_comparison = (
        dna_umap_root / "dna_100kb_umap_comparison_summary.tsv"
    )
    expected_tasks = 768 if tasks else 0
    expected_wells = 96 if tasks else 0
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
            "stage": "68-70_DNA_well_validation",
            "completed": validation_completed,
            "expected": expected_wells,
            "ready": int(
                expected_wells == 96
                and validation_completed == expected_wells
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
        {
            "stage": "85-86_DNA_IGV_species_BigWigs",
            "completed": int(is_validated(igv_bigwig_final)),
            "expected": 1,
            "ready": int(is_validated(igv_bigwig_final)),
        },
        {
            "stage": "87-89_DNA_100kb_well_matrices",
            "completed": dna_100kb_samples,
            "expected": expected_wells,
            "ready": int(
                expected_wells == 96
                and dna_100kb_samples == expected_wells
            ),
        },
        {
            "stage": "90-91_DNA_100kb_complete_matrix",
            "completed": int(is_validated(dna_100kb_final)),
            "expected": 1,
            "ready": int(is_validated(dna_100kb_final)),
        },
        {
            "stage": "92-93_96well_batch_species_summary",
            "completed": int(is_validated(cross_well_final)),
            "expected": 1,
            "ready": int(is_validated(cross_well_final)),
        },
        {
            "stage": "94-95_CB_overlap_across_96_wells",
            "completed": int(is_validated(cb_overlap_final)),
            "expected": 1,
            "ready": int(is_validated(cb_overlap_final)),
        },
        {
            "stage": "96-97_DNA_RNA_cell_QC_threshold_audit",
            "completed": int(is_validated(threshold_audit_final)),
            "expected": 1,
            "ready": int(is_validated(threshold_audit_final)),
        },
        {
            "stage": "98-100_RNA_species_alignment_audit",
            "completed": rna_species_alignment_completed,
            "expected": 96,
            "ready": int(rna_species_alignment_completed == 96),
        },
        {
            "stage": "101-102_RNA_species_residual_summary",
            "completed": int(is_validated(rna_species_audit_final)),
            "expected": 1,
            "ready": int(is_validated(rna_species_audit_final)),
        },
        {
            "stage": "103-106_DNA_100kb_UMAP_D20_D50",
            "completed": dna_umap_completed,
            "expected": 2,
            "ready": int(dna_umap_completed == 2),
        },
        {
            "stage": "107-108_DNA_100kb_UMAP_comparison",
            "completed": int(is_validated(dna_umap_comparison)),
            "expected": 1,
            "ready": int(is_validated(dna_umap_comparison)),
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
