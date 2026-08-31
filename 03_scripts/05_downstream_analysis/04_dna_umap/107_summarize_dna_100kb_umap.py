#!/usr/bin/env python3
"""Step 107: validate and compare the D20 and D50 DNA UMAP runs."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


RUN_SIGNATURE = "snmar_dna_100kb_umap_20260806"
PIPELINE_SIGNATURE = "snmar_dna_100kb_umap_comparison_20260806"
THRESHOLDS = ("D_total_20k", "D_total_50k")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def read_metrics(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    metrics = {row["metric"]: row["value"] for row in rows}
    if len(metrics) != len(rows):
        raise ValueError(f"duplicate metric: {path}")
    return metrics


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
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_root = args.input_root or (
        args.project / "08_qc_summary/32_dna_100kb_umap"
    )
    output_root = args.output_root or input_root
    comparison_rows = []
    for threshold in THRESHOLDS:
        summary_path = input_root / threshold / "dna_100kb_umap_summary.tsv"
        metrics = read_metrics(summary_path)
        if (
            metrics.get("status") != "validated"
            or metrics.get("pipeline_signature") != RUN_SIGNATURE
            or metrics.get("threshold_id") != threshold
        ):
            raise ValueError(f"invalid DNA UMAP run: {summary_path}")
        required_outputs = (
            "coordinates",
            "selected_features",
            "embeddings_npz",
            "overview_png",
            "overview_pdf",
        )
        for field in required_outputs:
            path = Path(metrics[field])
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"{threshold} missing {field}: {path}"
                )
        comparison_rows.append(
            {
                "threshold_id": threshold,
                "minimum_total_signal_covered_calls": metrics[
                    "minimum_total_WCG_plus_GCH_covered_calls"
                ],
                "retained_cells": metrics["retained_cells"],
                "physical_wells_represented": metrics[
                    "physical_wells_represented"
                ],
                "human_dominant_cells_posthoc": metrics[
                    "human_dominant_cells_posthoc"
                ],
                "pig_dominant_cells_posthoc": metrics[
                    "pig_dominant_cells_posthoc"
                ],
                "WCG_selected_bins": metrics["WCG_selected_bins"],
                "GCH_selected_bins": metrics["GCH_selected_bins"],
                "WCG_species_knn_agreement_posthoc": metrics[
                    "WCG_species_knn_agreement_posthoc"
                ],
                "GCH_species_knn_agreement_posthoc": metrics[
                    "GCH_species_knn_agreement_posthoc"
                ],
                "joint_species_knn_agreement_posthoc": metrics[
                    "joint_species_knn_agreement_posthoc"
                ],
                "WCG_species_silhouette_posthoc": metrics[
                    "WCG_species_silhouette_posthoc"
                ],
                "GCH_species_silhouette_posthoc": metrics[
                    "GCH_species_silhouette_posthoc"
                ],
                "joint_species_silhouette_posthoc": metrics[
                    "joint_species_silhouette_posthoc"
                ],
                "coordinates": metrics["coordinates"],
                "overview_png": metrics["overview_png"],
                "run_summary": summary_path,
                "status": "validated",
            }
        )
    comparison_path = (
        output_root / "dna_100kb_umap_threshold_comparison.tsv"
    )
    atomic_tsv(comparison_path, list(comparison_rows[0]), comparison_rows)
    summary_path = (
        output_root / "dna_100kb_umap_comparison_summary.tsv"
    )
    atomic_tsv(
        summary_path,
        ["metric", "value"],
        [
            {"metric": "status", "value": "validated"},
            {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE},
            {"metric": "primary_threshold", "value": "D_total_20k"},
            {"metric": "sensitivity_threshold", "value": "D_total_50k"},
            {"metric": "validated_runs", "value": len(comparison_rows)},
            {
                "metric": "species_role",
                "value": "post_hoc_validation_only",
            },
            {"metric": "comparison_table", "value": comparison_path},
        ],
    )
    print(f"[OK] DNA 100 kb UMAP comparison: {comparison_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
