#!/usr/bin/env python3
"""Step 149: summarize the fixed high-variable shared-promoter sensitivity grid."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


RUN_SIGNATURE = "snmar_dna_tss_umap_20260808"
PIPELINE_SIGNATURE = "snmar_dna_tss_hvf_comparison_20260827"
RUNS = (
    (20_000, 1_000, "HVF1000_F20_C100"),
    (20_000, 2_000, "HVF2000_F20_C100"),
    (50_000, 1_000, "HVF1000_F20_C100"),
    (50_000, 2_000, "HVF2000_F20_C100"),
)


def short_count(value: int) -> str:
    return f"{value // 1000}k" if value % 1000 == 0 else str(value)


def run_id(depth: int, tag: str) -> str:
    return f"ortholog_single_gene_D{short_count(depth)}_TSSctx100_{tag}"


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode + "t", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with open_text(path, "r") as handle:
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


def atomic_tsv(path: Path, fields: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        with open(temporary, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
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


def atomic_metrics(path: Path, rows: list[tuple[str, object]]) -> None:
    atomic_tsv(
        path,
        ["metric", "value"],
        ({"metric": metric, "value": value} for metric, value in rows),
    )


def neighbor_sets(coordinates, neighbors: int, NearestNeighbors):
    k = min(neighbors, len(coordinates) - 1)
    indices = NearestNeighbors(n_neighbors=k + 1).fit(coordinates).kneighbors(
        coordinates, return_distance=False
    )[:, 1:]
    return [set(map(int, row)) for row in indices], k


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.neighbors < 2:
        raise ValueError("neighbors must be at least 2")
    try:
        import numpy as np
        from sklearn.neighbors import NearestNeighbors
    except ImportError as error:
        raise RuntimeError("Step 149 requires numpy and scikit-learn") from error

    project = args.project
    qc_root = args.qc_root or project / "08_qc_summary/39_dna_tss_umap"
    summaries: dict[str, dict[str, str]] = {}
    coordinates: dict[str, dict[int, dict[str, str]]] = {}
    per_run = []
    for depth, maximum, tag in RUNS:
        name = run_id(depth, tag)
        summary_path = qc_root / name / "dna_tss_umap_summary.tsv"
        summary = read_metrics(summary_path)
        expected = {
            "status": "validated",
            "pipeline_signature": RUN_SIGNATURE,
            "run_id": name,
            "run_tag": tag,
            "feature_space": "ortholog_single_gene",
            "minimum_whole_genome_WCG_plus_GCH_calls": str(depth),
            "minimum_shared_promoter_calls_per_context": "100",
            "minimum_cells_per_feature": "20",
            "minimum_covered_calls_per_feature": "100",
            "maximum_features_per_context": str(maximum),
        }
        if any(summary.get(key) != value for key, value in expected.items()):
            raise ValueError(f"invalid high-variable TSS run: {summary_path}")
        rows = read_tsv(Path(summary["coordinates"]))
        keyed = {int(row["cell_index"]): row for row in rows}
        if len(keyed) != len(rows) or len(rows) != int(summary["retained_cells"]):
            raise ValueError(f"duplicate/missing coordinates for {name}")
        valid_labels = [
            row["DNA_mapping_derived_cell_line"]
            for row in rows
            if row["DNA_mapping_derived_cell_line"] in {"GM12878", "PK15"}
        ]
        fractions = [valid_labels.count(label) / len(valid_labels) for label in ("GM12878", "PK15")]
        random_baseline = sum(value * value for value in fractions)
        summaries[name] = summary
        coordinates[name] = keyed
        record = {
            "run_id": name,
            "minimum_whole_genome_calls": depth,
            "top_variable_features_requested": maximum,
            "retained_cells": summary["retained_cells"],
            "mapping_confident_cells": len(valid_labels),
            "random_same_species_neighbor_baseline": random_baseline,
            "WCG_selected_features": summary["WCG_selected_features"],
            "GCH_selected_features": summary["GCH_selected_features"],
        }
        for context in ("WCG", "GCH", "joint"):
            agreement = float(summary[f"{context}_cell_line_knn_agreement_posthoc"])
            record[f"{context}_knn_agreement"] = agreement
            record[f"{context}_knn_excess_over_random"] = agreement - random_baseline
            record[f"{context}_silhouette"] = summary[
                f"{context}_cell_line_silhouette_posthoc"
            ]
        per_run.append(record)

    comparisons = (
        (run_id(20_000, RUNS[0][2]), run_id(50_000, RUNS[2][2]), "D20_vs_D50_top1000"),
        (run_id(20_000, RUNS[1][2]), run_id(50_000, RUNS[3][2]), "D20_vs_D50_top2000"),
        (run_id(20_000, RUNS[0][2]), run_id(20_000, RUNS[1][2]), "top1000_vs_top2000_D20"),
        (run_id(50_000, RUNS[2][2]), run_id(50_000, RUNS[3][2]), "top1000_vs_top2000_D50"),
    )
    stability_rows = []
    for first, second, comparison in comparisons:
        shared = sorted(set(coordinates[first]) & set(coordinates[second]))
        if len(shared) < 20:
            raise ValueError(f"too few shared cells for {comparison}")
        for context in ("WCG", "GCH", "joint"):
            columns = (f"{context}_UMAP1", f"{context}_UMAP2")
            first_xy = np.asarray(
                [[float(coordinates[first][index][field]) for field in columns] for index in shared]
            )
            second_xy = np.asarray(
                [[float(coordinates[second][index][field]) for field in columns] for index in shared]
            )
            first_sets, k = neighbor_sets(first_xy, args.neighbors, NearestNeighbors)
            second_sets, _ = neighbor_sets(second_xy, args.neighbors, NearestNeighbors)
            overlaps = [len(a & b) / k for a, b in zip(first_sets, second_sets)]
            unions = [len(a | b) for a, b in zip(first_sets, second_sets)]
            jaccards = [len(a & b) / union for a, b, union in zip(first_sets, second_sets, unions)]
            stability_rows.append(
                {
                    "comparison": comparison,
                    "embedding": context,
                    "shared_cells": len(shared),
                    "neighbors": k,
                    "mean_neighbor_overlap_fraction": float(np.mean(overlaps)),
                    "median_neighbor_overlap_fraction": float(np.median(overlaps)),
                    "mean_neighbor_jaccard": float(np.mean(jaccards)),
                    "median_neighbor_jaccard": float(np.median(jaccards)),
                }
            )

    per_run_path = qc_root / "dna_tss_hvf_per_run.tsv"
    stability_path = qc_root / "dna_tss_hvf_neighbor_stability.tsv"
    atomic_tsv(per_run_path, list(per_run[0]), per_run)
    atomic_tsv(stability_path, list(stability_rows[0]), stability_rows)
    wcg_supported = any(
        float(row["WCG_silhouette"]) >= 0.10
        and float(row["WCG_knn_excess_over_random"]) >= 0.05
        for row in per_run
    )
    summary_path = qc_root / "dna_tss_hvf_comparison_summary.tsv"
    atomic_metrics(
        summary_path,
        [
            ("status", "validated"),
            ("pipeline_signature", PIPELINE_SIGNATURE),
            ("runs_compared", ",".join(name for name in summaries)),
            ("predeclared_minimum_WCG_silhouette", 0.10),
            ("predeclared_minimum_WCG_knn_excess_over_random", 0.05),
            ("any_run_meets_predeclared_WCG_separation", str(wcg_supported).lower()),
            ("labels_used_for_fitting_or_feature_selection", "no"),
            ("mapping_labels_used_posthoc", "yes"),
            ("per_run_metrics", per_run_path),
            ("neighbor_stability", stability_path),
        ],
    )
    print(f"[OK] high-variable TSS comparison: runs={len(per_run)} WCG_supported={wcg_supported}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
