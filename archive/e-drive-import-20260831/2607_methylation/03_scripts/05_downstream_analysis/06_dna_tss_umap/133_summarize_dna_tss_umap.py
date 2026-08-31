#!/usr/bin/env python3
"""Step 133: compare validated D20 and D50 shared-promoter TSS UMAPs."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


RUN_SIGNATURE = "snmar_dna_tss_umap_20260808"
PIPELINE_SIGNATURE = "snmar_dna_tss_umap_comparison_20260808"
RUNS = (
    "ortholog_single_gene_D20k_TSSctx25",
    "ortholog_single_gene_D50k_TSSctx25",
)


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
    result = {row["metric"]: row["value"] for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"duplicate metric: {path}")
    return result


def atomic_tsv(path: Path, fields: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        with open(temporary, "wt", encoding="utf-8", newline="") as handle:
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
        raise RuntimeError("Step 133 requires numpy and scikit-learn") from error
    project = args.project
    summaries = {}
    coordinates = {}
    per_run = []
    for run in RUNS:
        summary_path = (
            project
            / "08_qc_summary/39_dna_tss_umap"
            / run
            / "dna_tss_umap_summary.tsv"
        )
        summary = read_metrics(summary_path)
        if (
            summary.get("status") != "validated"
            or summary.get("pipeline_signature") != RUN_SIGNATURE
            or summary.get("run_id") != run
            or summary.get("feature_space") != "ortholog_single_gene"
        ):
            raise ValueError(f"invalid TSS UMAP summary: {summary_path}")
        rows = read_tsv(Path(summary["coordinates"]))
        keyed = {int(row["cell_index"]): row for row in rows}
        if len(keyed) != len(rows) or len(rows) != int(summary["retained_cells"]):
            raise ValueError(f"duplicate/missing coordinates for {run}")
        summaries[run] = summary
        coordinates[run] = keyed
        per_run.append(
            {
                "run_id": run,
                "minimum_whole_genome_calls": summary[
                    "minimum_whole_genome_WCG_plus_GCH_calls"
                ],
                "minimum_promoter_calls_per_context": summary[
                    "minimum_shared_promoter_calls_per_context"
                ],
                "retained_cells": summary["retained_cells"],
                "GM12878_cells": summary["GM12878_cells_mapping_derived"],
                "PK15_cells": summary["PK15_cells_mapping_derived"],
                "ambiguous_cells": summary["ambiguous_cells"],
                "WCG_knn_agreement": summary[
                    "WCG_cell_line_knn_agreement_posthoc"
                ],
                "GCH_knn_agreement": summary[
                    "GCH_cell_line_knn_agreement_posthoc"
                ],
                "joint_knn_agreement": summary[
                    "joint_cell_line_knn_agreement_posthoc"
                ],
                "WCG_silhouette": summary[
                    "WCG_cell_line_silhouette_posthoc"
                ],
                "GCH_silhouette": summary[
                    "GCH_cell_line_silhouette_posthoc"
                ],
                "joint_silhouette": summary[
                    "joint_cell_line_silhouette_posthoc"
                ],
            }
        )
    low, high = RUNS
    low_cells, high_cells = set(coordinates[low]), set(coordinates[high])
    if not high_cells.issubset(low_cells):
        raise ValueError("D50 retained cells are not a subset of D20")
    shared = sorted(high_cells)
    if len(shared) < 20:
        raise ValueError("too few shared D20/D50 cells")
    label_field = "DNA_mapping_derived_cell_line"
    if any(
        coordinates[low][index][label_field]
        != coordinates[high][index][label_field]
        for index in shared
    ):
        raise ValueError("mapping-derived labels differ across runs")
    stability_rows = []
    for embedding in ("WCG", "GCH", "joint"):
        columns = (f"{embedding}_UMAP1", f"{embedding}_UMAP2")
        low_xy = np.asarray(
            [[float(coordinates[low][index][field]) for field in columns] for index in shared]
        )
        high_xy = np.asarray(
            [[float(coordinates[high][index][field]) for field in columns] for index in shared]
        )
        low_sets, k = neighbor_sets(low_xy, args.neighbors, NearestNeighbors)
        high_sets, _ = neighbor_sets(high_xy, args.neighbors, NearestNeighbors)
        overlaps = [len(a & b) / k for a, b in zip(low_sets, high_sets)]
        jaccards = [len(a & b) / len(a | b) for a, b in zip(low_sets, high_sets)]
        stability_rows.append(
            {
                "embedding": embedding,
                "shared_cells": len(shared),
                "neighbors": k,
                "mean_neighbor_overlap_fraction": float(np.mean(overlaps)),
                "median_neighbor_overlap_fraction": float(np.median(overlaps)),
                "mean_neighbor_jaccard": float(np.mean(jaccards)),
                "median_neighbor_jaccard": float(np.median(jaccards)),
            }
        )
    qc_root = args.qc_root or (
        project / "08_qc_summary/39_dna_tss_umap"
    )
    per_run_path = qc_root / "dna_tss_umap_per_run.tsv"
    stability_path = qc_root / "dna_tss_umap_D20_D50_stability.tsv"
    atomic_tsv(per_run_path, list(per_run[0]), per_run)
    atomic_tsv(stability_path, list(stability_rows[0]), stability_rows)
    summary_path = qc_root / "dna_tss_umap_comparison_summary.tsv"
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("feature_space", "ortholog_single_gene"),
        ("runs_compared", ",".join(RUNS)),
        ("D20_cells", len(low_cells)),
        ("D50_cells", len(high_cells)),
        ("D50_subset_of_D20", "true"),
        ("shared_cells", len(shared)),
        ("stability_neighbors", min(args.neighbors, len(shared) - 1)),
        ("per_run_metrics", per_run_path),
        ("stability_metrics", stability_path),
        ("cell_line_labels", "mapping_derived_posthoc_not_independent_validation"),
    ]
    for row in stability_rows:
        embedding = row["embedding"]
        summary_rows.append(
            (f"{embedding}_mean_neighbor_overlap_fraction", row["mean_neighbor_overlap_fraction"])
        )
        summary_rows.append(
            (f"{embedding}_mean_neighbor_jaccard", row["mean_neighbor_jaccard"])
        )
    atomic_metrics(summary_path, summary_rows)
    print(f"[OK] TSS UMAP comparison: shared_cells={len(shared)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
