#!/usr/bin/env python3
"""Step 153: summarize D20/D50 reference-aware promoter UMAPs.

The mapping-derived cell-line label is post hoc. Metrics are reported both for
all retained cells and, critically, within the 49--96 mixed-well set so that
cell-line separation is not confused with the pure-well/mixed-well design.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from collections import Counter
from pathlib import Path


RUN_SIGNATURE = "snmar_dna_tss_umap_20260808"
PIPELINE_SIGNATURE = "snmar_dna_tss_reference_aware_audit_20260828"
RUNS = (
    "promoter_cluster_D20k_TSSctx25_reference_aware",
    "promoter_cluster_D50k_TSSctx25_reference_aware",
)
EMBEDDINGS = ("WCG", "GCH", "joint")


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


def class_baseline(labels) -> float:
    counts = Counter(labels)
    total = len(labels)
    if total < 2:
        return float("nan")
    return sum(value * (value - 1) for value in counts.values()) / (
        total * (total - 1)
    )


def separation_metrics(coordinates, labels, np, NearestNeighbors, silhouette_score):
    labels = np.asarray(labels)
    if len(coordinates) < 3 or len(set(labels)) < 2:
        return {
            "cells": len(coordinates),
            "class_1_cells": int(sum(labels == labels[0])) if len(labels) else 0,
            "class_2_cells": 0,
            "knn_agreement": float("nan"),
            "random_agreement": class_baseline(labels),
            "knn_excess_over_random": float("nan"),
            "silhouette": float("nan"),
        }
    k = min(15, len(coordinates) - 1)
    indices = NearestNeighbors(n_neighbors=k + 1).fit(coordinates).kneighbors(
        coordinates, return_distance=False
    )[:, 1:]
    agreement = float(np.mean(labels[indices] == labels[:, None]))
    baseline = class_baseline(labels)
    counts = sorted(Counter(labels).values(), reverse=True)
    return {
        "cells": len(coordinates),
        "class_1_cells": counts[0],
        "class_2_cells": counts[1],
        "knn_agreement": agreement,
        "random_agreement": baseline,
        "knn_excess_over_random": agreement - baseline,
        "silhouette": float(silhouette_score(coordinates, labels)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import numpy as np
        from sklearn.metrics import silhouette_score
        from sklearn.neighbors import NearestNeighbors
    except ImportError as error:
        raise RuntimeError("Step 153 requires numpy and scikit-learn") from error

    qc_root = args.qc_root or (
        args.project / "08_qc_summary/39_dna_tss_umap/reference_aware_audit"
    )
    output_rows = []
    for run_id in RUNS:
        summary_path = (
            args.project
            / "08_qc_summary/39_dna_tss_umap"
            / run_id
            / "dna_tss_umap_summary.tsv"
        )
        summary = read_metrics(summary_path)
        if (
            summary.get("status") != "validated"
            or summary.get("pipeline_signature") != RUN_SIGNATURE
            or summary.get("run_id") != run_id
            or summary.get("feature_space") != "promoter_cluster"
        ):
            raise ValueError(f"invalid reference-aware promoter run: {summary_path}")
        rows = read_tsv(Path(summary["coordinates"]))
        if len(rows) != int(summary["retained_cells"]):
            raise ValueError(f"coordinate count differs from summary: {run_id}")
        for embedding in EMBEDDINGS:
            fields = (f"{embedding}_UMAP1", f"{embedding}_UMAP2")
            for scope in ("all_confident_cells", "mixed_wells_only"):
                selected = [
                    row
                    for row in rows
                    if row["DNA_mapping_derived_cell_line"] in {"GM12878", "PK15"}
                    and (scope == "all_confident_cells" or row["well_set"] == "cellMix")
                ]
                xy = np.asarray(
                    [[float(row[field]) for field in fields] for row in selected],
                    dtype=float,
                )
                labels = [row["DNA_mapping_derived_cell_line"] for row in selected]
                metrics = separation_metrics(
                    xy, labels, np, NearestNeighbors, silhouette_score
                )
                counts = Counter(labels)
                output_rows.append(
                    {
                        "run_id": run_id,
                        "embedding": embedding,
                        "scope": scope,
                        "cells": metrics["cells"],
                        "GM12878_cells": counts["GM12878"],
                        "PK15_cells": counts["PK15"],
                        "knn_agreement": metrics["knn_agreement"],
                        "random_agreement": metrics["random_agreement"],
                        "knn_excess_over_random": metrics["knn_excess_over_random"],
                        "silhouette": metrics["silhouette"],
                    }
                )

    per_run_path = qc_root / "dna_tss_reference_aware_per_run_embedding.tsv"
    atomic_tsv(per_run_path, list(output_rows[0]), output_rows)
    summary_path = qc_root / "dna_tss_reference_aware_audit_summary.tsv"
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("feature_space", "species_specific_human_plus_pig_promoter_axes"),
        ("runs", ",".join(RUNS)),
        ("primary_context", "WCG"),
        ("primary_scope", "mixed_wells_only"),
        ("labels", "mapping_derived_posthoc_not_independent_validation"),
        ("comparison_warning", "reference_identity_and_missingness_remain_in_feature_space"),
        ("per_run_embedding_metrics", per_run_path),
    ]
    for row in output_rows:
        if row["embedding"] != "WCG" or row["scope"] != "mixed_wells_only":
            continue
        prefix = row["run_id"].split("_")[2]
        summary_rows.extend(
            [
                (f"{prefix}_mixed_WCG_cells", row["cells"]),
                (f"{prefix}_mixed_WCG_GM12878_cells", row["GM12878_cells"]),
                (f"{prefix}_mixed_WCG_PK15_cells", row["PK15_cells"]),
                (f"{prefix}_mixed_WCG_knn_excess", row["knn_excess_over_random"]),
                (f"{prefix}_mixed_WCG_silhouette", row["silhouette"]),
            ]
        )
    atomic_metrics(summary_path, summary_rows)
    print(f"[OK] reference-aware promoter audit: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
