#!/usr/bin/env python3
"""Step 120: exploratory UMAP for the low-depth exon-supported RNA matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_rna_exon_salvage_umap_20260807"


def opener_for(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with opener_for(path)(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        yield from reader


def atomic_tsv(
    path: Path, fields: list[str], rows: Iterable[dict[str, object]], gzip_output: bool = False
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        opener = gzip.open if gzip_output else open
        with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def balanced_variable_features(
    variances, species, detected_cells, maximum_per_species: int
):
    import numpy as np

    selected: list[int] = []
    for label in ("human", "pig"):
        candidates = np.flatnonzero((species == label) & (detected_cells >= 3))
        order = candidates[np.argsort(variances[candidates])[::-1]]
        selected.extend(order[:maximum_per_species].tolist())
    return np.asarray(sorted(set(selected)), dtype=np.int64)


def same_label_knn_fraction(embedding, labels, neighbors: int = 15) -> float | None:
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    labels = np.asarray(labels)
    keep = np.flatnonzero(np.isin(labels, ["human", "pig"]))
    if keep.size < 3:
        return None
    k = min(neighbors + 1, keep.size)
    indices = NearestNeighbors(n_neighbors=k).fit(embedding[keep]).kneighbors(return_distance=False)
    neighbor_indices = indices[:, 1:]
    if neighbor_indices.shape[1] == 0:
        return None
    return float(np.mean(labels[keep][:, None] == labels[keep][neighbor_indices]))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--matrix-h5", type=Path)
    parser.add_argument("--cell-metadata", type=Path)
    parser.add_argument("--feature-metadata", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--features-per-species", type=int, default=1000)
    parser.add_argument("--neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=20260807)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.features_per_species < 10 or args.neighbors < 2 or not 0 <= args.min_dist <= 1:
        raise ValueError("invalid UMAP parameter")
    project = args.project
    matrix_root = project / "10_downstream_analysis/06_rna_exon_salvage_matrix"
    matrix_h5 = args.matrix_h5 or (matrix_root / "rna_exon_salvage_complete_sparse.h5")
    cell_metadata_path = args.cell_metadata or (matrix_root / "rna_exon_salvage_cell_metadata.tsv.gz")
    feature_metadata_path = args.feature_metadata or (
        matrix_root / "raw_feature_bc_matrix/feature_metadata.tsv.gz"
    )
    output_root = args.output_root or (
        project / "10_downstream_analysis/07_rna_exon_salvage_umap"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/37_rna_exon_salvage_umap"
    )
    for path in (matrix_h5, cell_metadata_path, feature_metadata_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty input: {path}")

    import h5py
    import matplotlib
    import numpy as np
    import scipy.sparse as sp
    import umap
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import StandardScaler

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    all_metadata = {row["cell_id"]: row for row in iter_tsv(cell_metadata_path)}
    feature_rows = list(iter_tsv(feature_metadata_path))
    with h5py.File(matrix_h5, "r") as handle:
        if handle.attrs.get("pipeline_signature") != "snmar_rna_exon_salvage_matrix_20260807":
            raise ValueError("unexpected exon salvage matrix signature")
        feature_ids = [value.decode() for value in handle["feature_id"][:]]
        group = handle["filtered"]
        cell_ids = [value.decode() for value in group["cell_id"][:]]
        shape = tuple(int(value) for value in group.attrs["shape"])
        rows = group["row"][:]
        cols = group["col"][:]
        data = group["data"][:].astype(np.float64)
    if shape != (len(feature_ids), len(cell_ids)):
        raise ValueError("HDF5 feature/cell axes do not reconcile")
    if len(feature_rows) != len(feature_ids) or [row["gene_id"] for row in feature_rows] != feature_ids:
        raise ValueError("feature metadata does not match HDF5 feature axis")
    if any(cell not in all_metadata for cell in cell_ids):
        raise ValueError("filtered HDF5 contains cell absent from metadata")
    matrix = sp.coo_matrix((data, (rows, cols)), shape=shape).tocsr().T
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    genes = np.asarray((matrix > 0).sum(axis=1)).ravel()
    if np.any(totals < 5) or np.any(genes < 3):
        raise ValueError("filtered matrix violates UMI>=5/gene>=3")

    normalized = sp.diags(10_000.0 / totals) @ matrix
    normalized.data = np.log1p(normalized.data)
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    means = np.asarray(normalized.mean(axis=0)).ravel()
    squares = np.asarray(normalized.power(2).mean(axis=0)).ravel()
    variances = np.maximum(squares - means * means, 0.0)
    species = np.asarray([row["species"] for row in feature_rows])
    selected = balanced_variable_features(
        variances, species, detected, args.features_per_species
    )
    if selected.size < 20:
        raise ValueError(f"too few recurrent variable genes: {selected.size}")
    dense = normalized[:, selected].toarray()
    dense = StandardScaler().fit_transform(dense)
    n_components = min(30, dense.shape[0] - 2, dense.shape[1] - 1)
    if n_components < 2:
        raise ValueError("too few cells/features for PCA")
    pcs = PCA(n_components=n_components, random_state=args.seed).fit_transform(dense)
    effective_neighbors = min(args.neighbors, max(2, len(cell_ids) - 1))
    embedding = umap.UMAP(
        n_neighbors=effective_neighbors,
        min_dist=args.min_dist,
        metric="euclidean",
        random_state=args.seed,
        n_jobs=1,
    ).fit_transform(pcs)

    metadata = [all_metadata[cell] for cell in cell_ids]
    labels = np.asarray([row["dominant_species"] for row in metadata])
    well_groups = np.asarray([row["well_group"] for row in metadata])
    confident = np.flatnonzero(np.isin(labels, ["human", "pig"]))
    species_silhouette = (
        float(silhouette_score(embedding[confident], labels[confident]))
        if confident.size >= 3 and len(set(labels[confident])) == 2
        else None
    )
    knn = same_label_knn_fraction(embedding, labels, effective_neighbors)

    coordinate_rows: list[dict[str, object]] = []
    for index, (cell, row) in enumerate(zip(cell_ids, metadata)):
        coordinate_rows.append(
            {
                "cell_id": cell,
                "UMAP1": float(embedding[index, 0]),
                "UMAP2": float(embedding[index, 1]),
                "sample": row["sample"],
                "well_index": int(row["well_index"]),
                "well_group": row["well_group"],
                "CB": row["CB"],
                "exact_umis": int(row["exact_umis"]),
                "genes_detected": int(row["genes_detected"]),
                "human_exact_umis": int(row["human_exact_umis"]),
                "pig_exact_umis": int(row["pig_exact_umis"]),
                "dominant_species": row["dominant_species"],
                "major_species_fraction": row["major_species_fraction"],
            }
        )
    output_root.mkdir(parents=True, exist_ok=True)
    coordinates_path = output_root / "rna_exon_salvage_umap_coordinates.tsv.gz"
    atomic_tsv(coordinates_path, list(coordinate_rows[0].keys()), coordinate_rows, gzip_output=True)
    selected_rows = [
        {
            "feature_index": int(index),
            "gene_id": feature_rows[index]["gene_id"],
            "gene_name": feature_rows[index]["gene_name"],
            "species": feature_rows[index]["species"],
            "detected_cells": int(detected[index]),
            "log_normalized_variance": float(variances[index]),
        }
        for index in selected
    ]
    atomic_tsv(output_root / "rna_exon_salvage_selected_features.tsv.gz", list(selected_rows[0].keys()), selected_rows, gzip_output=True)

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    colors = {"human": "#4c8bf5", "pig": "#f05a5a", "tie": "#888888"}
    for label in ("human", "pig", "tie"):
        mask = labels == label
        if np.any(mask):
            axes[0, 0].scatter(embedding[mask, 0], embedding[mask, 1], s=14, alpha=0.75, color=colors[label], label=label, linewidths=0)
    axes[0, 0].legend(frameon=False)
    axes[0, 0].set_title("exon exact-UMI dominant species")
    for label, color in (("GM12878_only_wells", "#2ca58d"), ("GM12878_PK15_mix_wells", "#aaaaaa")):
        mask = well_groups == label
        axes[0, 1].scatter(embedding[mask, 0], embedding[mask, 1], s=14, alpha=0.7, color=color, label=label, linewidths=0)
    axes[0, 1].legend(frameon=False)
    axes[0, 1].set_title("physical-well design")
    scatter = axes[1, 0].scatter(embedding[:, 0], embedding[:, 1], c=np.log10(totals + 1), s=15, cmap="viridis", linewidths=0)
    axes[1, 0].set_title("log10(exon exact UMIs + 1)")
    fig.colorbar(scatter, ax=axes[1, 0], shrink=0.8)
    scatter = axes[1, 1].scatter(embedding[:, 0], embedding[:, 1], c=genes, s=15, cmap="magma", linewidths=0)
    axes[1, 1].set_title("exon genes detected")
    fig.colorbar(scatter, ax=axes[1, 1], shrink=0.8)
    for axis in axes.flat:
        axis.set_xlabel("UMAP1")
        axis.set_ylabel("UMAP2")
    fig.suptitle("snMAR exploratory RNA exon-supported UMAP (low depth)", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    figure_path = output_root / "rna_exon_salvage_umap.png"
    fig.savefig(figure_path, dpi=220)
    fig.savefig(figure_path.with_suffix(".pdf"))
    plt.close(fig)

    qc_root.mkdir(parents=True, exist_ok=True)
    label_counts = Counter(labels.tolist())
    summary_rows = [
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE, "interpretation": "validated implementation signature"},
        {"metric": "status", "value": "validated", "interpretation": "matrix, metadata, PCA and UMAP completed"},
        {"metric": "analysis_scope", "value": "exploratory_low_depth", "interpretation": "not formal high-quality RNA expression clustering"},
        {"metric": "cells", "value": len(cell_ids), "interpretation": "UMI>=5 and gene>=3 exon-supported cells"},
        {"metric": "human_cells", "value": label_counts["human"], "interpretation": "dominant human exon exact UMIs"},
        {"metric": "pig_cells", "value": label_counts["pig"], "interpretation": "dominant pig exon exact UMIs"},
        {"metric": "tie_cells", "value": label_counts["tie"], "interpretation": "equal human and pig exon exact UMIs"},
        {"metric": "selected_features", "value": len(selected), "interpretation": "balanced recurrent variable genes"},
        {"metric": "PCA_components", "value": n_components, "interpretation": "scaled log-normalized expression PCs"},
        {"metric": "UMAP_neighbors", "value": effective_neighbors, "interpretation": "exploratory embedding parameter"},
        {"metric": "species_knn_agreement", "value": knn if knn is not None else "NA", "interpretation": "internal, not independent validation"},
        {"metric": "species_silhouette", "value": species_silhouette if species_silhouette is not None else "NA", "interpretation": "internal, species-specific genes make this partly circular"},
        {"metric": "DNA_label_transfer_scope", "value": "matched_subset_only", "interpretation": "only cells with the same physical-well|CB may label DNA"},
    ]
    atomic_tsv(qc_root / "rna_exon_salvage_umap_summary.tsv", ["metric", "value", "interpretation"], summary_rows)
    print(f"[OK] exploratory RNA exon UMAP: cells={len(cell_ids)} features={len(selected)}")
    print(f"[OK] coordinates: {coordinates_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
