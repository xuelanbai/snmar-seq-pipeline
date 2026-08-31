#!/usr/bin/env python3
"""Step 103: build unsupervised DNA 100 kb WCG/GCH UMAP embeddings.

Cells are selected only by total WCG+GCH covered calls. Human/pig and
physical-well annotations are joined after the embedding and never enter
feature selection, SVD, or UMAP. Missing cell/bin entries mean no coverage;
they are represented by zero Pearson residual (the bin-level expectation),
not by a zero signal fraction.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import tempfile
from pathlib import Path
from typing import Iterable, Iterator


SOURCE_SIGNATURE = "snmar_dna_100kb_complete_20260804"
PIPELINE_SIGNATURE = "snmar_dna_100kb_umap_20260806"
CONTEXTS = {1: "WCG", 2: "GCH"}


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


def read_tsv(path: Path) -> list[dict[str, str]]:
    rows = list(iter_tsv(path))
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
    path: Path,
    fields: list[str],
    rows: Iterable[dict[str, object]],
    gzip_output: bool = False,
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


def threshold_id(min_total_calls: int) -> str:
    if min_total_calls < 1:
        raise ValueError("minimum total calls must be positive")
    if min_total_calls % 1000 == 0:
        return f"D_total_{min_total_calls // 1000}k"
    return f"D_total_{min_total_calls}"


def balanced_top_features(
    eligible_indices: Iterable[int],
    scores: Iterable[float],
    species_by_bin: list[str],
    maximum: int,
) -> list[int]:
    """Select high-variance bins while balancing reference species."""

    pairs = [
        (int(index), float(score))
        for index, score in zip(eligible_indices, scores)
        if math.isfinite(float(score))
    ]
    if maximum < 1 or len(pairs) <= maximum:
        return [index for index, _ in sorted(pairs, key=lambda x: -x[1])]
    ranked: dict[str, list[tuple[int, float]]] = {"human": [], "pig": []}
    for index, score in pairs:
        species = species_by_bin[index]
        if species in ranked:
            ranked[species].append((index, score))
    for values in ranked.values():
        values.sort(key=lambda x: -x[1])
    half = maximum // 2
    chosen = ranked["human"][:half] + ranked["pig"][:half]
    chosen_indices = {index for index, _ in chosen}
    leftovers = sorted(
        [pair for pair in pairs if pair[0] not in chosen_indices],
        key=lambda x: -x[1],
    )
    chosen.extend(leftovers[: maximum - len(chosen)])
    chosen.sort(key=lambda x: -x[1])
    return [index for index, _ in chosen]


def dominant_species(human_calls: int, pig_calls: int) -> str:
    if human_calls > pig_calls:
        return "human"
    if pig_calls > human_calls:
        return "pig"
    if human_calls:
        return "tie"
    return "no_signal"


def ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--min-total-calls", type=int, default=20_000)
    parser.add_argument("--min-bin-cells", type=int, default=20)
    parser.add_argument("--min-bin-covered-calls", type=int, default=100)
    parser.add_argument("--max-features-per-context", type=int, default=20_000)
    parser.add_argument("--svd-components", type=int, default=30)
    parser.add_argument("--neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=2607)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--chunk-size", type=int, default=2_000_000)
    parser.add_argument("--residual-clip", type=float, default=10.0)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def scan_retained_entries(
    matrix,
    local_lookup,
    chunk_size: int,
    np,
) -> dict[int, dict[str, object]]:
    fields = ("cell", "bin", "covered", "converted")
    parts: dict[int, dict[str, list[object]]] = {
        code: {field: [] for field in fields} for code in CONTEXTS
    }
    total = int(matrix["cell_index"].shape[0])
    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        global_cell = matrix["cell_index"][start:stop].astype(np.int64)
        if global_cell.size and (
            global_cell.min() < 0 or global_cell.max() >= local_lookup.size
        ):
            raise ValueError("matrix cell index exceeds canonical cell index")
        local_cell = local_lookup[global_cell]
        keep = local_cell >= 0
        if not np.any(keep):
            continue
        context = matrix["context_code"][start:stop]
        bin_index = matrix["bin_index"][start:stop]
        covered = matrix["covered_calls"][start:stop]
        converted = matrix["converted_calls"][start:stop]
        for code in CONTEXTS:
            take = keep & (context == code)
            if not np.any(take):
                continue
            parts[code]["cell"].append(local_cell[take].astype(np.int32))
            parts[code]["bin"].append(bin_index[take].astype(np.int32))
            parts[code]["covered"].append(covered[take].astype(np.float64))
            parts[code]["converted"].append(
                converted[take].astype(np.float64)
            )
    result: dict[int, dict[str, object]] = {}
    for code, values in parts.items():
        result[code] = {}
        for field in fields:
            if values[field]:
                result[code][field] = np.concatenate(values[field])
            else:
                dtype = np.int32 if field in {"cell", "bin"} else np.float64
                result[code][field] = np.asarray([], dtype=dtype)
    return result


def context_embedding(
    arrays: dict[str, object],
    context_name: str,
    n_cells: int,
    n_bins: int,
    species_by_bin: list[str],
    args: argparse.Namespace,
    np,
    sparse,
    TruncatedSVD,
    StandardScaler,
) -> tuple[object, list[dict[str, object]], dict[str, object]]:
    cells = arrays["cell"]
    bins = arrays["bin"]
    covered = arrays["covered"]
    converted = arrays["converted"]
    if len(cells) == 0:
        raise ValueError(f"no retained {context_name} matrix entries")
    observed_cells = np.bincount(bins, minlength=n_bins).astype(np.int64)
    total_covered = np.bincount(
        bins, weights=covered, minlength=n_bins
    ).astype(np.float64)
    total_converted = np.bincount(
        bins, weights=converted, minlength=n_bins
    ).astype(np.float64)
    pooled = np.divide(
        total_converted,
        total_covered,
        out=np.zeros(n_bins, dtype=np.float64),
        where=total_covered > 0,
    )
    entry_pooled = np.clip(pooled[bins], 1e-6, 1 - 1e-6)
    denominator = np.sqrt(
        covered * entry_pooled * (1.0 - entry_pooled) + 1.0
    )
    residual = (converted - covered * entry_pooled) / denominator
    residual = np.clip(residual, -args.residual_clip, args.residual_clip)
    score_sum = np.bincount(
        bins, weights=residual * residual, minlength=n_bins
    )
    variability = np.divide(
        score_sum,
        observed_cells,
        out=np.zeros(n_bins, dtype=np.float64),
        where=observed_cells > 0,
    )
    eligible = np.flatnonzero(
        (observed_cells >= args.min_bin_cells)
        & (total_covered >= args.min_bin_covered_calls)
    )
    if eligible.size < 3:
        raise ValueError(
            f"too few eligible {context_name} bins: {eligible.size}"
        )
    selected = balanced_top_features(
        eligible.tolist(),
        variability[eligible].tolist(),
        species_by_bin,
        args.max_features_per_context,
    )
    selected_array = np.asarray(selected, dtype=np.int64)
    feature_lookup = np.full(n_bins, -1, dtype=np.int32)
    feature_lookup[selected_array] = np.arange(
        selected_array.size, dtype=np.int32
    )
    columns = feature_lookup[bins]
    keep = columns >= 0
    matrix = sparse.coo_matrix(
        (residual[keep], (cells[keep], columns[keep])),
        shape=(n_cells, selected_array.size),
        dtype=np.float32,
    ).tocsr()
    row_norms = np.sqrt(matrix.multiply(matrix).sum(axis=1)).A1
    nonzero_rows = row_norms > 0
    if int(nonzero_rows.sum()) != n_cells:
        empty = n_cells - int(nonzero_rows.sum())
        raise ValueError(
            f"{context_name} selected features leave {empty} cells empty"
        )
    matrix = sparse.diags(1.0 / row_norms).dot(matrix).tocsr()
    components = min(
        args.svd_components, n_cells - 2, selected_array.size - 1
    )
    if components < 2:
        raise ValueError(f"too few dimensions for {context_name} SVD")
    svd = TruncatedSVD(
        n_components=components,
        n_iter=7,
        random_state=args.seed,
    )
    pcs = svd.fit_transform(matrix)
    coverage_by_cell = np.zeros(n_cells, dtype=np.float64)
    np.add.at(coverage_by_cell, cells, covered)
    design = np.column_stack(
        [np.ones(n_cells), np.log1p(coverage_by_cell)]
    )
    coefficients = np.linalg.lstsq(design, pcs, rcond=None)[0]
    pcs = pcs - design @ coefficients
    pcs = StandardScaler().fit_transform(pcs).astype(np.float32)
    feature_rows = []
    for rank, bin_index in enumerate(selected, 1):
        feature_rows.append(
            {
                "context": context_name,
                "feature_rank": rank,
                "bin_index": bin_index,
                "reference_species": species_by_bin[bin_index],
                "observed_cells": int(observed_cells[bin_index]),
                "total_covered_calls": int(total_covered[bin_index]),
                "total_converted_calls": int(total_converted[bin_index]),
                "pooled_signal_fraction": float(pooled[bin_index]),
                "pearson_residual_variability": float(
                    variability[bin_index]
                ),
            }
        )
    metrics = {
        "matrix_entries": len(cells),
        "eligible_bins": int(eligible.size),
        "selected_bins": len(selected),
        "selected_human_bins": sum(
            species_by_bin[index] == "human" for index in selected
        ),
        "selected_pig_bins": sum(
            species_by_bin[index] == "pig" for index in selected
        ),
        "svd_components": components,
        "svd_explained_variance": float(
            svd.explained_variance_ratio_.sum()
        ),
        "coverage_by_cell": coverage_by_cell,
    }
    return pcs, feature_rows, metrics


def run_umap(data, args: argparse.Namespace, umap_module, np):
    neighbors = min(args.neighbors, len(data) - 1)
    if neighbors < 2:
        raise ValueError("too few cells for UMAP")
    model = umap_module.UMAP(
        n_neighbors=neighbors,
        n_components=2,
        min_dist=args.min_dist,
        metric="euclidean",
        random_state=args.seed,
        transform_seed=args.seed,
        n_jobs=1,
        low_memory=False,
        verbose=False,
    )
    return np.asarray(model.fit_transform(data), dtype=np.float32)


def posthoc_embedding_metrics(
    coordinates,
    labels: list[str],
    neighbors: int,
    seed: int,
    np,
    NearestNeighbors,
    silhouette_score,
) -> dict[str, float]:
    valid = np.asarray([label in {"human", "pig"} for label in labels])
    encoded = np.asarray([1 if label == "human" else 0 for label in labels])
    if int(valid.sum()) < 3 or len(set(encoded[valid].tolist())) < 2:
        return {"species_knn_agreement": 0.0, "species_silhouette": 0.0}
    k = min(neighbors, len(coordinates) - 1)
    indices = NearestNeighbors(n_neighbors=k + 1).fit(coordinates).kneighbors(
        return_distance=False
    )[:, 1:]
    agreements = []
    for cell_index in np.flatnonzero(valid):
        eligible_neighbors = indices[cell_index][valid[indices[cell_index]]]
        if eligible_neighbors.size:
            agreements.append(
                float(
                    np.mean(
                        encoded[eligible_neighbors] == encoded[cell_index]
                    )
                )
            )
    valid_count = int(valid.sum())
    sample_size = min(5_000, valid_count)
    silhouette = silhouette_score(
        coordinates[valid],
        encoded[valid],
        sample_size=sample_size if sample_size < valid_count else None,
        random_state=seed,
    )
    return {
        "species_knn_agreement": float(np.mean(agreements)),
        "species_silhouette": float(silhouette),
    }


def save_plots(
    png_path: Path,
    pdf_path: Path,
    embeddings: dict[str, object],
    labels: list[str],
    groups: list[str],
    total_calls,
    threshold: str,
    plt,
    np,
) -> None:
    species_palette = {
        "human": "#2563EB",
        "pig": "#DC2626",
        "tie": "#6B7280",
        "no_signal": "#D1D5DB",
    }
    group_values = sorted(set(groups))
    group_colors = plt.cm.Set2(np.linspace(0.05, 0.95, len(group_values)))
    group_palette = dict(zip(group_values, group_colors))
    figure, axes = plt.subplots(
        3, 3, figsize=(17, 15), constrained_layout=True
    )
    log_calls = np.log10(np.asarray(total_calls, dtype=float) + 1.0)
    for row_index, context in enumerate(("WCG", "GCH", "joint")):
        xy = embeddings[context]
        ax = axes[row_index, 0]
        for label in ("human", "pig", "tie", "no_signal"):
            take = np.asarray([value == label for value in labels])
            if np.any(take):
                ax.scatter(
                    xy[take, 0],
                    xy[take, 1],
                    s=5,
                    alpha=0.70,
                    linewidths=0,
                    color=species_palette[label],
                    label=label,
                )
        ax.legend(markerscale=2.5, frameon=False, fontsize=8)
        ax.set_title(f"{context}: post-hoc DNA species")
        ax = axes[row_index, 1]
        for group in group_values:
            take = np.asarray([value == group for value in groups])
            ax.scatter(
                xy[take, 0],
                xy[take, 1],
                s=5,
                alpha=0.65,
                linewidths=0,
                color=group_palette[group],
                label=group,
            )
        ax.legend(markerscale=2.5, frameon=False, fontsize=8)
        ax.set_title(f"{context}: physical-well set")
        ax = axes[row_index, 2]
        points = ax.scatter(
            xy[:, 0],
            xy[:, 1],
            c=log_calls,
            cmap="viridis",
            s=5,
            alpha=0.75,
            linewidths=0,
        )
        figure.colorbar(
            points, ax=ax, label="log10(total covered calls + 1)"
        )
        ax.set_title(f"{context}: DNA signal depth")
        for column in range(3):
            axes[row_index, column].set_xlabel("UMAP1")
            axes[row_index, column].set_ylabel("UMAP2")
            axes[row_index, column].set_xticks([])
            axes[row_index, column].set_yticks([])
    figure.suptitle(
        f"snMAR DNA 100 kb UMAP ({threshold}; species used only post hoc)",
        fontsize=15,
    )
    for path in (png_path, pdf_path):
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(
            prefix=path.stem + ".", suffix=path.suffix, dir=path.parent
        )
        os.close(descriptor)
        try:
            figure.savefig(temporary, dpi=220, bbox_inches="tight")
            os.replace(temporary, path)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
    plt.close(figure)


def main() -> int:
    args = parse_args()
    if (
        args.min_total_calls < 1
        or args.min_bin_cells < 2
        or args.min_bin_covered_calls < 1
        or args.max_features_per_context < 10
        or args.svd_components < 2
        or args.neighbors < 2
        or not 0 <= args.min_dist <= 1
        or args.threads < 1
        or args.chunk_size < 1
        or args.residual_clip <= 0
    ):
        raise ValueError("invalid UMAP parameter")
    try:
        import h5py
        import matplotlib
        import numpy as np
        import scipy
        import scipy.sparse as sparse
        import sklearn
        import umap
        from sklearn.decomposition import TruncatedSVD
        from sklearn.metrics import silhouette_score
        from sklearn.neighbors import NearestNeighbors
        from sklearn.preprocessing import StandardScaler

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "Step 103 requires h5py, numpy, scipy, scikit-learn, "
            "umap-learn and matplotlib"
        ) from error

    project = args.project
    threshold = threshold_id(args.min_total_calls)
    output_root = args.output_root or (
        project / "10_downstream_analysis/02_dna_100kb_umap" / threshold
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/32_dna_100kb_umap" / threshold
    )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    summary_path = qc_root / "dna_100kb_umap_summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and existing.get("threshold_id") == threshold
        ):
            print(f"[SKIP] validated DNA UMAP exists: {summary_path}")
            return 0
        raise FileExistsError(
            "existing non-current summary; inspect or rerun with --force: "
            f"{summary_path}"
        )

    source_summary_path = (
        project
        / "08_qc_summary/27_dna_100kb/dna_100kb_complete_summary.tsv"
    )
    source = read_metrics(source_summary_path)
    if (
        source.get("status") != "validated"
        or source.get("pipeline_signature") != SOURCE_SIGNATURE
    ):
        raise ValueError("DNA 100 kb complete matrix is not current/validated")
    matrix_path = Path(source["matrix_h5_vds"])
    bin_path = Path(source["bin_index"])
    cell_index_path = Path(source["cell_index"])
    qc_path = (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc.tsv.gz"
    )
    qc_summary_path = (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc_summary.tsv"
    )
    if read_metrics(qc_summary_path).get("status") != "validated":
        raise ValueError("DNA single-cell QC is not validated")

    bin_rows = read_tsv(bin_path)
    if [int(row["bin_index"]) for row in bin_rows] != list(
        range(len(bin_rows))
    ):
        raise ValueError("100 kb bin index is not zero-based/contiguous")
    species_by_bin = [row["species"] for row in bin_rows]
    if not set(species_by_bin).issubset({"human", "pig"}):
        raise ValueError("unexpected species in 100 kb bin index")
    canonical_cells = read_tsv(cell_index_path)
    if [int(row["cell_index"]) for row in canonical_cells] != list(
        range(len(canonical_cells))
    ):
        raise ValueError("cell index is not zero-based/contiguous")
    canonical_by_index = {
        int(row["cell_index"]): row for row in canonical_cells
    }
    qc_by_index: dict[int, dict[str, str]] = {}
    for row in iter_tsv(qc_path):
        index = int(row["cell_index"])
        if index in qc_by_index:
            raise ValueError(f"duplicate DNA QC cell index: {index}")
        total = int(row["WCG_covered_calls"]) + int(row["GCH_covered_calls"])
        if total != int(row["total_signal_covered_calls"]):
            raise ValueError(f"DNA QC total calls do not reconcile: {index}")
        qc_by_index[index] = row
    if set(qc_by_index) != set(canonical_by_index):
        raise ValueError("DNA QC and canonical cell indices differ")
    retained_indices = sorted(
        index
        for index, row in qc_by_index.items()
        if int(row["total_signal_covered_calls"]) >= args.min_total_calls
    )
    if len(retained_indices) < 50:
        raise ValueError(f"too few retained DNA cells: {len(retained_indices)}")
    retained_rows = [qc_by_index[index] for index in retained_indices]
    local_lookup = np.full(len(canonical_cells), -1, dtype=np.int32)
    local_lookup[np.asarray(retained_indices, dtype=np.int64)] = np.arange(
        len(retained_indices), dtype=np.int32
    )

    with h5py.File(matrix_path, "r") as handle:
        if (
            str(handle.attrs.get("status")) != "validated"
            or str(handle.attrs.get("pipeline_signature")) != SOURCE_SIGNATURE
            or int(handle.attrs.get("total_bins", 0)) != len(bin_rows)
            or int(handle.attrs.get("indexed_cells", 0))
            != len(canonical_cells)
        ):
            raise ValueError("invalid DNA 100 kb VDS attributes")
        arrays = scan_retained_entries(
            handle["matrix"], local_lookup, args.chunk_size, np
        )

    expected_context_coverage = {
        1: np.asarray(
            [int(row["WCG_covered_calls"]) for row in retained_rows],
            dtype=np.int64,
        ),
        2: np.asarray(
            [int(row["GCH_covered_calls"]) for row in retained_rows],
            dtype=np.int64,
        ),
    }
    pcs: dict[str, object] = {}
    feature_rows: list[dict[str, object]] = []
    context_metrics: dict[str, dict[str, object]] = {}
    for code, context_name in CONTEXTS.items():
        context_pcs, rows, metrics = context_embedding(
            arrays[code],
            context_name,
            len(retained_indices),
            len(bin_rows),
            species_by_bin,
            args,
            np,
            sparse,
            TruncatedSVD,
            StandardScaler,
        )
        observed_coverage = np.rint(metrics["coverage_by_cell"]).astype(
            np.int64
        )
        if not np.array_equal(
            observed_coverage, expected_context_coverage[code]
        ):
            raise ValueError(
                f"retained {context_name} coverage differs from cell QC"
            )
        pcs[context_name] = context_pcs
        feature_rows.extend(rows)
        context_metrics[context_name] = metrics

    joint_input = np.concatenate(
        [
            pcs["WCG"] / math.sqrt(pcs["WCG"].shape[1]),
            pcs["GCH"] / math.sqrt(pcs["GCH"].shape[1]),
        ],
        axis=1,
    )
    embeddings = {
        "WCG": run_umap(pcs["WCG"], args, umap, np),
        "GCH": run_umap(pcs["GCH"], args, umap, np),
        "joint": run_umap(joint_input, args, umap, np),
    }

    labels = []
    groups = []
    total_calls = []
    for index, row in zip(retained_indices, retained_rows):
        canonical = canonical_by_index[index]
        if (
            row["sample_index"] != canonical["sample_index"]
            or row["CB"] != canonical["CB"]
        ):
            raise ValueError(f"DNA cell metadata differs: {index}")
        human_calls = int(row["human_signal_covered_calls"])
        pig_calls = int(row["pig_signal_covered_calls"])
        label = dominant_species(human_calls, pig_calls)
        if label != row["dominant_reference_signal"]:
            raise ValueError(f"DNA species label differs: {index}")
        labels.append(label)
        groups.append(canonical["group"])
        total_calls.append(int(row["total_signal_covered_calls"]))

    posthoc = {
        name: posthoc_embedding_metrics(
            coordinates,
            labels,
            min(15, args.neighbors),
            args.seed,
            np,
            NearestNeighbors,
            silhouette_score,
        )
        for name, coordinates in embeddings.items()
    }

    coordinates_path = (
        output_root / f"dna_100kb_{threshold}_umap_coordinates.tsv.gz"
    )
    coordinate_rows = []
    for local_index, (global_index, qc_row) in enumerate(
        zip(retained_indices, retained_rows)
    ):
        canonical = canonical_by_index[global_index]
        human_calls = int(qc_row["human_signal_covered_calls"])
        pig_calls = int(qc_row["pig_signal_covered_calls"])
        coordinate_rows.append(
            {
                "umap_cell_index": local_index,
                "cell_index": global_index,
                "well_index": canonical["sample_index"],
                "physical_well": canonical["sample"],
                "well_set": canonical["group"],
                "cell_id": canonical["cell_id"],
                "CB": canonical["CB"],
                "DNA_species_posthoc": labels[local_index],
                "human_signal_covered_calls": human_calls,
                "pig_signal_covered_calls": pig_calls,
                "major_reference_signal_fraction": ratio(
                    max(human_calls, pig_calls), human_calls + pig_calls
                ),
                "WCG_covered_calls": qc_row["WCG_covered_calls"],
                "WCG_methylation_fraction": qc_row[
                    "WCG_methylation_fraction"
                ],
                "GCH_covered_calls": qc_row["GCH_covered_calls"],
                "GCH_accessibility_fraction": qc_row[
                    "GCH_accessibility_fraction"
                ],
                "total_signal_covered_calls": qc_row[
                    "total_signal_covered_calls"
                ],
                "WCG_UMAP1": float(embeddings["WCG"][local_index, 0]),
                "WCG_UMAP2": float(embeddings["WCG"][local_index, 1]),
                "GCH_UMAP1": float(embeddings["GCH"][local_index, 0]),
                "GCH_UMAP2": float(embeddings["GCH"][local_index, 1]),
                "joint_UMAP1": float(embeddings["joint"][local_index, 0]),
                "joint_UMAP2": float(embeddings["joint"][local_index, 1]),
            }
        )
    atomic_tsv(
        coordinates_path,
        list(coordinate_rows[0]),
        coordinate_rows,
        gzip_output=True,
    )
    feature_path = (
        output_root / f"dna_100kb_{threshold}_selected_features.tsv.gz"
    )
    bin_metadata = {int(row["bin_index"]): row for row in bin_rows}
    enriched_features = []
    for row in feature_rows:
        feature = bin_metadata[int(row["bin_index"])]
        enriched_features.append(
            {
                **row,
                "chrom": feature["chrom"],
                "start0": feature["start0"],
                "end0": feature["end0"],
                "terminal_partial_bin": feature["terminal_partial_bin"],
            }
        )
    atomic_tsv(
        feature_path,
        list(enriched_features[0]),
        enriched_features,
        gzip_output=True,
    )

    embedding_path = output_root / f"dna_100kb_{threshold}_embeddings.npz"
    descriptor, temporary = tempfile.mkstemp(
        prefix=embedding_path.stem + ".",
        suffix=".npz",
        dir=output_root,
    )
    os.close(descriptor)
    try:
        np.savez_compressed(
            temporary,
            cell_index=np.asarray(retained_indices, dtype=np.int32),
            WCG_pcs=pcs["WCG"],
            GCH_pcs=pcs["GCH"],
            WCG_umap=embeddings["WCG"],
            GCH_umap=embeddings["GCH"],
            joint_umap=embeddings["joint"],
        )
        os.replace(temporary, embedding_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    png_path = output_root / f"dna_100kb_{threshold}_umap_overview.png"
    pdf_path = output_root / f"dna_100kb_{threshold}_umap_overview.pdf"
    save_plots(
        png_path,
        pdf_path,
        embeddings,
        labels,
        groups,
        total_calls,
        threshold,
        plt,
        np,
    )
    wells = {row["sample_index"] for row in retained_rows}
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("mode", "unsupervised_embedding_species_posthoc_only"),
        ("threshold_id", threshold),
        ("minimum_total_WCG_plus_GCH_covered_calls", args.min_total_calls),
        ("retained_cells", len(retained_indices)),
        ("physical_wells_represented", len(wells)),
        ("human_dominant_cells_posthoc", labels.count("human")),
        ("pig_dominant_cells_posthoc", labels.count("pig")),
        ("tie_cells_posthoc", labels.count("tie")),
        (
            "missing_value_policy",
            "absent_cell_bin_is_NA_encoded_as_zero_residual",
        ),
        (
            "feature_value",
            "coverage_weighted_binomial_Pearson_residual",
        ),
        (
            "depth_correction",
            "log_total_context_coverage_regressed_from_SVD",
        ),
        (
            "reference_feature_balance",
            "equal_human_pig_budget_when_available",
        ),
        ("minimum_cells_per_bin", args.min_bin_cells),
        ("minimum_covered_calls_per_bin", args.min_bin_covered_calls),
        ("maximum_features_per_context", args.max_features_per_context),
        ("requested_SVD_components", args.svd_components),
        ("UMAP_neighbors", min(args.neighbors, len(retained_indices) - 1)),
        ("UMAP_min_dist", args.min_dist),
        ("random_seed", args.seed),
        ("requested_threads", args.threads),
        ("UMAP_n_jobs", 1),
        (
            "UMAP_reproducibility",
            "seeded_single_thread_graph_optimization",
        ),
        ("numpy_version", np.__version__),
        ("h5py_version", h5py.__version__),
        ("scipy_version", scipy.__version__),
        ("sklearn_version", sklearn.__version__),
        ("umap_version", umap.__version__),
        ("source_100kb_summary", source_summary_path),
        ("source_matrix", matrix_path),
        ("source_cell_QC", qc_path),
        ("coordinates", coordinates_path),
        ("selected_features", feature_path),
        ("embeddings_npz", embedding_path),
        ("overview_png", png_path),
        ("overview_pdf", pdf_path),
    ]
    for context_name in ("WCG", "GCH"):
        metrics = context_metrics[context_name]
        for key in (
            "matrix_entries",
            "eligible_bins",
            "selected_bins",
            "selected_human_bins",
            "selected_pig_bins",
            "svd_components",
            "svd_explained_variance",
        ):
            summary_rows.append((f"{context_name}_{key}", metrics[key]))
    for embedding_name in ("WCG", "GCH", "joint"):
        for key, value in posthoc[embedding_name].items():
            summary_rows.append(
                (f"{embedding_name}_{key}_posthoc", value)
            )
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] DNA 100 kb UMAP {threshold}: cells={len(retained_indices)} "
        f"output={output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
