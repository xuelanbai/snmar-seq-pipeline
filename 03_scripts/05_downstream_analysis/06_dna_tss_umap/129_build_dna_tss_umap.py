#!/usr/bin/env python3
"""Step 129: build WCG, GCH and joint TSS-promoter DNA UMAPs.

The primary feature space maps conservative human/pig single-gene promoter
clusters onto shared, case-normalized gene-name axes. Mapping-derived
GM12878/PK15 labels are joined only after fitting and are not independent
validation. Missing cell-feature entries are encoded as zero Pearson residual
(the pooled expectation), never as zero methylation/accessibility.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator


SOURCE_SIGNATURE = "snmar_dna_tss_complete_matrix_20260808"
PIPELINE_SIGNATURE = "snmar_dna_tss_umap_20260808"
CONTEXTS = {1: "WCG", 2: "GCH"}


class EmptySelectedFeatureCells(ValueError):
    """Selected high-variable features leave label-blind cells unsupported."""

    def __init__(self, context_name: str, local_indices) -> None:
        self.context_name = context_name
        self.local_indices = tuple(int(index) for index in local_indices)
        super().__init__(
            f"{context_name} selected TSS features leave "
            f"{len(self.local_indices)} cells empty"
        )


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode + "t", encoding="utf-8", newline="")


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with open_text(path, "r") as handle:
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


def short_count(value: int) -> str:
    if value % 1000 == 0:
        return f"{value // 1000}k"
    return str(value)


def run_id(args: argparse.Namespace) -> str:
    base = (
        f"{args.feature_space}_D{short_count(args.min_total_calls)}_"
        f"TSSctx{short_count(args.min_promoter_calls_per_context)}"
    )
    return f"{base}_{args.run_tag}" if args.run_tag else base


def classify_cell_line(
    human_calls: int, pig_calls: int, minimum_major_fraction: float
) -> tuple[str, float]:
    total = human_calls + pig_calls
    if not total:
        return "ambiguous", 0.0
    major = max(human_calls, pig_calls) / total
    if major < minimum_major_fraction or human_calls == pig_calls:
        return "ambiguous", major
    return ("GM12878" if human_calls > pig_calls else "PK15"), major


def well_design(well_set: str) -> str:
    if well_set == "GM12878":
        return "GM12878_only_wells"
    if well_set == "cellMix":
        return "GM12878_PK15_mix_wells"
    raise ValueError(f"unexpected well_set: {well_set}")


def load_feature_space(
    promoter_rows: list[dict[str, str]],
    ortholog_rows: list[dict[str, str]],
    feature_space: str,
    np,
) -> tuple[object, list[dict[str, object]], list[str]]:
    if [int(row["feature_index"]) for row in promoter_rows] != list(
        range(len(promoter_rows))
    ):
        raise ValueError("promoter index is not contiguous")
    if feature_space == "promoter_cluster":
        mapping = np.arange(len(promoter_rows), dtype=np.int64)
        axes = [
            {
                "axis_index": int(row["feature_index"]),
                "axis_id": row["feature_id"],
                "axis_name": row["gene_names"] or row["feature_id"],
                "reference_species": row["species"],
                "human_feature_index": row["feature_index"]
                if row["species"] == "human"
                else "",
                "pig_feature_index": row["feature_index"]
                if row["species"] == "pig"
                else "",
                "pairing_rule": "none_species_specific_promoter_cluster",
            }
            for row in promoter_rows
        ]
        species = [str(row["species"]) for row in promoter_rows]
        return mapping, axes, species
    if feature_space != "ortholog_single_gene":
        raise ValueError(f"unsupported feature space: {feature_space}")
    ordered = sorted(ortholog_rows, key=lambda row: int(row["ortholog_index"]))
    if [int(row["ortholog_index"]) for row in ordered] != list(range(len(ordered))):
        raise ValueError("ortholog index is not contiguous")
    mapping = np.full(len(promoter_rows), -1, dtype=np.int64)
    axes: list[dict[str, object]] = []
    for row in ordered:
        index = int(row["ortholog_index"])
        human = int(row["human_feature_index"])
        pig = int(row["pig_feature_index"])
        if mapping[human] >= 0 or mapping[pig] >= 0:
            raise ValueError("promoter cluster belongs to multiple ortholog axes")
        mapping[human] = mapping[pig] = index
        axes.append(
            {
                "axis_index": index,
                "axis_id": f"ORTHO_{row['ortholog_key']}",
                "axis_name": row["ortholog_key"],
                "reference_species": "shared_human_pig",
                "human_feature_index": human,
                "human_gene_id": row["human_gene_id"],
                "pig_feature_index": pig,
                "pig_gene_id": row["pig_gene_id"],
                "pairing_rule": row["pairing_rule"],
            }
        )
    return mapping, axes, ["shared_human_pig"] * len(axes)


def collapse_arrays(cell, feature, covered, converted, n_features: int, np):
    if not len(cell):
        return {
            "cell": np.asarray([], dtype=np.int32),
            "feature": np.asarray([], dtype=np.int32),
            "covered": np.asarray([], dtype=np.float64),
            "converted": np.asarray([], dtype=np.float64),
        }
    keys = (
        np.asarray(cell, dtype=np.uint64) * np.uint64(n_features)
        + np.asarray(feature, dtype=np.uint64)
    )
    order = np.argsort(keys, kind="stable")
    keys = keys[order]
    covered = np.asarray(covered, dtype=np.float64)[order]
    converted = np.asarray(converted, dtype=np.float64)[order]
    starts = np.concatenate(
        (
            np.asarray([0], dtype=np.int64),
            np.flatnonzero(keys[1:] != keys[:-1]) + 1,
        )
    )
    unique = keys[starts]
    return {
        "cell": (unique // np.uint64(n_features)).astype(np.int32),
        "feature": (unique % np.uint64(n_features)).astype(np.int32),
        "covered": np.add.reduceat(covered, starts),
        "converted": np.add.reduceat(converted, starts),
    }


def scan_matrix(
    matrix,
    local_lookup,
    feature_map,
    n_features: int,
    chunk_size: int,
    np,
) -> dict[int, dict[str, object]]:
    fields = ("cell", "feature", "covered", "converted")
    parts = {code: {field: [] for field in fields} for code in CONTEXTS}
    total = int(matrix["cell_index"].shape[0])
    for start in range(0, total, chunk_size):
        stop = min(start + chunk_size, total)
        global_cell = matrix["cell_index"][start:stop].astype(np.int64)
        source_feature = matrix["feature_index"][start:stop].astype(np.int64)
        if global_cell.size and (
            global_cell.min() < 0 or global_cell.max() >= local_lookup.size
        ):
            raise ValueError("matrix cell index exceeds canonical index")
        if source_feature.size and (
            source_feature.min() < 0 or source_feature.max() >= feature_map.size
        ):
            raise ValueError("matrix feature index exceeds promoter index")
        local_cell = local_lookup[global_cell]
        axis = feature_map[source_feature]
        keep = (local_cell >= 0) & (axis >= 0)
        if not np.any(keep):
            continue
        context = matrix["context_code"][start:stop]
        covered = matrix["covered_calls"][start:stop]
        converted = matrix["converted_calls"][start:stop]
        for code in CONTEXTS:
            take = keep & (context == code)
            if not np.any(take):
                continue
            parts[code]["cell"].append(local_cell[take].astype(np.int32))
            parts[code]["feature"].append(axis[take].astype(np.int32))
            parts[code]["covered"].append(covered[take].astype(np.float64))
            parts[code]["converted"].append(converted[take].astype(np.float64))
    result: dict[int, dict[str, object]] = {}
    for code in CONTEXTS:
        values = parts[code]
        if values["cell"]:
            result[code] = collapse_arrays(
                np.concatenate(values["cell"]),
                np.concatenate(values["feature"]),
                np.concatenate(values["covered"]),
                np.concatenate(values["converted"]),
                n_features,
                np,
            )
        else:
            result[code] = collapse_arrays([], [], [], [], n_features, np)
    return result


def remap_cells(arrays, remap, np):
    old = arrays["cell"]
    new = remap[old]
    keep = new >= 0
    return {
        "cell": new[keep].astype(np.int32),
        "feature": arrays["feature"][keep],
        "covered": arrays["covered"][keep],
        "converted": arrays["converted"][keep],
    }


def supported_cell_mask(
    arrays,
    n_cells: int,
    n_features: int,
    min_feature_cells: int,
    min_feature_covered_calls: int,
    np,
):
    observed_cells = np.bincount(
        arrays["feature"], minlength=n_features
    ).astype(np.int64)
    total_covered = np.bincount(
        arrays["feature"],
        weights=arrays["covered"],
        minlength=n_features,
    )
    eligible = (
        (observed_cells >= min_feature_cells)
        & (total_covered >= min_feature_covered_calls)
    )
    result = np.zeros(n_cells, dtype=bool)
    take = eligible[arrays["feature"]]
    if np.any(take):
        result[np.unique(arrays["cell"][take])] = True
    return result


def choose_features(
    eligible,
    variability,
    species_by_axis: list[str],
    maximum: int,
    feature_space: str,
) -> list[int]:
    ranked = sorted(
        (int(index) for index in eligible),
        key=lambda index: (-float(variability[index]), index),
    )
    if len(ranked) <= maximum or feature_space == "ortholog_single_gene":
        return ranked[:maximum]
    by_species = {
        species: [index for index in ranked if species_by_axis[index] == species]
        for species in ("human", "pig")
    }
    half = maximum // 2
    selected = by_species["human"][:half] + by_species["pig"][:half]
    chosen = set(selected)
    selected.extend(index for index in ranked if index not in chosen)
    return selected[:maximum]


def context_embedding(
    arrays,
    context_name: str,
    n_cells: int,
    n_features: int,
    axes: list[dict[str, object]],
    species_by_axis: list[str],
    args: argparse.Namespace,
    np,
    sparse,
    TruncatedSVD,
    StandardScaler,
):
    cells = arrays["cell"]
    features = arrays["feature"]
    covered = arrays["covered"]
    converted = arrays["converted"]
    if not len(cells):
        raise ValueError(f"no retained {context_name} promoter entries")
    observed_cells = np.bincount(features, minlength=n_features).astype(np.int64)
    total_covered = np.bincount(
        features, weights=covered, minlength=n_features
    ).astype(np.float64)
    total_converted = np.bincount(
        features, weights=converted, minlength=n_features
    ).astype(np.float64)
    pooled = np.divide(
        total_converted,
        total_covered,
        out=np.zeros(n_features, dtype=np.float64),
        where=total_covered > 0,
    )
    expected = np.clip(pooled[features], 1e-6, 1 - 1e-6)
    residual = (converted - covered * expected) / np.sqrt(
        covered * expected * (1.0 - expected) + 1.0
    )
    residual = np.clip(residual, -args.residual_clip, args.residual_clip)
    score_sum = np.bincount(
        features, weights=residual * residual, minlength=n_features
    )
    variability = np.divide(
        score_sum,
        observed_cells,
        out=np.zeros(n_features, dtype=np.float64),
        where=observed_cells > 0,
    )
    eligible = np.flatnonzero(
        (observed_cells >= args.min_feature_cells)
        & (total_covered >= args.min_feature_covered_calls)
    )
    if eligible.size < 3:
        raise ValueError(f"too few eligible {context_name} TSS features")
    selected = choose_features(
        eligible,
        variability,
        species_by_axis,
        args.max_features_per_context,
        args.feature_space,
    )
    selected_array = np.asarray(selected, dtype=np.int64)
    lookup = np.full(n_features, -1, dtype=np.int32)
    lookup[selected_array] = np.arange(len(selected), dtype=np.int32)
    columns = lookup[features]
    take = columns >= 0
    matrix = sparse.coo_matrix(
        (residual[take], (cells[take], columns[take])),
        shape=(n_cells, len(selected)),
        dtype=np.float32,
    ).tocsr()
    row_norms = np.sqrt(matrix.multiply(matrix).sum(axis=1)).A1
    empty_cells = np.flatnonzero(row_norms <= 0)
    if empty_cells.size:
        if args.drop_empty_selected_feature_cells:
            raise EmptySelectedFeatureCells(context_name, empty_cells)
        raise ValueError(
            f"{context_name} selected TSS features leave "
            f"{empty_cells.size} cells empty; "
            "lower feature thresholds or raise promoter coverage threshold"
        )
    matrix = sparse.diags(1.0 / row_norms).dot(matrix).tocsr()
    components = min(args.svd_components, n_cells - 2, len(selected) - 1)
    if components < 2:
        raise ValueError(f"too few dimensions for {context_name} SVD")
    svd = TruncatedSVD(
        n_components=components, n_iter=7, random_state=args.seed
    )
    pcs = svd.fit_transform(matrix)
    coverage_by_cell = np.bincount(
        cells, weights=covered, minlength=n_cells
    ).astype(np.float64)
    design = np.column_stack([np.ones(n_cells), np.log1p(coverage_by_cell)])
    coefficients = np.linalg.lstsq(design, pcs, rcond=None)[0]
    pcs = pcs - design @ coefficients
    pcs = StandardScaler().fit_transform(pcs).astype(np.float32)
    selected_rows = []
    for rank, axis_index in enumerate(selected, 1):
        selected_rows.append(
            {
                "context": context_name,
                "feature_rank": rank,
                **axes[axis_index],
                "observed_cells": int(observed_cells[axis_index]),
                "total_covered_calls": int(total_covered[axis_index]),
                "total_converted_calls": int(total_converted[axis_index]),
                "pooled_signal_fraction": float(pooled[axis_index]),
                "pearson_residual_variability": float(variability[axis_index]),
            }
        )
    return pcs, selected_rows, {
        "matrix_entries": len(cells),
        "eligible_features": int(eligible.size),
        "selected_features": len(selected),
        "svd_components": components,
        "svd_explained_variance": float(svd.explained_variance_ratio_.sum()),
        "coverage_by_cell": coverage_by_cell,
    }


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


def embedding_metrics(
    coordinates,
    labels: list[str],
    neighbors: int,
    seed: int,
    np,
    NearestNeighbors,
    silhouette_score,
) -> dict[str, float]:
    valid = np.asarray([label in {"GM12878", "PK15"} for label in labels])
    if int(valid.sum()) < 3 or len(set(np.asarray(labels)[valid])) < 2:
        return {"cell_line_knn_agreement": 0.0, "cell_line_silhouette": 0.0}
    k = min(neighbors, len(coordinates) - 1)
    graph = NearestNeighbors(n_neighbors=k + 1).fit(coordinates)
    indices = graph.kneighbors(coordinates, return_distance=False)[:, 1:]
    agreements = []
    label_array = np.asarray(labels)
    for index in np.flatnonzero(valid):
        eligible = indices[index][valid[indices[index]]]
        if eligible.size:
            agreements.append(float(np.mean(label_array[eligible] == label_array[index])))
    valid_count = int(valid.sum())
    sample_size = min(5000, valid_count)
    silhouette = silhouette_score(
        coordinates[valid],
        label_array[valid],
        sample_size=sample_size if sample_size < valid_count else None,
        random_state=seed,
    )
    return {
        "cell_line_knn_agreement": float(np.mean(agreements)),
        "cell_line_silhouette": float(silhouette),
    }


def scatter_categories(ax, coordinates, labels, palette, order):
    for label in order:
        take = [value == label for value in labels]
        if any(take):
            selected = coordinates[take]
            ax.scatter(
                selected[:, 0],
                selected[:, 1],
                s=5,
                alpha=0.70,
                linewidths=0,
                color=palette[label],
                label=label,
            )
    ax.legend(markerscale=2.5, frameon=False, fontsize=8)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")


def save_plot(
    png_path: Path,
    pdf_path: Path,
    embeddings,
    labels,
    designs,
    confidence,
    run_name: str,
    feature_space: str,
    plt,
):
    label_palette = {
        "GM12878": "#377EB8",
        "PK15": "#E41A1C",
        "ambiguous": "#9E9E9E",
    }
    design_palette = {
        "GM12878_only_wells": "#1B9E77",
        "GM12878_PK15_mix_wells": "#BDBDBD",
    }
    figure, axes = plt.subplots(3, 3, figsize=(18, 16), constrained_layout=True)
    confidence_plot = None
    for row_index, name in enumerate(("WCG", "GCH", "joint")):
        coordinates = embeddings[name]
        scatter_categories(
            axes[row_index, 0],
            coordinates,
            labels,
            label_palette,
            ("GM12878", "PK15", "ambiguous"),
        )
        axes[row_index, 0].set_title(f"{name}: DNA mapping-derived cell line")
        scatter_categories(
            axes[row_index, 1],
            coordinates,
            designs,
            design_palette,
            ("GM12878_only_wells", "GM12878_PK15_mix_wells"),
        )
        axes[row_index, 1].set_title(f"{name}: physical-well design")
        confidence_plot = axes[row_index, 2].scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            c=confidence,
            cmap="viridis",
            vmin=0.5,
            vmax=1.0,
            s=5,
            alpha=0.75,
            linewidths=0,
        )
        axes[row_index, 2].set_title(f"{name}: mapping confidence")
        axes[row_index, 2].set_xlabel("UMAP1")
        axes[row_index, 2].set_ylabel("UMAP2")
    if confidence_plot is not None:
        figure.colorbar(
            confidence_plot,
            ax=axes[:, 2].tolist(),
            label="major human/pig reference fraction",
            shrink=0.78,
        )
    figure.suptitle(
        "snMAR DNA TSS±1.5 kb UMAP "
        f"({run_name}; {feature_space}; labels post hoc only)",
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument(
        "--feature-space",
        choices=("ortholog_single_gene", "promoter_cluster"),
        default="ortholog_single_gene",
    )
    parser.add_argument("--min-total-calls", type=int, default=20_000)
    parser.add_argument("--min-promoter-calls-per-context", type=int, default=25)
    parser.add_argument("--min-feature-cells", type=int, default=5)
    parser.add_argument("--min-feature-covered-calls", type=int, default=20)
    parser.add_argument("--max-features-per-context", type=int, default=20_000)
    parser.add_argument("--svd-components", type=int, default=30)
    parser.add_argument("--neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.30)
    parser.add_argument("--minimum-major-fraction", type=float, default=0.90)
    parser.add_argument("--residual-clip", type=float, default=10.0)
    parser.add_argument("--chunk-size", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=2607)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument(
        "--run-tag",
        default="",
        help=(
            "Optional filesystem-safe suffix used to keep sensitivity runs "
            "separate from the validated default analysis."
        ),
    )
    parser.add_argument(
        "--drop-empty-selected-feature-cells",
        action="store_true",
        help=(
            "Drop, without consulting mapping labels, cells that have no calls "
            "in the selected high-variable features."
        ),
    )
    parser.add_argument(
        "--max-empty-selected-feature-cell-fraction",
        type=float,
        default=0.01,
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.min_total_calls < 1
        or args.min_promoter_calls_per_context < 1
        or args.min_feature_cells < 2
        or args.min_feature_covered_calls < 1
        or args.max_features_per_context < 10
        or args.svd_components < 2
        or args.neighbors < 2
        or not 0 <= args.min_dist <= 1
        or not 0.5 < args.minimum_major_fraction <= 1
        or args.residual_clip <= 0
        or args.chunk_size < 1
        or args.threads < 1
        or (args.run_tag and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.run_tag) is None)
        or not 0 <= args.max_empty_selected_feature_cell_fraction < 1
    ):
        raise ValueError("invalid TSS UMAP parameter")
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
            "Step 129 requires h5py, numpy, scipy, scikit-learn, "
            "umap-learn and matplotlib"
        ) from error
    project = args.project
    name = run_id(args)
    output_root = args.output_root or (
        project / "10_downstream_analysis/09_dna_tss_umap" / name
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/39_dna_tss_umap" / name
    )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    summary_path = qc_root / "dna_tss_umap_summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and existing.get("run_id") == name
        ):
            print(f"[SKIP] validated TSS UMAP exists: {summary_path}")
            return 0
        raise FileExistsError(f"inspect or rerun with --force: {summary_path}")
    source_summary_path = (
        project / "08_qc_summary/38_dna_tss_promoters/dna_tss_complete_summary.tsv"
    )
    source = read_metrics(source_summary_path)
    if (
        source.get("status") != "validated"
        or source.get("pipeline_signature") != SOURCE_SIGNATURE
    ):
        raise ValueError("complete TSS matrix is not current/validated")
    matrix_path = Path(source["matrix_h5_vds"])
    promoter_rows = read_tsv(Path(source["promoter_index"]))
    ortholog_rows = read_tsv(Path(source["ortholog_index"]))
    canonical_cells = read_tsv(Path(source["cell_index"]))
    if [int(row["cell_index"]) for row in canonical_cells] != list(
        range(len(canonical_cells))
    ):
        raise ValueError("canonical cell index is not contiguous")
    feature_map, axes, species_by_axis = load_feature_space(
        promoter_rows, ortholog_rows, args.feature_space, np
    )
    n_features = len(axes)
    qc_path = (
        project / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc.tsv.gz"
    )
    qc_summary = (
        project / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc_summary.tsv"
    )
    if read_metrics(qc_summary).get("status") != "validated":
        raise ValueError("DNA single-cell QC is not validated")
    qc_by_index: dict[int, dict[str, str]] = {}
    for row in iter_tsv(qc_path):
        index = int(row["cell_index"])
        if index in qc_by_index:
            raise ValueError(f"duplicate DNA QC cell index: {index}")
        qc_by_index[index] = row
    if set(qc_by_index) != set(range(len(canonical_cells))):
        raise ValueError("DNA QC and canonical cell indices differ")
    initial_indices = [
        index
        for index in range(len(canonical_cells))
        if int(qc_by_index[index]["total_signal_covered_calls"])
        >= args.min_total_calls
    ]
    if len(initial_indices) < 50:
        raise ValueError(f"too few cells pass whole-genome DNA depth: {len(initial_indices)}")
    initial_lookup = np.full(len(canonical_cells), -1, dtype=np.int32)
    initial_lookup[np.asarray(initial_indices, dtype=np.int64)] = np.arange(
        len(initial_indices), dtype=np.int32
    )
    with h5py.File(matrix_path, "r") as handle:
        if (
            str(handle.attrs.get("status")) != "validated"
            or str(handle.attrs.get("pipeline_signature")) != SOURCE_SIGNATURE
            or int(handle.attrs.get("indexed_cells", 0)) != len(canonical_cells)
            or int(handle.attrs.get("total_features", 0)) != len(promoter_rows)
        ):
            raise ValueError("invalid complete TSS HDF5 attributes")
        arrays = scan_matrix(
            handle["matrix"],
            initial_lookup,
            feature_map,
            n_features,
            args.chunk_size,
            np,
        )
    initial_coverage = {
        code: np.bincount(
            arrays[code]["cell"],
            weights=arrays[code]["covered"],
            minlength=len(initial_indices),
        ).astype(np.int64)
        for code in CONTEXTS
    }
    depth_keep = (
        (initial_coverage[1] >= args.min_promoter_calls_per_context)
        & (initial_coverage[2] >= args.min_promoter_calls_per_context)
    )
    feature_support = {
        code: supported_cell_mask(
            arrays[code],
            len(initial_indices),
            n_features,
            args.min_feature_cells,
            args.min_feature_covered_calls,
            np,
        )
        for code in CONTEXTS
    }
    keep_initial = depth_keep & feature_support[1] & feature_support[2]
    retained_indices = [
        index for index, keep in zip(initial_indices, keep_initial) if keep
    ]
    if len(retained_indices) < 50:
        raise ValueError(
            f"too few cells pass shared-promoter context depth: {len(retained_indices)}"
        )
    remap = np.full(len(initial_indices), -1, dtype=np.int32)
    remap[np.flatnonzero(keep_initial)] = np.arange(
        len(retained_indices), dtype=np.int32
    )
    arrays = {code: remap_cells(arrays[code], remap, np) for code in CONTEXTS}
    cells_before_selected_feature_filter = len(retained_indices)
    selected_feature_drop_counts: Counter[str] = Counter()
    selected_feature_iterations = 0
    while True:
        selected_feature_iterations += 1
        if selected_feature_iterations > 20:
            raise RuntimeError("high-variable selected-feature cell filtering did not converge")
        pcs = {}
        selected_rows: list[dict[str, object]] = []
        context_metrics = {}
        try:
            for code, context_name in CONTEXTS.items():
                values, rows, metrics = context_embedding(
                    arrays[code],
                    context_name,
                    len(retained_indices),
                    n_features,
                    axes,
                    species_by_axis,
                    args,
                    np,
                    sparse,
                    TruncatedSVD,
                    StandardScaler,
                )
                pcs[context_name] = values
                selected_rows.extend(rows)
                context_metrics[context_name] = metrics
        except EmptySelectedFeatureCells as error:
            prospective_dropped = (
                cells_before_selected_feature_filter
                - len(retained_indices)
                + len(error.local_indices)
            )
            maximum_dropped = int(
                math.floor(
                    cells_before_selected_feature_filter
                    * args.max_empty_selected_feature_cell_fraction
                )
            )
            if prospective_dropped > maximum_dropped:
                raise ValueError(
                    "selected high-variable features would drop "
                    f"{prospective_dropped}/{cells_before_selected_feature_filter} cells, "
                    "exceeding the configured label-blind limit"
                ) from error
            keep_selected = np.ones(len(retained_indices), dtype=bool)
            keep_selected[np.asarray(error.local_indices, dtype=np.int64)] = False
            retained_indices = [
                index
                for index, keep in zip(retained_indices, keep_selected)
                if keep
            ]
            if len(retained_indices) < 50:
                raise ValueError("too few cells remain after selected-feature support filtering")
            selected_remap = np.full(len(keep_selected), -1, dtype=np.int32)
            selected_remap[np.flatnonzero(keep_selected)] = np.arange(
                len(retained_indices), dtype=np.int32
            )
            arrays = {
                code: remap_cells(values, selected_remap, np)
                for code, values in arrays.items()
            }
            selected_feature_drop_counts[error.context_name] += len(error.local_indices)
            continue
        break
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
    labels: list[str] = []
    confidence: list[float] = []
    designs: list[str] = []
    coordinate_rows: list[dict[str, object]] = []
    for local_index, global_index in enumerate(retained_indices):
        qc_row = qc_by_index[global_index]
        canonical = canonical_cells[global_index]
        if (
            qc_row["sample_index"] != canonical["sample_index"]
            or qc_row["CB"] != canonical["CB"]
        ):
            raise ValueError(f"DNA cell metadata differs: {global_index}")
        human = int(qc_row["human_signal_covered_calls"])
        pig = int(qc_row["pig_signal_covered_calls"])
        label, major = classify_cell_line(
            human, pig, args.minimum_major_fraction
        )
        design = well_design(canonical["group"])
        labels.append(label)
        confidence.append(major)
        designs.append(design)
        coordinate_rows.append(
            {
                "umap_cell_index": local_index,
                "cell_index": global_index,
                "well_index": canonical["sample_index"],
                "physical_well": canonical["sample"],
                "well_set": canonical["group"],
                "well_design": design,
                "cell_id": canonical["cell_id"],
                "CB": canonical["CB"],
                "DNA_mapping_derived_cell_line": label,
                "mapping_major_reference_fraction": major,
                "human_signal_covered_calls": human,
                "pig_signal_covered_calls": pig,
                "whole_genome_signal_covered_calls": qc_row["total_signal_covered_calls"],
                "shared_promoter_WCG_covered_calls": int(
                    context_metrics["WCG"]["coverage_by_cell"][local_index]
                ),
                "shared_promoter_GCH_covered_calls": int(
                    context_metrics["GCH"]["coverage_by_cell"][local_index]
                ),
                "WCG_UMAP1": float(embeddings["WCG"][local_index, 0]),
                "WCG_UMAP2": float(embeddings["WCG"][local_index, 1]),
                "GCH_UMAP1": float(embeddings["GCH"][local_index, 0]),
                "GCH_UMAP2": float(embeddings["GCH"][local_index, 1]),
                "joint_UMAP1": float(embeddings["joint"][local_index, 0]),
                "joint_UMAP2": float(embeddings["joint"][local_index, 1]),
            }
        )
    coordinates_path = output_root / "dna_tss_umap_coordinates.tsv.gz"
    selected_path = output_root / "dna_tss_selected_features.tsv.gz"
    atomic_tsv(
        coordinates_path,
        list(coordinate_rows[0]),
        coordinate_rows,
        gzip_output=True,
    )
    atomic_tsv(
        selected_path,
        list(selected_rows[0]),
        selected_rows,
        gzip_output=True,
    )
    embedding_path = output_root / "dna_tss_embeddings.npz"
    descriptor, temporary = tempfile.mkstemp(
        prefix=embedding_path.stem + ".", suffix=".npz", dir=output_root
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
    png_path = output_root / "dna_tss_umap_overview.png"
    pdf_path = output_root / "dna_tss_umap_overview.pdf"
    save_plot(
        png_path,
        pdf_path,
        embeddings,
        labels,
        designs,
        confidence,
        name,
        args.feature_space,
        plt,
    )
    posthoc = {
        embedding: embedding_metrics(
            coordinates,
            labels,
            min(15, args.neighbors),
            args.seed,
            np,
            NearestNeighbors,
            silhouette_score,
        )
        for embedding, coordinates in embeddings.items()
    }
    counts = Counter(labels)
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("run_id", name),
        ("run_tag", args.run_tag or "none"),
        ("mode", "unsupervised_TSS_embedding_mapping_labels_posthoc_only"),
        ("feature_space", args.feature_space),
        ("feature_space_axes", n_features),
        (
            "ortholog_scope",
            "gene_name_proxy_not_authoritative_orthology"
            if args.feature_space == "ortholog_single_gene"
            else "not_applicable",
        ),
        ("minimum_whole_genome_WCG_plus_GCH_calls", args.min_total_calls),
        (
            "minimum_shared_promoter_calls_per_context",
            args.min_promoter_calls_per_context,
        ),
        ("cells_passing_whole_genome_depth", len(initial_indices)),
        ("cells_passing_promoter_context_depth", int(depth_keep.sum())),
        (
            "cells_with_eligible_features_in_both_contexts",
            int((feature_support[1] & feature_support[2]).sum()),
        ),
        ("retained_cells", len(retained_indices)),
        ("cells_before_selected_feature_support_filter", cells_before_selected_feature_filter),
        (
            "cells_dropped_no_selected_feature_any_context",
            cells_before_selected_feature_filter - len(retained_indices),
        ),
        (
            "cells_dropped_no_selected_WCG_feature",
            selected_feature_drop_counts.get("WCG", 0),
        ),
        (
            "cells_dropped_no_selected_GCH_feature",
            selected_feature_drop_counts.get("GCH", 0),
        ),
        (
            "selected_feature_empty_cell_policy",
            "drop_label_blind" if args.drop_empty_selected_feature_cells else "error",
        ),
        (
            "maximum_empty_selected_feature_cell_fraction",
            args.max_empty_selected_feature_cell_fraction,
        ),
        ("physical_wells_represented", len({row["well_index"] for row in coordinate_rows})),
        ("GM12878_cells_mapping_derived", counts["GM12878"]),
        ("PK15_cells_mapping_derived", counts["PK15"]),
        ("ambiguous_cells", counts["ambiguous"]),
        ("minimum_mapping_major_fraction", args.minimum_major_fraction),
        ("independent_cell_line_validation", "false"),
        ("missing_value_policy", "absent_cell_feature_is_NA_encoded_as_zero_residual"),
        ("feature_value", "coverage_weighted_binomial_Pearson_residual"),
        ("depth_correction", "log_context_promoter_coverage_regressed_from_SVD"),
        ("minimum_cells_per_feature", args.min_feature_cells),
        ("minimum_covered_calls_per_feature", args.min_feature_covered_calls),
        ("maximum_features_per_context", args.max_features_per_context),
        ("requested_SVD_components", args.svd_components),
        ("UMAP_neighbors", min(args.neighbors, len(retained_indices) - 1)),
        ("UMAP_min_dist", args.min_dist),
        ("random_seed", args.seed),
        ("requested_threads", args.threads),
        ("UMAP_n_jobs", 1),
        ("source_TSS_summary", source_summary_path),
        ("source_matrix", matrix_path),
        ("source_cell_QC", qc_path),
        ("coordinates", coordinates_path),
        ("selected_features", selected_path),
        ("embeddings_npz", embedding_path),
        ("overview_png", png_path),
        ("overview_pdf", pdf_path),
        ("numpy_version", np.__version__),
        ("h5py_version", h5py.__version__),
        ("scipy_version", scipy.__version__),
        ("sklearn_version", sklearn.__version__),
        ("umap_version", umap.__version__),
    ]
    for context_name in ("WCG", "GCH"):
        metrics = context_metrics[context_name]
        for key in (
            "matrix_entries",
            "eligible_features",
            "selected_features",
            "svd_components",
            "svd_explained_variance",
        ):
            summary_rows.append((f"{context_name}_{key}", metrics[key]))
    for embedding in ("WCG", "GCH", "joint"):
        for key, value in posthoc[embedding].items():
            summary_rows.append((f"{embedding}_{key}_posthoc", value))
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] DNA TSS UMAP {name}: cells={len(retained_indices)} "
        f"features={n_features} output={output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
