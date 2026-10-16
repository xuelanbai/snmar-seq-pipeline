#!/usr/bin/env python3
"""Step 140: construct label-blind WCG matrices on shared 100 kb axes.

Four matrices with identical cells/features are exported for Seurat:

``ratio_paper``
    converted/covered WCG, with cell-bin coverage below the requested call
    threshold imputed to the pooled orthobin mean (published-style branch).
``beta_posterior``
    beta-binomial shrinkage divided by each cell's global WCG mean; an
    uncovered entry equals the neutral value 1, not zero methylation.
``pearson_residual``
    coverage-aware residual; an uncovered entry equals expectation (zero).
``coverage_control``
    log1p WCG covered calls, used only to detect coverage-driven structure.

Human/pig mapping labels are appended only after the matrices and variable
features have been selected.  They never enter filtering or feature ranking.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator


INDEX_SIGNATURE = "snmar_dna_ortholog_100kb_index_20260825"
WELL_SIGNATURE = "snmar_dna_ortholog_100kb_well_20260825"
PIPELINE_SIGNATURE = "snmar_dna_ortholog_100kb_seurat_input_20260825"
METHODS = (
    "ratio_paper",
    "beta_posterior",
    "pearson_residual",
    "coverage_control",
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
        ({"metric": key, "value": value} for key, value in rows),
    )


def short_count(value: int) -> str:
    if value % 1000 == 0:
        return f"{value // 1000}k"
    return str(value)


def run_id(args: argparse.Namespace) -> str:
    return (
        f"Dtotal{short_count(args.min_total_calls)}_"
        f"orthoWCG{short_count(args.min_orthobin_wcg_calls)}_"
        f"call{args.min_cell_bin_calls}_top{args.max_features}"
    )


def classify_cell_line(
    human_calls: int, pig_calls: int, minimum_major_fraction: float
) -> tuple[str, float]:
    total = human_calls + pig_calls
    if not total:
        return "ambiguous", 0.0
    confidence = max(human_calls, pig_calls) / total
    if confidence < minimum_major_fraction or human_calls == pig_calls:
        return "ambiguous", confidence
    return ("GM12878" if human_calls > pig_calls else "PK15"), confidence


def well_design(value: str) -> str:
    if value == "GM12878":
        return "GM12878_only_wells"
    if value == "cellMix":
        return "GM12878_PK15_mix_wells"
    return value


def collapse_arrays(cell, feature, covered, converted, n_features: int, np):
    if not len(cell):
        return {
            "cell": np.asarray([], dtype=np.int32),
            "feature": np.asarray([], dtype=np.int32),
            "covered": np.asarray([], dtype=np.float64),
            "converted": np.asarray([], dtype=np.float64),
        }
    keys = np.asarray(cell, dtype=np.uint64) * np.uint64(n_features) + np.asarray(
        feature, dtype=np.uint64
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


def load_feature_axis(feature_rows, axis_rows, np):
    if [int(row["feature_index"]) for row in feature_rows] != list(
        range(len(feature_rows))
    ):
        raise ValueError("physical feature index is not contiguous")
    if [int(row["orthobin_index"]) for row in axis_rows] != list(
        range(len(axis_rows))
    ):
        raise ValueError("orthobin axis is not contiguous")
    mapping = np.full(len(feature_rows), -1, dtype=np.int64)
    for row in feature_rows:
        feature = int(row["feature_index"])
        axis = int(row["orthobin_index"])
        if not 0 <= axis < len(axis_rows):
            raise ValueError("physical feature has an invalid orthobin index")
        mapping[feature] = axis
    if np.any(mapping < 0):
        raise ValueError("unmapped physical feature")
    for row in axis_rows:
        axis = int(row["orthobin_index"])
        human = int(row["human_feature_index"])
        pig = int(row["pig_feature_index"])
        if mapping[human] != axis or mapping[pig] != axis:
            raise ValueError("axis-to-physical feature mapping differs")
    return mapping


def load_well_arrays(
    project: Path,
    canonical_cells: list[dict[str, str]],
    initial_lookup,
    feature_map,
    feature_species,
    n_features: int,
    h5py,
    np,
):
    cell_parts = []
    feature_parts = []
    covered_parts = []
    converted_parts = []
    side_covered = {
        "human": np.zeros(n_features, dtype=np.float64),
        "pig": np.zeros(n_features, dtype=np.float64),
    }
    side_cells = {
        "human": np.zeros(n_features, dtype=np.int64),
        "pig": np.zeros(n_features, dtype=np.int64),
    }
    qc_root = project / "08_qc_summary/40_dna_ortholog_100kb/well_matrices"
    well_metadata: dict[int, tuple[str, str]] = {}
    for row in canonical_cells:
        index = int(row["sample_index"])
        value = (row["sample"], row["group"])
        previous = well_metadata.setdefault(index, value)
        if previous != value:
            raise ValueError("canonical metadata differs within a physical well")
    if sorted(well_metadata) != list(range(1, 97)):
        raise ValueError("canonical cell index does not contain 96 physical wells")
    for well_index in range(1, 97):
        sample, group = well_metadata[well_index]
        summary_path = qc_root / f"{sample}.dna_ortholog_100kb.summary.tsv"
        summary = read_metrics(summary_path)
        matrix_path = Path(summary.get("output_h5", ""))
        if (
            summary.get("status") != "validated"
            or summary.get("pipeline_signature") != WELL_SIGNATURE
            or int(summary.get("well_index", "0")) != well_index
            or summary.get("physical_well") != sample
            or summary.get("well_set") != group
            or not matrix_path.is_file()
        ):
            raise ValueError(f"invalid orthobin well summary: {summary_path}")
        with h5py.File(matrix_path, "r") as source:
            if (
                str(source.attrs.get("status")) != "validated"
                or str(source.attrs.get("pipeline_signature")) != WELL_SIGNATURE
                or int(source.attrs.get("well_index", 0)) != well_index
                or int(source.attrs.get("total_physical_features", 0))
                != feature_map.size
            ):
                raise ValueError(f"invalid orthobin well HDF5: {matrix_path}")
            matrix = source["matrix"]
            global_cell = matrix["cell_index"][:].astype(np.int64)
            source_feature = matrix["feature_index"][:].astype(np.int64)
            context = matrix["context_code"][:]
            covered = matrix["covered_calls"][:].astype(np.float64)
            converted = matrix["converted_calls"][:].astype(np.float64)
        if global_cell.size and (
            global_cell.min() < 0 or global_cell.max() >= initial_lookup.size
        ):
            raise ValueError("well matrix cell index exceeds canonical cells")
        if source_feature.size and (
            source_feature.min() < 0 or source_feature.max() >= feature_map.size
        ):
            raise ValueError("well matrix feature index exceeds index")
        local_cell = initial_lookup[global_cell]
        axis = feature_map[source_feature]
        keep = (local_cell >= 0) & (axis >= 0) & (context == 1)
        if np.any(keep):
            cell_parts.append(local_cell[keep].astype(np.int32))
            feature_parts.append(axis[keep].astype(np.int32))
            covered_parts.append(covered[keep])
            converted_parts.append(converted[keep])
            for species, code in (("human", 0), ("pig", 1)):
                side_take = keep & (feature_species[source_feature] == code)
                if np.any(side_take):
                    side_covered[species] += np.bincount(
                        axis[side_take],
                        weights=covered[side_take],
                        minlength=n_features,
                    )
                    side_cells[species] += np.bincount(
                        axis[side_take], minlength=n_features
                    ).astype(np.int64)
    if not cell_parts:
        raise ValueError("no WCG entries remain after initial cell filtering")
    arrays = collapse_arrays(
        np.concatenate(cell_parts),
        np.concatenate(feature_parts),
        np.concatenate(covered_parts),
        np.concatenate(converted_parts),
        n_features,
        np,
    )
    return arrays, {"covered": side_covered, "cells": side_cells}


def remap_cells(arrays, keep_mask, np):
    remap = np.full(len(keep_mask), -1, dtype=np.int32)
    remap[np.flatnonzero(keep_mask)] = np.arange(int(keep_mask.sum()), dtype=np.int32)
    new = remap[arrays["cell"]]
    keep = new >= 0
    return {
        "cell": new[keep],
        "feature": arrays["feature"][keep],
        "covered": arrays["covered"][keep],
        "converted": arrays["converted"][keep],
    }


def feature_statistics(arrays, side_support, n_features: int, args, np):
    observed_cells = np.bincount(
        arrays["feature"], minlength=n_features
    ).astype(np.int64)
    total_covered = np.bincount(
        arrays["feature"], weights=arrays["covered"], minlength=n_features
    )
    total_converted = np.bincount(
        arrays["feature"], weights=arrays["converted"], minlength=n_features
    )
    pooled = np.divide(
        total_converted,
        total_covered,
        out=np.zeros(n_features, dtype=np.float64),
        where=total_covered > 0,
    )
    entry_pooled = np.clip(pooled[arrays["feature"]], 1e-6, 1 - 1e-6)
    residual = (
        arrays["converted"] - arrays["covered"] * entry_pooled
    ) / np.sqrt(arrays["covered"] * entry_pooled * (1 - entry_pooled) + 1.0)
    residual = np.clip(residual, -args.residual_clip, args.residual_clip)
    variability = np.divide(
        np.bincount(
            arrays["feature"], weights=residual * residual, minlength=n_features
        ),
        observed_cells,
        out=np.zeros(n_features, dtype=np.float64),
        where=observed_cells > 0,
    )
    eligible = np.flatnonzero(
        (observed_cells >= args.min_feature_cells)
        & (total_covered >= args.min_feature_covered_calls)
        & (
            side_support["cells"]["human"]
            >= args.min_feature_side_cells
        )
        & (
            side_support["cells"]["pig"]
            >= args.min_feature_side_cells
        )
        & (
            side_support["covered"]["human"]
            >= args.min_feature_side_covered_calls
        )
        & (
            side_support["covered"]["pig"]
            >= args.min_feature_side_covered_calls
        )
    )
    selected = sorted(
        map(int, eligible), key=lambda index: (-float(variability[index]), index)
    )[: args.max_features]
    if len(selected) < 3:
        raise ValueError(f"too few eligible shared WCG orthobins: {len(selected)}")
    return {
        "observed_cells": observed_cells,
        "total_covered": total_covered,
        "total_converted": total_converted,
        "pooled": pooled,
        "entry_residual": residual,
        "variability": variability,
        "eligible": eligible,
        "selected": selected,
        "side_support": side_support,
    }


def estimate_cell_prior(raw_ratios, global_mean: float, fallback: float, np):
    mean = min(max(global_mean, 1e-4), 1 - 1e-4)
    if len(raw_ratios) >= 10:
        variance = float(np.var(raw_ratios, ddof=1))
        concentration = mean * (1 - mean) / variance - 1 if variance > 0 else fallback
    else:
        concentration = fallback
    if not math.isfinite(concentration):
        concentration = fallback
    concentration = min(max(concentration, 2.0), 200.0)
    return mean * concentration, (1 - mean) * concentration, concentration


def build_matrices(
    arrays,
    selected: list[int],
    pooled,
    residual,
    global_means,
    args,
    np,
    sparse,
):
    n_cells = len(global_means)
    n_features = len(selected)
    lookup = np.full(len(pooled), -1, dtype=np.int32)
    lookup[np.asarray(selected, dtype=np.int64)] = np.arange(n_features, dtype=np.int32)
    columns = lookup[arrays["feature"]]
    take = columns >= 0
    cells = arrays["cell"][take]
    columns = columns[take]
    covered = arrays["covered"][take]
    converted = arrays["converted"][take]
    residual_values = residual[take]
    residual_matrix = sparse.coo_matrix(
        (residual_values, (cells, columns)), shape=(n_cells, n_features)
    ).tocsr()
    coverage_matrix = sparse.coo_matrix(
        (np.log1p(covered), (cells, columns)), shape=(n_cells, n_features)
    ).tocsr()
    pooled_selected = pooled[np.asarray(selected, dtype=np.int64)]
    ratio = np.broadcast_to(pooled_selected, (n_cells, n_features)).copy()
    sufficient = covered >= args.min_cell_bin_calls
    ratio[cells[sufficient], columns[sufficient]] = (
        converted[sufficient] / covered[sufficient]
    )
    posterior = np.ones((n_cells, n_features), dtype=np.float32)
    prior_strengths = np.zeros(n_cells, dtype=np.float32)
    for cell in range(n_cells):
        cell_take = cells == cell
        raw_take = cell_take & sufficient
        raw_ratios = converted[raw_take] / covered[raw_take]
        alpha, beta, strength = estimate_cell_prior(
            raw_ratios, float(global_means[cell]), args.prior_strength, np
        )
        prior_strengths[cell] = strength
        if np.any(cell_take):
            values = (converted[cell_take] + alpha) / (
                covered[cell_take] + alpha + beta
            )
            posterior[cell, columns[cell_take]] = values / (alpha / (alpha + beta))
    return {
        "ratio_paper": sparse.csr_matrix(ratio.astype(np.float32)),
        "beta_posterior": sparse.csr_matrix(posterior),
        "pearson_residual": residual_matrix.astype(np.float32),
        "coverage_control": coverage_matrix.astype(np.float32),
        "prior_strengths": prior_strengths,
        "selected_coverage_by_cell": np.bincount(
            cells, weights=covered, minlength=n_cells
        ).astype(np.int64),
    }


def atomic_matrix_market(path: Path, matrix, scipy_io) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(tempfile.mkdtemp(prefix=path.name + ".", dir=path.parent))
    try:
        plain = temporary_dir / "matrix.mtx"
        scipy_io.mmwrite(plain, matrix)
        descriptor, compressed = tempfile.mkstemp(
            prefix=path.name + ".", suffix=".tmp", dir=path.parent
        )
        os.close(descriptor)
        with plain.open("rb") as source, gzip.open(compressed, "wb", compresslevel=1) as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
        os.replace(compressed, path)
    finally:
        shutil.rmtree(temporary_dir, ignore_errors=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--min-total-calls", type=int, default=20_000)
    parser.add_argument("--min-orthobin-wcg-calls", type=int, default=100)
    parser.add_argument("--min-cell-bin-calls", type=int, default=5)
    parser.add_argument("--min-feature-cells", type=int, default=20)
    parser.add_argument("--min-feature-covered-calls", type=int, default=100)
    parser.add_argument("--min-feature-side-cells", type=int, default=2)
    parser.add_argument("--min-feature-side-covered-calls", type=int, default=20)
    parser.add_argument("--max-features", type=int, default=2_000)
    parser.add_argument("--prior-strength", type=float, default=20.0)
    parser.add_argument("--residual-clip", type=float, default=10.0)
    parser.add_argument("--minimum-major-fraction", type=float, default=0.90)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.min_total_calls < 1
        or args.min_orthobin_wcg_calls < 1
        or args.min_cell_bin_calls < 1
        or args.min_feature_cells < 2
        or args.min_feature_covered_calls < 1
        or args.min_feature_side_cells < 1
        or args.min_feature_side_covered_calls < 1
        or args.max_features < 3
        or args.prior_strength <= 0
        or args.residual_clip <= 0
        or not 0.5 < args.minimum_major_fraction <= 1
    ):
        raise ValueError("invalid shared-orthobin matrix parameter")
    try:
        import h5py
        import numpy as np
        import scipy.io
        import scipy.sparse as sparse
    except ImportError as error:
        raise RuntimeError("Step 140 requires h5py, numpy and scipy") from error
    project = args.project
    name = run_id(args)
    output_root = args.output_root or (
        project / "10_downstream_analysis/11_dna_ortholog_100kb_seurat_input" / name
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/41_dna_ortholog_100kb_seurat_input" / name
    )
    summary_path = qc_root / "dna_ortholog_100kb_seurat_input_summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and existing.get("run_id") == name
        ):
            print(f"[SKIP] validated Seurat input exists: {summary_path}")
            return 0
        raise FileExistsError(f"inspect or rerun with --force: {summary_path}")
    index_summary_path = (
        project
        / "08_qc_summary/40_dna_ortholog_100kb/"
        "dna_ortholog_100kb_index_summary.tsv"
    )
    index_summary = read_metrics(index_summary_path)
    if (
        index_summary.get("status") != "validated"
        or index_summary.get("pipeline_signature") != INDEX_SIGNATURE
    ):
        raise ValueError("shared ortholog 100 kb index is not validated")
    feature_rows = read_tsv(Path(index_summary["physical_feature_index"]))
    axis_rows = read_tsv(Path(index_summary["orthobin_axis"]))
    feature_map = load_feature_axis(feature_rows, axis_rows, np)
    feature_species = np.asarray(
        [0 if row["species"] == "human" else 1 if row["species"] == "pig" else -1
         for row in feature_rows],
        dtype=np.int8,
    )
    if np.any(feature_species < 0):
        raise ValueError("unexpected species in physical feature index")
    n_features = len(axis_rows)
    cell_index_path = (
        project
        / "08_qc_summary/20_single_cytosine_inputs/"
        "dna_single_cytosine_cell_index.tsv.gz"
    )
    canonical_cells = read_tsv(cell_index_path)
    if [int(row["cell_index"]) for row in canonical_cells] != list(
        range(len(canonical_cells))
    ):
        raise ValueError("canonical cell index is not contiguous")
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
        raise ValueError(f"too few cells pass whole-genome depth: {len(initial_indices)}")
    initial_lookup = np.full(len(canonical_cells), -1, dtype=np.int32)
    initial_lookup[np.asarray(initial_indices, dtype=np.int64)] = np.arange(
        len(initial_indices), dtype=np.int32
    )
    arrays, side_support = load_well_arrays(
        project,
        canonical_cells,
        initial_lookup,
        feature_map,
        feature_species,
        n_features,
        h5py,
        np,
    )
    orthobin_coverage = np.bincount(
        arrays["cell"], weights=arrays["covered"], minlength=len(initial_indices)
    ).astype(np.int64)
    keep = orthobin_coverage >= args.min_orthobin_wcg_calls
    retained_indices = [
        index for index, retain in zip(initial_indices, keep) if retain
    ]
    if len(retained_indices) < 50:
        raise ValueError(
            f"too few cells pass shared-orthobin WCG depth: {len(retained_indices)}"
        )
    arrays = remap_cells(arrays, keep, np)
    statistics = feature_statistics(arrays, side_support, n_features, args, np)
    selected = statistics["selected"]
    global_means = np.asarray(
        [float(qc_by_index[index]["WCG_methylation_fraction"]) for index in retained_indices],
        dtype=np.float64,
    )
    matrices = build_matrices(
        arrays,
        selected,
        statistics["pooled"],
        statistics["entry_residual"],
        global_means,
        args,
        np,
        sparse,
    )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    matrix_paths: dict[str, Path] = {}
    for method in METHODS:
        path = output_root / f"{method}.features_by_cells.mtx.gz"
        atomic_matrix_market(path, matrices[method].transpose().tocoo(), scipy.io)
        matrix_paths[method] = path
    selected_rows = []
    for rank, axis_index in enumerate(selected, 1):
        axis = axis_rows[axis_index]
        selected_rows.append(
            {
                "matrix_feature_index": rank - 1,
                "variability_rank": rank,
                **axis,
                "observed_cells": int(statistics["observed_cells"][axis_index]),
                "total_covered_calls": int(statistics["total_covered"][axis_index]),
                "total_converted_calls": int(statistics["total_converted"][axis_index]),
                "human_observed_cells": int(
                    statistics["side_support"]["cells"]["human"][axis_index]
                ),
                "pig_observed_cells": int(
                    statistics["side_support"]["cells"]["pig"][axis_index]
                ),
                "human_covered_calls": int(
                    statistics["side_support"]["covered"]["human"][axis_index]
                ),
                "pig_covered_calls": int(
                    statistics["side_support"]["covered"]["pig"][axis_index]
                ),
                "pooled_WCG_methylation_fraction": float(statistics["pooled"][axis_index]),
                "pearson_residual_variability": float(statistics["variability"][axis_index]),
            }
        )
    feature_path = output_root / "features.tsv.gz"
    atomic_tsv(feature_path, list(selected_rows[0]), selected_rows, gzip_output=True)
    cell_rows = []
    label_counts = Counter()
    prior_strengths = matrices["prior_strengths"]
    selected_coverage = matrices["selected_coverage_by_cell"]
    all_shared_coverage = np.bincount(
        arrays["cell"],
        weights=arrays["covered"],
        minlength=len(retained_indices),
    ).astype(np.int64)
    for matrix_index, global_index in enumerate(retained_indices):
        canonical = canonical_cells[global_index]
        qc_row = qc_by_index[global_index]
        if canonical["CB"] != qc_row["CB"]:
            raise ValueError(f"cell metadata differs at {global_index}")
        human = int(qc_row["human_signal_covered_calls"])
        pig = int(qc_row["pig_signal_covered_calls"])
        label, confidence = classify_cell_line(
            human, pig, args.minimum_major_fraction
        )
        label_counts[label] += 1
        cell_rows.append(
            {
                "matrix_cell_index": matrix_index,
                "cell_index": global_index,
                "cell_id": canonical["cell_id"],
                "physical_well": canonical["sample"],
                "well_index": canonical["sample_index"],
                "well_set": canonical["group"],
                "well_design": well_design(canonical["group"]),
                "CB": canonical["CB"],
                "DNA_mapping_derived_cell_line_posthoc": label,
                "mapping_major_reference_fraction": confidence,
                "human_signal_covered_calls": human,
                "pig_signal_covered_calls": pig,
                "whole_genome_WCG_covered_calls": qc_row["WCG_covered_calls"],
                "whole_genome_WCG_methylation_fraction": qc_row[
                    "WCG_methylation_fraction"
                ],
                "whole_genome_signal_covered_calls": qc_row[
                    "total_signal_covered_calls"
                ],
                "all_shared_orthobin_WCG_covered_calls": int(
                    all_shared_coverage[matrix_index]
                ),
                "selected_shared_orthobin_WCG_covered_calls": int(
                    selected_coverage[matrix_index]
                ),
                "beta_prior_strength": float(prior_strengths[matrix_index]),
            }
        )
    cells_path = output_root / "cells.tsv.gz"
    atomic_tsv(cells_path, list(cell_rows[0]), cell_rows, gzip_output=True)
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("run_id", name),
        ("mode", "label_blind_WCG_only_shared_human_pig_orthobins"),
        ("initial_cells_after_whole_genome_depth", len(initial_indices)),
        ("retained_cells", len(retained_indices)),
        ("orthobin_axes_in_index", n_features),
        ("eligible_orthobins", len(statistics["eligible"])),
        ("selected_variable_orthobins", len(selected)),
        ("minimum_total_WCG_plus_GCH_calls", args.min_total_calls),
        ("minimum_shared_orthobin_WCG_calls", args.min_orthobin_wcg_calls),
        ("minimum_cell_bin_WCG_calls_for_raw_ratio", args.min_cell_bin_calls),
        ("minimum_feature_cells", args.min_feature_cells),
        ("minimum_feature_covered_calls", args.min_feature_covered_calls),
        ("minimum_feature_cells_each_reference_side", args.min_feature_side_cells),
        (
            "minimum_feature_covered_calls_each_reference_side",
            args.min_feature_side_covered_calls,
        ),
        ("fallback_beta_prior_strength", args.prior_strength),
        ("GM12878_posthoc_cells", label_counts["GM12878"]),
        ("PK15_posthoc_cells", label_counts["PK15"]),
        ("ambiguous_posthoc_cells", label_counts["ambiguous"]),
        ("posthoc_label_policy", f"major_mapping_reference_fraction>={args.minimum_major_fraction}"),
        ("posthoc_labels_used_for_matrix_or_feature_selection", "no"),
        ("missing_ratio_paper", "pooled_orthobin_mean"),
        ("missing_beta_posterior", "cell_prior_normalized_to_1"),
        ("missing_pearson_residual", "feature_expectation_equal_to_0"),
        ("cells", cells_path),
        ("features", feature_path),
        ("source_index_summary", index_summary_path),
    ]
    for method, path in matrix_paths.items():
        summary_rows.append((f"matrix_{method}", path))
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] Seurat matrices {name}: cells={len(retained_indices)} "
        f"features={len(selected)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
