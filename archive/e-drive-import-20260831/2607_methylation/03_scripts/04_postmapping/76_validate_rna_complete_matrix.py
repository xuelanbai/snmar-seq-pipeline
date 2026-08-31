#!/usr/bin/env python3
"""Step 76: validate and package complete raw/filtered RNA matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import tempfile
from pathlib import Path


EXPECTED_SAMPLES = 96


def read_tsv(path: Path) -> list[dict[str, str]]:
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


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def inspect_matrix(path: Path) -> dict[str, int]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        first = handle.readline().strip()
        if first != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"unsupported MatrixMarket header: {path}")
        dimensions = None
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"invalid MatrixMarket dimensions: {path}")
            dimensions = tuple(map(int, fields))
            break
        if dimensions is None:
            raise ValueError(f"missing MatrixMarket dimensions: {path}")
        features, cells, expected_entries = dimensions
        entries = 0
        umi_sum = 0
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"invalid MatrixMarket entry: {path}")
            feature, cell, count = map(int, fields)
            if not (
                1 <= feature <= features
                and 1 <= cell <= cells
                and count > 0
            ):
                raise ValueError(f"out-of-range MatrixMarket entry: {path}")
            entries += 1
            umi_sum += count
    if entries != expected_entries:
        raise ValueError(f"MatrixMarket nnz does not reconcile: {path}")
    return {
        "features": features,
        "cells": cells,
        "nonzero_entries": entries,
        "exact_umis": umi_sum,
    }


def read_matrix_coo(path: Path):
    """Read validated MatrixMarket coordinates as zero-based NumPy arrays."""
    try:
        import numpy as np
    except ImportError as error:
        raise RuntimeError("Step 76 requires numpy") from error
    feature_indices = []
    cell_indices = []
    counts = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        first = handle.readline().strip()
        if first != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"unsupported MatrixMarket header: {path}")
        dimensions = None
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            dimensions = tuple(map(int, line.split()))
            break
        if dimensions is None or len(dimensions) != 3:
            raise ValueError(f"missing MatrixMarket dimensions: {path}")
        features, cells, expected_entries = dimensions
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            feature, cell, count = map(int, line.split())
            feature_indices.append(feature - 1)
            cell_indices.append(cell - 1)
            counts.append(count)
    if len(counts) != expected_entries:
        raise ValueError(f"MatrixMarket nnz differs: {path}")
    return {
        "features": features,
        "cells": cells,
        "feature_index": np.asarray(feature_indices, dtype=np.uint32),
        "cell_index": np.asarray(cell_indices, dtype=np.uint32),
        "exact_umis": np.asarray(counts, dtype=np.uint32),
    }


def write_rna_h5(
    path: Path,
    matrices: dict[str, dict[str, object]],
    matrix_roots: dict[str, Path],
) -> None:
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("Step 76 requires h5py") from error
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        with h5py.File(temporary, "w") as output:
            output.attrs["status"] = "validated"
            output.attrs["modality"] = "RNA_expression"
            output.attrs["matrix_format"] = "sparse_COO_gene_by_cell"
            output.attrs["cell_axes"] = "raw_and_filtered_are_independent"
            for role in ("raw", "filtered"):
                matrix = matrices[role]
                group = output.create_group(role)
                for field in ("feature_index", "cell_index", "exact_umis"):
                    values = matrix[field]
                    options = {}
                    if len(values):
                        options = {
                            "compression": "gzip",
                            "compression_opts": 1,
                            "shuffle": True,
                            "chunks": True,
                        }
                    group.create_dataset(field, data=values, **options)
                group.attrs["features"] = int(matrix["features"])
                group.attrs["cells"] = int(matrix["cells"])
                group.attrs["nonzero_entries"] = len(matrix["exact_umis"])
                group.attrs["matrix_market_path"] = str(
                    matrix_roots[role] / "matrix.mtx.gz"
                )
                group.attrs["barcodes_path"] = str(
                    matrix_roots[role] / "barcodes.tsv.gz"
                )
                group.attrs["features_path"] = str(
                    matrix_roots[role] / "features.tsv.gz"
                )
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
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--minimum-exact-umis", type=int, default=10)
    parser.add_argument("--minimum-genes", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.minimum_exact_umis < 1 or args.minimum_genes < 1:
        raise ValueError("RNA thresholds must be positive")
    project = args.project
    suffix = f"minUMI{args.minimum_exact_umis}_minGene{args.minimum_genes}"
    source_root = project / "07_postmapping/rna_reverse_genefull_full"
    raw_root = source_root / "raw_feature_bc_matrix"
    filtered_root = (
        source_root / f"filtered_feature_bc_matrix_{suffix}"
    )
    raw_qc = project / "08_qc_summary/16_rna_reverse_genefull_full"
    filtered_qc = (
        project / "08_qc_summary/19_rna_filtered_matrix_full" / suffix
    )
    output_root = (
        args.output_root
        or project / "07_postmapping/complete_matrices/rna_expression"
    )
    qc_root = (
        args.qc_root
        or project / "08_qc_summary/20_complete_matrices/rna_expression"
    )
    raw_summary_path = raw_qc / "rna_reverse_genefull_full_summary.tsv"
    raw_cell_qc_path = raw_qc / "rna_reverse_cell_qc.tsv.gz"
    filtered_summary_path = (
        filtered_qc / "rna_filtered_matrix_full_summary.tsv"
    )
    filtered_per_sample_path = (
        filtered_qc / "rna_filtered_matrix_full_per_sample.tsv"
    )
    filtered_cell_qc_path = (
        filtered_qc / "rna_filtered_cell_qc.tsv.gz"
    )
    matrix_roots = {"raw": raw_root, "filtered": filtered_root}
    prerequisites = [
        raw_summary_path,
        raw_cell_qc_path,
        filtered_summary_path,
        filtered_per_sample_path,
        filtered_cell_qc_path,
    ]
    for root in matrix_roots.values():
        prerequisites.extend(
            [
                root / "matrix.mtx.gz",
                root / "barcodes.tsv.gz",
                root / "features.tsv.gz",
                root / "feature_metadata.tsv.gz",
            ]
        )
    for path in prerequisites:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty prerequisite: {path}")

    raw_summary = read_metrics(raw_summary_path)
    filtered_summary = read_metrics(filtered_summary_path)
    if (
        raw_summary.get("status") != "validated"
        or raw_summary.get("cell_filtering")
        != "not_applied_raw_matrix"
        or filtered_summary.get("status") != "validated"
        or filtered_summary.get("raw_matrix_retained") != "1"
    ):
        raise ValueError("RNA matrix summaries are not validated")

    matrix_stats = {
        role: inspect_matrix(root / "matrix.mtx.gz")
        for role, root in matrix_roots.items()
    }
    expected_raw = {
        "features": int(raw_summary["raw_matrix_features"]),
        "cells": int(raw_summary["raw_matrix_cells"]),
        "nonzero_entries": int(raw_summary["raw_matrix_nonzero_entries"]),
        "exact_umis": int(raw_summary["raw_matrix_umi_sum"]),
    }
    expected_filtered = {
        "features": int(filtered_summary["features"]),
        "cells": int(filtered_summary["filtered_cells"]),
        "nonzero_entries": int(
            filtered_summary["filtered_nonzero_entries"]
        ),
        "exact_umis": int(filtered_summary["filtered_exact_umis"]),
    }
    if (
        matrix_stats["raw"] != expected_raw
        or matrix_stats["filtered"] != expected_filtered
    ):
        raise ValueError("RNA MatrixMarket statistics do not reconcile")

    raw_barcodes = read_gzip_lines(raw_root / "barcodes.tsv.gz")
    filtered_barcodes = read_gzip_lines(
        filtered_root / "barcodes.tsv.gz"
    )
    raw_features = read_gzip_lines(raw_root / "features.tsv.gz")
    filtered_features = read_gzip_lines(
        filtered_root / "features.tsv.gz"
    )
    if (
        len(raw_barcodes) != expected_raw["cells"]
        or len(set(raw_barcodes)) != len(raw_barcodes)
        or len(filtered_barcodes) != expected_filtered["cells"]
        or len(set(filtered_barcodes)) != len(filtered_barcodes)
        or not set(filtered_barcodes).issubset(raw_barcodes)
    ):
        raise ValueError("RNA barcode lists do not reconcile")
    if (
        len(raw_features) != expected_raw["features"]
        or raw_features != filtered_features
    ):
        raise ValueError("RNA feature lists differ")

    filtered_cell_rows = []
    with gzip.open(
        filtered_cell_qc_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        filtered_cell_rows = list(csv.DictReader(handle, delimiter="\t"))
    if (
        len(filtered_cell_rows) != len(filtered_barcodes)
        or {row["cell_id"] for row in filtered_cell_rows}
        != set(filtered_barcodes)
        or sum(int(row["exact_umis"]) for row in filtered_cell_rows)
        != expected_filtered["exact_umis"]
    ):
        raise ValueError("RNA filtered cell QC does not match matrix")

    raw_cell_rows = []
    with gzip.open(
        raw_cell_qc_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        raw_cell_rows = list(csv.DictReader(handle, delimiter="\t"))
    if (
        len(raw_cell_rows) != len(raw_barcodes)
        or {row["cell_id"] for row in raw_cell_rows}
        != set(raw_barcodes)
        or sum(int(row["exact_umis"]) for row in raw_cell_rows)
        != expected_raw["exact_umis"]
    ):
        raise ValueError("RNA raw cell QC does not match matrix")

    per_sample = read_tsv(filtered_per_sample_path)
    if (
        len(per_sample) != EXPECTED_SAMPLES
        or sum(int(row["filtered_cells"]) for row in per_sample)
        != expected_filtered["cells"]
        or sum(int(row["filtered_exact_umis"]) for row in per_sample)
        != expected_filtered["exact_umis"]
    ):
        raise ValueError("RNA filtered per-sample table does not reconcile")

    output_root.mkdir(parents=True, exist_ok=True)
    h5_path = output_root / "rna_expression_complete_sparse.h5"
    coo_matrices = {
        role: read_matrix_coo(root / "matrix.mtx.gz")
        for role, root in matrix_roots.items()
    }
    for role in ("raw", "filtered"):
        coo = coo_matrices[role]
        stats = matrix_stats[role]
        if (
            int(coo["features"]) != stats["features"]
            or int(coo["cells"]) != stats["cells"]
            or len(coo["exact_umis"]) != stats["nonzero_entries"]
            or int(coo["exact_umis"].sum()) != stats["exact_umis"]
        ):
            raise ValueError(f"RNA {role} COO conversion differs")
    write_rna_h5(h5_path, coo_matrices, matrix_roots)

    manifest_rows = []
    for role in ("raw", "filtered"):
        root = matrix_roots[role]
        stats = matrix_stats[role]
        matrix_path = root / "matrix.mtx.gz"
        manifest_rows.append(
            {
                "matrix_role": role,
                "matrix_path": matrix_path,
                "barcodes_path": root / "barcodes.tsv.gz",
                "features_path": root / "features.tsv.gz",
                "feature_metadata_path": root
                / "feature_metadata.tsv.gz",
                "cell_metadata_path": (
                    raw_cell_qc_path
                    if role == "raw"
                    else filtered_cell_qc_path
                ),
                "sparse_h5_path": h5_path,
                "sparse_h5_group": f"/{role}",
                "features": stats["features"],
                "cells": stats["cells"],
                "nonzero_entries": stats["nonzero_entries"],
                "exact_umis": stats["exact_umis"],
                "matrix_sha256": sha256(matrix_path),
                "status": "validated",
            }
        )
    atomic_tsv(
        output_root / "rna_expression_matrix_manifest.tsv",
        list(manifest_rows[0]),
        manifest_rows,
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "modality", "value": "RNA_expression"},
        {"metric": "matrix_format", "value": "MatrixMarket_gene_by_cell"},
        {"metric": "sparse_h5", "value": h5_path},
        {
            "metric": "sparse_h5_format",
            "value": "zero_based_COO_gene_by_cell",
        },
        {"metric": "raw_matrix_path", "value": raw_root},
        {"metric": "raw_cell_metadata", "value": raw_cell_qc_path},
        {"metric": "raw_cells", "value": expected_raw["cells"]},
        {"metric": "raw_features", "value": expected_raw["features"]},
        {"metric": "raw_exact_umis", "value": expected_raw["exact_umis"]},
        {"metric": "filtered_matrix_path", "value": filtered_root},
        {
            "metric": "filtered_cell_metadata",
            "value": filtered_cell_qc_path,
        },
        {"metric": "filtered_cells", "value": expected_filtered["cells"]},
        {
            "metric": "filtered_features",
            "value": expected_filtered["features"],
        },
        {
            "metric": "filtered_exact_umis",
            "value": expected_filtered["exact_umis"],
        },
        {
            "metric": "filter_rule",
            "value": (
                f"exact_umis>={args.minimum_exact_umis}_and_"
                f"genes>={args.minimum_genes}"
            ),
        },
        {
            "metric": "cross_modality_cell_pairing",
            "value": "not_assumed",
        },
        {"metric": "rna_umap_status", "value": "not_run_separate_step"},
    ]
    atomic_tsv(
        qc_root / "rna_complete_expression_matrix_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] RNA complete matrices: raw={expected_raw['cells']} cells, "
        f"filtered={expected_filtered['cells']} cells"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
