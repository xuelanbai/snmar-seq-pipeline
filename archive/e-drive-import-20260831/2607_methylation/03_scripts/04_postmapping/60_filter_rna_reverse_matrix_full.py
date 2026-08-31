#!/usr/bin/env python3
"""Step 60: create the full filtered RNA matrix from the raw matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable


def read_metrics(path: Path) -> dict[str, str]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    metrics = {row["metric"]: row["value"] for row in rows}
    if not rows or len(metrics) != len(rows):
        raise ValueError(f"invalid metric TSV: {path}")
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


def atomic_gzip_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as out:
            for line in lines:
                out.write(line)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def read_gzip_lines(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def parse_matrix_dimensions(handle) -> tuple[int, int, int]:
    first = handle.readline()
    if not first.startswith("%%MatrixMarket matrix coordinate integer"):
        raise ValueError("unsupported MatrixMarket header")
    for line in handle:
        if line.startswith("%"):
            continue
        fields = line.split()
        if len(fields) != 3:
            raise ValueError("invalid MatrixMarket dimensions")
        return int(fields[0]), int(fields[1]), int(fields[2])
    raise ValueError("missing MatrixMarket dimensions")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--raw-matrix-root", type=Path)
    parser.add_argument("--cell-qc", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--min-exact-umis", type=int, default=10)
    parser.add_argument("--min-genes", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_exact_umis < 1 or args.min_genes < 1:
        raise ValueError("RNA filtering thresholds must be positive")
    project = args.project
    source_qc = project / "08_qc_summary/16_rna_reverse_genefull_full"
    raw_root = (
        args.raw_matrix_root
        or project
        / "07_postmapping/rna_reverse_genefull_full/raw_feature_bc_matrix"
    )
    cell_qc_path = (
        args.cell_qc or source_qc / "rna_reverse_cell_qc.tsv.gz"
    )
    summary_path = (
        args.summary or source_qc / "rna_reverse_genefull_full_summary.tsv"
    )
    suffix = (
        f"minUMI{args.min_exact_umis}_minGene{args.min_genes}"
    )
    output_root = (
        args.output_root
        or project
        / "07_postmapping/rna_reverse_genefull_full"
        / f"filtered_feature_bc_matrix_{suffix}"
    )
    qc_root = (
        args.qc_root
        or project / "08_qc_summary/19_rna_filtered_matrix_full" / suffix
    )
    matrix_path = raw_root / "matrix.mtx.gz"
    barcodes_path = raw_root / "barcodes.tsv.gz"
    features_path = raw_root / "features.tsv.gz"
    metadata_path = raw_root / "feature_metadata.tsv.gz"
    for path in (
        summary_path,
        cell_qc_path,
        matrix_path,
        barcodes_path,
        features_path,
        metadata_path,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty prerequisite: {path}")

    summary = read_metrics(summary_path)
    if (
        summary.get("status") != "validated"
        or summary.get("cell_filtering") != "not_applied_raw_matrix"
    ):
        raise ValueError("RNA source is not a validated raw matrix")

    all_qc_rows: dict[str, dict[str, str]] = {}
    cell_qc_fields: list[str] | None = None
    with gzip.open(
        cell_qc_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        cell_qc_fields = reader.fieldnames
        if not cell_qc_fields:
            raise ValueError("RNA cell-QC table has no header")
        for row in reader:
            cell_id = row["cell_id"]
            if cell_id in all_qc_rows:
                raise ValueError(f"duplicate RNA cell_id: {cell_id}")
            all_qc_rows[cell_id] = row

    raw_barcodes = read_gzip_lines(barcodes_path)
    if (
        len(raw_barcodes) != int(summary["raw_matrix_cells"])
        or len(raw_barcodes) != len(set(raw_barcodes))
    ):
        raise ValueError("raw barcode list does not match matrix summary")
    missing_qc = [cell for cell in raw_barcodes if cell not in all_qc_rows]
    if missing_qc:
        raise ValueError(f"raw barcode missing cell QC: {missing_qc[0]}")

    selected_barcodes = []
    selected_qc_rows = []
    for cell_id in raw_barcodes:
        row = all_qc_rows[cell_id]
        if int(row["in_raw_matrix"]) != 1:
            raise ValueError(f"raw barcode not marked in matrix: {cell_id}")
        if (
            int(row["exact_umis"]) >= args.min_exact_umis
            and int(row["genes_detected"]) >= args.min_genes
        ):
            selected_barcodes.append(cell_id)
            selected_qc_rows.append(row)
    if not selected_barcodes:
        raise ValueError("RNA thresholds retained zero cells")
    selected_set = set(selected_barcodes)
    old_to_new = {
        old_index: new_index
        for new_index, (old_index, cell_id) in enumerate(
            (
                (index, cell)
                for index, cell in enumerate(raw_barcodes, 1)
                if cell in selected_set
            ),
            1,
        )
    }
    if len(old_to_new) != len(selected_barcodes):
        raise ValueError("RNA selected-column map is inconsistent")

    output_root.mkdir(parents=True, exist_ok=True)
    matrix_output = output_root / "matrix.mtx.gz"
    descriptor, temporary = tempfile.mkstemp(
        prefix=matrix_output.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    raw_umi_sum = 0
    raw_entries = 0
    filtered_umi_sum = 0
    filtered_entries = 0
    try:
        with gzip.open(
            matrix_path, "rt", encoding="utf-8"
        ) as source, gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as out:
            features, cells, expected_entries = parse_matrix_dimensions(source)
            if (
                features != int(summary["raw_matrix_features"])
                or cells != len(raw_barcodes)
                or expected_entries != int(summary["raw_matrix_nonzero_entries"])
            ):
                raise ValueError("raw MatrixMarket dimensions do not reconcile")
            out.write("%%MatrixMarket matrix coordinate integer general\n")
            out.write(
                f"% snMAR reverse GeneFull filtered: "
                f"exact_umis>={args.min_exact_umis}, "
                f"genes>={args.min_genes}\n"
            )
            # Write entries to a temporary body because dimensions require nnz.
            body_descriptor, body_path_text = tempfile.mkstemp(
                prefix="matrix.body.", suffix=".tmp", dir=output_root
            )
            try:
                with os.fdopen(
                    body_descriptor, "w", encoding="utf-8", newline=""
                ) as body:
                    for line in source:
                        if not line.strip() or line.startswith("%"):
                            continue
                        fields = line.split()
                        if len(fields) != 3:
                            raise ValueError("invalid MatrixMarket entry")
                        feature_index, old_cell, count = map(int, fields)
                        if (
                            feature_index < 1
                            or feature_index > features
                            or old_cell < 1
                            or old_cell > cells
                            or count <= 0
                        ):
                            raise ValueError("out-of-range MatrixMarket entry")
                        raw_entries += 1
                        raw_umi_sum += count
                        new_cell = old_to_new.get(old_cell)
                        if new_cell is not None:
                            filtered_entries += 1
                            filtered_umi_sum += count
                            body.write(
                                f"{feature_index} {new_cell} {count}\n"
                            )
                if raw_entries != expected_entries:
                    raise ValueError("raw MatrixMarket nnz does not reconcile")
                out.write(
                    f"{features} {len(selected_barcodes)} "
                    f"{filtered_entries}\n"
                )
                with open(
                    body_path_text, "rt", encoding="utf-8"
                ) as body:
                    shutil_copyfileobj(body, out)
            finally:
                try:
                    os.unlink(body_path_text)
                except FileNotFoundError:
                    pass
        if raw_umi_sum != int(summary["raw_matrix_umi_sum"]):
            raise ValueError("raw MatrixMarket UMI sum does not reconcile")
        expected_filtered_umis = sum(
            int(row["exact_umis"]) for row in selected_qc_rows
        )
        if filtered_umi_sum != expected_filtered_umis:
            raise ValueError("filtered matrix and selected-cell UMI sums differ")
        os.replace(temporary, matrix_output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    atomic_gzip_lines(
        output_root / "barcodes.tsv.gz",
        (f"{cell_id}\n" for cell_id in selected_barcodes),
    )
    for source, destination in (
        (features_path, output_root / "features.tsv.gz"),
        (metadata_path, output_root / "feature_metadata.tsv.gz"),
    ):
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            atomic_gzip_lines(destination, iter(handle))

    assert cell_qc_fields is not None
    selected_qc_output = qc_root / "rna_filtered_cell_qc.tsv.gz"
    selected_qc_output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=selected_qc_output.name + ".",
        suffix=".tmp",
        dir=selected_qc_output.parent,
    )
    os.close(descriptor)
    try:
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as out:
            writer = csv.DictWriter(
                out,
                fieldnames=cell_qc_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(selected_qc_rows)
        os.replace(temporary, selected_qc_output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    raw_by_sample: Counter[str] = Counter()
    filtered_by_sample: Counter[str] = Counter()
    raw_umis_by_sample: Counter[str] = Counter()
    filtered_umis_by_sample: Counter[str] = Counter()
    for cell_id in raw_barcodes:
        row = all_qc_rows[cell_id]
        sample = row["sample"]
        raw_by_sample[sample] += 1
        raw_umis_by_sample[sample] += int(row["exact_umis"])
    for row in selected_qc_rows:
        sample = row["sample"]
        filtered_by_sample[sample] += 1
        filtered_umis_by_sample[sample] += int(row["exact_umis"])
    per_sample_rows = []
    for sample in sorted(raw_by_sample):
        raw_cells = raw_by_sample[sample]
        kept_cells = filtered_by_sample[sample]
        raw_umis = raw_umis_by_sample[sample]
        kept_umis = filtered_umis_by_sample[sample]
        per_sample_rows.append(
            {
                "sample": sample,
                "raw_cells": raw_cells,
                "filtered_cells": kept_cells,
                "cell_fraction": kept_cells / raw_cells,
                "raw_exact_umis": raw_umis,
                "filtered_exact_umis": kept_umis,
                "umi_fraction": kept_umis / raw_umis if raw_umis else 0.0,
                "minimum_exact_umis": args.min_exact_umis,
                "minimum_genes": args.min_genes,
                "status": "filtered",
            }
        )
    qc_root.mkdir(parents=True, exist_ok=True)
    atomic_tsv(
        qc_root / "rna_filtered_matrix_full_per_sample.tsv",
        list(per_sample_rows[0]),
        per_sample_rows,
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "source_matrix", "value": raw_root},
        {"metric": "output_matrix", "value": output_root},
        {
            "metric": "filter_rule",
            "value": "exact_umis>=minimum AND genes_detected>=minimum",
        },
        {"metric": "minimum_exact_umis", "value": args.min_exact_umis},
        {"metric": "minimum_genes", "value": args.min_genes},
        {"metric": "raw_cells", "value": len(raw_barcodes)},
        {"metric": "filtered_cells", "value": len(selected_barcodes)},
        {
            "metric": "cell_fraction",
            "value": len(selected_barcodes) / len(raw_barcodes),
        },
        {"metric": "raw_exact_umis", "value": raw_umi_sum},
        {"metric": "filtered_exact_umis", "value": filtered_umi_sum},
        {
            "metric": "umi_fraction",
            "value": filtered_umi_sum / raw_umi_sum,
        },
        {"metric": "features", "value": int(summary["raw_matrix_features"])},
        {"metric": "filtered_nonzero_entries", "value": filtered_entries},
        {"metric": "samples", "value": len(per_sample_rows)},
        {"metric": "raw_matrix_retained", "value": 1},
    ]
    atomic_tsv(
        qc_root / "rna_filtered_matrix_full_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] RNA filtered matrix: cells={len(selected_barcodes)}/"
        f"{len(raw_barcodes)}, UMIs={filtered_umi_sum}/{raw_umi_sum}"
    )
    return 0


def shutil_copyfileobj(source, destination, length: int = 1024 * 1024) -> None:
    while True:
        chunk = source.read(length)
        if not chunk:
            break
        destination.write(chunk)


if __name__ == "__main__":
    raise SystemExit(main())
