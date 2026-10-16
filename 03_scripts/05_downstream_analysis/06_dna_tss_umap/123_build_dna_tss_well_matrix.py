#!/usr/bin/env python3
"""Step 123: aggregate one physical well into TSS promoter matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable


EXPECTED_WELLS = 96
EXPECTED_SHARDS = 8
SOURCE_SIGNATURE = "snmar_single_cytosine_validation_20260801"
INDEX_SIGNATURE = "snmar_dna_tss_promoter_index_20260808"
PIPELINE_SIGNATURE = "snmar_dna_tss_well_matrix_20260808"
CONTEXTS = {1: "WCG", 2: "GCH"}
MATRIX_FIELDS = (
    "cell_index",
    "feature_index",
    "context_code",
    "observed_sites",
    "covered_calls",
    "converted_calls",
    "signal_fraction",
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
    metrics = {row["metric"]: row["value"] for row in rows}
    if len(metrics) != len(rows):
        raise ValueError(f"duplicate metric: {path}")
    return metrics


def atomic_tsv(
    path: Path,
    fields: list[str],
    rows: Iterable[dict[str, object]],
) -> None:
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


def build_interval_lookup(
    feature_rows: list[dict[str, str]], np
) -> dict[int, tuple[object, object, object]]:
    indices = [int(row["feature_index"]) for row in feature_rows]
    if indices != list(range(len(feature_rows))):
        raise ValueError("feature index must be zero-based and contiguous")
    grouped: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    for row in feature_rows:
        grouped[int(row["contig_index"])].append(
            (int(row["start0"]), int(row["end0"]), int(row["feature_index"]))
        )
    result: dict[int, tuple[object, object, object]] = {}
    for contig_index, intervals in grouped.items():
        intervals.sort()
        previous_end = -1
        for start0, end0, _feature_index in intervals:
            if start0 < 0 or end0 <= start0 or start0 < previous_end:
                raise ValueError("promoter intervals overlap or are invalid")
            previous_end = end0
        result[contig_index] = (
            np.asarray([value[0] for value in intervals], dtype=np.int64),
            np.asarray([value[1] for value in intervals], dtype=np.int64),
            np.asarray([value[2] for value in intervals], dtype=np.int64),
        )
    return result


def assign_features(contig_index, position0, lookup, np):
    contigs = np.asarray(contig_index, dtype=np.int64)
    positions = np.asarray(position0, dtype=np.int64)
    if contigs.shape != positions.shape:
        raise ValueError("contig and position arrays differ")
    feature = np.full(contigs.size, -1, dtype=np.int64)
    for contig in np.unique(contigs):
        intervals = lookup.get(int(contig))
        if intervals is None:
            continue
        starts, ends, indices = intervals
        selected_rows = np.flatnonzero(contigs == contig)
        selected_positions = positions[selected_rows]
        candidate = np.searchsorted(starts, selected_positions, side="right") - 1
        valid = candidate >= 0
        if np.any(valid):
            candidate_rows = selected_rows[valid]
            candidate_indices = candidate[valid]
            inside = selected_positions[valid] < ends[candidate_indices]
            feature[candidate_rows[inside]] = indices[candidate_indices[inside]]
    return feature


def aggregate_source_arrays(
    cell_index,
    contig_index,
    position0,
    context_code,
    covered_calls,
    converted_calls,
    lookup,
    total_features: int,
    np,
) -> tuple[dict[str, object], dict[str, int]]:
    arrays = (
        cell_index,
        contig_index,
        position0,
        context_code,
        covered_calls,
        converted_calls,
    )
    lengths = {len(values) for values in arrays}
    if len(lengths) != 1:
        raise ValueError("source arrays have different lengths")
    entries = lengths.pop()
    if entries == 0:
        return empty_matrix(np), {
            "source_entries": 0,
            "selected_entries": 0,
            "WCG_covered": 0,
            "WCG_converted": 0,
            "GCH_covered": 0,
            "GCH_converted": 0,
        }
    cell = np.asarray(cell_index, dtype=np.uint64)
    context = np.asarray(context_code, dtype=np.int64)
    covered = np.asarray(covered_calls, dtype=np.uint64)
    converted = np.asarray(converted_calls, dtype=np.uint64)
    if (
        np.any((context != 1) & (context != 2))
        or np.any(covered < 1)
        or np.any(converted > covered)
    ):
        raise ValueError("invalid source count/context")
    feature = assign_features(contig_index, position0, lookup, np)
    keep = feature >= 0
    selected = int(keep.sum())
    metrics = {"source_entries": entries, "selected_entries": selected}
    for code, name in CONTEXTS.items():
        take = keep & (context == code)
        metrics[f"{name}_covered"] = int(covered[take].sum(dtype=np.uint64))
        metrics[f"{name}_converted"] = int(converted[take].sum(dtype=np.uint64))
    if not selected:
        return empty_matrix(np), metrics
    cell = cell[keep]
    feature = feature[keep].astype(np.uint64)
    context = context[keep]
    covered = covered[keep]
    converted = converted[keep]
    keys = (cell * np.uint64(total_features) + feature) * 2 + (
        context.astype(np.uint64) - 1
    )
    return reduce_entries(keys, covered, converted, total_features, np), metrics


def reduce_entries(keys, covered, converted, total_features: int, np):
    if not len(keys):
        return empty_matrix(np)
    order = np.argsort(keys, kind="stable")
    keys = keys[order]
    covered = np.asarray(covered, dtype=np.uint64)[order]
    converted = np.asarray(converted, dtype=np.uint64)[order]
    starts = np.concatenate(
        (
            np.asarray([0], dtype=np.int64),
            np.flatnonzero(keys[1:] != keys[:-1]) + 1,
        )
    )
    unique = keys[starts]
    observed = np.diff(np.append(starts, len(keys))).astype(np.uint64)
    covered_sum = np.add.reduceat(covered, starts)
    converted_sum = np.add.reduceat(converted, starts)
    if (
        np.any(converted_sum > covered_sum)
        or max(
            int(observed.max()), int(covered_sum.max()), int(converted_sum.max())
        )
        > np.iinfo(np.uint32).max
    ):
        raise ValueError("aggregated count exceeds uint32")
    context = (unique % 2 + 1).astype(np.uint8)
    cell_feature = unique // 2
    feature = (cell_feature % np.uint64(total_features)).astype(np.uint32)
    cell = (cell_feature // np.uint64(total_features)).astype(np.uint32)
    return {
        "cell_index": cell,
        "feature_index": feature,
        "context_code": context,
        "observed_sites": observed.astype(np.uint32),
        "covered_calls": covered_sum.astype(np.uint32),
        "converted_calls": converted_sum.astype(np.uint32),
        "signal_fraction": (converted_sum / covered_sum).astype(np.float32),
    }


def combine_parts(parts: list[dict[str, object]], total_features: int, np):
    if not parts or not any(len(part["cell_index"]) for part in parts):
        return empty_matrix(np)
    cell = np.concatenate([part["cell_index"] for part in parts])
    feature = np.concatenate([part["feature_index"] for part in parts])
    context = np.concatenate([part["context_code"] for part in parts])
    observed = np.concatenate([part["observed_sites"] for part in parts])
    covered = np.concatenate([part["covered_calls"] for part in parts])
    converted = np.concatenate([part["converted_calls"] for part in parts])
    keys = (
        (cell.astype(np.uint64) * total_features + feature.astype(np.uint64)) * 2
        + context.astype(np.uint64)
        - 1
    )
    result = reduce_entries(keys, covered, converted, total_features, np)
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    starts = np.concatenate(
        (
            np.asarray([0], dtype=np.int64),
            np.flatnonzero(sorted_keys[1:] != sorted_keys[:-1]) + 1,
        )
    )
    observed_sum = np.add.reduceat(observed[order].astype(np.uint64), starts)
    if int(observed_sum.max()) > np.iinfo(np.uint32).max:
        raise ValueError("observed-site count exceeds uint32")
    result["observed_sites"] = observed_sum.astype(np.uint32)
    return result


def empty_matrix(np):
    return {
        "cell_index": np.asarray([], dtype=np.uint32),
        "feature_index": np.asarray([], dtype=np.uint32),
        "context_code": np.asarray([], dtype=np.uint8),
        "observed_sites": np.asarray([], dtype=np.uint32),
        "covered_calls": np.asarray([], dtype=np.uint32),
        "converted_calls": np.asarray([], dtype=np.uint32),
        "signal_fraction": np.asarray([], dtype=np.float32),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--well-index", "--sample-index", dest="well_index", type=int, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import h5py
        import numpy as np
    except ImportError as error:
        raise RuntimeError("Step 123 requires h5py and numpy") from error
    if not 1 <= args.well_index <= EXPECTED_WELLS:
        raise ValueError("well index must be in 1..96")
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
    cell_index_path = input_root / "dna_single_cytosine_cell_index.tsv.gz"
    promoter_summary_path = (
        project / "08_qc_summary/38_dna_tss_promoters/dna_tss_promoter_index_summary.tsv"
    )
    promoter_summary = read_metrics(promoter_summary_path)
    if (
        promoter_summary.get("status") != "validated"
        or promoter_summary.get("pipeline_signature") != INDEX_SIGNATURE
    ):
        raise ValueError("TSS promoter index is not current/validated")
    promoter_path = Path(promoter_summary["promoter_index"])
    feature_rows = read_tsv(promoter_path)
    total_features = len(feature_rows)
    lookup = build_interval_lookup(feature_rows, np)
    tasks = read_tsv(tasks_path)
    well_tasks = [
        row for row in tasks if int(row["sample_index"]) == args.well_index
    ]
    if (
        len(tasks) != EXPECTED_WELLS * EXPECTED_SHARDS
        or len(well_tasks) != EXPECTED_SHARDS
        or sorted(int(row["shard"]) for row in well_tasks)
        != list(range(1, EXPECTED_SHARDS + 1))
    ):
        raise ValueError("single-cytosine task table is incomplete")
    sample = well_tasks[0]["sample"]
    group = well_tasks[0]["group"]
    if any(row["sample"] != sample or row["group"] != group for row in well_tasks):
        raise ValueError("well metadata differs across shards")
    validation_path = (
        project
        / "08_qc_summary/21_dna_single_cytosine_validation"
        / sample
        / f"{sample}.single_cytosine.validation_summary.tsv"
    )
    validation = read_metrics(validation_path)
    if (
        validation.get("status") != "validated"
        or validation.get("pipeline_signature") != SOURCE_SIGNATURE
        or int(validation.get("sample_index", "0")) != args.well_index
        or validation.get("sample") != sample
        or int(validation.get("validated_shards", "0")) != EXPECTED_SHARDS
    ):
        raise ValueError(f"invalid well validation: {validation_path}")
    sample_cells = {
        int(row["cell_index"])
        for row in read_tsv(cell_index_path)
        if int(row["sample_index"]) == args.well_index
    }
    if not sample_cells:
        raise ValueError("well has no indexed cells")
    output_root = args.output_root or (
        project / "10_downstream_analysis/08_dna_tss_promoters/well_matrices"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/38_dna_tss_promoters/well_matrices"
    )
    output_path = output_root / f"{sample}.dna_tss_promoters.h5"
    summary_path = qc_root / f"{sample}.dna_tss_promoters.summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and int(existing.get("well_index", "0")) == args.well_index
            and Path(existing.get("output_h5", "")).is_file()
        ):
            print(f"[SKIP] validated well matrix exists: {summary_path}")
            return 0
        raise FileExistsError(f"inspect or rerun with --force: {summary_path}")

    parts: list[dict[str, object]] = []
    totals = {
        "source_entries": 0,
        "selected_entries": 0,
        "WCG_covered": 0,
        "WCG_converted": 0,
        "GCH_covered": 0,
        "GCH_converted": 0,
    }
    for task in sorted(well_tasks, key=lambda row: int(row["shard"])):
        shard = int(task["shard"])
        source_path = Path(task["output_dir"]) / (
            f"{sample}.shard_{shard:02d}.single_cytosine.h5"
        )
        if not source_path.is_file() or source_path.stat().st_size == 0:
            raise FileNotFoundError(source_path)
        with h5py.File(source_path, "r") as source:
            if (
                str(source.attrs.get("status")) != "validated"
                or str(source.attrs.get("sample")) != sample
                or int(source.attrs.get("shard", 0)) != shard
            ):
                raise ValueError(f"invalid source HDF5: {source_path}")
            matrix = source["matrix"]
            source_arrays = {
                "cell_index": matrix["cell_index"][:],
                "contig_index": matrix["contig_index"][:],
                "position0": matrix["cytosine_position0"][:],
                "context_code": matrix["context_code"][:],
                "covered_calls": matrix["covered"][:],
                "converted_calls": matrix["converted"][:],
            }
        if len(source_arrays["cell_index"]):
            observed_cells = set(map(int, np.unique(source_arrays["cell_index"])))
            if not observed_cells.issubset(sample_cells):
                raise ValueError(f"foreign cell index in {source_path}")
        part, metrics = aggregate_source_arrays(
            **source_arrays,
            lookup=lookup,
            total_features=total_features,
            np=np,
        )
        parts.append(part)
        for key in totals:
            totals[key] += int(metrics[key])
    combined = combine_parts(parts, total_features, np)
    if (
        totals["source_entries"]
        != int(validation["observed_cell_cytosine_entries"])
        or len(sample_cells) != int(validation["indexed_cells"])
    ):
        raise ValueError("complete well source totals do not reconcile")
    if int(combined["observed_sites"].sum(dtype=np.uint64)) != totals["selected_entries"]:
        raise ValueError("selected source entries do not reconcile")
    for code, name in CONTEXTS.items():
        take = combined["context_code"] == code
        if (
            int(combined["covered_calls"][take].sum(dtype=np.uint64))
            != totals[f"{name}_covered"]
            or int(combined["converted_calls"][take].sum(dtype=np.uint64))
            != totals[f"{name}_converted"]
            or totals[f"{name}_covered"] > int(validation[f"{name}_covered_calls"])
            or totals[f"{name}_converted"] > int(validation[f"{name}_converted_calls"])
        ):
            raise ValueError(f"{name} promoter totals do not reconcile")

    output_root.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=output_path.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    try:
        with h5py.File(temporary, "w") as output:
            matrix = output.create_group("matrix")
            for field in MATRIX_FIELDS:
                values = combined[field]
                options = {}
                if len(values):
                    options = {
                        "compression": "gzip",
                        "compression_opts": 1,
                        "shuffle": True,
                        "chunks": True,
                    }
                matrix.create_dataset(field, data=values, **options)
            output.attrs["status"] = "validated"
            output.attrs["pipeline_signature"] = PIPELINE_SIGNATURE
            output.attrs["well_index"] = args.well_index
            output.attrs["sample_index"] = args.well_index
            output.attrs["physical_well"] = sample
            output.attrs["sample"] = sample
            output.attrs["well_set"] = group
            output.attrs["total_features"] = total_features
            output.attrs["indexed_well_cells"] = len(sample_cells)
            output.attrs["matrix_entries"] = len(combined["cell_index"])
            output.attrs["missing_semantics"] = "absent_is_NA_not_zero"
            output.attrs["feature_index_path"] = str(promoter_path)
            output.attrs["context_code_1"] = "WCG_endogenous_methylation"
            output.attrs["context_code_2"] = "GCH_GpC_accessibility"
        os.replace(temporary, output_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("well_index", args.well_index),
        ("sample_index", args.well_index),
        ("physical_well", sample),
        ("sample", sample),
        ("well_set", group),
        ("promoter_index_signature", INDEX_SIGNATURE),
        ("promoter_index", promoter_path),
        ("total_reference_features", total_features),
        ("indexed_well_cells", len(sample_cells)),
        ("source_h5_shards", EXPECTED_SHARDS),
        ("source_cell_cytosine_entries", totals["source_entries"]),
        ("promoter_cell_cytosine_entries", totals["selected_entries"]),
        ("matrix_cell_feature_context_entries", len(combined["cell_index"])),
        ("WCG_covered_calls", totals["WCG_covered"]),
        ("WCG_converted_calls", totals["WCG_converted"]),
        ("GCH_covered_calls", totals["GCH_covered"]),
        ("GCH_converted_calls", totals["GCH_converted"]),
        ("missing_semantics", "absent_is_NA_not_zero"),
        ("output_h5", output_path),
        ("source_validation", validation_path),
    ]
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] well {args.well_index:02d} TSS matrix: "
        f"cells={len(sample_cells)} entries={len(combined['cell_index'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
