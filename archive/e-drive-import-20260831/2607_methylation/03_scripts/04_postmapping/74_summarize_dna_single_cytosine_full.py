#!/usr/bin/env python3
"""Step 74: validate all WCG/GCH outputs and create one HDF5 VDS."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


EXPECTED_SAMPLES = 96
EXPECTED_SHARDS = 8
EXPECTED_TASKS = EXPECTED_SAMPLES * EXPECTED_SHARDS
VALIDATION_SIGNATURE = "snmar_single_cytosine_validation_20260801"
MATRIX_FIELDS = (
    "cell_index",
    "contig_index",
    "cytosine_position0",
    "strand",
    "context_code",
    "covered",
    "converted",
    "signal_fraction",
)
CONTEXTS = ("WCG", "GCH")


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


def atomic_concat_gzip_tsv(
    path: Path,
    sources: list[tuple[str, Path]],
    required_fields: set[str],
) -> int:
    """Concatenate validated per-sample cell tables without loading them."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    rows_written = 0
    expected_fields: list[str] | None = None
    try:
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as output:
            writer: csv.DictWriter | None = None
            for sample, source_path in sources:
                if (
                    not source_path.is_file()
                    or source_path.stat().st_size == 0
                ):
                    raise FileNotFoundError(source_path)
                sample_rows = 0
                with gzip.open(
                    source_path, "rt", encoding="utf-8", newline=""
                ) as source:
                    reader = csv.DictReader(source, delimiter="\t")
                    if reader.fieldnames is None:
                        raise ValueError(f"missing fields: {source_path}")
                    if expected_fields is None:
                        expected_fields = reader.fieldnames
                        missing = required_fields - set(expected_fields)
                        if missing:
                            raise ValueError(
                                "per-cell table lacks current QC fields: "
                                + ", ".join(sorted(missing))
                            )
                        writer = csv.DictWriter(
                            output,
                            fieldnames=expected_fields,
                            delimiter="\t",
                            lineterminator="\n",
                        )
                        writer.writeheader()
                    elif reader.fieldnames != expected_fields:
                        raise ValueError(
                            f"per-cell fields differ: {source_path}"
                        )
                    assert writer is not None
                    for row in reader:
                        if row["sample"] != sample:
                            raise ValueError(
                                f"mixed sample in {source_path}: {row['sample']}"
                            )
                        writer.writerow(row)
                        rows_written += 1
                        sample_rows += 1
                if sample_rows == 0:
                    raise ValueError(f"empty per-cell table: {source_path}")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return rows_written


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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("Step 74 requires h5py in snTAPS_align") from error
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
    cell_index_path = (
        input_root / "dna_single_cytosine_cell_index.tsv.gz"
    )
    contig_index_path = (
        input_root / "dna_single_cytosine_contig_index.tsv"
    )
    prep_summary_path = input_root / "single_cytosine_prepare_summary.tsv"
    for path in (
        tasks_path,
        cell_index_path,
        contig_index_path,
        prep_summary_path,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    prep = read_metrics(prep_summary_path)
    tasks = read_tsv(tasks_path)
    if (
        prep.get("status") != "validated"
        or len(tasks) != EXPECTED_TASKS
        or [int(row["task_index"]) for row in tasks]
        != list(range(1, EXPECTED_TASKS + 1))
    ):
        raise ValueError("single-cytosine preparation/task table is invalid")
    output_root = (
        args.output_root
        or project
        / "07_postmapping/complete_matrices/dna_single_cytosine"
    )
    qc_root = (
        args.qc_root
        or project / "08_qc_summary/23_dna_single_cytosine_complete"
    )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)

    validation_preflight: list[dict[str, object]] = []
    for sample_index in range(1, EXPECTED_SAMPLES + 1):
        sample_tasks = [
            row
            for row in tasks
            if int(row["sample_index"]) == sample_index
        ]
        sample = sample_tasks[0]["sample"] if sample_tasks else "NA"
        validation = (
            project
            / "08_qc_summary/21_dna_single_cytosine_validation"
            / sample
            / f"{sample}.single_cytosine.validation_summary.tsv"
        )
        per_cell = validation.with_name(
            f"{sample}.single_cytosine.per_cell.tsv.gz"
        )
        state = "validated"
        detail = "ok"
        if len(sample_tasks) != EXPECTED_SHARDS:
            state = "invalid_task_count"
            detail = f"expected_8_tasks_observed_{len(sample_tasks)}"
        elif not validation.is_file() or validation.stat().st_size == 0:
            state = "missing_validation_summary"
            detail = "run_step_70"
        else:
            try:
                metrics = read_metrics(validation)
                current = (
                    metrics.get("status") == "validated"
                    and metrics.get("pipeline_signature")
                    == VALIDATION_SIGNATURE
                    and int(metrics.get("sample_index", "0"))
                    == sample_index
                    and metrics.get("sample") == sample
                    and int(metrics.get("validated_shards", "0"))
                    == EXPECTED_SHARDS
                )
            except (OSError, KeyError, TypeError, ValueError) as error:
                state = "unreadable_validation_summary"
                detail = type(error).__name__
            else:
                if not current:
                    state = "stale_or_invalid_validation_summary"
                    detail = "rerun_step_70_with_FORCE_1"
                elif not per_cell.is_file() or per_cell.stat().st_size == 0:
                    state = "missing_per_cell_table"
                    detail = "rerun_step_70_with_FORCE_1"
        validation_preflight.append(
            {
                "sample_index": sample_index,
                "sample": sample,
                "state": state,
                "detail": detail,
                "validation_summary": validation,
                "per_cell_table": per_cell,
            }
        )
    validation_preflight_path = (
        qc_root / "dna_single_cytosine_validation_preflight.tsv"
    )
    atomic_tsv(
        validation_preflight_path,
        list(validation_preflight[0]),
        validation_preflight,
    )
    incomplete = [
        row for row in validation_preflight if row["state"] != "validated"
    ]
    if incomplete:
        raise RuntimeError(
            "DNA sample validation is incomplete or stale: "
            f"{len(incomplete)}/{EXPECTED_SAMPLES} samples failed preflight. "
            f"Inspect {validation_preflight_path}. Rerun step 70 with "
            "FORCE=1, wait for 96/96 validations, run step 73, and after "
            "8/8 tracks submit steps 75 and 82 (they may run in parallel)."
        )

    sample_rows: list[dict[str, object]] = []
    per_cell_sources: list[tuple[str, Path]] = []
    totals = {
        "entries": 0,
        "WCG_covered": 0,
        "WCG_converted": 0,
        "GCH_covered": 0,
        "GCH_converted": 0,
        "GCG_ambiguous": 0,
        "GCG_ambiguous_converted": 0,
        "legacy_CpG_covered": 0,
        "cells_WCG": 0,
        "cells_GCH": 0,
        "cells_either": 0,
    }
    for sample_index in range(1, EXPECTED_SAMPLES + 1):
        sample_tasks = [
            row
            for row in tasks
            if int(row["sample_index"]) == sample_index
        ]
        if len(sample_tasks) != EXPECTED_SHARDS:
            raise ValueError(f"sample {sample_index} does not have 8 tasks")
        sample = sample_tasks[0]["sample"]
        validation = (
            project
            / "08_qc_summary/21_dna_single_cytosine_validation"
            / sample
            / f"{sample}.single_cytosine.validation_summary.tsv"
        )
        if not validation.is_file() or validation.stat().st_size == 0:
            raise FileNotFoundError(validation)
        metrics = read_metrics(validation)
        if (
            metrics.get("status") != "validated"
            or metrics.get("pipeline_signature") != VALIDATION_SIGNATURE
            or int(metrics["sample_index"]) != sample_index
            or metrics["sample"] != sample
            or int(metrics["validated_shards"]) != EXPECTED_SHARDS
        ):
            raise ValueError(f"invalid sample validation: {validation}")
        values = {
            "entries": int(metrics["observed_cell_cytosine_entries"]),
            "WCG_covered": int(metrics["WCG_covered_calls"]),
            "WCG_converted": int(metrics["WCG_converted_calls"]),
            "GCH_covered": int(metrics["GCH_covered_calls"]),
            "GCH_converted": int(metrics["GCH_converted_calls"]),
            "GCG_ambiguous": int(
                metrics["GCG_ambiguous_template_calls_excluded"]
            ),
            "GCG_ambiguous_converted": int(
                metrics[
                    "GCG_ambiguous_converted_template_calls_excluded"
                ]
            ),
            "legacy_CpG_covered": int(
                metrics["legacy_CpG_covered_calls_without_end_exclusion"]
            ),
            "cells_WCG": int(metrics["cells_with_WCG_calls"]),
            "cells_GCH": int(metrics["cells_with_GCH_calls"]),
            "cells_either": int(metrics["cells_with_either_signal"]),
        }
        if (
            values["WCG_converted"] > values["WCG_covered"]
            or values["GCH_converted"] > values["GCH_covered"]
            or values["GCG_ambiguous_converted"]
            > values["GCG_ambiguous"]
            or values["entries"] < 0
        ):
            raise ValueError(f"invalid signal counts: {validation}")
        sample_rows.append(
            {
                "sample_index": sample_index,
                "sample": sample,
                "group": sample_tasks[0]["group"],
                "indexed_cells": int(metrics["indexed_cells"]),
                "cells_with_WCG_calls": values["cells_WCG"],
                "cells_with_GCH_calls": values["cells_GCH"],
                "cells_with_either_signal": values["cells_either"],
                "observed_cell_cytosine_entries": values["entries"],
                "WCG_covered_calls": values["WCG_covered"],
                "WCG_converted_calls": values["WCG_converted"],
                "WCG_methylation_fraction": metrics[
                    "WCG_methylation_fraction"
                ],
                "GCH_covered_calls": values["GCH_covered"],
                "GCH_converted_calls": values["GCH_converted"],
                "GCH_accessibility_fraction": metrics[
                    "GCH_accessibility_fraction"
                ],
                "GCG_ambiguous_template_calls_excluded": values[
                    "GCG_ambiguous"
                ],
                "GCG_ambiguous_converted_template_calls_excluded": values[
                    "GCG_ambiguous_converted"
                ],
                "legacy_CpG_covered_calls_without_end_exclusion": values[
                    "legacy_CpG_covered"
                ],
                "validation_summary": validation,
                "status": "validated",
            }
        )
        per_cell_sources.append(
            (
                sample,
                project
                / "08_qc_summary/21_dna_single_cytosine_validation"
                / sample
                / f"{sample}.single_cytosine.per_cell.tsv.gz",
            )
        )
        for key, value in values.items():
            totals[key] += value

    track_totals = {
        "WCG_covered": 0,
        "WCG_converted": 0,
        "WCG_sites": 0,
        "GCH_covered": 0,
        "GCH_converted": 0,
        "GCH_sites": 0,
    }
    track_rows: list[dict[str, object]] = []
    for shard in range(1, EXPECTED_SHARDS + 1):
        summary = (
            project
            / "08_qc_summary/22_dna_single_cytosine_tracks"
            / f"dna_single_cytosine.track_shard_{shard:02d}.summary.tsv"
        )
        if not summary.is_file() or summary.stat().st_size == 0:
            raise FileNotFoundError(summary)
        metrics = read_metrics(summary)
        if (
            metrics.get("status") != "validated"
            or int(metrics["shard"]) != shard
            or int(metrics["input_samples"]) != EXPECTED_SAMPLES
        ):
            raise ValueError(f"invalid track summary: {summary}")
        paths = {
            "detailed_track": Path(metrics["detailed_track"]),
            "WCG_bedgraph": Path(metrics["WCG_bedgraph"]),
            "GCH_bedgraph": Path(metrics["GCH_bedgraph"]),
        }
        for path in paths.values():
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        row: dict[str, object] = {
            "shard": shard,
            **paths,
            "summary": summary,
            "status": "validated",
        }
        for context in CONTEXTS:
            sites = int(metrics[f"{context}_observed_sites"])
            covered = int(metrics[f"{context}_covered_calls"])
            converted = int(metrics[f"{context}_converted_calls"])
            if converted > covered:
                raise ValueError(f"invalid {context} track counts: {summary}")
            row[f"{context}_observed_sites"] = sites
            row[f"{context}_covered_calls"] = covered
            row[f"{context}_converted_calls"] = converted
            track_totals[f"{context}_sites"] += sites
            track_totals[f"{context}_covered"] += covered
            track_totals[f"{context}_converted"] += converted
        track_rows.append(row)
    for context in CONTEXTS:
        if (
            track_totals[f"{context}_covered"]
            != totals[f"{context}_covered"]
            or track_totals[f"{context}_converted"]
            != totals[f"{context}_converted"]
        ):
            raise ValueError(
                f"{context} track and sample-validation totals differ"
            )

    vds_path = (
        output_root / "dna_single_cytosine_complete_matrix.vds.h5"
    )
    descriptor, vds_temp = tempfile.mkstemp(
        prefix=vds_path.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    shard_rows: list[dict[str, object]] = []
    source_specs: list[tuple[Path, int, dict[str, object]]] = []
    offset = 0
    try:
        for task in tasks:
            sample = task["sample"]
            shard = int(task["shard"])
            prefix = f"{sample}.shard_{shard:02d}"
            h5_path = (
                Path(task["output_dir"])
                / f"{prefix}.single_cytosine.h5"
            )
            if not h5_path.is_file() or h5_path.stat().st_size == 0:
                raise FileNotFoundError(h5_path)
            with h5py.File(h5_path, "r") as source:
                if (
                    str(source.attrs.get("status")) != "validated"
                    or str(source.attrs.get("sample")) != sample
                    or int(source.attrs.get("shard", 0)) != shard
                ):
                    raise ValueError(f"invalid shard HDF5: {h5_path}")
                shapes = {
                    field: source[f"matrix/{field}"].shape
                    for field in MATRIX_FIELDS
                }
                dtypes = {
                    field: source[f"matrix/{field}"].dtype
                    for field in MATRIX_FIELDS
                }
            if len(set(shapes.values())) != 1:
                raise ValueError(f"shard HDF5 shapes differ: {h5_path}")
            length = shapes["cell_index"][0]
            source_specs.append((h5_path, length, dtypes))
            shard_rows.append(
                {
                    "task_index": int(task["task_index"]),
                    "sample_index": int(task["sample_index"]),
                    "sample": sample,
                    "group": task["group"],
                    "shard": shard,
                    "matrix_h5": h5_path,
                    "row_start0": offset,
                    "row_end0_exclusive": offset + length,
                    "observed_cell_cytosine_entries": length,
                    "status": "validated",
                }
            )
            offset += length
        if offset != totals["entries"]:
            raise ValueError("HDF5 shard rows and validation totals differ")
        with h5py.File(vds_temp, "w", libver="latest") as output:
            group = output.create_group("matrix")
            for field in MATRIX_FIELDS:
                dtype = source_specs[0][2][field]
                if any(spec[2][field] != dtype for spec in source_specs):
                    raise ValueError(f"source dtype differs for {field}")
                layout = h5py.VirtualLayout(shape=(offset,), dtype=dtype)
                start = 0
                for source_path, length, _dtypes in source_specs:
                    if length:
                        relative = os.path.relpath(source_path, output_root)
                        source = h5py.VirtualSource(
                            relative,
                            f"/matrix/{field}",
                            shape=(length,),
                        )
                        layout[start : start + length] = source
                    start += length
                group.create_virtual_dataset(field, layout)
            output.attrs["status"] = "validated"
            output.attrs["matrix_format"] = "sparse_COO_virtual_dataset"
            output.attrs["coordinate_system"] = (
                "0-based_reference_C_plus_reference_G_minus"
            )
            output.attrs["context_code_1"] = (
                "WCG_endogenous_methylation"
            )
            output.attrs["context_code_2"] = "GCH_GpC_accessibility"
            output.attrs["ambiguous_context_policy"] = (
                "GCG_excluded_from_both_signals"
            )
            output.attrs["indexed_cells"] = int(prep["indexed_cells"])
            output.attrs["indexed_contigs"] = int(prep["indexed_contigs"])
            output.attrs["observed_cell_cytosine_entries"] = offset
            output.attrs["WCG_covered_calls"] = totals["WCG_covered"]
            output.attrs["GCH_covered_calls"] = totals["GCH_covered"]
            output.attrs["missing_semantics"] = (
                "absent_coordinate_is_NA_not_zero"
            )
            output.attrs["cell_index_path"] = os.path.relpath(
                cell_index_path, output_root
            )
            output.attrs["contig_index_path"] = os.path.relpath(
                contig_index_path, output_root
            )
        os.replace(vds_temp, vds_path)
    except BaseException:
        try:
            os.unlink(vds_temp)
        except FileNotFoundError:
            pass
        raise

    atomic_tsv(
        output_root / "dna_single_cytosine_shard_manifest.tsv",
        list(shard_rows[0]),
        shard_rows,
    )
    atomic_tsv(
        qc_root / "dna_single_cytosine_complete_per_sample.tsv",
        list(sample_rows[0]),
        sample_rows,
    )
    complete_per_cell = (
        qc_root / "dna_single_cytosine_complete_per_cell.tsv.gz"
    )
    complete_cells = atomic_concat_gzip_tsv(
        complete_per_cell,
        per_cell_sources,
        {
            "cell_index",
            "sample",
            "cell_id",
            "WCG_observed_sites",
            "WCG_covered_calls",
            "GCH_observed_sites",
            "GCH_covered_calls",
            "human_WCG_observed_sites",
            "human_GCH_observed_sites",
            "pig_WCG_observed_sites",
            "pig_GCH_observed_sites",
            "human_signal_call_fraction",
        },
    )
    expected_cells = sum(int(row["indexed_cells"]) for row in sample_rows)
    if (
        complete_cells != expected_cells
        or complete_cells != int(prep["indexed_cells"])
    ):
        raise ValueError(
            "complete per-cell row count does not match indexed cells"
        )
    atomic_tsv(
        qc_root / "dna_single_cytosine_track_manifest.tsv",
        list(track_rows[0]),
        track_rows,
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "modality", "value": "DNA_WCG_and_GCH"},
        {"metric": "matrix_h5_vds", "value": vds_path},
        {"metric": "matrix_format", "value": "sparse_COO_HDF5_VDS"},
        {"metric": "source_h5_shards", "value": EXPECTED_TASKS},
        {"metric": "validated_samples", "value": EXPECTED_SAMPLES},
        {"metric": "track_shards", "value": EXPECTED_SHARDS},
        {"metric": "indexed_cells", "value": int(prep["indexed_cells"])},
        {
            "metric": "complete_per_cell_qc",
            "value": complete_per_cell,
        },
        {"metric": "indexed_contigs", "value": int(prep["indexed_contigs"])},
        {
            "metric": "cells_with_WCG_calls",
            "value": totals["cells_WCG"],
        },
        {
            "metric": "cells_with_GCH_calls",
            "value": totals["cells_GCH"],
        },
        {
            "metric": "cells_with_either_signal",
            "value": totals["cells_either"],
        },
        {
            "metric": "observed_cell_cytosine_entries",
            "value": totals["entries"],
        },
        {
            "metric": "WCG_observed_sites_sum_across_shards",
            "value": track_totals["WCG_sites"],
        },
        {"metric": "WCG_covered_calls", "value": totals["WCG_covered"]},
        {
            "metric": "WCG_converted_calls",
            "value": totals["WCG_converted"],
        },
        {
            "metric": "WCG_methylation_fraction",
            "value": totals["WCG_converted"] / totals["WCG_covered"]
            if totals["WCG_covered"]
            else "NA",
        },
        {
            "metric": "GCH_observed_sites_sum_across_shards",
            "value": track_totals["GCH_sites"],
        },
        {"metric": "GCH_covered_calls", "value": totals["GCH_covered"]},
        {
            "metric": "GCH_converted_calls",
            "value": totals["GCH_converted"],
        },
        {
            "metric": "GCH_accessibility_fraction",
            "value": totals["GCH_converted"] / totals["GCH_covered"]
            if totals["GCH_covered"]
            else "NA",
        },
        {
            "metric": "GCG_ambiguous_template_calls_excluded",
            "value": totals["GCG_ambiguous"],
        },
        {
            "metric": "GCG_ambiguous_converted_template_calls_excluded",
            "value": totals["GCG_ambiguous_converted"],
        },
        {
            "metric": "legacy_CpG_covered_calls_without_end_exclusion",
            "value": totals["legacy_CpG_covered"],
        },
        {
            "metric": "legacy_comparison",
            "value": "reference_only_context_and_strand_model_changed",
        },
        {
            "metric": "context_model",
            "value": "WCG=methylation;GCH=accessibility;GCG=excluded",
        },
        {
            "metric": "missing_cell_cytosine_semantics",
            "value": "NA_not_zero_signal",
        },
        {
            "metric": "read_end_exclusion",
            "value": prep["excluded_bases_each_read_end"],
        },
        {"metric": "cell_index", "value": cell_index_path},
        {"metric": "contig_index", "value": contig_index_path},
    ]
    atomic_tsv(
        qc_root / "dna_single_cytosine_complete_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        "[OK] DNA context-aware VDS: "
        f"entries={totals['entries']} WCG={totals['WCG_covered']} "
        f"GCH={totals['GCH_covered']} "
        f"GCG_excluded={totals['GCG_ambiguous']} "
        f"GCG_converted={totals['GCG_ambiguous_converted']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
