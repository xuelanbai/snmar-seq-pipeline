#!/usr/bin/env python3
"""Step 136: aggregate one physical well into ortholog-100 kb matrices.

The two physical intervals for each orthobin remain separate in this HDF5.
They are collapsed onto one shared axis only after cell filtering, without a
human/pig label entering matrix construction.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import tempfile
from pathlib import Path


EXPECTED_WELLS = 96
EXPECTED_SHARDS = 8
SOURCE_SIGNATURE = "snmar_single_cytosine_validation_20260801"
INDEX_SIGNATURE = "snmar_dna_ortholog_100kb_index_20260825"
PIPELINE_SIGNATURE = "snmar_dna_ortholog_100kb_well_20260825"
MATRIX_FIELDS = (
    "cell_index",
    "feature_index",
    "context_code",
    "observed_sites",
    "covered_calls",
    "converted_calls",
    "signal_fraction",
)
CONTEXTS = {1: "WCG", 2: "GCH"}


def load_interval_helpers():
    path = (
        Path(__file__).resolve().parents[1]
        / "06_dna_tss_umap/123_build_dna_tss_well_matrix.py"
    )
    if not path.is_file():
        raise FileNotFoundError(f"missing shared interval helper: {path}")
    spec = importlib.util.spec_from_file_location("snmar_interval_helpers", path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise ImportError(path)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument(
        "--well-index", "--sample-index", dest="well_index", type=int, required=True
    )
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
        raise RuntimeError("Step 136 requires h5py and numpy") from error
    if not 1 <= args.well_index <= EXPECTED_WELLS:
        raise ValueError("well index must be in 1..96")
    helper = load_interval_helpers()
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
    cell_index_path = input_root / "dna_single_cytosine_cell_index.tsv.gz"
    index_summary_path = (
        project
        / "08_qc_summary/40_dna_ortholog_100kb/"
        "dna_ortholog_100kb_index_summary.tsv"
    )
    index_summary = helper.read_metrics(index_summary_path)
    if (
        index_summary.get("status") != "validated"
        or index_summary.get("pipeline_signature") != INDEX_SIGNATURE
    ):
        raise ValueError("ortholog 100 kb index is not current/validated")
    feature_path = Path(index_summary["physical_feature_index"])
    feature_rows = helper.read_tsv(feature_path)
    total_features = len(feature_rows)
    lookup = helper.build_interval_lookup(feature_rows, np)
    tasks = helper.read_tsv(tasks_path)
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
        raise ValueError("physical-well metadata differs across shards")
    validation_path = (
        project
        / "08_qc_summary/21_dna_single_cytosine_validation"
        / sample
        / f"{sample}.single_cytosine.validation_summary.tsv"
    )
    validation = helper.read_metrics(validation_path)
    if (
        validation.get("status") != "validated"
        or validation.get("pipeline_signature") != SOURCE_SIGNATURE
        or int(validation.get("sample_index", "0")) != args.well_index
        or validation.get("sample") != sample
        or int(validation.get("validated_shards", "0")) != EXPECTED_SHARDS
    ):
        raise ValueError(f"invalid physical-well validation: {validation_path}")
    sample_cells = {
        int(row["cell_index"])
        for row in helper.read_tsv(cell_index_path)
        if int(row["sample_index"]) == args.well_index
    }
    if not sample_cells:
        raise ValueError("physical well has no indexed cells")
    output_root = args.output_root or (
        project / "10_downstream_analysis/10_dna_ortholog_100kb/well_matrices"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/40_dna_ortholog_100kb/well_matrices"
    )
    output_path = output_root / f"{sample}.dna_ortholog_100kb.h5"
    summary_path = qc_root / f"{sample}.dna_ortholog_100kb.summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = helper.read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and int(existing.get("well_index", "0")) == args.well_index
            and Path(existing.get("output_h5", "")).is_file()
        ):
            print(f"[SKIP] validated orthobin well matrix exists: {summary_path}")
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
            arrays = {
                "cell_index": matrix["cell_index"][:],
                "contig_index": matrix["contig_index"][:],
                "position0": matrix["cytosine_position0"][:],
                "context_code": matrix["context_code"][:],
                "covered_calls": matrix["covered"][:],
                "converted_calls": matrix["converted"][:],
            }
        if len(arrays["cell_index"]):
            observed = set(map(int, np.unique(arrays["cell_index"])))
            if not observed.issubset(sample_cells):
                raise ValueError(f"foreign cell index in {source_path}")
        part, metrics = helper.aggregate_source_arrays(
            **arrays,
            lookup=lookup,
            total_features=total_features,
            np=np,
        )
        parts.append(part)
        for key in totals:
            totals[key] += int(metrics[key])
    combined = helper.combine_parts(parts, total_features, np)
    if (
        totals["source_entries"]
        != int(validation["observed_cell_cytosine_entries"])
        or len(sample_cells) != int(validation["indexed_cells"])
    ):
        raise ValueError("complete well source totals do not reconcile")
    if int(combined["observed_sites"].sum(dtype=np.uint64)) != totals["selected_entries"]:
        raise ValueError("selected orthobin entries do not reconcile")
    for code, name in CONTEXTS.items():
        take = combined["context_code"] == code
        if (
            int(combined["covered_calls"][take].sum(dtype=np.uint64))
            != totals[f"{name}_covered"]
            or int(combined["converted_calls"][take].sum(dtype=np.uint64))
            != totals[f"{name}_converted"]
            or totals[f"{name}_covered"] > int(validation[f"{name}_covered_calls"])
        ):
            raise ValueError(f"{name} orthobin totals do not reconcile")

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
            output.attrs["physical_well"] = sample
            output.attrs["well_set"] = group
            output.attrs["total_physical_features"] = total_features
            output.attrs["orthobin_axes"] = int(index_summary["orthobin_axes"])
            output.attrs["indexed_well_cells"] = len(sample_cells)
            output.attrs["matrix_entries"] = len(combined["cell_index"])
            output.attrs["missing_semantics"] = "absent_is_NA_not_zero"
            output.attrs["feature_index_path"] = str(feature_path)
            output.attrs["context_code_1"] = "WCG_endogenous_methylation"
            output.attrs["context_code_2"] = "GCH_GpC_accessibility_not_used_in_primary_mC"
        os.replace(temporary, output_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    helper.atomic_metrics(
        summary_path,
        [
            ("status", "validated"),
            ("pipeline_signature", PIPELINE_SIGNATURE),
            ("well_index", args.well_index),
            ("physical_well", sample),
            ("sample", sample),
            ("well_set", group),
            ("orthobin_index_signature", INDEX_SIGNATURE),
            ("physical_feature_index", feature_path),
            ("total_physical_features", total_features),
            ("orthobin_axes", index_summary["orthobin_axes"]),
            ("indexed_well_cells", len(sample_cells)),
            ("source_h5_shards", EXPECTED_SHARDS),
            ("source_cell_cytosine_entries", totals["source_entries"]),
            ("orthobin_cell_cytosine_entries", totals["selected_entries"]),
            ("matrix_cell_feature_context_entries", len(combined["cell_index"])),
            ("WCG_covered_calls", totals["WCG_covered"]),
            ("WCG_converted_calls", totals["WCG_converted"]),
            ("GCH_covered_calls", totals["GCH_covered"]),
            ("GCH_converted_calls", totals["GCH_converted"]),
            ("missing_semantics", "absent_is_NA_not_zero"),
            ("output_h5", output_path),
            ("source_validation", validation_path),
        ],
    )
    print(
        f"[OK] well {args.well_index:02d} ortholog 100 kb matrix: "
        f"cells={len(sample_cells)} entries={len(combined['cell_index'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
