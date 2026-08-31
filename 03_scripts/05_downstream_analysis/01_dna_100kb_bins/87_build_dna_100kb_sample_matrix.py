#!/usr/bin/env python3
"""Step 87: aggregate one DNA sample into cell-by-100 kb matrices.

WCG and GCH observations remain distinct.  A missing cell/bin/context entry
means no coverage and must not be interpreted as zero signal.  The final bin
on each contig is retained and clipped to the true contig length.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


EXPECTED_SAMPLES = 96
EXPECTED_SHARDS = 8
CONTEXTS = {1: "WCG", 2: "GCH"}
SOURCE_SIGNATURE = "snmar_single_cytosine_validation_20260801"
PIPELINE_SIGNATURE = "snmar_dna_100kb_sample_20260804"
MATRIX_FIELDS = (
    "cell_index",
    "bin_index",
    "context_code",
    "observed_sites",
    "covered_calls",
    "converted_calls",
    "signal_fraction",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
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


def atomic_metrics(path: Path, rows: list[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as handle:
            handle.write("metric\tvalue\n")
            for metric, value in rows:
                handle.write(f"{metric}\t{value}\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def build_bins(
    contig_rows: list[dict[str, str]], bin_size: int
) -> tuple[list[dict[str, object]], list[int], list[int]]:
    """Return bins, first-bin offsets and contig lengths by contig index."""
    if bin_size < 1:
        raise ValueError("bin size must be positive")
    ordered = sorted(contig_rows, key=lambda row: int(row["contig_index"]))
    observed_indices = [int(row["contig_index"]) for row in ordered]
    if observed_indices != list(range(len(ordered))):
        raise ValueError("contig indices must be zero-based and contiguous")
    bins: list[dict[str, object]] = []
    offsets: list[int] = []
    lengths: list[int] = []
    seen_contigs: set[str] = set()
    for row in ordered:
        contig = row["contig"]
        species = row["species"]
        length = int(row["length"])
        contig_index = int(row["contig_index"])
        if (
            not contig
            or contig in seen_contigs
            or species not in {"human", "pig"}
            or length < 1
        ):
            raise ValueError(f"invalid contig index row: {row}")
        seen_contigs.add(contig)
        offsets.append(len(bins))
        lengths.append(length)
        for start0 in range(0, length, bin_size):
            end0 = min(start0 + bin_size, length)
            bins.append(
                {
                    "bin_index": len(bins),
                    "species": species,
                    "contig_index": contig_index,
                    "chrom": contig,
                    "start0": start0,
                    "end0": end0,
                    "width": end0 - start0,
                    "terminal_partial_bin": int(end0 - start0 < bin_size),
                }
            )
    if not bins:
        raise ValueError("no bins generated")
    return bins, offsets, lengths


def aggregate_arrays(
    cell_index,
    contig_index,
    position0,
    context_code,
    covered_calls,
    converted_calls,
    bin_offsets: list[int],
    contig_lengths: list[int],
    bin_size: int,
    total_bins: int,
):
    """Aggregate cell-cytosine COO arrays to sorted cell-bin-context COO."""
    import numpy as np

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
        raise ValueError("source matrix arrays have different lengths")
    entries = lengths.pop()
    if entries == 0:
        return {
            "cell_index": np.asarray([], dtype=np.uint32),
            "bin_index": np.asarray([], dtype=np.uint32),
            "context_code": np.asarray([], dtype=np.uint8),
            "observed_sites": np.asarray([], dtype=np.uint32),
            "covered_calls": np.asarray([], dtype=np.uint32),
            "converted_calls": np.asarray([], dtype=np.uint32),
            "signal_fraction": np.asarray([], dtype=np.float32),
        }
    cell = np.asarray(cell_index, dtype=np.uint64)
    contig = np.asarray(contig_index, dtype=np.int64)
    position = np.asarray(position0, dtype=np.int64)
    context = np.asarray(context_code, dtype=np.int64)
    covered = np.asarray(covered_calls, dtype=np.uint64)
    converted = np.asarray(converted_calls, dtype=np.uint64)
    if (
        np.any(contig < 0)
        or np.any(contig >= len(bin_offsets))
        or np.any((context != 1) & (context != 2))
        or np.any(covered < 1)
        or np.any(converted > covered)
    ):
        raise ValueError("invalid source matrix value")
    contig_length_array = np.asarray(contig_lengths, dtype=np.int64)
    if np.any(position < 0) or np.any(position >= contig_length_array[contig]):
        raise ValueError("cytosine position lies outside its contig")
    offset_array = np.asarray(bin_offsets, dtype=np.uint64)
    bin_index = offset_array[contig] + (
        position.astype(np.uint64) // np.uint64(bin_size)
    )
    if np.any(bin_index >= total_bins):
        raise ValueError("calculated bin index is out of range")
    maximum_cell = int(cell.max())
    maximum_key = ((maximum_cell * total_bins + total_bins - 1) * 2) + 1
    if maximum_key > np.iinfo(np.uint64).max:
        raise OverflowError("cell/bin key exceeds uint64")
    keys = (cell * np.uint64(total_bins) + bin_index) * 2 + (
        context.astype(np.uint64) - 1
    )
    order = np.argsort(keys, kind="stable")
    sorted_keys = keys[order]
    starts = np.concatenate(
        (
            np.asarray([0], dtype=np.int64),
            np.flatnonzero(sorted_keys[1:] != sorted_keys[:-1]) + 1,
        )
    )
    unique_keys = sorted_keys[starts]
    observed = np.diff(
        np.append(starts, np.asarray([entries], dtype=np.int64))
    ).astype(np.uint64)
    covered_sum = np.add.reduceat(covered[order], starts)
    converted_sum = np.add.reduceat(converted[order], starts)
    if (
        np.any(converted_sum > covered_sum)
        or max(
            int(observed.max()),
            int(covered_sum.max()),
            int(converted_sum.max()),
        )
        > np.iinfo(np.uint32).max
    ):
        raise ValueError("aggregated count exceeds valid uint32 range")
    context_out = (unique_keys % 2 + 1).astype(np.uint8)
    cell_bin = unique_keys // 2
    bin_out = (cell_bin % np.uint64(total_bins)).astype(np.uint32)
    cell_out = (cell_bin // np.uint64(total_bins)).astype(np.uint32)
    fractions = (converted_sum / covered_sum).astype(np.float32)
    return {
        "cell_index": cell_out,
        "bin_index": bin_out,
        "context_code": context_out,
        "observed_sites": observed.astype(np.uint32),
        "covered_calls": covered_sum.astype(np.uint32),
        "converted_calls": converted_sum.astype(np.uint32),
        "signal_fraction": fractions,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--sample-index", type=int, required=True)
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import h5py
        import numpy as np
    except ImportError as error:
        raise RuntimeError("Step 87 requires numpy and h5py") from error
    if not 1 <= args.sample_index <= EXPECTED_SAMPLES:
        raise ValueError("sample index must be in 1..96")
    if args.bin_size < 1:
        raise ValueError("bin size must be positive")
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
    cell_index_path = input_root / "dna_single_cytosine_cell_index.tsv.gz"
    contig_index_path = input_root / "dna_single_cytosine_contig_index.tsv"
    complete_summary_path = (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete/"
        "dna_single_cytosine_complete_summary.tsv"
    )
    tasks = read_tsv(tasks_path)
    sample_tasks = [
        row
        for row in tasks
        if int(row["sample_index"]) == args.sample_index
    ]
    if (
        len(tasks) != EXPECTED_SAMPLES * EXPECTED_SHARDS
        or len(sample_tasks) != EXPECTED_SHARDS
        or sorted(int(row["shard"]) for row in sample_tasks)
        != list(range(1, EXPECTED_SHARDS + 1))
    ):
        raise ValueError("single-cytosine task table is incomplete")
    sample = sample_tasks[0]["sample"]
    group = sample_tasks[0]["group"]
    if any(
        row["sample"] != sample or row["group"] != group
        for row in sample_tasks
    ):
        raise ValueError("sample task metadata differs across shards")
    complete_metrics = read_metrics(complete_summary_path)
    if complete_metrics.get("status") != "validated":
        raise ValueError("complete single-cytosine matrix is not validated")
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
        or int(validation.get("sample_index", "0")) != args.sample_index
        or validation.get("sample") != sample
        or int(validation.get("validated_shards", "0")) != EXPECTED_SHARDS
    ):
        raise ValueError(f"invalid sample validation: {validation_path}")
    contig_rows = read_tsv(contig_index_path)
    bins, bin_offsets, contig_lengths = build_bins(
        contig_rows, args.bin_size
    )
    total_bins = len(bins)
    sample_cell_rows = [
        row
        for row in read_tsv(cell_index_path)
        if int(row["sample_index"]) == args.sample_index
    ]
    sample_cells = {int(row["cell_index"]) for row in sample_cell_rows}
    if not sample_cells or len(sample_cells) != len(sample_cell_rows):
        raise ValueError("sample cell index is empty or duplicated")

    parts: dict[str, list[object]] = {field: [] for field in MATRIX_FIELDS}
    source_entries = 0
    source_covered = {1: 0, 2: 0}
    source_converted = {1: 0, 2: 0}
    source_paths: list[Path] = []
    for task in sorted(sample_tasks, key=lambda row: int(row["shard"])):
        shard = int(task["shard"])
        prefix = f"{sample}.shard_{shard:02d}"
        source_path = (
            Path(task["output_dir"]) / f"{prefix}.single_cytosine.h5"
        )
        if not source_path.is_file() or source_path.stat().st_size == 0:
            raise FileNotFoundError(source_path)
        source_paths.append(source_path)
        with h5py.File(source_path, "r") as source:
            if (
                str(source.attrs.get("status")) != "validated"
                or str(source.attrs.get("sample")) != sample
                or int(source.attrs.get("shard", 0)) != shard
            ):
                raise ValueError(f"invalid source HDF5: {source_path}")
            source_arrays = {
                "cell_index": source["matrix/cell_index"][:],
                "contig_index": source["matrix/contig_index"][:],
                "position0": source["matrix/cytosine_position0"][:],
                "context_code": source["matrix/context_code"][:],
                "covered_calls": source["matrix/covered"][:],
                "converted_calls": source["matrix/converted"][:],
            }
        if len(source_arrays["cell_index"]):
            observed_cells = set(
                map(int, np.unique(source_arrays["cell_index"]))
            )
            if not observed_cells.issubset(sample_cells):
                raise ValueError(f"foreign cell index in {source_path}")
        aggregated = aggregate_arrays(
            **source_arrays,
            bin_offsets=bin_offsets,
            contig_lengths=contig_lengths,
            bin_size=args.bin_size,
            total_bins=total_bins,
        )
        for field in MATRIX_FIELDS:
            parts[field].append(aggregated[field])
        source_entries += len(source_arrays["cell_index"])
        for context_code in CONTEXTS:
            selected = source_arrays["context_code"] == context_code
            source_covered[context_code] += int(
                source_arrays["covered_calls"][selected].sum(
                    dtype=np.uint64
                )
            )
            source_converted[context_code] += int(
                source_arrays["converted_calls"][selected].sum(
                    dtype=np.uint64
                )
            )

    combined = {
        field: np.concatenate(values)
        if values
        else np.asarray([], dtype=np.float32)
        for field, values in parts.items()
    }
    if len(combined["cell_index"]):
        keys = (
            combined["cell_index"].astype(np.uint64) * total_bins
            + combined["bin_index"].astype(np.uint64)
        ) * 2 + (combined["context_code"].astype(np.uint64) - 1)
        order = np.argsort(keys, kind="stable")
        for field in MATRIX_FIELDS:
            combined[field] = combined[field][order]
        sorted_keys = keys[order]
        if np.any(sorted_keys[1:] == sorted_keys[:-1]):
            raise ValueError("duplicate cell/bin/context across source shards")
    observed_total = int(
        combined["observed_sites"].sum(dtype=np.uint64)
    )
    if observed_total != source_entries:
        raise ValueError("100 kb observed-site totals do not reconcile")
    for context_code, context_name in CONTEXTS.items():
        selected = combined["context_code"] == context_code
        if (
            int(
                combined["covered_calls"][selected].sum(dtype=np.uint64)
            )
            != source_covered[context_code]
            or int(
                combined["converted_calls"][selected].sum(
                    dtype=np.uint64
                )
            )
            != source_converted[context_code]
            or source_covered[context_code]
            != int(validation[f"{context_name}_covered_calls"])
            or source_converted[context_code]
            != int(validation[f"{context_name}_converted_calls"])
        ):
            raise ValueError(f"{context_name} totals do not reconcile")

    output_root = args.output_root or (
        project
        / "10_downstream_analysis/01_dna_100kb/sample_matrices"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/27_dna_100kb/sample_matrices"
    )
    output_root.mkdir(parents=True, exist_ok=True)
    output_path = output_root / f"{sample}.dna_100kb.h5"
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
            output.attrs["sample_index"] = args.sample_index
            output.attrs["sample"] = sample
            output.attrs["group"] = group
            output.attrs["bin_size"] = args.bin_size
            output.attrs["total_bins"] = total_bins
            output.attrs["indexed_sample_cells"] = len(sample_cells)
            output.attrs["matrix_entries"] = len(combined["cell_index"])
            output.attrs["missing_semantics"] = "absent_is_NA_not_zero"
            output.attrs["context_code_1"] = "WCG_endogenous_methylation"
            output.attrs["context_code_2"] = "GCH_GpC_accessibility"
            output.attrs["terminal_bin_policy"] = (
                "retain_and_clip_to_contig_length"
            )
            output.attrs["cell_index_path"] = str(cell_index_path)
            output.attrs["contig_index_path"] = str(contig_index_path)
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
        ("sample_index", args.sample_index),
        ("sample", sample),
        ("group", group),
        ("bin_size", args.bin_size),
        ("total_reference_bins", total_bins),
        (
            "terminal_partial_bins",
            sum(int(row["terminal_partial_bin"]) for row in bins),
        ),
        ("indexed_sample_cells", len(sample_cells)),
        ("source_h5_shards", len(source_paths)),
        ("source_cell_cytosine_entries", source_entries),
        ("matrix_cell_bin_context_entries", len(combined["cell_index"])),
        ("observed_sites_sum", observed_total),
    ]
    for context_code, context_name in CONTEXTS.items():
        selected = combined["context_code"] == context_code
        summary_rows.extend(
            [
                (
                    f"{context_name}_matrix_entries",
                    int(selected.sum()),
                ),
                (
                    f"{context_name}_observed_sites",
                    int(
                        combined["observed_sites"][selected].sum(
                            dtype=np.uint64
                        )
                    ),
                ),
                (
                    f"{context_name}_covered_calls",
                    source_covered[context_code],
                ),
                (
                    f"{context_name}_converted_calls",
                    source_converted[context_code],
                ),
            ]
        )
    summary_rows.extend(
        [
            ("output_h5", output_path),
            ("source_validation", validation_path),
            ("missing_semantics", "absent_is_NA_not_zero"),
            (
                "terminal_bin_policy",
                "retain_and_clip_to_contig_length",
            ),
        ]
    )
    atomic_metrics(
        qc_root / f"{sample}.dna_100kb.summary.tsv", summary_rows
    )
    print(
        "[OK] DNA 100 kb sample matrix: "
        f"sample={sample} entries={len(combined['cell_index'])} "
        f"output={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
