#!/usr/bin/env python3
"""Step 90: combine one batch's 96 DNA-well 100 kb matrices as a VDS."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


EXPECTED_SAMPLES = 96
SAMPLE_SIGNATURE = "snmar_dna_100kb_sample_20260804"
PIPELINE_SIGNATURE = "snmar_dna_100kb_complete_20260804"
CONTEXTS = {1: "WCG", 2: "GCH"}
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


def atomic_tsv(
    path: Path, fields: list[str], rows: list[dict[str, object]], gzip_output=False
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        opener = gzip.open if gzip_output else open
        with opener(
            temporary, "wt", encoding="utf-8", newline=""
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


def atomic_metrics(path: Path, rows: list[tuple[str, object]]) -> None:
    atomic_tsv(
        path,
        ["metric", "value"],
        [{"metric": metric, "value": value} for metric, value in rows],
    )


def build_bins(
    contig_rows: list[dict[str, str]], bin_size: int
) -> list[dict[str, object]]:
    ordered = sorted(contig_rows, key=lambda row: int(row["contig_index"]))
    if [int(row["contig_index"]) for row in ordered] != list(
        range(len(ordered))
    ):
        raise ValueError("contig indices must be zero-based and contiguous")
    bins: list[dict[str, object]] = []
    for row in ordered:
        length = int(row["length"])
        if length < 1 or row["species"] not in {"human", "pig"}:
            raise ValueError(f"invalid contig row: {row}")
        for start0 in range(0, length, bin_size):
            end0 = min(start0 + bin_size, length)
            bins.append(
                {
                    "bin_index": len(bins),
                    "species": row["species"],
                    "contig_index": int(row["contig_index"]),
                    "chrom": row["contig"],
                    "start0": start0,
                    "end0": end0,
                    "width": end0 - start0,
                    "terminal_partial_bin": int(end0 - start0 < bin_size),
                }
            )
    if not bins:
        raise ValueError("no bins generated")
    return bins


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--sample-root", type=Path)
    parser.add_argument("--sample-qc-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import h5py
        import numpy as np
    except ImportError as error:
        raise RuntimeError("Step 90 requires numpy and h5py") from error
    if args.bin_size < 1:
        raise ValueError("bin size must be positive")
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    cell_index_path = input_root / "dna_single_cytosine_cell_index.tsv.gz"
    contig_index_path = input_root / "dna_single_cytosine_contig_index.tsv"
    complete_summary_path = (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete/"
        "dna_single_cytosine_complete_summary.tsv"
    )
    sample_root = args.sample_root or (
        project
        / "10_downstream_analysis/01_dna_100kb/sample_matrices"
    )
    sample_qc_root = args.sample_qc_root or (
        project / "08_qc_summary/27_dna_100kb/sample_matrices"
    )
    output_root = args.output_root or (
        project / "10_downstream_analysis/01_dna_100kb"
    )
    qc_root = args.qc_root or project / "08_qc_summary/27_dna_100kb"
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    complete = read_metrics(complete_summary_path)
    if complete.get("status") != "validated":
        raise ValueError("complete single-cytosine matrix is not validated")
    cell_rows = read_tsv(cell_index_path)
    cell_indices = sorted(int(row["cell_index"]) for row in cell_rows)
    if cell_indices != list(range(len(cell_rows))):
        raise ValueError("global cell index is not zero-based/contiguous")
    sample_metadata: dict[int, tuple[str, str]] = {}
    for row in cell_rows:
        sample_index = int(row["sample_index"])
        value = (row["sample"], row["group"])
        previous = sample_metadata.setdefault(sample_index, value)
        if previous != value:
            raise ValueError("sample metadata differs within cell index")
    if sorted(sample_metadata) != list(range(1, EXPECTED_SAMPLES + 1)):
        raise ValueError("cell index does not contain exactly 96 physical wells")
    bins = build_bins(read_tsv(contig_index_path), args.bin_size)
    total_bins = len(bins)

    manifest: list[dict[str, object]] = []
    totals = {
        "matrix_entries": 0,
        "observed_sites": 0,
        "WCG_covered": 0,
        "WCG_converted": 0,
        "GCH_covered": 0,
        "GCH_converted": 0,
    }
    source_specs: list[tuple[Path, int, dict[str, object]]] = []
    offset = 0
    for sample_index in range(1, EXPECTED_SAMPLES + 1):
        sample, group = sample_metadata[sample_index]
        summary_path = (
            sample_qc_root / f"{sample}.dna_100kb.summary.tsv"
        )
        metrics = read_metrics(summary_path)
        matrix_path = Path(metrics.get("output_h5", ""))
        if (
            metrics.get("status") != "validated"
            or metrics.get("pipeline_signature") != SAMPLE_SIGNATURE
            or int(metrics.get("sample_index", "0")) != sample_index
            or metrics.get("sample") != sample
            or metrics.get("group") != group
            or int(metrics.get("bin_size", "0")) != args.bin_size
            or int(metrics.get("total_reference_bins", "0")) != total_bins
            or not matrix_path.is_file()
            or matrix_path.stat().st_size == 0
        ):
                raise ValueError(f"invalid well 100 kb summary: {summary_path}")
        with h5py.File(matrix_path, "r") as source:
            if (
                str(source.attrs.get("status")) != "validated"
                or str(source.attrs.get("pipeline_signature"))
                != SAMPLE_SIGNATURE
                or int(source.attrs.get("sample_index", 0)) != sample_index
                or str(source.attrs.get("sample")) != sample
                or int(source.attrs.get("bin_size", 0)) != args.bin_size
                or int(source.attrs.get("total_bins", 0)) != total_bins
            ):
                raise ValueError(f"invalid well HDF5: {matrix_path}")
            shapes = {
                field: source[f"matrix/{field}"].shape
                for field in MATRIX_FIELDS
            }
            dtypes = {
                field: source[f"matrix/{field}"].dtype
                for field in MATRIX_FIELDS
            }
            if len(set(shapes.values())) != 1:
                raise ValueError(f"well matrix shapes differ: {matrix_path}")
            length = shapes["cell_index"][0]
            if length != int(metrics["matrix_cell_bin_context_entries"]):
                raise ValueError(f"well matrix length differs: {matrix_path}")
            if length:
                cells = source["matrix/cell_index"][:]
                bins_observed = source["matrix/bin_index"][:]
                contexts = source["matrix/context_code"][:]
                covered = source["matrix/covered_calls"][:]
                converted = source["matrix/converted_calls"][:]
                observed = source["matrix/observed_sites"][:]
                expected_cells = {
                    int(row["cell_index"])
                    for row in cell_rows
                    if int(row["sample_index"]) == sample_index
                }
                if (
                    not set(map(int, np.unique(cells))).issubset(
                        expected_cells
                    )
                    or int(bins_observed.max()) >= total_bins
                    or not set(map(int, np.unique(contexts))).issubset({1, 2})
                    or np.any(observed < 1)
                    or np.any(covered < 1)
                    or np.any(converted > covered)
                ):
                    raise ValueError(f"invalid matrix values: {matrix_path}")
        source_specs.append((matrix_path, length, dtypes))
        row = {
            "sample_index": sample_index,
            "sample": sample,
            "group": group,
            "matrix_h5": matrix_path,
            "row_start0": offset,
            "row_end0_exclusive": offset + length,
            "matrix_entries": length,
            "observed_sites": int(metrics["observed_sites_sum"]),
            "WCG_covered_calls": int(metrics["WCG_covered_calls"]),
            "WCG_converted_calls": int(metrics["WCG_converted_calls"]),
            "GCH_covered_calls": int(metrics["GCH_covered_calls"]),
            "GCH_converted_calls": int(metrics["GCH_converted_calls"]),
            "summary": summary_path,
            "status": "validated",
        }
        manifest.append(row)
        offset += length
        totals["matrix_entries"] += length
        totals["observed_sites"] += int(row["observed_sites"])
        for context in ("WCG", "GCH"):
            totals[f"{context}_covered"] += int(
                row[f"{context}_covered_calls"]
            )
            totals[f"{context}_converted"] += int(
                row[f"{context}_converted_calls"]
            )

    expected = {
        "observed_sites": int(complete["observed_cell_cytosine_entries"]),
        "WCG_covered": int(complete["WCG_covered_calls"]),
        "WCG_converted": int(complete["WCG_converted_calls"]),
        "GCH_covered": int(complete["GCH_covered_calls"]),
        "GCH_converted": int(complete["GCH_converted_calls"]),
    }
    for key, value in expected.items():
        if totals[key] != value:
            raise ValueError(
                f"complete 100 kb {key} differs: {totals[key]} != {value}"
            )

    output_path = output_root / "dna_100kb_complete_matrix.vds.h5"
    descriptor, temporary = tempfile.mkstemp(
        prefix=output_path.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    try:
        with h5py.File(temporary, "w", libver="latest") as output:
            matrix = output.create_group("matrix")
            for field in MATRIX_FIELDS:
                dtype = source_specs[0][2][field]
                if any(spec[2][field] != dtype for spec in source_specs):
                    raise ValueError(f"well dtype differs for {field}")
                layout = h5py.VirtualLayout(
                    shape=(totals["matrix_entries"],), dtype=dtype
                )
                start = 0
                for path, length, _dtypes in source_specs:
                    if length:
                        source = h5py.VirtualSource(
                            str(path), f"matrix/{field}", shape=(length,)
                        )
                        layout[start : start + length] = source
                    start += length
                matrix.create_virtual_dataset(field, layout)
            output.attrs["status"] = "validated"
            output.attrs["pipeline_signature"] = PIPELINE_SIGNATURE
            output.attrs["biological_sample_scope"] = "not_inferred"
            output.attrs["aggregation_scope"] = "one_batch_96_physical_wells"
            output.attrs["physical_wells"] = EXPECTED_SAMPLES
            output.attrs["legacy_sample_field_semantics"] = (
                "physical_well_id"
            )
            output.attrs["matrix_format"] = "sparse_COO_HDF5_VDS"
            output.attrs["bin_size"] = args.bin_size
            output.attrs["total_bins"] = total_bins
            output.attrs["indexed_cells"] = len(cell_rows)
            output.attrs["matrix_entries"] = totals["matrix_entries"]
            output.attrs["missing_semantics"] = "absent_is_NA_not_zero"
            output.attrs["terminal_bin_policy"] = (
                "retain_and_clip_to_contig_length"
            )
            output.attrs["context_code_1"] = "WCG_endogenous_methylation"
            output.attrs["context_code_2"] = "GCH_GpC_accessibility"
            output.attrs["cell_index_path"] = str(cell_index_path)
            output.attrs["bin_index_path"] = str(
                output_root / "dna_100kb_bin_index.tsv.gz"
            )
        os.replace(temporary, output_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    bin_path = output_root / "dna_100kb_bin_index.tsv.gz"
    atomic_tsv(bin_path, list(bins[0]), bins, gzip_output=True)
    manifest_path = qc_root / "dna_100kb_well_manifest.tsv"
    atomic_tsv(manifest_path, list(manifest[0]), manifest)
    summary_path = qc_root / "dna_100kb_complete_summary.tsv"
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("biological_samples", "not_inferred"),
        ("aggregation_scope", "one_batch_96_physical_wells"),
        ("physical_wells", EXPECTED_SAMPLES),
        ("legacy_sample_field_semantics", "physical_well_id"),
        ("modality", "DNA_WCG_methylation_and_GCH_accessibility"),
        ("matrix_h5_vds", output_path),
        ("matrix_format", "sparse_COO_HDF5_VDS"),
        ("bin_index", bin_path),
        ("cell_index", cell_index_path),
        ("bin_size", args.bin_size),
        ("total_reference_bins", total_bins),
        (
            "terminal_partial_bins",
            sum(int(row["terminal_partial_bin"]) for row in bins),
        ),
        ("indexed_cells", len(cell_rows)),
        ("validated_wells", len(manifest)),
        ("matrix_cell_bin_context_entries", totals["matrix_entries"]),
        ("observed_sites_sum", totals["observed_sites"]),
        ("WCG_covered_calls", totals["WCG_covered"]),
        ("WCG_converted_calls", totals["WCG_converted"]),
        ("GCH_covered_calls", totals["GCH_covered"]),
        ("GCH_converted_calls", totals["GCH_converted"]),
        ("well_manifest", manifest_path),
        ("missing_semantics", "absent_is_NA_not_zero"),
        ("terminal_bin_policy", "retain_and_clip_to_contig_length"),
    ]
    atomic_metrics(summary_path, summary_rows)
    print(
        "[OK] complete DNA 100 kb VDS: "
        f"samples={len(manifest)} bins={total_bins} "
        f"entries={totals['matrix_entries']} output={output_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
