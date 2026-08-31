#!/usr/bin/env python3
"""Step 127: validate 96 TSS well matrices and expose one HDF5 VDS."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


EXPECTED_WELLS = 96
INDEX_SIGNATURE = "snmar_dna_tss_promoter_index_20260808"
WELL_SIGNATURE = "snmar_dna_tss_well_matrix_20260808"
PIPELINE_SIGNATURE = "snmar_dna_tss_complete_matrix_20260808"
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


def atomic_tsv(path: Path, fields: list[str], rows) -> None:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--well-root", type=Path)
    parser.add_argument("--well-qc-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import h5py
        import numpy as np
    except ImportError as error:
        raise RuntimeError("Step 127 requires h5py and numpy") from error
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    cell_index_path = input_root / "dna_single_cytosine_cell_index.tsv.gz"
    index_summary_path = (
        project / "08_qc_summary/38_dna_tss_promoters/"
        "dna_tss_promoter_index_summary.tsv"
    )
    index_summary = read_metrics(index_summary_path)
    if (
        index_summary.get("status") != "validated"
        or index_summary.get("pipeline_signature") != INDEX_SIGNATURE
    ):
        raise ValueError("TSS promoter index is not current/validated")
    feature_path = Path(index_summary["promoter_index"])
    ortholog_path = Path(index_summary["ortholog_index"])
    feature_rows = read_tsv(feature_path)
    total_features = len(feature_rows)
    if [int(row["feature_index"]) for row in feature_rows] != list(
        range(total_features)
    ):
        raise ValueError("promoter feature index is not contiguous")
    cell_rows = read_tsv(cell_index_path)
    if [int(row["cell_index"]) for row in cell_rows] != list(
        range(len(cell_rows))
    ):
        raise ValueError("canonical cell index is not contiguous")
    well_metadata: dict[int, tuple[str, str]] = {}
    cells_by_well: dict[int, set[int]] = {}
    for row in cell_rows:
        well_index = int(row["sample_index"])
        value = (row["sample"], row["group"])
        if well_metadata.setdefault(well_index, value) != value:
            raise ValueError("well metadata differs within cell index")
        cells_by_well.setdefault(well_index, set()).add(int(row["cell_index"]))
    if sorted(well_metadata) != list(range(1, EXPECTED_WELLS + 1)):
        raise ValueError("canonical index does not contain exactly 96 wells")
    well_root = args.well_root or (
        project / "10_downstream_analysis/08_dna_tss_promoters/well_matrices"
    )
    well_qc_root = args.well_qc_root or (
        project / "08_qc_summary/38_dna_tss_promoters/well_matrices"
    )
    output_root = args.output_root or (
        project / "10_downstream_analysis/08_dna_tss_promoters"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/38_dna_tss_promoters"
    )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    sources: list[tuple[Path, int, dict[str, object]]] = []
    totals = {
        "matrix_entries": 0,
        "observed_sites": 0,
        "WCG_covered": 0,
        "WCG_converted": 0,
        "GCH_covered": 0,
        "GCH_converted": 0,
    }
    offset = 0
    for well_index in range(1, EXPECTED_WELLS + 1):
        sample, group = well_metadata[well_index]
        summary_path = well_qc_root / f"{sample}.dna_tss_promoters.summary.tsv"
        metrics = read_metrics(summary_path)
        matrix_path = Path(metrics.get("output_h5", ""))
        if (
            metrics.get("status") != "validated"
            or metrics.get("pipeline_signature") != WELL_SIGNATURE
            or int(metrics.get("well_index", "0")) != well_index
            or metrics.get("physical_well") != sample
            or metrics.get("well_set") != group
            or int(metrics.get("total_reference_features", "0"))
            != total_features
            or not matrix_path.is_file()
            or matrix_path.stat().st_size == 0
        ):
            raise ValueError(f"invalid TSS well summary: {summary_path}")
        with h5py.File(matrix_path, "r") as source:
            if (
                str(source.attrs.get("status")) != "validated"
                or str(source.attrs.get("pipeline_signature"))
                != WELL_SIGNATURE
                or int(source.attrs.get("well_index", 0)) != well_index
                or str(source.attrs.get("physical_well")) != sample
                or int(source.attrs.get("total_features", 0))
                != total_features
            ):
                raise ValueError(f"invalid TSS well HDF5: {matrix_path}")
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
            if length != int(metrics["matrix_cell_feature_context_entries"]):
                raise ValueError(f"well matrix length differs: {matrix_path}")
            if length:
                cells = source["matrix/cell_index"][:]
                features = source["matrix/feature_index"][:]
                contexts = source["matrix/context_code"][:]
                observed = source["matrix/observed_sites"][:]
                covered = source["matrix/covered_calls"][:]
                converted = source["matrix/converted_calls"][:]
                if (
                    not set(map(int, np.unique(cells))).issubset(
                        cells_by_well[well_index]
                    )
                    or int(features.max()) >= total_features
                    or not set(map(int, np.unique(contexts))).issubset({1, 2})
                    or np.any(observed < 1)
                    or np.any(covered < 1)
                    or np.any(converted > covered)
                ):
                    raise ValueError(f"invalid TSS values: {matrix_path}")
        sources.append((matrix_path, length, dtypes))
        row = {
            "well_index": well_index,
            "physical_well": sample,
            "well_set": group,
            "matrix_h5": matrix_path,
            "row_start0": offset,
            "row_end0_exclusive": offset + length,
            "matrix_entries": length,
            "observed_sites": int(metrics["promoter_cell_cytosine_entries"]),
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

    output_path = output_root / "dna_tss_complete_matrix.vds.h5"
    descriptor, temporary = tempfile.mkstemp(
        prefix=output_path.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    try:
        with h5py.File(temporary, "w", libver="latest") as output:
            matrix = output.create_group("matrix")
            for field in MATRIX_FIELDS:
                dtype = sources[0][2][field]
                if any(spec[2][field] != dtype for spec in sources):
                    raise ValueError(f"well dtype differs for {field}")
                layout = h5py.VirtualLayout(
                    shape=(totals["matrix_entries"],), dtype=dtype
                )
                start = 0
                for path, length, _dtypes in sources:
                    if length:
                        virtual = h5py.VirtualSource(
                            str(path), f"matrix/{field}", shape=(length,)
                        )
                        layout[start : start + length] = virtual
                    start += length
                matrix.create_virtual_dataset(field, layout)
            output.attrs["status"] = "validated"
            output.attrs["pipeline_signature"] = PIPELINE_SIGNATURE
            output.attrs["aggregation_scope"] = (
                "one_batch_96_physical_wells"
            )
            output.attrs["physical_wells"] = EXPECTED_WELLS
            output.attrs["matrix_format"] = "sparse_COO_HDF5_VDS"
            output.attrs["indexed_cells"] = len(cell_rows)
            output.attrs["total_features"] = total_features
            output.attrs["matrix_entries"] = totals["matrix_entries"]
            output.attrs["missing_semantics"] = "absent_is_NA_not_zero"
            output.attrs["feature_index_path"] = str(feature_path)
            output.attrs["ortholog_index_path"] = str(ortholog_path)
            output.attrs["cell_index_path"] = str(cell_index_path)
            output.attrs["context_code_1"] = "WCG_endogenous_methylation"
            output.attrs["context_code_2"] = "GCH_GpC_accessibility"
        os.replace(temporary, output_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    manifest_path = qc_root / "dna_tss_well_manifest.tsv"
    atomic_tsv(manifest_path, list(manifest[0]), manifest)
    summary_path = qc_root / "dna_tss_complete_summary.tsv"
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("aggregation_scope", "one_batch_96_physical_wells"),
        ("physical_wells", EXPECTED_WELLS),
        ("modality", "DNA_TSS_WCG_methylation_and_GCH_accessibility"),
        ("matrix_h5_vds", output_path),
        ("matrix_format", "sparse_COO_HDF5_VDS"),
        ("promoter_index", feature_path),
        ("ortholog_index", ortholog_path),
        ("cell_index", cell_index_path),
        ("total_reference_features", total_features),
        (
            "conservative_ortholog_features",
            index_summary["conservative_ortholog_features"],
        ),
        ("indexed_cells", len(cell_rows)),
        ("validated_wells", len(manifest)),
        ("matrix_cell_feature_context_entries", totals["matrix_entries"]),
        ("promoter_cell_cytosine_entries", totals["observed_sites"]),
        ("WCG_covered_calls", totals["WCG_covered"]),
        ("WCG_converted_calls", totals["WCG_converted"]),
        ("GCH_covered_calls", totals["GCH_covered"]),
        ("GCH_converted_calls", totals["GCH_converted"]),
        ("well_manifest", manifest_path),
        ("missing_semantics", "absent_is_NA_not_zero"),
    ]
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] complete TSS VDS: wells={len(manifest)} "
        f"features={total_features} entries={totals['matrix_entries']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
