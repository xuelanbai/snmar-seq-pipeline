#!/usr/bin/env python3
"""Step 78: final integrity summary for complete DNA and RNA matrices."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        import h5py
    except ImportError as error:
        raise RuntimeError("Step 78 requires h5py in snTAPS_align") from error
    project = args.project
    dna_summary_path = (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete"
        / "dna_single_cytosine_complete_summary.tsv"
    )
    rna_summary_path = (
        project
        / "08_qc_summary/20_complete_matrices/rna_expression"
        / "rna_complete_expression_matrix_summary.tsv"
    )
    for path in (dna_summary_path, rna_summary_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    dna = read_metrics(dna_summary_path)
    rna = read_metrics(rna_summary_path)
    if dna.get("status") != "validated" or rna.get("status") != "validated":
        raise ValueError("DNA or RNA complete matrix is not validated")
    dna_h5 = Path(dna["matrix_h5_vds"])
    rna_h5 = Path(rna["sparse_h5"])
    for path in (dna_h5, rna_h5):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    with h5py.File(dna_h5, "r") as source:
        if (
            str(source.attrs.get("status")) != "validated"
            or int(source.attrs["indexed_cells"])
            != int(dna["indexed_cells"])
            or int(source.attrs["observed_cell_cytosine_entries"])
            != int(dna["observed_cell_cytosine_entries"])
        ):
            raise ValueError("DNA VDS attributes differ from summary")
        dna_lengths = {
            len(source[f"matrix/{field}"])
            for field in (
                "cell_index",
                "contig_index",
                "cytosine_position0",
                "strand",
                "context_code",
                "covered",
                "converted",
                "signal_fraction",
            )
        }
        if dna_lengths != {int(dna["observed_cell_cytosine_entries"])}:
            raise ValueError("DNA VDS arrays differ in length")
        if (
            str(source.attrs.get("context_code_1"))
            != "WCG_endogenous_methylation"
            or str(source.attrs.get("context_code_2"))
            != "GCH_GpC_accessibility"
            or str(source.attrs.get("ambiguous_context_policy"))
            != "GCG_excluded_from_both_signals"
        ):
            raise ValueError("DNA VDS context metadata is invalid")
    with h5py.File(rna_h5, "r") as source:
        if str(source.attrs.get("status")) != "validated":
            raise ValueError("RNA HDF5 is not validated")
        for role, expected_cells, expected_features in (
            ("raw", int(rna["raw_cells"]), int(rna["raw_features"])),
            (
                "filtered",
                int(rna["filtered_cells"]),
                int(rna["filtered_features"]),
            ),
        ):
            group = source[role]
            lengths = {
                len(group[field])
                for field in ("feature_index", "cell_index", "exact_umis")
            }
            if (
                len(lengths) != 1
                or int(group.attrs["cells"]) != expected_cells
                or int(group.attrs["features"]) != expected_features
            ):
                raise ValueError(f"RNA {role} HDF5 differs from summary")

    output_root = (
        args.output_root
        or project / "08_qc_summary/24_complete_modalities"
    )
    modality_rows: list[dict[str, object]] = [
        {
            "modality": "DNA",
            "matrix_role": "WCG_endogenous_methylation_cell_by_site",
            "primary_matrix": dna_h5,
            "cells": int(dna["indexed_cells"]),
            "features": int(dna["WCG_observed_sites_sum_across_shards"]),
            "nonzero_or_observed_entries": int(
                dna["observed_cell_cytosine_entries"]
            ),
            "value": "signal_fraction_where_context_code=1",
            "missing_value": "NA",
            "status": "downstream_ready",
        },
        {
            "modality": "DNA",
            "matrix_role": "GCH_GpC_accessibility_cell_by_site",
            "primary_matrix": dna_h5,
            "cells": int(dna["indexed_cells"]),
            "features": int(dna["GCH_observed_sites_sum_across_shards"]),
            "nonzero_or_observed_entries": int(
                dna["observed_cell_cytosine_entries"]
            ),
            "value": "signal_fraction_where_context_code=2",
            "missing_value": "NA",
            "status": "downstream_ready",
        },
        {
            "modality": "RNA",
            "matrix_role": "raw_gene_by_cell",
            "primary_matrix": rna_h5,
            "cells": int(rna["raw_cells"]),
            "features": int(rna["raw_features"]),
            "nonzero_or_observed_entries": "see_RNA_HDF5",
            "value": "exact_UMI_count",
            "missing_value": "zero_expression",
            "status": "downstream_ready",
        },
        {
            "modality": "RNA",
            "matrix_role": "filtered_gene_by_cell",
            "primary_matrix": rna_h5,
            "cells": int(rna["filtered_cells"]),
            "features": int(rna["filtered_features"]),
            "nonzero_or_observed_entries": "see_RNA_HDF5",
            "value": "exact_UMI_count",
            "missing_value": "zero_expression",
            "status": "downstream_ready",
        },
    ]
    atomic_tsv(
        output_root / "complete_modality_matrix_manifest.tsv",
        list(modality_rows[0]),
        modality_rows,
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "DNA_WCG_matrix_status", "value": "downstream_ready"},
        {"metric": "DNA_GCH_matrix_status", "value": "downstream_ready"},
        {
            "metric": "DNA_primary_resolution",
            "value": "single_cytosine_context_and_strand_aware",
        },
        {"metric": "DNA_matrix", "value": dna_h5},
        {
            "metric": "DNA_context_semantics",
            "value": "WCG=methylation;GCH=GpC_accessibility;GCG=excluded",
        },
        {
            "metric": "DNA_converted_semantics",
            "value": "forward_T_or_reverse_A_after_TAPS",
        },
        {"metric": "DNA_missing_semantics", "value": "NA_not_zero"},
        {"metric": "RNA_raw_matrix_status", "value": "downstream_ready"},
        {
            "metric": "RNA_filtered_matrix_status",
            "value": "downstream_ready",
        },
        {"metric": "RNA_matrix", "value": rna_h5},
        {
            "metric": "cross_modality_cell_pairing",
            "value": "not_assumed",
        },
        {
            "metric": "downstream_derivations",
            "value": "100kb_bins,TSS,gene_body_and_other_analysis_steps",
        },
        {"metric": "biological_interpretation", "value": "not_performed"},
    ]
    atomic_tsv(
        output_root / "complete_modalities_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print("[OK] complete DNA and RNA matrices are downstream-ready")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
