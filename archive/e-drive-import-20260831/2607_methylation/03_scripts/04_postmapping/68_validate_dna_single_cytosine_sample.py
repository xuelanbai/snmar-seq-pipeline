#!/usr/bin/env python3
"""Step 68: validate eight context-aware cytosine shards for one sample."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path
from typing import Iterable


EXPECTED_SHARDS = 8
PIPELINE_SIGNATURE = "snmar_single_cytosine_validation_20260801"
CONTEXT_WCG = 1
CONTEXT_GCH = 2
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


def add_by_cell(target, cells, values, first_cell: int, cell_count: int) -> None:
    import numpy as np

    if len(cells):
        target += np.bincount(
            cells - first_cell,
            weights=values,
            minlength=cell_count,
        ).astype(np.uint64)


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
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.sample_index < 1 or args.sample_index > 96:
        raise ValueError("sample-index must be 1..96")
    try:
        import h5py
        import numpy as np
    except ImportError as error:
        raise RuntimeError(
            "Step 68 requires numpy and h5py in snTAPS_align"
        ) from error
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
    if prep.get("status") != "validated":
        raise ValueError("single-cytosine preparation is not validated")
    task_rows = [
        row
        for row in read_tsv(tasks_path)
        if int(row["sample_index"]) == args.sample_index
    ]
    task_rows.sort(key=lambda row: int(row["shard"]))
    if (
        len(task_rows) != EXPECTED_SHARDS
        or [int(row["shard"]) for row in task_rows]
        != list(range(1, EXPECTED_SHARDS + 1))
    ):
        raise ValueError("sample does not have eight ordered tasks")
    sample = task_rows[0]["sample"]
    if len({row["sample"] for row in task_rows}) != 1:
        raise ValueError("sample tasks contain mixed sample names")

    cell_rows = []
    with gzip.open(
        cell_index_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] == sample:
                cell_rows.append(row)
    if not cell_rows:
        raise ValueError(f"sample has no indexed cells: {sample}")
    cell_rows.sort(key=lambda row: int(row["cell_index"]))
    cell_indices = np.asarray(
        [int(row["cell_index"]) for row in cell_rows], dtype=np.int64
    )
    if not np.array_equal(
        cell_indices,
        np.arange(cell_indices[0], cell_indices[-1] + 1),
    ):
        raise ValueError("sample cell indices are not contiguous")
    first_cell = int(cell_indices[0])
    cell_count = len(cell_rows)
    legacy_covered = np.asarray(
        [int(row["legacy_covered_template_cpg_calls"]) for row in cell_rows],
        dtype=np.uint64,
    )
    wcg_covered = np.zeros(cell_count, dtype=np.uint64)
    wcg_converted = np.zeros(cell_count, dtype=np.uint64)
    wcg_sites = np.zeros(cell_count, dtype=np.uint64)
    gch_covered = np.zeros(cell_count, dtype=np.uint64)
    gch_converted = np.zeros(cell_count, dtype=np.uint64)
    gch_sites = np.zeros(cell_count, dtype=np.uint64)
    species_context_stats = {
        (species, context, metric): np.zeros(
            cell_count, dtype=np.uint64
        )
        for species in ("human", "pig")
        for context in ("WCG", "GCH")
        for metric in ("observed_sites", "covered_calls", "converted_calls")
    }
    shard_contigs: dict[int, set[int]] = {}
    contig_lengths: dict[int, int] = {}
    contig_species: dict[int, str] = {}
    for row in read_tsv(contig_index_path):
        contig_index = int(row["contig_index"])
        species = row["species"]
        if species not in {"human", "pig"}:
            raise ValueError(f"invalid contig species: {row}")
        shard_contigs.setdefault(int(row["shard"]), set()).add(contig_index)
        contig_lengths[contig_index] = int(row["length"])
        contig_species[contig_index] = species
    contig_length_lookup = np.zeros(
        max(contig_lengths) + 1, dtype=np.int64
    )
    contig_species_lookup = np.zeros(
        max(contig_lengths) + 1, dtype=np.uint8
    )
    for contig_index, length in contig_lengths.items():
        contig_length_lookup[contig_index] = length
        contig_species_lookup[contig_index] = (
            1 if contig_species[contig_index] == "human" else 2
        )

    totals = {
        "entries": 0,
        "WCG_covered": 0,
        "WCG_converted": 0,
        "GCH_covered": 0,
        "GCH_converted": 0,
        "GCG_ambiguous": 0,
        "GCG_ambiguous_converted": 0,
    }
    shard_manifest = []
    for task in task_rows:
        shard = int(task["shard"])
        output_dir = Path(task["output_dir"])
        prefix = f"{sample}.shard_{shard:02d}"
        h5_path = output_dir / f"{prefix}.single_cytosine.h5"
        site_path = output_dir / f"{prefix}.site_summary.tsv.gz"
        summary_path = output_dir / f"{prefix}.summary.tsv"
        for path in (h5_path, site_path, summary_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        summary = read_metrics(summary_path)
        expected_mapq = int(task["minimum_mapq"])
        expected_baseq = int(task["minimum_baseq"])
        expected_ends = int(task["excluded_bases_each_read_end"])
        if (
            summary.get("status") != "validated"
            or summary.get("sample") != sample
            or int(summary.get("shard", "0")) != shard
            or int(summary["minimum_mapq"]) != expected_mapq
            or int(summary["minimum_baseq"]) != expected_baseq
            or int(summary["excluded_bases_each_read_end"]) != expected_ends
        ):
            raise ValueError(f"invalid caller summary: {summary_path}")
        with h5py.File(h5_path, "r") as source:
            if (
                str(source.attrs.get("status")) != "validated"
                or str(source.attrs.get("sample")) != sample
                or int(source.attrs.get("shard", 0)) != shard
                or int(source.attrs.get("minimum_mapq", -1))
                != expected_mapq
                or int(source.attrs.get("minimum_baseq", -1))
                != expected_baseq
                or int(
                    source.attrs.get("excluded_bases_each_read_end", -1)
                )
                != expected_ends
                or str(source.attrs.get("ambiguous_context_policy"))
                != "GCG_excluded_from_both_signals"
            ):
                raise ValueError(f"invalid HDF5 attributes: {h5_path}")
            group = source["matrix"]
            arrays = {field: group[field][:] for field in MATRIX_FIELDS}
        lengths = {len(values) for values in arrays.values()}
        if len(lengths) != 1:
            raise ValueError(f"HDF5 arrays differ in length: {h5_path}")
        entries = lengths.pop()
        cells = arrays["cell_index"].astype(np.int64, copy=False)
        contig_indices = arrays["contig_index"].astype(
            np.int64, copy=False
        )
        positions = arrays["cytosine_position0"].astype(
            np.int64, copy=False
        )
        strands = arrays["strand"].astype(np.int8, copy=False)
        contexts = arrays["context_code"].astype(np.uint8, copy=False)
        covered = arrays["covered"].astype(np.uint64, copy=False)
        converted = arrays["converted"].astype(np.uint64, copy=False)
        fractions = arrays["signal_fraction"]
        if entries:
            if (
                cells.min() < first_cell
                or cells.max() >= first_cell + cell_count
                or not set(np.unique(contig_indices)).issubset(
                    shard_contigs[shard]
                )
                or np.any(positions < 0)
                or np.any(
                    positions >= contig_length_lookup[contig_indices]
                )
                or not set(np.unique(strands)).issubset({-1, 1})
                or not set(np.unique(contexts)).issubset(
                    {CONTEXT_WCG, CONTEXT_GCH}
                )
                or np.any(covered == 0)
                or np.any(converted > covered)
                or not np.allclose(
                    fractions,
                    converted / covered,
                    rtol=1e-6,
                    atol=1e-7,
                )
            ):
                raise ValueError(f"invalid HDF5 values: {h5_path}")
            observed_order = np.lexsort(
                (cells, contexts, strands, positions, contig_indices)
            )
            if not np.array_equal(
                observed_order, np.arange(entries, dtype=observed_order.dtype)
            ):
                raise ValueError(
                    f"HDF5 coordinates are not sorted: {h5_path}"
                )
            duplicate = (
                (contig_indices[1:] == contig_indices[:-1])
                & (positions[1:] == positions[:-1])
                & (strands[1:] == strands[:-1])
                & (contexts[1:] == contexts[:-1])
                & (cells[1:] == cells[:-1])
            )
            if np.any(duplicate):
                raise ValueError(
                    f"HDF5 contains duplicate cell-site rows: {h5_path}"
                )
        wcg = contexts == CONTEXT_WCG
        gch = contexts == CONTEXT_GCH
        ones = np.ones(entries, dtype=np.uint64)
        add_by_cell(
            wcg_sites,
            cells[wcg],
            ones[wcg],
            first_cell,
            cell_count,
        )
        add_by_cell(
            wcg_covered,
            cells[wcg],
            covered[wcg],
            first_cell,
            cell_count,
        )
        add_by_cell(
            wcg_converted,
            cells[wcg],
            converted[wcg],
            first_cell,
            cell_count,
        )
        add_by_cell(
            gch_sites,
            cells[gch],
            ones[gch],
            first_cell,
            cell_count,
        )
        add_by_cell(
            gch_covered,
            cells[gch],
            covered[gch],
            first_cell,
            cell_count,
        )
        species_codes = contig_species_lookup[contig_indices]
        for species, species_code in (("human", 1), ("pig", 2)):
            for context_name, context_mask in (
                ("WCG", wcg),
                ("GCH", gch),
            ):
                selected = context_mask & (species_codes == species_code)
                add_by_cell(
                    species_context_stats[
                        (species, context_name, "observed_sites")
                    ],
                    cells[selected],
                    ones[selected],
                    first_cell,
                    cell_count,
                )
                add_by_cell(
                    species_context_stats[
                        (species, context_name, "covered_calls")
                    ],
                    cells[selected],
                    covered[selected],
                    first_cell,
                    cell_count,
                )
                add_by_cell(
                    species_context_stats[
                        (species, context_name, "converted_calls")
                    ],
                    cells[selected],
                    converted[selected],
                    first_cell,
                    cell_count,
                )
        add_by_cell(
            gch_converted,
            cells[gch],
            converted[gch],
            first_cell,
            cell_count,
        )
        shard_wcg_covered = int(covered[wcg].sum(dtype=np.uint64))
        shard_wcg_converted = int(converted[wcg].sum(dtype=np.uint64))
        shard_gch_covered = int(covered[gch].sum(dtype=np.uint64))
        shard_gch_converted = int(converted[gch].sum(dtype=np.uint64))
        if (
            entries != int(summary["observed_cell_cytosine_entries"])
            or shard_wcg_covered != int(summary["WCG_covered_calls"])
            or shard_wcg_converted != int(summary["WCG_converted_calls"])
            or shard_gch_covered != int(summary["GCH_covered_calls"])
            or shard_gch_converted != int(summary["GCH_converted_calls"])
        ):
            raise ValueError(f"HDF5 and summary differ: {h5_path}")
        ambiguous = int(
            summary.get("ambiguous_GCG_template_calls_excluded", "0")
        )
        ambiguous_converted = int(
            summary.get(
                "ambiguous_GCG_converted_template_calls_excluded", "0"
            )
        )
        if ambiguous_converted > ambiguous:
            raise ValueError(f"invalid GCG audit counts: {summary_path}")
        totals["entries"] += entries
        totals["WCG_covered"] += shard_wcg_covered
        totals["WCG_converted"] += shard_wcg_converted
        totals["GCH_covered"] += shard_gch_covered
        totals["GCH_converted"] += shard_gch_converted
        totals["GCG_ambiguous"] += ambiguous
        totals["GCG_ambiguous_converted"] += ambiguous_converted
        shard_manifest.append(
            {
                "sample_index": args.sample_index,
                "sample": sample,
                "shard": shard,
                "matrix_h5": h5_path,
                "site_summary": site_path,
                "summary": summary_path,
                "observed_cell_cytosine_entries": entries,
                "WCG_covered_calls": shard_wcg_covered,
                "WCG_converted_calls": shard_wcg_converted,
                "GCH_covered_calls": shard_gch_covered,
                "GCH_converted_calls": shard_gch_converted,
                "GCG_ambiguous_calls_excluded": ambiguous,
                "GCG_ambiguous_converted_calls_excluded": (
                    ambiguous_converted
                ),
                "status": "validated",
            }
        )
    if np.any(wcg_converted > wcg_covered) or np.any(
        gch_converted > gch_covered
    ):
        raise ValueError("converted calls exceed coverage")
    for context_name, sites, covered_calls, converted_calls in (
        ("WCG", wcg_sites, wcg_covered, wcg_converted),
        ("GCH", gch_sites, gch_covered, gch_converted),
    ):
        for metric, expected in (
            ("observed_sites", sites),
            ("covered_calls", covered_calls),
            ("converted_calls", converted_calls),
        ):
            separated = sum(
                (
                    species_context_stats[(species, context_name, metric)]
                    for species in ("human", "pig")
                ),
                np.zeros(cell_count, dtype=np.uint64),
            )
            if not np.array_equal(separated, expected):
                raise ValueError(
                    f"{context_name} species-separated {metric} do not "
                    "reconcile"
                )

    output_root = (
        args.output_root
        or project / "08_qc_summary/21_dna_single_cytosine_validation"
    )
    sample_root = output_root / sample
    per_cell_rows = []
    for offset, row in enumerate(cell_rows):
        wcg_cov = int(wcg_covered[offset])
        wcg_conv = int(wcg_converted[offset])
        gch_cov = int(gch_covered[offset])
        gch_conv = int(gch_converted[offset])
        per_cell_rows.append(
            {
                "sample_index": args.sample_index,
                "sample": sample,
                "group": row["group"],
                "cell_index": int(row["cell_index"]),
                "cell_id": row["cell_id"],
                "CB": row["CB"],
                "WCG_observed_sites": int(wcg_sites[offset]),
                "WCG_covered_calls": wcg_cov,
                "WCG_converted_calls": wcg_conv,
                "WCG_methylation_fraction": (
                    wcg_conv / wcg_cov if wcg_cov else "NA"
                ),
                "GCH_observed_sites": int(gch_sites[offset]),
                "GCH_covered_calls": gch_cov,
                "GCH_converted_calls": gch_conv,
                "GCH_accessibility_fraction": (
                    gch_conv / gch_cov if gch_cov else "NA"
                ),
                "legacy_CpG_covered_calls_without_end_exclusion": int(
                    legacy_covered[offset]
                ),
                "legacy_comparison": "reference_only_not_directly_equal",
            }
        )
        output_row = per_cell_rows[-1]
        for species in ("human", "pig"):
            for context_name, fraction_label in (
                ("WCG", "methylation_fraction"),
                ("GCH", "accessibility_fraction"),
            ):
                sites = int(
                    species_context_stats[
                        (species, context_name, "observed_sites")
                    ][offset]
                )
                covered_calls = int(
                    species_context_stats[
                        (species, context_name, "covered_calls")
                    ][offset]
                )
                converted_calls = int(
                    species_context_stats[
                        (species, context_name, "converted_calls")
                    ][offset]
                )
                prefix = f"{species}_{context_name}"
                output_row[f"{prefix}_observed_sites"] = sites
                output_row[f"{prefix}_covered_calls"] = covered_calls
                output_row[f"{prefix}_converted_calls"] = converted_calls
                output_row[f"{prefix}_{fraction_label}"] = (
                    converted_calls / covered_calls
                    if covered_calls
                    else "NA"
                )
        human_calls = sum(
            int(
                species_context_stats[
                    ("human", context_name, "covered_calls")
                ][offset]
            )
            for context_name in ("WCG", "GCH")
        )
        pig_calls = sum(
            int(
                species_context_stats[
                    ("pig", context_name, "covered_calls")
                ][offset]
            )
            for context_name in ("WCG", "GCH")
        )
        signal_calls = human_calls + pig_calls
        output_row["human_signal_covered_calls"] = human_calls
        output_row["pig_signal_covered_calls"] = pig_calls
        output_row["human_signal_call_fraction"] = (
            human_calls / signal_calls if signal_calls else "NA"
        )
        output_row["dominant_reference_signal"] = (
            "human"
            if human_calls > pig_calls
            else "pig"
            if pig_calls > human_calls
            else "tie"
            if signal_calls
            else "no_signal"
        )
    atomic_tsv(
        sample_root / f"{sample}.single_cytosine.per_cell.tsv.gz",
        list(per_cell_rows[0]),
        per_cell_rows,
        gzip_output=True,
    )
    atomic_tsv(
        sample_root / f"{sample}.single_cytosine.shard_manifest.tsv",
        list(shard_manifest[0]),
        shard_manifest,
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE},
        {"metric": "sample_index", "value": args.sample_index},
        {"metric": "sample", "value": sample},
        {"metric": "validated_shards", "value": EXPECTED_SHARDS},
        {"metric": "indexed_cells", "value": cell_count},
        {
            "metric": "cells_with_WCG_calls",
            "value": int(np.count_nonzero(wcg_covered)),
        },
        {
            "metric": "cells_with_GCH_calls",
            "value": int(np.count_nonzero(gch_covered)),
        },
        {
            "metric": "cells_with_either_signal",
            "value": int(np.count_nonzero(wcg_covered + gch_covered)),
        },
        {"metric": "observed_cell_cytosine_entries", "value": totals["entries"]},
        {"metric": "WCG_covered_calls", "value": totals["WCG_covered"]},
        {"metric": "WCG_converted_calls", "value": totals["WCG_converted"]},
        {
            "metric": "WCG_methylation_fraction",
            "value": totals["WCG_converted"] / totals["WCG_covered"]
            if totals["WCG_covered"]
            else "NA",
        },
        {"metric": "GCH_covered_calls", "value": totals["GCH_covered"]},
        {"metric": "GCH_converted_calls", "value": totals["GCH_converted"]},
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
            "value": int(legacy_covered.sum(dtype=np.uint64)),
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
    ]
    atomic_tsv(
        sample_root / f"{sample}.single_cytosine.validation_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] {sample}: WCG={totals['WCG_covered']} "
        f"GCH={totals['GCH_covered']} GCG_excluded={totals['GCG_ambiguous']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
