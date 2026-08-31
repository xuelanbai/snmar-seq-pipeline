#!/usr/bin/env python3
"""Step 94: audit reuse of exact CB sequences across 96 physical wells.

The outer well is an additional experimental partition and is not encoded in
the T7+T5+I7+I5 CB sequence.  This audit therefore reports cross-well barcode
reuse without collapsing or filtering any cell.  It evaluates four scopes:
all indexed DNA CBs, DNA CBs with WCG/GCH signal, raw-matrix RNA CBs, and
filtered RNA CBs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import statistics
import tempfile
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterator


EXPECTED_WELLS = 96
T7_COUNT = 96
T5_COUNT = 96
I7_COUNT = 384
I5_COUNT = 8
THEORETICAL_CB_SPACE = T7_COUNT * T5_COUNT * I7_COUNT * I5_COUNT
PIPELINE_SIGNATURE = "snmar_cb_overlap_96wells_20260804"


def opener_for(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with opener_for(path)(
        path, "rt", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        yield from reader


def read_metrics(path: Path) -> dict[str, str]:
    rows = list(iter_tsv(path))
    metrics = {row["metric"]: row["value"] for row in rows}
    if not rows or len(metrics) != len(rows):
        raise ValueError(f"empty or duplicate metrics: {path}")
    return metrics


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"identifier lacks trailing well number: {value}")
    return int(match.group(1))


def validate_cb(cb: str) -> None:
    if len(cb) != 36 or re.fullmatch(r"[ACGT]+", cb) is None:
        raise ValueError(f"CB is not an exact 36 bp T7+T5+I7+I5 sequence: {cb}")


def atomic_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        if path.suffix == ".gz":
            os.close(descriptor)
            context = gzip.open(
                temporary,
                "wt",
                encoding="utf-8",
                newline="",
                compresslevel=1,
            )
        else:
            context = os.fdopen(
                descriptor, "w", encoding="utf-8", newline=""
            )
        with context as handle:
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


def make_well_sets() -> dict[int, set[str]]:
    return {index: set() for index in range(1, EXPECTED_WELLS + 1)}


def add_cb(
    well_sets: dict[int, set[str]], well_index: int, cb: str, scope: str
) -> None:
    if not 1 <= well_index <= EXPECTED_WELLS:
        raise ValueError(f"{scope}: invalid well index {well_index}")
    validate_cb(cb)
    if cb in well_sets[well_index]:
        raise ValueError(f"{scope}: duplicate CB within well {well_index}: {cb}")
    well_sets[well_index].add(cb)


def audit_scope(
    scope: str,
    well_sets: dict[int, set[str]],
    well_metadata: dict[int, tuple[str, str]],
) -> tuple[
    dict[str, object],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    """Return summary, per-well, pairwise, reused-CB, and distribution rows."""
    if sorted(well_sets) != list(range(1, EXPECTED_WELLS + 1)):
        raise ValueError(f"{scope}: well set is incomplete")
    barcode_wells: defaultdict[str, list[int]] = defaultdict(list)
    for well_index, barcodes in well_sets.items():
        for cb in barcodes:
            barcode_wells[cb].append(well_index)
    reused = {
        cb: sorted(indices)
        for cb, indices in barcode_wells.items()
        if len(indices) > 1
    }
    reused_cb = set(reused)
    per_well_rows: list[dict[str, object]] = []
    for well_index in range(1, EXPECTED_WELLS + 1):
        barcodes = well_sets[well_index]
        shared = len(barcodes & reused_cb)
        well_id, well_set = well_metadata[well_index]
        per_well_rows.append(
            {
                "scope": scope,
                "well_index": well_index,
                "well_id": well_id,
                "well_set": well_set,
                "CB_count": len(barcodes),
                "CB_shared_with_any_other_well": shared,
                "CB_unique_to_this_well": len(barcodes) - shared,
                "shared_CB_fraction": shared / len(barcodes) if barcodes else "NA",
            }
        )

    pairwise_rows: list[dict[str, object]] = []
    intersections: list[int] = []
    for well_a, well_b in combinations(range(1, EXPECTED_WELLS + 1), 2):
        a = well_sets[well_a]
        b = well_sets[well_b]
        intersection = len(a & b)
        union = len(a | b)
        minimum = min(len(a), len(b))
        expected = len(a) * len(b) / THEORETICAL_CB_SPACE
        intersections.append(intersection)
        pairwise_rows.append(
            {
                "scope": scope,
                "well_a": well_a,
                "well_b": well_b,
                "well_set_a": well_metadata[well_a][1],
                "well_set_b": well_metadata[well_b][1],
                "same_well_set": int(
                    well_metadata[well_a][1] == well_metadata[well_b][1]
                ),
                "CB_a": len(a),
                "CB_b": len(b),
                "intersection": intersection,
                "union": union,
                "jaccard": intersection / union if union else "NA",
                "overlap_coefficient": (
                    intersection / minimum if minimum else "NA"
                ),
                "expected_intersection_uniform_full_CB_space": expected,
                "observed_over_uniform_expectation": (
                    intersection / expected if expected else "NA"
                ),
            }
        )

    reused_rows: list[dict[str, object]] = []
    for cb, well_indices in sorted(
        reused.items(), key=lambda item: (-len(item[1]), item[0])
    ):
        well_sets_seen = sorted(
            {well_metadata[index][1] for index in well_indices}
        )
        reused_rows.append(
            {
                "scope": scope,
                "CB": cb,
                "well_count": len(well_indices),
                "well_indices": ",".join(map(str, well_indices)),
                "well_sets": ",".join(well_sets_seen),
                "crosses_well_sets": int(len(well_sets_seen) > 1),
            }
        )

    total_cell_keys = sum(len(barcodes) for barcodes in well_sets.values())
    reused_cell_keys = sum(len(indices) for indices in reused.values())
    well_count_distribution = Counter(map(len, barcode_wells.values()))
    summary = {
        "scope": scope,
        "physical_wells": EXPECTED_WELLS,
        "well_CB_cell_keys": total_cell_keys,
        "unique_CB_sequences": len(barcode_wells),
        "CB_sequences_in_multiple_wells": len(reused),
        "well_CB_cell_keys_with_reused_CB": reused_cell_keys,
        "fraction_cell_keys_with_reused_CB": (
            reused_cell_keys / total_cell_keys if total_cell_keys else "NA"
        ),
        "maximum_wells_per_CB": max(well_count_distribution, default=0),
        "pairwise_comparisons": len(intersections),
        "well_pairs_with_nonzero_overlap": sum(value > 0 for value in intersections),
        "total_pairwise_CB_intersections": sum(intersections),
        "median_pairwise_intersection": (
            statistics.median(intersections) if intersections else "NA"
        ),
        "maximum_pairwise_intersection": max(intersections, default=0),
        "theoretical_T7_T5_I7_I5_CB_space": THEORETICAL_CB_SPACE,
        "automatic_filtering": "none_audit_only",
    }
    distribution_rows = [
        {
            "scope": scope,
            "number_of_wells_per_CB": number_of_wells,
            "CB_sequences": count,
            "well_CB_cell_keys": number_of_wells * count,
        }
        for number_of_wells, count in sorted(well_count_distribution.items())
    ]
    return (
        summary,
        per_well_rows,
        pairwise_rows,
        reused_rows,
        distribution_rows,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    dna_root = project / "08_qc_summary/23_dna_single_cytosine_complete"
    dna_summary = read_metrics(
        dna_root / "dna_single_cytosine_complete_summary.tsv"
    )
    dna_per_well = list(
        iter_tsv(dna_root / "dna_single_cytosine_complete_per_sample.tsv")
    )
    dna_cells_path = (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc.tsv.gz"
    )
    rna_raw_root = project / "08_qc_summary/16_rna_reverse_genefull_full"
    rna_raw_summary = read_metrics(
        rna_raw_root / "rna_reverse_genefull_full_summary.tsv"
    )
    rna_raw_cells_path = rna_raw_root / "rna_reverse_cell_qc.tsv.gz"
    rna_filtered_root = (
        project
        / "08_qc_summary/19_rna_filtered_matrix_full/minUMI10_minGene5"
    )
    rna_filtered_summary = read_metrics(
        rna_filtered_root / "rna_filtered_matrix_full_summary.tsv"
    )
    rna_filtered_cells_path = rna_filtered_root / "rna_filtered_cell_qc.tsv.gz"
    output_root = args.output_root or (
        project / "08_qc_summary/29_cb_overlap_96wells"
    )

    if (
        dna_summary.get("status") != "validated"
        or rna_raw_summary.get("status") != "validated"
        or rna_filtered_summary.get("status") != "validated"
    ):
        raise ValueError("upstream DNA/RNA summaries are not validated")
    if len(dna_per_well) != EXPECTED_WELLS:
        raise ValueError("DNA per-well summary does not contain 96 wells")
    well_metadata: dict[int, tuple[str, str]] = {}
    for row in dna_per_well:
        well_index = int(row["sample_index"])
        value = (row["sample"], row["group"])
        if well_index in well_metadata or not 1 <= well_index <= EXPECTED_WELLS:
            raise ValueError(f"invalid/duplicate DNA well: {row}")
        well_metadata[well_index] = value
    if sorted(well_metadata) != list(range(1, EXPECTED_WELLS + 1)):
        raise ValueError("DNA well metadata is incomplete")
    well_set_counts = Counter(value[1] for value in well_metadata.values())
    if sorted(well_set_counts.values()) != [48, 48]:
        raise ValueError(f"expected two 48-well sets: {well_set_counts}")

    scopes = {
        "DNA_indexed": make_well_sets(),
        "DNA_WCG_or_GCH_signal": make_well_sets(),
        "RNA_raw_matrix": make_well_sets(),
        "RNA_filtered_minUMI10_minGene5": make_well_sets(),
    }
    for row in iter_tsv(dna_cells_path):
        well_index = int(row["sample_index"])
        cb = row["CB"]
        add_cb(scopes["DNA_indexed"], well_index, cb, "DNA_indexed")
        signal_calls = int(row["WCG_covered_calls"]) + int(
            row["GCH_covered_calls"]
        )
        if signal_calls:
            add_cb(
                scopes["DNA_WCG_or_GCH_signal"],
                well_index,
                cb,
                "DNA_WCG_or_GCH_signal",
            )
    for row in iter_tsv(rna_raw_cells_path):
        if row["in_raw_matrix"] != "1":
            continue
        well_index = trailing_number(row["sample"])
        add_cb(scopes["RNA_raw_matrix"], well_index, row["CB"], "RNA_raw_matrix")
    for row in iter_tsv(rna_filtered_cells_path):
        well_index = trailing_number(row["sample"])
        add_cb(
            scopes["RNA_filtered_minUMI10_minGene5"],
            well_index,
            row["CB"],
            "RNA_filtered_minUMI10_minGene5",
        )

    if not all(
        scopes["RNA_filtered_minUMI10_minGene5"][index]
        <= scopes["RNA_raw_matrix"][index]
        for index in range(1, EXPECTED_WELLS + 1)
    ):
        raise ValueError("filtered RNA CBs are not a subset of raw RNA CBs")
    if not all(
        scopes["DNA_WCG_or_GCH_signal"][index]
        <= scopes["DNA_indexed"][index]
        for index in range(1, EXPECTED_WELLS + 1)
    ):
        raise ValueError("signal DNA CBs are not a subset of indexed DNA CBs")
    if sum(map(len, scopes["DNA_indexed"].values())) != int(
        dna_summary["indexed_cells"]
    ):
        raise ValueError("DNA indexed CB total differs from complete summary")
    if sum(map(len, scopes["DNA_WCG_or_GCH_signal"].values())) != int(
        dna_summary["cells_with_either_signal"]
    ):
        raise ValueError("DNA signal CB total differs from complete summary")
    if sum(map(len, scopes["RNA_raw_matrix"].values())) != int(
        rna_raw_summary["raw_matrix_cells"]
    ):
        raise ValueError("RNA raw CB total differs from matrix summary")
    if sum(map(len, scopes["RNA_filtered_minUMI10_minGene5"].values())) != int(
        rna_filtered_summary["filtered_cells"]
    ):
        raise ValueError("RNA filtered CB total differs from filter summary")

    summary_rows: list[dict[str, object]] = []
    per_well_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []
    reused_rows: list[dict[str, object]] = []
    distribution_rows: list[dict[str, object]] = []
    for scope, well_sets in scopes.items():
        summary, per_well, pairwise, reused, distribution = audit_scope(
            scope, well_sets, well_metadata
        )
        summary_rows.append(summary)
        per_well_rows.extend(per_well)
        pairwise_rows.extend(pairwise)
        reused_rows.extend(reused)
        distribution_rows.extend(distribution)

    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "cb_overlap_scope_summary.tsv"
    per_well_path = output_root / "cb_overlap_per_well.tsv"
    pairwise_path = output_root / "cb_overlap_pairwise_wells.tsv.gz"
    reused_path = output_root / "cb_reused_across_wells.tsv.gz"
    distribution_path = output_root / "cb_well_count_distribution.tsv"
    atomic_tsv(summary_path, list(summary_rows[0]), summary_rows)
    atomic_tsv(per_well_path, list(per_well_rows[0]), per_well_rows)
    atomic_tsv(pairwise_path, list(pairwise_rows[0]), pairwise_rows)
    atomic_tsv(
        reused_path,
        [
            "scope",
            "CB",
            "well_count",
            "well_indices",
            "well_sets",
            "crosses_well_sets",
        ],
        reused_rows,
    )
    atomic_tsv(
        distribution_path,
        [
            "scope",
            "number_of_wells_per_CB",
            "CB_sequences",
            "well_CB_cell_keys",
        ],
        distribution_rows,
    )
    run_summary_rows = [
        {"metric": "status", "value": "validated"},
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE},
        {"metric": "biological_sample_interpretation", "value": "not_inferred"},
        {
            "metric": "technical_repeat_status",
            "value": "not_resolved_by_CB_overlap_alone",
        },
        {"metric": "physical_wells", "value": EXPECTED_WELLS},
        {"metric": "cell_key", "value": "well_index|CB"},
        {"metric": "CB_sequence", "value": "T7+T5+I7+I5_36bp"},
        {"metric": "automatic_filtering", "value": "none_audit_only"},
        {
            "metric": "cross_well_reused_CB_interpretation",
            "value": "report_only_not_a_within_cell_doublet_definition",
        },
        {"metric": "scope_summary", "value": summary_path},
        {"metric": "per_well_table", "value": per_well_path},
        {"metric": "pairwise_table", "value": pairwise_path},
        {"metric": "reused_CB_table", "value": reused_path},
        {"metric": "well_count_distribution", "value": distribution_path},
    ]
    run_summary_path = output_root / "cb_overlap_audit_summary.tsv"
    atomic_tsv(run_summary_path, ["metric", "value"], run_summary_rows)
    print(f"[OK] 96-well CB overlap audit: {run_summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
