#!/usr/bin/env python3
"""Step 92: validate and summarize the 96 physical wells across DNA/RNA.

The upstream global summary files are already merged across all wells.  This
step does not merge matrices or recompute BigWigs.  It reconciles count totals
and emits one overall table, one well-set table, one physical-well table and one
species/context BigWig table.  Ratios are always recomputed from summed
numerators and denominators; per-well percentages are never averaged.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED_WELLS = 96
EXPECTED_GROUP_WELLS = 48
EXPECTED_TRACK_SHARDS = 8
PIPELINE_SIGNATURE = "snmar_96well_batch_species_overview_20260804"
CONTEXTS = ("WCG", "GCH")
SPECIES = ("human", "pig")


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


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"identifier lacks trailing well number: {value}")
    return int(match.group(1))


def ratio(numerator: int, denominator: int) -> float:
    if denominator < 1:
        raise ValueError("ratio denominator must be positive")
    if numerator < 0 or numerator > denominator:
        raise ValueError("ratio numerator is outside 0..denominator")
    return numerator / denominator


def ratio_or_na(numerator: int, denominator: int) -> float | str:
    if denominator == 0:
        if numerator != 0:
            raise ValueError("nonzero numerator with zero denominator")
        return "NA"
    return ratio(numerator, denominator)


def dominant_species(
    human_count: int, pig_count: int, zero_label: str
) -> str:
    """Classify one well|CB without ever collapsing CBs across wells."""
    if human_count < 0 or pig_count < 0:
        raise ValueError("species counts must be nonnegative")
    if human_count > pig_count:
        return "human"
    if pig_count > human_count:
        return "pig"
    return "tie" if human_count + pig_count else zero_label


def atomic_tsv(
    path: Path,
    fields: list[str],
    rows: list[dict[str, object]],
    allow_empty: bool = False,
) -> None:
    if not rows and not allow_empty:
        raise ValueError(f"refusing to write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        if path.suffix == ".gz":
            os.close(descriptor)
            handle_context = gzip.open(
                temporary,
                "wt",
                encoding="utf-8",
                newline="",
                compresslevel=1,
            )
        else:
            handle_context = os.fdopen(
                descriptor, "w", encoding="utf-8", newline=""
            )
        with handle_context as handle:
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


def require_integer_sum(
    rows: list[dict[str, str]], field: str, expected: int, label: str
) -> None:
    observed = sum(int(row[field]) for row in rows)
    if observed != expected:
        raise ValueError(f"{label} differs: {observed} != {expected}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--approx-nuclei-per-well", type=int, default=40)
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.approx_nuclei_per_well < 1:
        raise ValueError("approximate nuclei per well must be positive")
    project = args.project
    dna_root = project / "08_qc_summary/23_dna_single_cytosine_complete"
    dna_summary_path = dna_root / "dna_single_cytosine_complete_summary.tsv"
    dna_per_well_path = dna_root / "dna_single_cytosine_complete_per_sample.tsv"
    dna_cell_qc_path = (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc_summary.tsv"
    )
    dna_cell_table_path = (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc.tsv.gz"
    )
    bigwig_path = (
        project
        / "08_qc_summary/24_dna_single_cytosine_bigwig/"
        "dna_single_base_bigwig_summary.tsv"
    )
    track_root = project / "08_qc_summary/22_dna_single_cytosine_tracks"
    rna_raw_root = project / "08_qc_summary/16_rna_reverse_genefull_full"
    rna_raw_summary_path = rna_raw_root / "rna_reverse_genefull_full_summary.tsv"
    rna_per_file_path = rna_raw_root / "rna_reverse_genefull_full_per_file.tsv"
    rna_filter_root = (
        project
        / "08_qc_summary/19_rna_filtered_matrix_full/"
        "minUMI10_minGene5"
    )
    rna_filter_summary_path = (
        rna_filter_root / "rna_filtered_matrix_full_summary.tsv"
    )
    rna_filter_per_well_path = (
        rna_filter_root / "rna_filtered_matrix_full_per_sample.tsv"
    )
    rna_filtered_cell_path = (
        rna_filter_root / "rna_filtered_cell_qc.tsv.gz"
    )
    output_root = args.output_root or (
        project / "08_qc_summary/28_cross_well_overview"
    )

    dna = read_metrics(dna_summary_path)
    dna_wells = read_tsv(dna_per_well_path)
    dna_cell_qc = read_metrics(dna_cell_qc_path)
    dna_cells = read_tsv(dna_cell_table_path)
    bigwig = read_metrics(bigwig_path)
    rna_raw = read_metrics(rna_raw_summary_path)
    rna_files = read_tsv(rna_per_file_path)
    rna_filter = read_metrics(rna_filter_summary_path)
    rna_wells = read_tsv(rna_filter_per_well_path)
    rna_filtered_cells = read_tsv(rna_filtered_cell_path)
    if (
        dna.get("status") != "validated"
        or dna_cell_qc.get("status") != "validated"
        or bigwig.get("status") != "validated"
        or bigwig.get("aggregation")
        != "pseudobulk_all_CBs_by_reference_species"
        or rna_raw.get("status") != "validated"
        or rna_raw.get("cell_filtering") != "not_applied_raw_matrix"
        or rna_filter.get("status") != "validated"
        or rna_filter.get("minimum_exact_umis") != "10"
        or rna_filter.get("minimum_genes") != "5"
    ):
        raise ValueError("one or more upstream summaries are not validated")

    if len(dna_wells) != EXPECTED_WELLS:
        raise ValueError("DNA per-well table does not have 96 rows")
    dna_by_index: dict[int, dict[str, str]] = {}
    group_counts: Counter[str] = Counter()
    for row in dna_wells:
        index = int(row["sample_index"])
        if (
            index in dna_by_index
            or not 1 <= index <= EXPECTED_WELLS
            or trailing_number(row["sample"]) != index
            or row["status"] != "validated"
        ):
            raise ValueError(f"invalid DNA physical-well row: {row}")
        dna_by_index[index] = row
        group_counts[row["group"]] += 1
    if sorted(dna_by_index) != list(range(1, EXPECTED_WELLS + 1)):
        raise ValueError("DNA physical-well indices are incomplete")
    if (
        len(group_counts) != 2
        or sorted(group_counts.values())
        != [EXPECTED_GROUP_WELLS, EXPECTED_GROUP_WELLS]
    ):
        raise ValueError(
            f"expected two 48-well experimental sets: {group_counts}"
        )

    dna_reconciliation = {
        "indexed_cells": "indexed_cells",
        "cells_with_WCG_calls": "cells_with_WCG_calls",
        "cells_with_GCH_calls": "cells_with_GCH_calls",
        "cells_with_either_signal": "cells_with_either_signal",
        "observed_cell_cytosine_entries": "observed_cell_cytosine_entries",
        "WCG_covered_calls": "WCG_covered_calls",
        "WCG_converted_calls": "WCG_converted_calls",
        "GCH_covered_calls": "GCH_covered_calls",
        "GCH_converted_calls": "GCH_converted_calls",
        "GCG_ambiguous_template_calls_excluded": (
            "GCG_ambiguous_template_calls_excluded"
        ),
    }
    for field, metric in dna_reconciliation.items():
        require_integer_sum(
            dna_wells, field, int(dna[metric]), f"DNA {metric}"
        )
    if (
        int(dna_cell_qc["cells"]) != int(dna["indexed_cells"])
        or int(dna_cell_qc["cells_with_WCG_observed_sites"])
        != int(dna["cells_with_WCG_calls"])
        or int(dna_cell_qc["cells_with_GCH_observed_sites"])
        != int(dna["cells_with_GCH_calls"])
    ):
        raise ValueError("DNA complete and single-cell QC summaries differ")
    dominant_total = sum(
        int(dna_cell_qc[f"cells_dominant_reference_{category}"])
        for category in ("human", "pig", "tie", "no_signal")
    )
    if dominant_total != int(dna["indexed_cells"]):
        raise ValueError("DNA dominant-reference cell counts do not reconcile")

    track_summaries = []
    for shard in range(1, EXPECTED_TRACK_SHARDS + 1):
        path = track_root / f"dna_single_cytosine.track_shard_{shard:02d}.summary.tsv"
        metrics = read_metrics(path)
        if (
            metrics.get("status") != "validated"
            or int(metrics.get("shard", "0")) != shard
            or int(metrics.get("input_samples", "0")) != EXPECTED_WELLS
        ):
            raise ValueError(f"track shard is not a 96-well aggregate: {path}")
        track_summaries.append(path)

    species_rows: list[dict[str, object]] = []
    for species in SPECIES:
        for context in CONTEXTS:
            prefix = f"{species}_{context}"
            reference_sites = int(bigwig[f"{prefix}_reference_assayable_sites"])
            observed_sites = int(bigwig[f"{prefix}_observed_sites"])
            covered_entries = int(bigwig[f"{prefix}_covered_cell_site_entries"])
            covered_calls = int(bigwig[f"{prefix}_covered_calls"])
            converted_calls = int(bigwig[f"{prefix}_converted_calls"])
            observed_fraction = float(
                bigwig[f"{prefix}_observed_site_coverage_fraction"]
            )
            global_fraction = float(bigwig[f"{prefix}_global_fraction"])
            if (
                not 0 <= observed_sites <= reference_sites
                or covered_entries < observed_sites
                or not 0 <= converted_calls <= covered_calls
                or abs(observed_fraction - observed_sites / reference_sites)
                > 1e-12
                or abs(global_fraction - converted_calls / covered_calls)
                > 1e-12
            ):
                raise ValueError(f"invalid BigWig species/context totals: {prefix}")
            species_rows.append(
                {
                    "species": species,
                    "context": context,
                    "aggregation_scope": "all_96_wells_all_CBs",
                    "reference_assayable_sites": reference_sites,
                    "observed_sites": observed_sites,
                    "observed_site_coverage_fraction": observed_fraction,
                    "covered_cell_site_entries": covered_entries,
                    "covered_calls": covered_calls,
                    "converted_calls": converted_calls,
                    "global_fraction": global_fraction,
                    "fraction_bigwig": bigwig[f"{prefix}_fraction_bigwig"],
                    "coverage_bigwig": bigwig[f"{prefix}_coverage_bigwig"],
                    "covered_cells_bigwig": bigwig[
                        f"{prefix}_covered_cells_bigwig"
                    ],
                }
            )
    for context in CONTEXTS:
        context_rows = [row for row in species_rows if row["context"] == context]
        if (
            sum(int(row["covered_calls"]) for row in context_rows)
            != int(dna[f"{context}_covered_calls"])
            or sum(int(row["converted_calls"]) for row in context_rows)
            != int(dna[f"{context}_converted_calls"])
        ):
            raise ValueError(f"BigWig and DNA complete {context} totals differ")
    if (
        sum(int(row["covered_cell_site_entries"]) for row in species_rows)
        != int(dna["observed_cell_cytosine_entries"])
    ):
        raise ValueError("BigWig cell-site entries and DNA matrix rows differ")

    if (
        len(rna_files) != EXPECTED_WELLS
        or int(rna_raw.get("completed_files", "0")) != EXPECTED_WELLS
        or len(rna_wells) != EXPECTED_WELLS
        or int(rna_filter.get("samples", "0")) != EXPECTED_WELLS
    ):
        raise ValueError("RNA tables do not contain 96 physical wells")
    rna_files_by_sample = {row["sample"]: row for row in rna_files}
    rna_filter_by_sample = {row["sample"]: row for row in rna_wells}
    if (
        len(rna_files_by_sample) != EXPECTED_WELLS
        or len(rna_filter_by_sample) != EXPECTED_WELLS
        or set(rna_files_by_sample) != set(rna_filter_by_sample)
    ):
        raise ValueError("RNA raw and filtered physical-well identities differ")
    rna_by_index: dict[int, tuple[dict[str, str], dict[str, str]]] = {}
    for sample, raw_row in rna_files_by_sample.items():
        index = trailing_number(sample)
        filtered_row = rna_filter_by_sample[sample]
        if (
            index in rna_by_index
            or not 1 <= index <= EXPECTED_WELLS
            or raw_row["status"] != "validated"
            or filtered_row["status"] != "filtered"
        ):
            raise ValueError(f"invalid RNA physical-well row: {sample}")
        rna_by_index[index] = (raw_row, filtered_row)
    if sorted(rna_by_index) != list(range(1, EXPECTED_WELLS + 1)):
        raise ValueError("RNA physical-well indices are incomplete")

    raw_sum_fields = (
        "star_input_reads",
        "primary_mapped_reads",
        "unique_primary_reads",
        "tag_eligible_unique_primary",
        "gene_assigned_reads",
        "gene_ambiguous_reads",
        "gene_no_feature_reads",
        "exact_umis",
        "human_exact_umis",
        "pig_exact_umis",
    )
    for field in raw_sum_fields:
        require_integer_sum(
            rna_files, field, int(rna_raw[field]), f"RNA raw {field}"
        )
    filtered_reconciliation = {
        "raw_cells": "raw_cells",
        "filtered_cells": "filtered_cells",
        "raw_exact_umis": "raw_exact_umis",
        "filtered_exact_umis": "filtered_exact_umis",
    }
    for field, metric in filtered_reconciliation.items():
        require_integer_sum(
            rna_wells, field, int(rna_filter[metric]), f"RNA {metric}"
        )
    if (
        int(rna_filter["raw_cells"]) != int(rna_raw["raw_matrix_cells"])
        or int(rna_filter["raw_exact_umis"])
        != int(rna_raw["raw_matrix_umi_sum"])
    ):
        raise ValueError("RNA raw and filtered summaries differ")

    # The batch contains 96 physical wells.  Their biological-sample
    # relationship is metadata-dependent and is not inferred here.
    # Pool cells by dominant reference species without collapsing identical CB
    # sequences across wells; cell_id already includes the well identifier.
    if len(dna_cells) != int(dna["indexed_cells"]):
        raise ValueError("DNA per-cell table does not match indexed well|CB")
    dna_required = {
        "sample_index",
        "sample",
        "group",
        "cell_index",
        "cell_id",
        "CB",
        "dominant_reference_signal",
        "human_signal_covered_calls",
        "pig_signal_covered_calls",
        "WCG_covered_calls",
        "WCG_converted_calls",
        "GCH_covered_calls",
        "GCH_converted_calls",
    }
    for species in SPECIES:
        for context in CONTEXTS:
            dna_required.update(
                {
                    f"{species}_{context}_covered_calls",
                    f"{species}_{context}_converted_calls",
                }
            )
    if dna_required - set(dna_cells[0]):
        raise ValueError(
            "DNA per-cell table lacks species fields: "
            + ", ".join(sorted(dna_required - set(dna_cells[0])))
        )
    dna_species_totals: defaultdict[str, Counter[str]] = defaultdict(Counter)
    dna_label_rows: list[dict[str, object]] = []
    dna_seen_cells: set[str] = set()
    for row in dna_cells:
        label = row["dominant_reference_signal"]
        if label not in {"human", "pig", "tie", "no_signal"}:
            raise ValueError(f"invalid DNA species label: {label}")
        cell_id = row["cell_id"]
        well_index = int(row["sample_index"])
        if (
            cell_id in dna_seen_cells
            or not 1 <= well_index <= EXPECTED_WELLS
            or dna_by_index[well_index]["sample"] != row["sample"]
        ):
            raise ValueError(f"invalid/duplicate DNA cell: {cell_id}")
        dna_seen_cells.add(cell_id)
        human_calls = int(row["human_signal_covered_calls"])
        pig_calls = int(row["pig_signal_covered_calls"])
        total_calls = human_calls + pig_calls
        expected_label = dominant_species(human_calls, pig_calls, "no_signal")
        if label != expected_label:
            raise ValueError(f"DNA species label does not reconcile: {cell_id}")
        values = dna_species_totals[label]
        values["cells"] += 1
        values["human_signal_calls"] += human_calls
        values["pig_signal_calls"] += pig_calls
        for context in CONTEXTS:
            values[f"{context}_covered_all_references"] += int(
                row[f"{context}_covered_calls"]
            )
            values[f"{context}_converted_all_references"] += int(
                row[f"{context}_converted_calls"]
            )
            if label in SPECIES:
                values[f"matching_{context}_covered"] += int(
                    row[f"{label}_{context}_covered_calls"]
                )
                values[f"matching_{context}_converted"] += int(
                    row[f"{label}_{context}_converted_calls"]
                )
        matching_calls = (
            human_calls
            if label == "human"
            else pig_calls
            if label == "pig"
            else 0
        )
        other_calls = (
            pig_calls
            if label == "human"
            else human_calls
            if label == "pig"
            else total_calls
        )
        dna_label_rows.append(
            {
                "cell_index": int(row["cell_index"]),
                "cell_id": cell_id,
                "well_index": well_index,
                "well_id": row["sample"],
                "well_set": row["group"],
                "CB": row["CB"],
                "species_label": label,
                "human_signal_covered_calls": human_calls,
                "pig_signal_covered_calls": pig_calls,
                "total_signal_covered_calls": total_calls,
                "matching_species_calls": matching_calls,
                "other_species_calls": other_calls,
                "major_species_fraction": ratio_or_na(
                    max(human_calls, pig_calls), total_calls
                ),
                "label_method": "dominant_reference_covered_calls",
            }
        )
    for category in ("human", "pig", "tie", "no_signal"):
        if dna_species_totals[category]["cells"] != int(
            dna_cell_qc[f"cells_dominant_reference_{category}"]
        ):
            raise ValueError(f"DNA {category} cell count does not reconcile")

    dna_species_rows: list[dict[str, object]] = []
    for label in ("human", "pig", "tie", "no_signal"):
        values = dna_species_totals[label]
        matching_calls = (
            values["human_signal_calls"]
            if label == "human"
            else values["pig_signal_calls"]
            if label == "pig"
            else 0
        )
        other_calls = (
            values["pig_signal_calls"]
            if label == "human"
            else values["human_signal_calls"]
            if label == "pig"
            else values["human_signal_calls"] + values["pig_signal_calls"]
        )
        output: dict[str, object] = {
            "species_class": label,
            "cells": values["cells"],
            "human_reference_signal_calls": values["human_signal_calls"],
            "pig_reference_signal_calls": values["pig_signal_calls"],
            "matching_species_calls": matching_calls,
            "other_species_calls": other_calls,
            "matching_species_fraction": ratio_or_na(
                matching_calls, matching_calls + other_calls
            ),
            "cell_identity": "unique_well_id|CB",
            "pooling_scope": "one_batch_all_96_wells",
        }
        for context in CONTEXTS:
            all_covered = values[f"{context}_covered_all_references"]
            all_converted = values[f"{context}_converted_all_references"]
            matching_covered = values[f"matching_{context}_covered"]
            matching_converted = values[f"matching_{context}_converted"]
            output[f"{context}_covered_all_references"] = all_covered
            output[f"{context}_converted_all_references"] = all_converted
            output[f"{context}_fraction_all_references"] = ratio_or_na(
                all_converted, all_covered
            )
            output[f"matching_{context}_covered_calls"] = matching_covered
            output[f"matching_{context}_converted_calls"] = matching_converted
            output[f"matching_{context}_fraction"] = ratio_or_na(
                matching_converted, matching_covered
            )
        dna_species_rows.append(output)

    if len(rna_filtered_cells) != int(rna_filter["filtered_cells"]):
        raise ValueError("RNA filtered per-cell table does not reconcile")
    rna_required = {
        "sample",
        "cell_id",
        "CB",
        "exact_umis",
        "human_exact_umis",
        "pig_exact_umis",
    }
    if rna_required - set(rna_filtered_cells[0]):
        raise ValueError(
            "RNA filtered cell table lacks fields: "
            + ", ".join(sorted(rna_required - set(rna_filtered_cells[0])))
        )
    rna_species_totals: defaultdict[str, Counter[str]] = defaultdict(Counter)
    rna_label_rows: list[dict[str, object]] = []
    rna_seen_cells: set[str] = set()
    for row in rna_filtered_cells:
        cell_id = row["cell_id"]
        well_index = trailing_number(row["sample"])
        human_umis = int(row["human_exact_umis"])
        pig_umis = int(row["pig_exact_umis"])
        exact_umis = int(row["exact_umis"])
        if (
            cell_id in rna_seen_cells
            or human_umis + pig_umis != exact_umis
            or rna_by_index[well_index][0]["sample"] != row["sample"]
        ):
            raise ValueError(f"invalid/duplicate RNA cell: {cell_id}")
        rna_seen_cells.add(cell_id)
        label = dominant_species(human_umis, pig_umis, "tie")
        values = rna_species_totals[label]
        values["cells"] += 1
        values["human_exact_umis"] += human_umis
        values["pig_exact_umis"] += pig_umis
        values["exact_umis"] += exact_umis
        matching_umis = (
            human_umis
            if label == "human"
            else pig_umis
            if label == "pig"
            else 0
        )
        other_umis = exact_umis - matching_umis
        rna_label_rows.append(
            {
                "cell_id": cell_id,
                "well_index": well_index,
                "well_id": row["sample"],
                "CB": row["CB"],
                "species_label": label,
                "human_exact_umis": human_umis,
                "pig_exact_umis": pig_umis,
                "total_exact_umis": exact_umis,
                "matching_species_umis": matching_umis,
                "other_species_umis": other_umis,
                "major_species_fraction": ratio(
                    max(human_umis, pig_umis), exact_umis
                ),
                "label_method": "dominant_species_exact_UMIs",
            }
        )
    if sum(
        values["exact_umis"] for values in rna_species_totals.values()
    ) != int(rna_filter["filtered_exact_umis"]):
        raise ValueError("RNA species-pooled filtered UMI total differs")
    rna_species_rows: list[dict[str, object]] = []
    for label in ("human", "pig", "tie"):
        values = rna_species_totals[label]
        matching_umis = (
            values["human_exact_umis"]
            if label == "human"
            else values["pig_exact_umis"]
            if label == "pig"
            else 0
        )
        other_umis = values["exact_umis"] - matching_umis
        rna_species_rows.append(
            {
                "species_class": label,
                "filtered_cells": values["cells"],
                "human_exact_umis": values["human_exact_umis"],
                "pig_exact_umis": values["pig_exact_umis"],
                "total_exact_umis": values["exact_umis"],
                "matching_species_umis": matching_umis,
                "other_species_umis": other_umis,
                "matching_species_fraction": ratio_or_na(
                    matching_umis, values["exact_umis"]
                ),
                "cell_identity": "unique_well_id|CB",
                "pooling_scope": "one_batch_all_96_wells",
            }
        )

    # DNA and RNA use different filename prefixes, so the cross-modal cell key
    # is normalized to well_index|CB.  Species is annotation, never part of the
    # key and never used to pre-split the DNA UMAP input.
    dna_by_cross_key: dict[tuple[int, str], dict[str, object]] = {}
    for row in dna_label_rows:
        key = (int(row["well_index"]), str(row["CB"]))
        if key in dna_by_cross_key:
            raise ValueError(f"duplicate DNA well_index|CB: {key}")
        dna_by_cross_key[key] = row
    rna_by_cross_key: dict[tuple[int, str], dict[str, object]] = {}
    for row in rna_label_rows:
        key = (int(row["well_index"]), str(row["CB"]))
        if key in rna_by_cross_key:
            raise ValueError(f"duplicate RNA well_index|CB: {key}")
        rna_by_cross_key[key] = row
    shared_cross_keys = sorted(set(dna_by_cross_key) & set(rna_by_cross_key))
    cross_modal_rows: list[dict[str, object]] = []
    comparable_cells = 0
    concordant_cells = 0
    for well_index, cb in shared_cross_keys:
        dna_row = dna_by_cross_key[(well_index, cb)]
        rna_row = rna_by_cross_key[(well_index, cb)]
        dna_label = str(dna_row["species_label"])
        rna_label = str(rna_row["species_label"])
        comparable = dna_label in SPECIES and rna_label in SPECIES
        concordant: int | str = "NA"
        if comparable:
            comparable_cells += 1
            concordant = int(dna_label == rna_label)
            concordant_cells += int(concordant)
        cross_modal_rows.append(
            {
                "cross_modal_cell_id": f"{well_index}|{cb}",
                "well_index": well_index,
                "CB": cb,
                "dna_cell_id": dna_row["cell_id"],
                "rna_cell_id": rna_row["cell_id"],
                "dna_species_label": dna_label,
                "rna_species_label": rna_label,
                "species_comparable": int(comparable),
                "species_concordant": concordant,
                "dna_human_signal_covered_calls": dna_row[
                    "human_signal_covered_calls"
                ],
                "dna_pig_signal_covered_calls": dna_row[
                    "pig_signal_covered_calls"
                ],
                "rna_human_exact_umis": rna_row["human_exact_umis"],
                "rna_pig_exact_umis": rna_row["pig_exact_umis"],
                "join_key_definition": "well_index|CB",
                "species_role": "post_hoc_UMAP_validation_only",
            }
        )

    per_well_rows: list[dict[str, object]] = []
    for index in range(1, EXPECTED_WELLS + 1):
        dna_row = dna_by_index[index]
        rna_raw_row, rna_filter_row = rna_by_index[index]
        dna_indexed = int(dna_row["indexed_cells"])
        dna_signal = int(dna_row["cells_with_either_signal"])
        wcg_covered = int(dna_row["WCG_covered_calls"])
        wcg_converted = int(dna_row["WCG_converted_calls"])
        gch_covered = int(dna_row["GCH_covered_calls"])
        gch_converted = int(dna_row["GCH_converted_calls"])
        rna_raw_cells = int(rna_filter_row["raw_cells"])
        rna_filtered_cells = int(rna_filter_row["filtered_cells"])
        rna_raw_umis = int(rna_filter_row["raw_exact_umis"])
        rna_filtered_umis = int(rna_filter_row["filtered_exact_umis"])
        if (
            int(rna_raw_row["exact_umis"]) != rna_raw_umis
            or dna_signal > dna_indexed
        ):
            raise ValueError(f"physical well {index} totals differ")
        per_well_rows.append(
            {
                "well_index": index,
                "experimental_unit": "physical_well",
                "well_set": dna_row["group"],
                "dna_well_id": dna_row["sample"],
                "rna_well_id": rna_raw_row["sample"],
                "design_approx_loaded_nuclei": args.approx_nuclei_per_well,
                "dna_indexed_well_CB": dna_indexed,
                "dna_CB_with_WCG": int(dna_row["cells_with_WCG_calls"]),
                "dna_CB_with_GCH": int(dna_row["cells_with_GCH_calls"]),
                "dna_CB_with_either_signal": dna_signal,
                "dna_signal_CB_fraction": ratio(dna_signal, dna_indexed),
                "dna_observed_cell_cytosine_entries": int(
                    dna_row["observed_cell_cytosine_entries"]
                ),
                "dna_WCG_covered_calls": wcg_covered,
                "dna_WCG_converted_calls": wcg_converted,
                "dna_WCG_weighted_fraction": ratio(
                    wcg_converted, wcg_covered
                ),
                "dna_GCH_covered_calls": gch_covered,
                "dna_GCH_converted_calls": gch_converted,
                "dna_GCH_weighted_fraction": ratio(
                    gch_converted, gch_covered
                ),
                "dna_GCG_ambiguous_calls_excluded": int(
                    dna_row["GCG_ambiguous_template_calls_excluded"]
                ),
                "rna_STAR_input_reads": int(rna_raw_row["star_input_reads"]),
                "rna_unique_primary_reads": int(
                    rna_raw_row["unique_primary_reads"]
                ),
                "rna_gene_assigned_reads": int(
                    rna_raw_row["gene_assigned_reads"]
                ),
                "rna_gene_assignment_fraction": float(
                    rna_raw_row["gene_assignment_fraction"]
                ),
                "rna_raw_cells": rna_raw_cells,
                "rna_filtered_cells": rna_filtered_cells,
                "rna_cell_retention_fraction": ratio(
                    rna_filtered_cells, rna_raw_cells
                ),
                "rna_raw_exact_umis": rna_raw_umis,
                "rna_filtered_exact_umis": rna_filtered_umis,
                "rna_umi_retention_fraction": ratio(
                    rna_filtered_umis, rna_raw_umis
                ),
                "rna_human_exact_umis": int(
                    rna_raw_row["human_exact_umis"]
                ),
                "rna_pig_exact_umis": int(rna_raw_row["pig_exact_umis"]),
            }
        )

    grouped: defaultdict[str, list[dict[str, object]]] = defaultdict(list)
    for row in per_well_rows:
        grouped[str(row["well_set"])].append(row)
    per_well_set_rows: list[dict[str, object]] = []
    for well_set in sorted(grouped):
        rows = grouped[well_set]

        def total(field: str) -> int:
            return sum(int(row[field]) for row in rows)

        dna_indexed = total("dna_indexed_well_CB")
        dna_signal = total("dna_CB_with_either_signal")
        wcg_covered = total("dna_WCG_covered_calls")
        wcg_converted = total("dna_WCG_converted_calls")
        gch_covered = total("dna_GCH_covered_calls")
        gch_converted = total("dna_GCH_converted_calls")
        rna_raw_cells = total("rna_raw_cells")
        rna_filtered_cells = total("rna_filtered_cells")
        rna_raw_umis = total("rna_raw_exact_umis")
        rna_filtered_umis = total("rna_filtered_exact_umis")
        per_well_set_rows.append(
            {
                "well_set": well_set,
                "physical_wells": len(rows),
                "design_approx_loaded_nuclei": (
                    len(rows) * args.approx_nuclei_per_well
                ),
                "dna_indexed_well_CB": dna_indexed,
                "dna_CB_with_either_signal": dna_signal,
                "dna_signal_CB_fraction": ratio(dna_signal, dna_indexed),
                "dna_WCG_covered_calls": wcg_covered,
                "dna_WCG_converted_calls": wcg_converted,
                "dna_WCG_weighted_fraction": ratio(
                    wcg_converted, wcg_covered
                ),
                "dna_GCH_covered_calls": gch_covered,
                "dna_GCH_converted_calls": gch_converted,
                "dna_GCH_weighted_fraction": ratio(
                    gch_converted, gch_covered
                ),
                "rna_raw_cells": rna_raw_cells,
                "rna_filtered_cells": rna_filtered_cells,
                "rna_cell_retention_fraction": ratio(
                    rna_filtered_cells, rna_raw_cells
                ),
                "rna_raw_exact_umis": rna_raw_umis,
                "rna_filtered_exact_umis": rna_filtered_umis,
                "rna_umi_retention_fraction": ratio(
                    rna_filtered_umis, rna_raw_umis
                ),
                "rna_human_exact_umis": total("rna_human_exact_umis"),
                "rna_pig_exact_umis": total("rna_pig_exact_umis"),
                "ratio_method": "sum_numerator_divided_by_sum_denominator",
            }
        )

    per_well_path = output_root / "snmar_96well_batch_per_well.tsv"
    per_well_set_path = (
        output_root / "snmar_96well_batch_per_well_set.tsv"
    )
    species_path = (
        output_root / "snmar_96well_batch_bigwig_by_species_context.tsv"
    )
    dna_species_path = (
        output_root / "snmar_96well_batch_DNA_cells_by_species.tsv"
    )
    rna_species_path = (
        output_root / "snmar_96well_batch_RNA_filtered_cells_by_species.tsv"
    )
    dna_labels_path = (
        output_root / "snmar_96well_batch_DNA_cell_species_labels.tsv.gz"
    )
    rna_labels_path = (
        output_root / "snmar_96well_batch_RNA_filtered_cell_species_labels.tsv.gz"
    )
    cross_modal_path = (
        output_root / "snmar_96well_batch_DNA_RNA_CB_match.tsv.gz"
    )
    overall_path = output_root / "snmar_96well_batch_overall_summary.tsv"
    atomic_tsv(per_well_path, list(per_well_rows[0]), per_well_rows)
    atomic_tsv(
        per_well_set_path,
        list(per_well_set_rows[0]),
        per_well_set_rows,
    )
    atomic_tsv(species_path, list(species_rows[0]), species_rows)
    atomic_tsv(dna_species_path, list(dna_species_rows[0]), dna_species_rows)
    atomic_tsv(rna_species_path, list(rna_species_rows[0]), rna_species_rows)
    atomic_tsv(dna_labels_path, list(dna_label_rows[0]), dna_label_rows)
    atomic_tsv(rna_labels_path, list(rna_label_rows[0]), rna_label_rows)
    cross_modal_fields = [
        "cross_modal_cell_id",
        "well_index",
        "CB",
        "dna_cell_id",
        "rna_cell_id",
        "dna_species_label",
        "rna_species_label",
        "species_comparable",
        "species_concordant",
        "dna_human_signal_covered_calls",
        "dna_pig_signal_covered_calls",
        "rna_human_exact_umis",
        "rna_pig_exact_umis",
        "join_key_definition",
        "species_role",
    ]
    atomic_tsv(
        cross_modal_path,
        cross_modal_fields,
        cross_modal_rows,
        allow_empty=True,
    )

    overall_rows: list[dict[str, object]] = []

    def add(metric: str, value: object, definition: str) -> None:
        overall_rows.append(
            {"metric": metric, "value": value, "definition": definition}
        )

    add("status", "validated", "all input counts and ratios reconciled")
    add("pipeline_signature", PIPELINE_SIGNATURE, "summary code version")
    add(
        "biological_samples",
        "not_inferred",
        "requires experimental metadata; not inferred from file or CB structure",
    )
    add(
        "physical_wells",
        EXPECTED_WELLS,
        "96 physical wells in this batch; biological relationship not inferred",
    )
    add(
        "legacy_upstream_sample_field_semantics",
        "physical_well_id",
        "upstream field name sample is retained for path compatibility only",
    )
    add(
        "design_approx_nuclei_per_well",
        args.approx_nuclei_per_well,
        "experimental-design annotation; not inferred from CB counts",
    )
    add(
        "design_approx_loaded_nuclei_total",
        EXPECTED_WELLS * args.approx_nuclei_per_well,
        "approximate design total; not an observed cell count",
    )
    add(
        "dna_indexed_well_CB",
        int(dna["indexed_cells"]),
        "well-level CB records across 96 wells; exact cross-well CB overlap is audited by step 94",
    )
    add(
        "dna_CB_with_either_WCG_or_GCH",
        int(dna["cells_with_either_signal"]),
        "well_id|CB with at least one WCG or GCH observation",
    )
    add(
        "dna_signal_CB_fraction",
        ratio(
            int(dna["cells_with_either_signal"]), int(dna["indexed_cells"])
        ),
        "signal CB divided by indexed well_id|CB",
    )
    for context, fraction_name in (
        ("WCG", "endogenous_methylation"),
        ("GCH", "GpC_accessibility"),
    ):
        covered = int(dna[f"{context}_covered_calls"])
        converted = int(dna[f"{context}_converted_calls"])
        add(
            f"dna_{context}_covered_calls",
            covered,
            "PCR-deduplicated template calls across all 96 wells",
        )
        add(
            f"dna_{context}_converted_calls",
            converted,
            "converted template calls across all 96 wells",
        )
        add(
            f"dna_{context}_weighted_fraction",
            ratio(converted, covered),
            f"{fraction_name}: total converted / total covered",
        )
    add(
        "dna_GCG_ambiguous_calls_excluded",
        int(dna["GCG_ambiguous_template_calls_excluded"]),
        "excluded from both WCG and GCH",
    )
    for metric in (
        "WCG_methylation_fraction_median",
        "GCH_accessibility_fraction_median",
        "cells_WCG_observed_sites_ge_1000",
        "cells_GCH_observed_sites_ge_1000",
        "cells_WCG_observed_sites_ge_10000",
        "cells_GCH_observed_sites_ge_10000",
        "cells_dominant_reference_human",
        "cells_dominant_reference_pig",
        "cells_dominant_reference_tie",
    ):
        add(
            f"dna_cell_qc_{metric}",
            dna_cell_qc[metric],
            "global single-cell QC across 96 wells",
        )
    for metric in (
        "star_input_reads",
        "primary_mapped_reads",
        "unique_primary_reads",
        "gene_assigned_reads",
        "gene_assignment_fraction",
        "exact_umis",
        "human_exact_umis",
        "pig_exact_umis",
        "raw_matrix_cells",
    ):
        add(
            f"rna_raw_{metric}",
            rna_raw[metric],
            "raw RNA aggregate across 96 wells",
        )
    for metric in (
        "raw_cells",
        "filtered_cells",
        "cell_fraction",
        "raw_exact_umis",
        "filtered_exact_umis",
        "umi_fraction",
        "features",
        "filtered_nonzero_entries",
    ):
        add(
            f"rna_filtered_{metric}",
            rna_filter[metric],
            "RNA minUMI10/minGene5 aggregate across 96 wells",
        )
    add("per_well_table", per_well_path, "96 physical-well rows")
    add(
        "per_well_set_table",
        per_well_set_path,
        "two 48-well experimental sets; not biological samples",
    )
    add(
        "bigwig_species_context_table",
        species_path,
        "human/pig WCG/GCH global BigWig counts and fractions",
    )
    add(
        "bigwig_aggregation_scope",
        "one_batch_all_96_wells_all_CBs_then_split_by_reference_species",
        "validated through eight track shards with input_samples=96",
    )
    add(
        "DNA_species_pooled_table",
        dna_species_path,
        "human/pig/tie/no-signal well|CB pools across all 96 wells",
    )
    add(
        "RNA_species_pooled_table",
        rna_species_path,
        "human/pig/tie filtered RNA cell pools across all 96 wells",
    )
    add(
        "DNA_species_labels",
        dna_labels_path,
        "cell-level dominant-reference labels for UMAP validation",
    )
    add(
        "RNA_species_labels",
        rna_labels_path,
        "cell-level dominant-species UMI labels for UMAP validation",
    )
    add(
        "DNA_100kb_UMAP_input_policy",
        "all_DNA_cells_together_no_species_pre_split",
        "species is not used to construct or partition the UMAP input",
    )
    add(
        "cross_modal_cell_key",
        "well_index|CB",
        "DNA and RNA filename prefixes differ; normalized well index and CB are joined",
    )
    add(
        "DNA_RNA_filtered_CB_matches",
        len(cross_modal_rows),
        "intersection of DNA indexed cells and filtered RNA cells",
    )
    add(
        "DNA_RNA_species_comparable_matches",
        comparable_cells,
        "matched cells with non-tie human/pig labels in both modalities",
    )
    add(
        "DNA_RNA_species_concordant_matches",
        concordant_cells,
        "comparable matched cells with the same dominant species",
    )
    add(
        "DNA_RNA_species_concordance_fraction",
        ratio_or_na(concordant_cells, comparable_cells),
        "post-hoc validation only; not an input to DNA UMAP",
    )
    add(
        "DNA_RNA_CB_match_table",
        cross_modal_path,
        "one row per shared well_index|CB; header-only if no matches; no cross-well CB collapse",
    )
    add(
        "ratio_policy",
        "sum_numerator_divided_by_sum_denominator",
        "never arithmetic mean of per-well percentages",
    )
    add(
        "DNA_RNA_cell_axis",
        "joinable_by_normalized_well_index_and_CB",
        "RNA validates DNA UMAP labels after embedding; it does not define the embedding",
    )
    add("dna_summary_input", dna_summary_path, "validated upstream input")
    add("dna_per_well_input", dna_per_well_path, "96 DNA wells")
    add("dna_cell_qc_input", dna_cell_qc_path, "global DNA cell QC")
    add(
        "dna_cell_table_input",
        dna_cell_table_path,
        "DNA well|CB species calls across all wells",
    )
    add("bigwig_summary_input", bigwig_path, "global species BigWigs")
    add("rna_raw_summary_input", rna_raw_summary_path, "global raw RNA")
    add("rna_per_file_input", rna_per_file_path, "96 RNA wells")
    add(
        "rna_filter_summary_input",
        rna_filter_summary_path,
        "global filtered RNA",
    )
    add(
        "rna_filter_per_well_input",
        rna_filter_per_well_path,
        "96 filtered RNA wells",
    )
    add(
        "rna_filtered_cell_input",
        rna_filtered_cell_path,
        "filtered RNA well|CB species UMIs",
    )
    atomic_tsv(
        overall_path, ["metric", "value", "definition"], overall_rows
    )
    print(
        "[OK] 96-well batch species-aware DNA/RNA overview: "
        f"DNA_CB={dna['indexed_cells']} "
        f"DNA_signal_CB={dna['cells_with_either_signal']} "
        f"RNA_filtered_cells={rna_filter['filtered_cells']} "
        f"output={overall_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
