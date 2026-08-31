#!/usr/bin/env python3
"""Step 96: audit DNA/RNA cell-QC thresholds without filtering matrices.

The assay contains 96 physical wells with approximately 40 loaded nuclei per
well.  A whitelist-matched CB record is not automatically interpreted as an
intact cell.  This script evaluates a fixed, documented threshold grid against
the complete per-cell QC tables and reports retention globally, per well, and
across the paired DNA/RNA modalities.  It never edits upstream matrices or
chooses a final threshold automatically.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import statistics
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_cell_qc_threshold_audit_20260805"
DEFAULT_EXPECTED_WELLS = 96
DEFAULT_APPROX_NUCLEI_PER_WELL = 40


@dataclass(frozen=True)
class DnaRule:
    threshold_id: str
    min_total_calls: int
    min_wcg_sites: int = 0
    min_gch_sites: int = 0
    min_major_fraction: float = 0.0


@dataclass(frozen=True)
class RnaRule:
    threshold_id: str
    min_exact_umis: int
    min_genes: int


@dataclass(frozen=True)
class DnaCell:
    well_index: int
    cb: str
    wcg_sites: int
    gch_sites: int
    wcg_calls: int
    gch_calls: int
    human_calls: int
    pig_calls: int

    @property
    def total_calls(self) -> int:
        return self.wcg_calls + self.gch_calls

    @property
    def species_label(self) -> str:
        return dominant_label(self.human_calls, self.pig_calls)

    @property
    def major_fraction(self) -> float:
        return major_fraction(self.human_calls, self.pig_calls)


@dataclass(frozen=True)
class RnaCell:
    well_index: int
    cb: str
    exact_umis: int
    genes: int
    human_umis: int
    pig_umis: int

    @property
    def species_label(self) -> str:
        return dominant_label(self.human_umis, self.pig_umis)


DNA_RULES = (
    DnaRule("D_signal_any", 1),
    DnaRule("D_total_1k", 1_000),
    DnaRule("D_total_2k", 2_000),
    DnaRule("D_total_5k", 5_000),
    DnaRule("D_total_10k", 10_000),
    DnaRule("D_total_20k", 20_000),
    DnaRule("D_total_50k", 50_000),
    DnaRule("D_total_100k", 100_000),
    DnaRule("D_WCG_1k", 1, min_wcg_sites=1_000),
    DnaRule("D_GCH_1k", 1, min_gch_sites=1_000),
    DnaRule(
        "D_WCG_1k_GCH_1k", 1, min_wcg_sites=1_000, min_gch_sites=1_000
    ),
    DnaRule(
        "D_WCG_1k_GCH_10k", 1, min_wcg_sites=1_000, min_gch_sites=10_000
    ),
    DnaRule("D_total_5k_purity80", 5_000, min_major_fraction=0.80),
    DnaRule("D_total_5k_purity90", 5_000, min_major_fraction=0.90),
    DnaRule("D_total_10k_purity90", 10_000, min_major_fraction=0.90),
    DnaRule("D_total_10k_purity95", 10_000, min_major_fraction=0.95),
)

RNA_RULES = (
    RnaRule("R_UMI10_gene5", 10, 5),
    RnaRule("R_UMI20_gene10", 20, 10),
    RnaRule("R_UMI25_gene10", 25, 10),
    RnaRule("R_UMI30_gene10", 30, 10),
    RnaRule("R_UMI50_gene20", 50, 20),
    RnaRule("R_UMI100_gene50", 100, 50),
    RnaRule("R_UMI200_gene100", 200, 100),
    RnaRule("R_UMI500_gene200", 500, 200),
)


def opener_for(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with opener_for(path)(path, "rt", encoding="utf-8", newline="") as handle:
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


def atomic_tsv(
    path: Path, fields: list[str], rows: Iterable[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
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


def parse_nonnegative(row: dict[str, str], field: str) -> int:
    value = int(row[field])
    if value < 0:
        raise ValueError(f"negative {field}: {row}")
    return value


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"identifier lacks trailing well number: {value}")
    return int(match.group(1))


def dominant_label(first: int, second: int) -> str:
    if first == 0 and second == 0:
        return "no_signal"
    if first == second:
        return "tie"
    return "human" if first > second else "pig"


def major_fraction(first: int, second: int) -> float:
    total = first + second
    return max(first, second) / total if total else 0.0


def ratio(numerator: int | float, denominator: int | float) -> float | str:
    return numerator / denominator if denominator else "NA"


def median_or_na(values: list[int]) -> float | str:
    return statistics.median(values) if values else "NA"


def dna_passes(cell: DnaCell, rule: DnaRule) -> bool:
    return (
        cell.total_calls >= rule.min_total_calls
        and cell.wcg_sites >= rule.min_wcg_sites
        and cell.gch_sites >= rule.min_gch_sites
        and cell.major_fraction >= rule.min_major_fraction
    )


def rna_passes(cell: RnaCell, rule: RnaRule) -> bool:
    return cell.exact_umis >= rule.min_exact_umis and cell.genes >= rule.min_genes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--dna-cell-qc", type=Path)
    parser.add_argument("--dna-summary", type=Path)
    parser.add_argument("--rna-cell-qc", type=Path)
    parser.add_argument("--rna-summary", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--expected-wells", type=int, default=DEFAULT_EXPECTED_WELLS)
    parser.add_argument(
        "--approx-nuclei-per-well",
        type=int,
        default=DEFAULT_APPROX_NUCLEI_PER_WELL,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.expected_wells < 1 or args.approx_nuclei_per_well < 1:
        raise ValueError("expected wells and nuclei per well must be positive")
    project = args.project
    dna_path = args.dna_cell_qc or (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc.tsv.gz"
    )
    dna_summary_path = args.dna_summary or (
        project
        / "08_qc_summary/25_dna_single_cell_qc/"
        "dna_single_cell_coverage_qc_summary.tsv"
    )
    rna_path = args.rna_cell_qc or (
        project
        / "08_qc_summary/16_rna_reverse_genefull_full/"
        "rna_reverse_cell_qc.tsv.gz"
    )
    rna_summary_path = args.rna_summary or (
        project
        / "08_qc_summary/16_rna_reverse_genefull_full/"
        "rna_reverse_genefull_full_summary.tsv"
    )
    output_root = args.output_root or (
        project / "08_qc_summary/30_cell_qc_threshold_audit"
    )

    dna_summary = read_metrics(dna_summary_path)
    rna_summary = read_metrics(rna_summary_path)
    if dna_summary.get("status") != "validated":
        raise ValueError("DNA single-cell QC is not validated")
    if rna_summary.get("status") != "validated":
        raise ValueError("RNA raw matrix is not validated")

    well_metadata: dict[int, tuple[str, str]] = {}
    dna_cells: dict[tuple[int, str], DnaCell] = {}
    for row in iter_tsv(dna_path):
        well_index = int(row["sample_index"])
        if not 1 <= well_index <= args.expected_wells:
            raise ValueError(f"DNA well outside expected range: {well_index}")
        metadata = (row["sample"], row["group"])
        if well_index in well_metadata and well_metadata[well_index] != metadata:
            raise ValueError(f"inconsistent DNA well metadata: {well_index}")
        well_metadata[well_index] = metadata
        cell = DnaCell(
            well_index=well_index,
            cb=row["CB"],
            wcg_sites=parse_nonnegative(row, "WCG_observed_sites"),
            gch_sites=parse_nonnegative(row, "GCH_observed_sites"),
            wcg_calls=parse_nonnegative(row, "WCG_covered_calls"),
            gch_calls=parse_nonnegative(row, "GCH_covered_calls"),
            human_calls=parse_nonnegative(row, "human_signal_covered_calls"),
            pig_calls=parse_nonnegative(row, "pig_signal_covered_calls"),
        )
        if cell.human_calls + cell.pig_calls != cell.total_calls:
            raise ValueError(f"DNA species calls do not reconcile: {well_index}|{cell.cb}")
        key = (well_index, cell.cb)
        if not cell.cb or key in dna_cells:
            raise ValueError(f"duplicate/empty DNA cell key: {key}")
        dna_cells[key] = cell
    if sorted(well_metadata) != list(range(1, args.expected_wells + 1)):
        raise ValueError("DNA cell QC does not cover every expected physical well")
    if len(dna_cells) != int(dna_summary["cells"]):
        raise ValueError("DNA cell count differs from validated summary")

    rna_cells: dict[tuple[int, str], RnaCell] = {}
    rna_umi_sum = 0
    for row in iter_tsv(rna_path):
        if row["in_raw_matrix"] != "1":
            continue
        well_index = trailing_number(row["sample"])
        if not 1 <= well_index <= args.expected_wells:
            raise ValueError(f"RNA well outside expected range: {well_index}")
        cell = RnaCell(
            well_index=well_index,
            cb=row["CB"],
            exact_umis=parse_nonnegative(row, "exact_umis"),
            genes=parse_nonnegative(row, "genes_detected"),
            human_umis=parse_nonnegative(row, "human_exact_umis"),
            pig_umis=parse_nonnegative(row, "pig_exact_umis"),
        )
        if cell.human_umis + cell.pig_umis != cell.exact_umis:
            raise ValueError(f"RNA species UMIs do not reconcile: {well_index}|{cell.cb}")
        key = (well_index, cell.cb)
        if not cell.cb or key in rna_cells:
            raise ValueError(f"duplicate/empty RNA cell key: {key}")
        rna_cells[key] = cell
        rna_umi_sum += cell.exact_umis
    if len(rna_cells) != int(rna_summary["raw_matrix_cells"]):
        raise ValueError("RNA raw cell count differs from validated summary")
    if rna_umi_sum != int(rna_summary["raw_matrix_umi_sum"]):
        raise ValueError("RNA UMI sum differs from validated summary")

    expected_design_total = args.expected_wells * args.approx_nuclei_per_well
    dna_signal_cells = sum(cell.total_calls > 0 for cell in dna_cells.values())
    dna_total_signal_calls = sum(cell.total_calls for cell in dna_cells.values())
    dna_rule_sets: dict[str, set[tuple[int, str]]] = {}
    rna_rule_sets: dict[str, set[tuple[int, str]]] = {}
    dna_grid_rows: list[dict[str, object]] = []
    rna_grid_rows: list[dict[str, object]] = []
    per_well_rows: list[dict[str, object]] = []

    dna_input_by_well = Counter(cell.well_index for cell in dna_cells.values())
    rna_input_by_well = Counter(cell.well_index for cell in rna_cells.values())

    for rule in DNA_RULES:
        passed = {key for key, cell in dna_cells.items() if dna_passes(cell, rule)}
        dna_rule_sets[rule.threshold_id] = passed
        by_well = Counter(key[0] for key in passed)
        wcg_sites_sum = sum(dna_cells[key].wcg_sites for key in passed)
        gch_sites_sum = sum(dna_cells[key].gch_sites for key in passed)
        wcg_calls_sum = sum(dna_cells[key].wcg_calls for key in passed)
        gch_calls_sum = sum(dna_cells[key].gch_calls for key in passed)
        signal_sum = wcg_calls_sum + gch_calls_sum
        signal_by_well = Counter()
        for key in passed:
            signal_by_well[key[0]] += dna_cells[key].total_calls
        call_values = [dna_cells[key].total_calls for key in passed]
        labels = Counter(dna_cells[key].species_label for key in passed)
        well_counts = [by_well[index] for index in range(1, args.expected_wells + 1)]
        dna_grid_rows.append(
            {
                "threshold_id": rule.threshold_id,
                "min_total_signal_covered_calls": rule.min_total_calls,
                "min_WCG_observed_sites": rule.min_wcg_sites,
                "min_GCH_observed_sites": rule.min_gch_sites,
                "min_major_reference_signal_fraction": rule.min_major_fraction,
                "retained_cells": len(passed),
                "retained_fraction_of_indexed_CB": ratio(len(passed), len(dna_cells)),
                "retained_fraction_of_signal_CB": ratio(len(passed), dna_signal_cells),
                "retained_WCG_observed_sites": wcg_sites_sum,
                "retained_GCH_observed_sites": gch_sites_sum,
                "retained_WCG_covered_calls": wcg_calls_sum,
                "retained_GCH_covered_calls": gch_calls_sum,
                "retained_total_signal_covered_calls": signal_sum,
                "retained_fraction_of_total_signal_covered_calls": ratio(
                    signal_sum, dna_total_signal_calls
                ),
                "median_total_signal_covered_calls": median_or_na(call_values),
                "human_dominant_cells": labels["human"],
                "pig_dominant_cells": labels["pig"],
                "tie_cells": labels["tie"],
                "no_signal_cells": labels["no_signal"],
                "median_cells_per_well": statistics.median(well_counts),
                "minimum_cells_per_well": min(well_counts),
                "maximum_cells_per_well": max(well_counts),
                "wells_with_20_to_60_cells": sum(20 <= value <= 60 for value in well_counts),
                "absolute_difference_from_approx_design_total": abs(
                    len(passed) - expected_design_total
                ),
            }
        )
        for well_index in range(1, args.expected_wells + 1):
            well_id, well_set = well_metadata[well_index]
            retained = by_well[well_index]
            per_well_rows.append(
                {
                    "modality": "DNA",
                    "threshold_id": rule.threshold_id,
                    "well_index": well_index,
                    "well_id": well_id,
                    "well_set": well_set,
                    "input_CB": dna_input_by_well[well_index],
                    "retained_cells": retained,
                    "retained_fraction": ratio(retained, dna_input_by_well[well_index]),
                    "retained_signal_units": signal_by_well[well_index],
                    "signal_unit": "WCG_plus_GCH_covered_calls",
                }
            )

    for rule in RNA_RULES:
        passed = {key for key, cell in rna_cells.items() if rna_passes(cell, rule)}
        rna_rule_sets[rule.threshold_id] = passed
        by_well = Counter(key[0] for key in passed)
        umi_sum = sum(rna_cells[key].exact_umis for key in passed)
        umi_by_well = Counter()
        for key in passed:
            umi_by_well[key[0]] += rna_cells[key].exact_umis
        umi_values = [rna_cells[key].exact_umis for key in passed]
        gene_values = [rna_cells[key].genes for key in passed]
        labels = Counter(rna_cells[key].species_label for key in passed)
        well_counts = [by_well[index] for index in range(1, args.expected_wells + 1)]
        rna_grid_rows.append(
            {
                "threshold_id": rule.threshold_id,
                "min_exact_umis": rule.min_exact_umis,
                "min_genes": rule.min_genes,
                "retained_cells": len(passed),
                "retained_fraction_of_raw_matrix": ratio(len(passed), len(rna_cells)),
                "retained_exact_umis": umi_sum,
                "retained_fraction_of_raw_exact_umis": ratio(umi_sum, rna_umi_sum),
                "median_exact_umis": median_or_na(umi_values),
                "median_genes_detected": median_or_na(gene_values),
                "human_dominant_cells": labels["human"],
                "pig_dominant_cells": labels["pig"],
                "tie_cells": labels["tie"],
                "median_cells_per_well": statistics.median(well_counts),
                "minimum_cells_per_well": min(well_counts),
                "maximum_cells_per_well": max(well_counts),
                "wells_with_20_to_60_cells": sum(20 <= value <= 60 for value in well_counts),
                "absolute_difference_from_approx_design_total": abs(
                    len(passed) - expected_design_total
                ),
            }
        )
        for well_index in range(1, args.expected_wells + 1):
            well_id, well_set = well_metadata[well_index]
            retained = by_well[well_index]
            per_well_rows.append(
                {
                    "modality": "RNA",
                    "threshold_id": rule.threshold_id,
                    "well_index": well_index,
                    "well_id": well_id,
                    "well_set": well_set,
                    "input_CB": rna_input_by_well[well_index],
                    "retained_cells": retained,
                    "retained_fraction": ratio(retained, rna_input_by_well[well_index]),
                    "retained_signal_units": umi_by_well[well_index],
                    "signal_unit": "exact_UMIs",
                }
            )

    intersection_rows: list[dict[str, object]] = []
    intersection_per_well_rows: list[dict[str, object]] = []
    for dna_rule in DNA_RULES:
        dna_set = dna_rule_sets[dna_rule.threshold_id]
        for rna_rule in RNA_RULES:
            rna_set = rna_rule_sets[rna_rule.threshold_id]
            shared = dna_set & rna_set
            pairs = Counter(
                (dna_cells[key].species_label, rna_cells[key].species_label)
                for key in shared
            )
            comparable = sum(
                count
                for (dna_label, rna_label), count in pairs.items()
                if dna_label in {"human", "pig"} and rna_label in {"human", "pig"}
            )
            concordant = pairs[("human", "human")] + pairs[("pig", "pig")]
            intersection_rows.append(
                {
                    "DNA_threshold_id": dna_rule.threshold_id,
                    "RNA_threshold_id": rna_rule.threshold_id,
                    "DNA_retained_cells": len(dna_set),
                    "RNA_retained_cells": len(rna_set),
                    "shared_well_CB": len(shared),
                    "fraction_DNA_with_RNA_match": ratio(len(shared), len(dna_set)),
                    "fraction_RNA_with_DNA_match": ratio(len(shared), len(rna_set)),
                    "species_comparable_shared": comparable,
                    "species_concordant_shared": concordant,
                    "species_concordance_fraction": ratio(concordant, comparable),
                    "human_human": pairs[("human", "human")],
                    "pig_pig": pairs[("pig", "pig")],
                    "human_pig": pairs[("human", "pig")],
                    "pig_human": pairs[("pig", "human")],
                    "tie_or_no_signal_shared": len(shared) - comparable,
                }
            )
            dna_by_well = Counter(key[0] for key in dna_set)
            rna_by_well = Counter(key[0] for key in rna_set)
            shared_by_well = Counter(key[0] for key in shared)
            for well_index in range(1, args.expected_wells + 1):
                well_id, well_set = well_metadata[well_index]
                intersection_per_well_rows.append(
                    {
                        "DNA_threshold_id": dna_rule.threshold_id,
                        "RNA_threshold_id": rna_rule.threshold_id,
                        "well_index": well_index,
                        "well_id": well_id,
                        "well_set": well_set,
                        "DNA_retained_cells": dna_by_well[well_index],
                        "RNA_retained_cells": rna_by_well[well_index],
                        "shared_well_CB": shared_by_well[well_index],
                    }
                )

    output_root.mkdir(parents=True, exist_ok=True)
    dna_grid_path = output_root / "dna_threshold_grid.tsv"
    rna_grid_path = output_root / "rna_threshold_grid.tsv"
    per_well_path = output_root / "cell_qc_threshold_per_well.tsv"
    intersection_path = output_root / "dna_rna_threshold_intersection.tsv"
    intersection_per_well_path = (
        output_root / "dna_rna_threshold_intersection_per_well.tsv"
    )
    atomic_tsv(dna_grid_path, list(dna_grid_rows[0]), dna_grid_rows)
    atomic_tsv(rna_grid_path, list(rna_grid_rows[0]), rna_grid_rows)
    atomic_tsv(per_well_path, list(per_well_rows[0]), per_well_rows)
    atomic_tsv(intersection_path, list(intersection_rows[0]), intersection_rows)
    atomic_tsv(
        intersection_per_well_path,
        list(intersection_per_well_rows[0]),
        intersection_per_well_rows,
    )
    summary_rows = [
        {"metric": "status", "value": "validated"},
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE},
        {"metric": "mode", "value": "audit_only_no_matrix_filtering"},
        {"metric": "physical_wells", "value": args.expected_wells},
        {
            "metric": "approx_nuclei_per_well_design_annotation",
            "value": args.approx_nuclei_per_well,
        },
        {"metric": "approx_loaded_nuclei_total", "value": expected_design_total},
        {"metric": "cell_key", "value": "well_index|CB"},
        {"metric": "DNA_indexed_CB", "value": len(dna_cells)},
        {"metric": "DNA_CB_with_any_WCG_or_GCH_call", "value": dna_signal_cells},
        {"metric": "DNA_total_WCG_plus_GCH_covered_calls", "value": dna_total_signal_calls},
        {"metric": "RNA_raw_matrix_CB", "value": len(rna_cells)},
        {"metric": "RNA_raw_exact_UMIs", "value": rna_umi_sum},
        {"metric": "DNA_threshold_rules", "value": len(DNA_RULES)},
        {"metric": "RNA_threshold_rules", "value": len(RNA_RULES)},
        {"metric": "automatic_threshold_selection", "value": "none"},
        {"metric": "automatic_cell_filtering", "value": "none"},
        {"metric": "dna_threshold_grid", "value": dna_grid_path},
        {"metric": "rna_threshold_grid", "value": rna_grid_path},
        {"metric": "threshold_per_well", "value": per_well_path},
        {"metric": "DNA_RNA_intersection_grid", "value": intersection_path},
        {
            "metric": "DNA_RNA_intersection_per_well",
            "value": intersection_per_well_path,
        },
        {"metric": "DNA_cell_QC_input", "value": dna_path},
        {"metric": "RNA_cell_QC_input", "value": rna_path},
    ]
    summary_path = output_root / "cell_qc_threshold_audit_summary.tsv"
    atomic_tsv(summary_path, ["metric", "value"], summary_rows)
    print(
        "[OK] cell-QC threshold audit: "
        f"DNA={len(dna_cells)} RNA={len(rna_cells)} output={summary_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
