#!/usr/bin/env python3
"""Step 101: summarize RNA species and residual-DNA diagnostics across wells.

The output is diagnostic only.  A validated status means that counts reconcile;
it does not mean that the biological species labels have passed QC.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_rna_species_residual_audit_20260806"
ALIGNMENT_SIGNATURE = "snmar_rna_species_alignment_audit_20260806"
CELL_RULES = (
    ("R_UMI20_gene10", 20, 10),
    ("R_UMI25_gene10", 25, 10),
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


def read_tsv(path: Path) -> list[dict[str, str]]:
    return list(iter_tsv(path))


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"sample lacks trailing well number: {value}")
    return int(match.group(1))


def as_nonnegative(row: dict[str, str], field: str) -> int:
    value = int(row[field])
    if value < 0:
        raise ValueError(f"negative {field}: {row}")
    return value


def ratio(numerator: int | float, denominator: int | float) -> float | str:
    return numerator / denominator if denominator else "NA"


def dominant_label(human: int, pig: int) -> str:
    if human == pig:
        return "tie"
    return "human" if human > pig else "pig"


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


def pivot_demux(path: Path) -> dict[str, dict[str, str]]:
    by_pair: dict[str, dict[str, str]] = defaultdict(dict)
    samples: dict[str, str] = {}
    for row in iter_tsv(path):
        pair_id = row["pair_id"]
        sample = row["sample"]
        if pair_id in samples and samples[pair_id] != sample:
            raise ValueError(f"inconsistent demultiplex sample: {pair_id}")
        samples[pair_id] = sample
        metric = row["metric"]
        if metric in by_pair[pair_id]:
            raise ValueError(f"duplicate demultiplex metric: {pair_id} {metric}")
        by_pair[pair_id][metric] = row["count"]
    required = {
        "total_pairs",
        "strict_assigned",
        "DNA_residual_GTCA_linker_given_T7",
    }
    for pair_id, metrics in by_pair.items():
        missing = sorted(required - metrics.keys())
        if missing:
            raise ValueError(f"{pair_id}: missing demultiplex metrics: {missing}")
    return dict(by_pair)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--alignment-audit-root", type=Path)
    parser.add_argument("--rna-per-file", type=Path)
    parser.add_argument("--rna-cell-qc", type=Path)
    parser.add_argument("--demux-per-pair", type=Path)
    parser.add_argument("--genefull-by-file-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--expected-wells", type=int, default=96)
    parser.add_argument("--human-only-wells", type=int, default=48)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.expected_wells < 2
        or args.human_only_wells < 1
        or args.human_only_wells >= args.expected_wells
    ):
        raise ValueError("invalid well design")
    project = args.project
    manifest_path = args.manifest or (
        project
        / "08_qc_summary/14_postmapping_full_inputs/"
        "rna_reverse_genefull_full_manifest.tsv"
    )
    alignment_root = args.alignment_audit_root or (
        project / "08_qc_summary/31_rna_species_residual_audit/by_file"
    )
    per_file_path = args.rna_per_file or (
        project
        / "08_qc_summary/16_rna_reverse_genefull_full/"
        "rna_reverse_genefull_full_per_file.tsv"
    )
    cell_qc_path = args.rna_cell_qc or (
        project
        / "08_qc_summary/16_rna_reverse_genefull_full/"
        "rna_reverse_cell_qc.tsv.gz"
    )
    demux_path = args.demux_per_pair or (
        project / "08_qc_summary/02_rna_demultiplex/full_strict_per_pair.tsv"
    )
    genefull_root = args.genefull_by_file_root or (
        project / "07_postmapping/rna_reverse_genefull_full/by_file"
    )
    output_root = args.output_root or (
        project / "08_qc_summary/31_rna_species_residual_audit"
    )

    manifest_rows = read_tsv(manifest_path)
    if len(manifest_rows) != args.expected_wells:
        raise ValueError(
            f"expected {args.expected_wells} manifest rows, found {len(manifest_rows)}"
        )
    per_file_rows = read_tsv(per_file_path)
    if len(per_file_rows) != args.expected_wells:
        raise ValueError(
            f"expected {args.expected_wells} RNA per-file rows, found {len(per_file_rows)}"
        )
    per_file = {row["pair_id"]: row for row in per_file_rows}
    if len(per_file) != len(per_file_rows):
        raise ValueError("duplicate RNA per-file pair_id")
    demux = pivot_demux(demux_path)
    if len(demux) != args.expected_wells:
        raise ValueError(
            f"expected {args.expected_wells} demultiplex pairs, found {len(demux)}"
        )

    cell_counts: Counter[tuple[str, int, str]] = Counter()
    cell_umis: Counter[tuple[str, int, str]] = Counter()
    raw_cells_by_well: Counter[int] = Counter()
    raw_umis_by_well: Counter[int] = Counter()
    for row in iter_tsv(cell_qc_path):
        if row["in_raw_matrix"] != "1":
            continue
        well_index = trailing_number(row["sample"])
        if not 1 <= well_index <= args.expected_wells:
            raise ValueError(f"RNA cell outside expected wells: {well_index}")
        exact_umis = as_nonnegative(row, "exact_umis")
        genes = as_nonnegative(row, "genes_detected")
        human_umis = as_nonnegative(row, "human_exact_umis")
        pig_umis = as_nonnegative(row, "pig_exact_umis")
        if human_umis + pig_umis != exact_umis:
            raise ValueError(f"RNA species UMI counts do not reconcile: {row}")
        raw_cells_by_well[well_index] += 1
        raw_umis_by_well[well_index] += exact_umis
        label = dominant_label(human_umis, pig_umis)
        for rule_id, min_umis, min_genes in CELL_RULES:
            if exact_umis >= min_umis and genes >= min_genes:
                cell_counts[(rule_id, well_index, label)] += 1
                cell_umis[(rule_id, well_index, "human")] += human_umis
                cell_umis[(rule_id, well_index, "pig")] += pig_umis

    seen_wells: set[int] = set()
    per_well_rows: list[dict[str, object]] = []
    for manifest in manifest_rows:
        pair_id = manifest["pair_id"]
        sample = manifest["sample"]
        well_index = trailing_number(sample)
        if not 1 <= well_index <= args.expected_wells or well_index in seen_wells:
            raise ValueError(f"invalid/duplicate RNA well: {well_index}")
        seen_wells.add(well_index)
        expected_set = (
            "GM12878" if well_index <= args.human_only_wells else "cellMix"
        )
        if pair_id not in per_file or pair_id not in demux:
            raise ValueError(f"missing per-file or demultiplex row: {pair_id}")
        pf = per_file[pair_id]
        if pf["sample"] != sample:
            raise ValueError(f"RNA per-file sample mismatch: {pair_id}")
        audit_path = alignment_root / f"{pair_id}.rna_species_alignment_audit.tsv"
        audit_rows = read_tsv(audit_path)
        if len(audit_rows) != 1:
            raise ValueError(f"expected one alignment audit row: {audit_path}")
        audit = audit_rows[0]
        if (
            audit.get("status") != "validated"
            or audit.get("pipeline_signature") != ALIGNMENT_SIGNATURE
            or audit.get("pair_id") != pair_id
            or audit.get("sample") != sample
        ):
            raise ValueError(f"invalid alignment audit: {audit_path}")
        summary_path = genefull_root / pair_id / f"{pair_id}.summary.json"
        if not summary_path.is_file() or summary_path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty GeneFull summary: {summary_path}")
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        if payload.get("status") != "validated" or payload.get("pair_id") != pair_id:
            raise ValueError(f"invalid GeneFull summary: {summary_path}")
        counters = payload["counters"]

        primary_mapped = as_nonnegative(audit, "primary_mapped")
        unique_primary = as_nonnegative(audit, "unique_primary")
        eligible = as_nonnegative(audit, "tag_eligible_unique_primary")
        if primary_mapped != int(manifest["primary_mapped_reads"]):
            raise ValueError(f"{pair_id}: primary mapping total differs")
        if unique_primary != int(pf["unique_primary_reads"]):
            raise ValueError(f"{pair_id}: unique mapping total differs")
        if eligible != int(pf["tag_eligible_unique_primary"]):
            raise ValueError(f"{pair_id}: eligible unique total differs")

        human_assigned = int(counters.get("human_assigned_reads", 0))
        pig_assigned = int(counters.get("pig_assigned_reads", 0))
        other_assigned = int(counters.get("other_assigned_reads", 0))
        gene_assigned = int(counters.get("gene_assigned", 0))
        if human_assigned + pig_assigned + other_assigned != gene_assigned:
            raise ValueError(f"{pair_id}: species gene-assignment totals differ")
        human_umis = int(pf["human_exact_umis"])
        pig_umis = int(pf["pig_exact_umis"])
        exact_umis = int(pf["exact_umis"])
        if human_umis + pig_umis != exact_umis:
            raise ValueError(f"{pair_id}: species exact-UMI totals differ")

        dm = demux[pair_id]
        total_pairs = int(dm["total_pairs"])
        strict_assigned = int(dm["strict_assigned"])
        dna_residual = int(dm["DNA_residual_GTCA_linker_given_T7"])
        if min(total_pairs, strict_assigned, dna_residual) < 0:
            raise ValueError(f"{pair_id}: negative demultiplex count")

        row_out: dict[str, object] = {
            "well_index": well_index,
            "expected_well_set": expected_set,
            "pair_id": pair_id,
            "sample": sample,
            "demux_total_pairs": total_pairs,
            "demux_strict_assigned": strict_assigned,
            "demux_strict_fraction": ratio(strict_assigned, total_pairs),
            "demux_DNA_residual_GTCA": dna_residual,
            "demux_DNA_residual_fraction": ratio(dna_residual, total_pairs),
            "STAR_input_reads": int(pf["star_input_reads"]),
            "primary_mapped": primary_mapped,
            "primary_human": int(audit["primary_human"]),
            "primary_pig": int(audit["primary_pig"]),
            "primary_human_fraction": audit["primary_human_fraction"],
            "unique_primary": unique_primary,
            "unique_human": int(audit["unique_human"]),
            "unique_pig": int(audit["unique_pig"]),
            "unique_human_fraction": audit["unique_human_fraction"],
            "eligible_unique_primary": eligible,
            "eligible_human": int(audit["eligible_human"]),
            "eligible_pig": int(audit["eligible_pig"]),
            "eligible_spliced_human": int(audit["eligible_spliced_human"]),
            "eligible_spliced_pig": int(audit["eligible_spliced_pig"]),
            "eligible_human_spliced_fraction": audit[
                "eligible_human_spliced_fraction"
            ],
            "eligible_pig_spliced_fraction": audit[
                "eligible_pig_spliced_fraction"
            ],
            "gene_assigned_reads": gene_assigned,
            "gene_assigned_human": human_assigned,
            "gene_assigned_pig": pig_assigned,
            "gene_assigned_human_fraction": ratio(
                human_assigned, human_assigned + pig_assigned
            ),
            "exact_umis": exact_umis,
            "human_exact_umis": human_umis,
            "pig_exact_umis": pig_umis,
            "human_exact_umi_fraction": ratio(human_umis, human_umis + pig_umis),
            "raw_matrix_cells": raw_cells_by_well[well_index],
            "raw_matrix_exact_umis": raw_umis_by_well[well_index],
        }
        for rule_id, _, _ in CELL_RULES:
            human_cells = cell_counts[(rule_id, well_index, "human")]
            pig_cells = cell_counts[(rule_id, well_index, "pig")]
            tie_cells = cell_counts[(rule_id, well_index, "tie")]
            prefix = rule_id.removeprefix("R_")
            row_out[f"{prefix}_cells"] = human_cells + pig_cells + tie_cells
            row_out[f"{prefix}_human_cells"] = human_cells
            row_out[f"{prefix}_pig_cells"] = pig_cells
            row_out[f"{prefix}_tie_cells"] = tie_cells
            row_out[f"{prefix}_human_cell_fraction"] = ratio(
                human_cells, human_cells + pig_cells
            )
            row_out[f"{prefix}_human_umis"] = cell_umis[
                (rule_id, well_index, "human")
            ]
            row_out[f"{prefix}_pig_umis"] = cell_umis[
                (rule_id, well_index, "pig")
            ]
        unique_human_fraction = ratio(
            int(audit["unique_human"]),
            int(audit["unique_human"]) + int(audit["unique_pig"]),
        )
        row_out["unexpected_pig_unique_signal_in_human_only_well"] = int(
            expected_set == "GM12878"
            and unique_human_fraction != "NA"
            and float(unique_human_fraction) < 0.95
        )
        per_well_rows.append(row_out)

    if seen_wells != set(range(1, args.expected_wells + 1)):
        raise ValueError("RNA manifest does not cover every expected physical well")
    per_well_rows.sort(key=lambda row: int(row["well_index"]))

    sum_fields = [
        "demux_total_pairs",
        "demux_strict_assigned",
        "demux_DNA_residual_GTCA",
        "STAR_input_reads",
        "primary_mapped",
        "primary_human",
        "primary_pig",
        "unique_primary",
        "unique_human",
        "unique_pig",
        "eligible_unique_primary",
        "eligible_human",
        "eligible_pig",
        "eligible_spliced_human",
        "eligible_spliced_pig",
        "gene_assigned_reads",
        "gene_assigned_human",
        "gene_assigned_pig",
        "exact_umis",
        "human_exact_umis",
        "pig_exact_umis",
        "raw_matrix_cells",
        "raw_matrix_exact_umis",
    ]
    for rule_id, _, _ in CELL_RULES:
        prefix = rule_id.removeprefix("R_")
        sum_fields.extend(
            [
                f"{prefix}_cells",
                f"{prefix}_human_cells",
                f"{prefix}_pig_cells",
                f"{prefix}_tie_cells",
                f"{prefix}_human_umis",
                f"{prefix}_pig_umis",
            ]
        )

    per_set_rows: list[dict[str, object]] = []
    by_set = {
        name: [row for row in per_well_rows if row["expected_well_set"] == name]
        for name in ("GM12878", "cellMix")
    }
    for set_name, rows in by_set.items():
        if not rows:
            raise ValueError(f"empty expected well set: {set_name}")
        totals = {
            field: sum(int(row[field]) for row in rows) for field in sum_fields
        }
        set_row: dict[str, object] = {
            "expected_well_set": set_name,
            "wells": len(rows),
            **totals,
            "demux_strict_fraction": ratio(
                totals["demux_strict_assigned"], totals["demux_total_pairs"]
            ),
            "demux_DNA_residual_fraction": ratio(
                totals["demux_DNA_residual_GTCA"], totals["demux_total_pairs"]
            ),
            "primary_human_fraction": ratio(
                totals["primary_human"],
                totals["primary_human"] + totals["primary_pig"],
            ),
            "unique_human_fraction": ratio(
                totals["unique_human"],
                totals["unique_human"] + totals["unique_pig"],
            ),
            "eligible_human_spliced_fraction": ratio(
                totals["eligible_spliced_human"], totals["eligible_human"]
            ),
            "eligible_pig_spliced_fraction": ratio(
                totals["eligible_spliced_pig"], totals["eligible_pig"]
            ),
            "gene_assigned_human_fraction": ratio(
                totals["gene_assigned_human"],
                totals["gene_assigned_human"] + totals["gene_assigned_pig"],
            ),
            "human_exact_umi_fraction": ratio(
                totals["human_exact_umis"],
                totals["human_exact_umis"] + totals["pig_exact_umis"],
            ),
            "wells_with_unexpected_pig_unique_signal": sum(
                int(row["unexpected_pig_unique_signal_in_human_only_well"])
                for row in rows
            ),
        }
        for rule_id, _, _ in CELL_RULES:
            prefix = rule_id.removeprefix("R_")
            set_row[f"{prefix}_human_cell_fraction"] = ratio(
                totals[f"{prefix}_human_cells"],
                totals[f"{prefix}_human_cells"] + totals[f"{prefix}_pig_cells"],
            )
        per_set_rows.append(set_row)

    per_well_path = output_root / "rna_species_residual_audit_per_well.tsv"
    per_set_path = output_root / "rna_species_residual_audit_per_set.tsv"
    summary_path = output_root / "rna_species_residual_audit_summary.tsv"
    atomic_tsv(per_well_path, list(per_well_rows[0]), per_well_rows)
    atomic_tsv(per_set_path, list(per_set_rows[0]), per_set_rows)

    gm = next(row for row in per_set_rows if row["expected_well_set"] == "GM12878")
    all_totals = {
        field: sum(int(row[field]) for row in per_set_rows) for field in sum_fields
    }
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE},
        {"metric": "mode", "value": "audit_only_no_matrix_modification"},
        {"metric": "physical_wells", "value": args.expected_wells},
        {"metric": "human_only_wells_1_to_N", "value": args.human_only_wells},
        {
            "metric": "known_DNA_residual_GTCA_linker_pairs",
            "value": all_totals["demux_DNA_residual_GTCA"],
        },
        {
            "metric": "known_DNA_residual_GTCA_linker_fraction_of_raw_pairs",
            "value": ratio(
                all_totals["demux_DNA_residual_GTCA"],
                all_totals["demux_total_pairs"],
            ),
        },
        {
            "metric": "strict_RNA_assignment_fraction_of_raw_pairs",
            "value": ratio(
                all_totals["demux_strict_assigned"], all_totals["demux_total_pairs"]
            ),
        },
        {
            "metric": "GM12878_primary_human_fraction",
            "value": gm["primary_human_fraction"],
        },
        {
            "metric": "GM12878_unique_human_fraction",
            "value": gm["unique_human_fraction"],
        },
        {
            "metric": "GM12878_eligible_human_spliced_fraction",
            "value": gm["eligible_human_spliced_fraction"],
        },
        {
            "metric": "GM12878_eligible_pig_spliced_fraction",
            "value": gm["eligible_pig_spliced_fraction"],
        },
        {
            "metric": "GM12878_gene_assigned_human_fraction",
            "value": gm["gene_assigned_human_fraction"],
        },
        {
            "metric": "GM12878_human_exact_UMI_fraction",
            "value": gm["human_exact_umi_fraction"],
        },
        {
            "metric": "GM12878_R_UMI20_gene10_human_cell_fraction",
            "value": gm["UMI20_gene10_human_cell_fraction"],
        },
        {
            "metric": "GM12878_R_UMI25_gene10_human_cell_fraction",
            "value": gm["UMI25_gene10_human_cell_fraction"],
        },
        {
            "metric": "GM12878_wells_with_unique_human_fraction_below_0.95",
            "value": gm["wells_with_unexpected_pig_unique_signal"],
        },
        {
            "metric": "biological_species_validation",
            "value": (
                "passed_expected_human_only_reference"
                if gm["unique_human_fraction"] != "NA"
                and float(gm["unique_human_fraction"]) >= 0.95
                else "requires_review"
            ),
        },
        {
            "metric": "DNA_residual_interpretation",
            "value": (
                "GTCA-linker residual is counted before strict RNA assignment; "
                "causal leakage must be assessed with unique and splice-layer metrics"
            ),
        },
        {"metric": "per_well_output", "value": per_well_path},
        {"metric": "per_set_output", "value": per_set_path},
    ]
    atomic_tsv(summary_path, ["metric", "value"], summary_rows)
    print(
        f"[OK] RNA species/residual audit: GM unique human fraction="
        f"{gm['unique_human_fraction']} output={summary_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
