#!/usr/bin/env python3
"""Step 116: summarize full RNA GeneFull/exon/junction evidence recount."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import statistics
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_rna_expression_evidence_summary_20260806"
RECOUNT_SIGNATURE = "snmar_rna_expression_evidence_recount_20260806"
MODES = ("genefull", "exon_reverse", "junction_genefull")
RULES = (
    ("UMI5_gene3", 5, 3),
    ("UMI10_gene5", 10, 5),
    ("UMI20_gene10", 20, 10),
    ("UMI25_gene10", 25, 10),
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
        raise ValueError(f"sample lacks trailing physical-well number: {value}")
    return int(match.group(1))


def ratio(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def dominant_label(human: int, pig: int) -> str:
    if human == pig:
        return "tie"
    return "human" if human > pig else "pig"


def format_value(value: object) -> object:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.10g}" if math.isfinite(value) else "NA"
    return value


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
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: format_value(row.get(field)) for field in fields})
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def summarize_cells(rows: list[dict[str, object]]) -> dict[str, object]:
    labels = Counter(str(row["dominant_species"]) for row in rows)
    human_reads = sum(int(row["human_assigned_reads"]) for row in rows)
    pig_reads = sum(int(row["pig_assigned_reads"]) for row in rows)
    human_umis = sum(int(row["human_exact_umis"]) for row in rows)
    pig_umis = sum(int(row["pig_exact_umis"]) for row in rows)
    human_fractions = [
        float(row["human_UMI_fraction"])
        for row in rows
        if row["human_UMI_fraction"] is not None
    ]
    return {
        "cells": len(rows),
        "exact_umis": human_umis + pig_umis,
        "human_cells": labels["human"],
        "pig_cells": labels["pig"],
        "tie_cells": labels["tie"],
        "human_cell_fraction": ratio(labels["human"], len(rows)),
        "human_assigned_read_fraction": ratio(human_reads, human_reads + pig_reads),
        "human_exact_UMI_fraction": ratio(human_umis, human_umis + pig_umis),
        "median_cell_human_UMI_fraction": statistics.median(human_fractions) if human_fractions else None,
    }


def build_plot(per_group_rows: list[dict[str, object]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gm = [row for row in per_group_rows if row["well_group"] == "GM12878_only_wells"]
    all_wells = [row for row in per_group_rows if row["well_group"] == "all_wells"]
    rule_names = [rule[0] for rule in RULES]
    colors = {"genefull": "#777777", "exon_reverse": "#2ca58d", "junction_genefull": "#7a5195"}
    labels = {"genefull": "GeneFull", "exon_reverse": "exon reverse", "junction_genefull": "junction GeneFull"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    x = list(range(len(rule_names)))
    for mode in MODES:
        gm_mode = {str(row["rule"]): row for row in gm if row["mode"] == mode}
        all_mode = {str(row["rule"]): row for row in all_wells if row["mode"] == mode}
        axes[0, 0].plot(x, [all_mode[rule]["cells"] for rule in rule_names], marker="o", color=colors[mode], label=labels[mode])
        axes[0, 1].plot(x, [gm_mode[rule]["human_cell_fraction"] for rule in rule_names], marker="o", color=colors[mode], label=labels[mode])
        axes[1, 0].plot(x, [gm_mode[rule]["human_exact_UMI_fraction"] for rule in rule_names], marker="o", color=colors[mode], label=labels[mode])
        axes[1, 1].plot(x, [gm_mode[rule]["human_assigned_read_fraction"] for rule in rule_names], marker="o", color=colors[mode], label=labels[mode])
    for axis in axes.flat:
        axis.set_xticks(x, rule_names, rotation=25, ha="right")
        axis.grid(axis="y", alpha=0.2)
    axes[0, 0].set(title="Cells retained across 96 wells", ylabel="cells")
    axes[0, 1].set(title="GM-only wells: human-dominant cells", ylabel="fraction", ylim=(0, 1.02))
    axes[1, 0].set(title="GM-only wells: human exact UMIs", ylabel="fraction", ylim=(0, 1.02))
    axes[1, 1].set(title="GM-only wells: human assigned reads", ylabel="fraction", ylim=(0, 1.02))
    axes[0, 0].legend(frameon=False)
    fig.suptitle("snMAR RNA expression-evidence salvage audit", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--by-file-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--expected-wells", type=int, default=96)
    parser.add_argument("--human-only-wells", type=int, default=48)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.expected_wells != 96 or args.human_only_wells != 48:
        raise ValueError("expected confirmed 48 GM12878-only + 48 mix-well design")
    project = args.project
    manifest_path = args.manifest or (
        project / "08_qc_summary/14_postmapping_full_inputs/rna_reverse_genefull_full_manifest.tsv"
    )
    by_file_root = args.by_file_root or (
        project / "10_downstream_analysis/05_rna_expression_evidence_audit/by_file"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/35_rna_expression_evidence_audit"
    )
    output_root = args.output_root or (
        project / "10_downstream_analysis/05_rna_expression_evidence_audit"
    )

    manifests = read_tsv(manifest_path)
    if len(manifests) != args.expected_wells:
        raise ValueError(f"expected 96 manifest rows, found {len(manifests)}")
    all_cells: list[dict[str, object]] = []
    seen_wells: set[int] = set()
    per_file_validation: list[dict[str, object]] = []
    for manifest in manifests:
        pair_id = manifest["pair_id"]
        sample = manifest["sample"]
        well_index = trailing_number(sample)
        if well_index in seen_wells or not 1 <= well_index <= args.expected_wells:
            raise ValueError(f"invalid/duplicate well: {well_index}")
        seen_wells.add(well_index)
        root = by_file_root / pair_id
        summary_path = root / f"{pair_id}.evidence_mode_summary.json"
        if not summary_path.is_file() or summary_path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty summary: {summary_path}")
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        if (
            payload.get("status") != "validated"
            or payload.get("pipeline_signature") != RECOUNT_SIGNATURE
            or payload.get("pair_id") != pair_id
            or payload.get("sample") != sample
            or tuple(payload.get("modes", ())) != MODES
        ):
            raise ValueError(f"invalid evidence recount: {summary_path}")
        counters = payload["counters"]
        if int(counters["primary_records"]) != int(manifest["star_input_reads"]):
            raise ValueError(f"{pair_id}: STAR input reconciliation failed")
        if int(counters["primary_mapped"]) != int(manifest["primary_mapped_reads"]):
            raise ValueError(f"{pair_id}: primary mapping reconciliation failed")
        cell_qc_path = Path(payload["cell_qc"])
        mode_seen: Counter[str] = Counter()
        for source in iter_tsv(cell_qc_path):
            mode = source["mode"]
            if mode not in MODES:
                raise ValueError(f"unknown mode: {mode}")
            if source["pair_id"] != pair_id or source["sample"] != sample:
                raise ValueError(f"cell-QC identity mismatch: {pair_id}")
            exact_umis = int(source["exact_umis"])
            genes = int(source["genes_detected"])
            human_reads = int(source["human_assigned_reads"])
            pig_reads = int(source["pig_assigned_reads"])
            human_umis = int(source["human_exact_umis"])
            pig_umis = int(source["pig_exact_umis"])
            if min(exact_umis, genes, human_reads, pig_reads, human_umis, pig_umis) < 0:
                raise ValueError(f"negative RNA evidence count: {source}")
            if human_umis + pig_umis != exact_umis:
                raise ValueError(f"species UMI counts do not reconcile: {source}")
            if human_umis > human_reads or pig_umis > pig_reads:
                raise ValueError(f"species UMI count exceeds reads: {source}")
            mode_seen[mode] += 1
            all_cells.append(
                {
                    "well_index": well_index,
                    "well_group": "GM12878_only_wells" if well_index <= args.human_only_wells else "GM12878_PK15_mix_wells",
                    "pair_id": pair_id,
                    "sample": sample,
                    "cell_id": source["cell_id"],
                    "CB": source["CB"],
                    "mode": mode,
                    "assigned_reads": human_reads + pig_reads,
                    "exact_umis": exact_umis,
                    "genes_detected": genes,
                    "human_assigned_reads": human_reads,
                    "pig_assigned_reads": pig_reads,
                    "human_exact_umis": human_umis,
                    "pig_exact_umis": pig_umis,
                    "human_UMI_fraction": ratio(human_umis, exact_umis),
                    "dominant_species": dominant_label(human_umis, pig_umis),
                }
            )
        if set(mode_seen) != set(MODES) or len(set(mode_seen.values())) != 1:
            raise ValueError(f"{pair_id}: incomplete/unbalanced cell-mode rows: {mode_seen}")
        per_file_validation.append(
            {
                "well_index": well_index,
                "pair_id": pair_id,
                "sample": sample,
                **{f"{mode}_exact_umis": int(payload["mode_totals"][mode]["exact_umis"]) for mode in MODES},
                "status": "validated",
            }
        )

    per_group_rule: list[dict[str, object]] = []
    per_well_rule: list[dict[str, object]] = []
    groups = ("all_wells", "GM12878_only_wells", "GM12878_PK15_mix_wells")
    for mode in MODES:
        mode_rows = [row for row in all_cells if row["mode"] == mode]
        for rule, min_umis, min_genes in RULES:
            selected = [row for row in mode_rows if int(row["exact_umis"]) >= min_umis and int(row["genes_detected"]) >= min_genes]
            for group in groups:
                subset = selected if group == "all_wells" else [row for row in selected if row["well_group"] == group]
                per_group_rule.append(
                    {
                        "mode": mode,
                        "rule": rule,
                        "min_umis": min_umis,
                        "min_genes": min_genes,
                        "well_group": group,
                        **summarize_cells(subset),
                    }
                )
            for well_index in range(1, args.expected_wells + 1):
                subset = [row for row in selected if row["well_index"] == well_index]
                per_well_rule.append(
                    {
                        "well_index": well_index,
                        "well_group": "GM12878_only_wells" if well_index <= args.human_only_wells else "GM12878_PK15_mix_wells",
                        "mode": mode,
                        "rule": rule,
                        **summarize_cells(subset),
                    }
                )

    candidate_rule = "UMI10_gene5"
    candidate_rows = [
        row
        for row in per_group_rule
        if row["rule"] == candidate_rule
        and row["well_group"] == "GM12878_only_wells"
        and row["mode"] in {"exon_reverse", "junction_genefull"}
    ]
    viable = [row for row in candidate_rows if int(row["cells"]) >= 100]
    candidate = max(
        viable,
        key=lambda row: (
            float(row["human_cell_fraction"] or 0.0),
            int(row["cells"]),
        ),
        default=None,
    )
    candidate_mode = str(candidate["mode"]) if candidate else "none_insufficient_cells"
    candidate_human_fraction = candidate["human_cell_fraction"] if candidate else None
    ready = int(
        candidate is not None
        and float(candidate["human_cell_fraction"] or 0.0) >= 0.80
        and float(candidate["human_exact_UMI_fraction"] or 0.0) >= 0.80
    )
    summary_rows = [
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE, "interpretation": "validated implementation signature"},
        {"metric": "status", "value": "validated", "interpretation": "all 96 BAM recounts and count invariants passed"},
        {"metric": "completed_wells", "value": len(seen_wells), "interpretation": "physical wells"},
        {"metric": "candidate_filter_rule", "value": candidate_rule, "interpretation": "low-depth salvage comparison rule"},
        {"metric": "candidate_expression_mode", "value": candidate_mode, "interpretation": "strict mode with best GM-only human-cell recovery and at least 100 cells"},
        {"metric": "candidate_GM_only_human_cell_fraction", "value": candidate_human_fraction, "interpretation": "must be interpreted with retained-cell count"},
        {"metric": "RNA_expression_UMAP_ready", "value": ready, "interpretation": "1 requires >=80% human cells and human UMIs in GM-only wells"},
        {"metric": "old_genefull_matrix_replaced", "value": 0, "interpretation": "this stage is audit-only until a strict mode is selected"},
    ]
    qc_root.mkdir(parents=True, exist_ok=True)
    atomic_tsv(qc_root / "rna_expression_evidence_audit_summary.tsv", ["metric", "value", "interpretation"], summary_rows)
    atomic_tsv(qc_root / "rna_expression_evidence_per_group_rule.tsv", list(per_group_rule[0].keys()), per_group_rule)
    atomic_tsv(qc_root / "rna_expression_evidence_per_well_rule.tsv", list(per_well_rule[0].keys()), per_well_rule)
    per_file_validation.sort(key=lambda row: int(row["well_index"]))
    atomic_tsv(qc_root / "rna_expression_evidence_per_well_validation.tsv", list(per_file_validation[0].keys()), per_file_validation)
    build_plot(per_group_rule, output_root / "rna_expression_evidence_audit.png")
    print(f"[OK] RNA evidence audit: candidate={candidate_mode} ready={ready}")
    print(f"[OK] output: {qc_root / 'rna_expression_evidence_audit_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
