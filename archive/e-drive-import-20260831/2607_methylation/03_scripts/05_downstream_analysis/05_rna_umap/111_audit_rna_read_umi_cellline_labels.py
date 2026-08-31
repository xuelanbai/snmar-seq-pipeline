#!/usr/bin/env python3
"""Step 111: audit RNA human/pig labels before expression UMAP.

The same GeneFull counting run reports species at two layers:

* assigned reads, before exact-UMI collapsing;
* exact UMIs, after collapsing within (physical well, CB, gene, UMI).

For this batch the two layers disagree strongly in the GM12878-only wells.
This script keeps the layers separate, annotates every raw-matrix cell, and
produces the diagnostic needed before RNA labels are used to validate DNA UMAP.
It does not modify the RNA matrix and does not silently choose one label layer.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import statistics
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_rna_read_umi_cellline_audit_20260806"
GROUP_ORDER = ("all_wells", "GM12878_only_wells", "GM12878_PK15_mix_wells")


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


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"sample lacks trailing physical-well number: {value}")
    return int(match.group(1))


def as_nonnegative_int(row: dict[str, str], field: str) -> int:
    try:
        value = int(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid integer field {field}: {row}") from exc
    if value < 0:
        raise ValueError(f"negative {field}: {row}")
    return value


def fraction(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def format_value(value: object) -> object:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.10g}"
    return value


def classify_species(
    human: int, pig: int, minimum_major_fraction: float = 0.90
) -> tuple[str, float | None]:
    """Return a high-confidence species label and its major fraction."""
    if not 0.5 < minimum_major_fraction <= 1.0:
        raise ValueError("minimum_major_fraction must be in (0.5, 1]")
    total = human + pig
    if total == 0:
        return "no_signal", None
    major_fraction = max(human, pig) / total
    if major_fraction < minimum_major_fraction:
        return "ambiguous", major_fraction
    return ("GM12878" if human > pig else "PK15"), major_fraction


def comparison_status(read_label: str, umi_label: str) -> str:
    confident = {"GM12878", "PK15"}
    read_confident = read_label in confident
    umi_confident = umi_label in confident
    if read_confident and umi_confident:
        return "concordant" if read_label == umi_label else "discordant"
    if read_confident:
        return "read_only_confident"
    if umi_confident:
        return "umi_only_confident"
    return "neither_confident"


def reads_per_umi(reads: int, umis: int) -> float | None:
    return reads / umis if umis else None


def median_or_none(values: Iterable[float | None]) -> float | None:
    present = [value for value in values if value is not None and math.isfinite(value)]
    return statistics.median(present) if present else None


def atomic_tsv(
    path: Path, fields: list[str], rows: Iterable[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as raw:
            if path.suffix == ".gz":
                text_handle = gzip.open(raw, "wt", encoding="utf-8", newline="")
            else:
                text_handle = os.fdopen(os.dup(raw.fileno()), "w", encoding="utf-8", newline="")
            with text_handle as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
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


def well_design(well_index: int, human_only_wells: int) -> str:
    return (
        "GM12878_only_wells"
        if well_index <= human_only_wells
        else "GM12878_PK15_mix_wells"
    )


def recommendation_from_gm_group(metrics: dict[str, object]) -> str:
    """Conservative recommendation; never hides read/UMI disagreement."""
    read_human = metrics.get("read_GM12878_fraction")
    umi_human = metrics.get("umi_GM12878_fraction")
    discordant = metrics.get("confident_discordant_fraction")
    if not all(isinstance(value, float) for value in (read_human, umi_human, discordant)):
        return "insufficient_GM12878_only_cells"
    if read_human >= 0.80 and umi_human >= 0.80 and discordant <= 0.10:
        return "read_and_UMI_species_labels_supported"
    if read_human >= 0.80 and umi_human < 0.50:
        return "read_layer_provisional_only_do_not_use_UMI_species_labels"
    return "do_not_use_RNA_species_labels_until_reviewed"


def summarize_group(rows: list[dict[str, object]]) -> dict[str, object]:
    labels_read = Counter(str(row["read_species_label"]) for row in rows)
    labels_umi = Counter(str(row["umi_species_label"]) for row in rows)
    statuses = Counter(str(row["read_umi_status"]) for row in rows)
    confident_read = labels_read["GM12878"] + labels_read["PK15"]
    confident_umi = labels_umi["GM12878"] + labels_umi["PK15"]
    both_confident = statuses["concordant"] + statuses["discordant"]
    human_reads = sum(int(row["human_assigned_reads"]) for row in rows)
    pig_reads = sum(int(row["pig_assigned_reads"]) for row in rows)
    human_umis = sum(int(row["human_exact_umis"]) for row in rows)
    pig_umis = sum(int(row["pig_exact_umis"]) for row in rows)
    return {
        "retained_cells": len(rows),
        "read_GM12878_cells": labels_read["GM12878"],
        "read_PK15_cells": labels_read["PK15"],
        "read_ambiguous_or_no_signal_cells": len(rows) - confident_read,
        "read_GM12878_fraction": fraction(labels_read["GM12878"], len(rows)),
        "umi_GM12878_cells": labels_umi["GM12878"],
        "umi_PK15_cells": labels_umi["PK15"],
        "umi_ambiguous_or_no_signal_cells": len(rows) - confident_umi,
        "umi_GM12878_fraction": fraction(labels_umi["GM12878"], len(rows)),
        "both_layers_confident_cells": both_confident,
        "confident_concordant_cells": statuses["concordant"],
        "confident_discordant_cells": statuses["discordant"],
        "confident_discordant_fraction": fraction(statuses["discordant"], both_confident),
        "aggregate_human_assigned_read_fraction": fraction(human_reads, human_reads + pig_reads),
        "aggregate_human_exact_UMI_fraction": fraction(human_umis, human_umis + pig_umis),
        "median_read_human_fraction": median_or_none(
            row["read_human_fraction"] for row in rows  # type: ignore[arg-type]
        ),
        "median_UMI_human_fraction": median_or_none(
            row["umi_human_fraction"] for row in rows  # type: ignore[arg-type]
        ),
        "median_human_reads_per_human_UMI": median_or_none(
            row["human_reads_per_human_UMI"] for row in rows  # type: ignore[arg-type]
        ),
        "median_pig_reads_per_pig_UMI": median_or_none(
            row["pig_reads_per_pig_UMI"] for row in rows  # type: ignore[arg-type]
        ),
    }


def build_plot(rows: list[dict[str, object]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    retained = [row for row in rows if row["passes_expression_QC"] == 1]
    gm = [row for row in retained if row["well_design"] == "GM12878_only_wells"]
    mix = [row for row in retained if row["well_design"] == "GM12878_PK15_mix_wells"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    for subset, color, label in ((gm, "#2ca58d", "GM12878-only wells"), (mix, "#888888", "mix wells")):
        x = [row["read_human_fraction"] for row in subset if row["read_human_fraction"] is not None and row["umi_human_fraction"] is not None]
        y = [row["umi_human_fraction"] for row in subset if row["read_human_fraction"] is not None and row["umi_human_fraction"] is not None]
        axes[0, 0].scatter(x, y, s=7, alpha=0.35, color=color, label=label, linewidths=0)
    axes[0, 0].plot([0, 1], [0, 1], color="black", linewidth=1, linestyle="--")
    axes[0, 0].set(xlabel="human fraction of assigned reads", ylabel="human fraction of exact UMIs", title="Per-cell read vs exact-UMI species fraction")
    axes[0, 0].legend(frameon=False, markerscale=2)

    axes[0, 1].hist(
        [row["read_human_fraction"] for row in gm if row["read_human_fraction"] is not None],
        bins=40,
        range=(0, 1),
        histtype="step",
        linewidth=2,
        color="#1f77b4",
        label="assigned reads",
    )
    axes[0, 1].hist(
        [row["umi_human_fraction"] for row in gm if row["umi_human_fraction"] is not None],
        bins=40,
        range=(0, 1),
        histtype="step",
        linewidth=2,
        color="#d62728",
        label="exact UMIs",
    )
    axes[0, 1].set(xlabel="human fraction", ylabel="cells", title="GM12878-only wells")
    axes[0, 1].legend(frameon=False)

    dup_rows = [
        row
        for row in gm
        if row["human_reads_per_human_UMI"] is not None
        and row["pig_reads_per_pig_UMI"] is not None
        and float(row["human_reads_per_human_UMI"]) > 0
        and float(row["pig_reads_per_pig_UMI"]) > 0
    ]
    axes[1, 0].scatter(
        [math.log10(float(row["human_reads_per_human_UMI"])) for row in dup_rows],
        [math.log10(float(row["pig_reads_per_pig_UMI"])) for row in dup_rows],
        s=8,
        alpha=0.4,
        color="#6a51a3",
        linewidths=0,
    )
    axes[1, 0].plot([0, 4], [0, 4], color="black", linewidth=1, linestyle="--")
    axes[1, 0].set(xlabel="log10(human assigned reads / human exact UMI)", ylabel="log10(pig assigned reads / pig exact UMI)", title="Species-specific collapsing in GM12878-only wells")

    by_well: Counter[int] = Counter(int(row["well_index"]) for row in retained)
    axes[1, 1].bar(range(1, 97), [by_well[index] for index in range(1, 97)], color=["#2ca58d" if index <= 48 else "#888888" for index in range(1, 97)], width=0.9)
    axes[1, 1].axvline(48.5, color="black", linewidth=1)
    axes[1, 1].set(xlabel="physical well", ylabel="retained cells", title="Cells passing RNA expression QC")

    fig.suptitle("snMAR RNA species audit: assigned reads versus exact UMIs", fontsize=15)
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
    parser.add_argument("--rna-cell-qc", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--expected-wells", type=int, default=96)
    parser.add_argument("--human-only-wells", type=int, default=48)
    parser.add_argument("--min-umis", type=int, default=25)
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--minimum-major-fraction", type=float, default=0.90)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.expected_wells != 96 or args.human_only_wells != 48:
        raise ValueError("this audit expects the confirmed 48 GM12878-only + 48 mix-well design")
    if args.min_umis < 1 or args.min_genes < 1:
        raise ValueError("RNA expression thresholds must be positive")
    project = args.project
    cell_qc_path = args.rna_cell_qc or (
        project / "08_qc_summary/16_rna_reverse_genefull_full/rna_reverse_cell_qc.tsv.gz"
    )
    output_root = args.output_root or (
        project / "10_downstream_analysis/04_rna_cell_line_audit"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/34_rna_cell_line_audit"
    )

    annotated: list[dict[str, object]] = []
    seen_cells: set[str] = set()
    seen_wells: set[int] = set()
    for source in iter_tsv(cell_qc_path):
        if source.get("in_raw_matrix") != "1":
            continue
        well_index = trailing_number(source["sample"])
        if not 1 <= well_index <= args.expected_wells:
            raise ValueError(f"cell outside expected wells: {source}")
        seen_wells.add(well_index)
        cell_id = source["cell_id"]
        if cell_id in seen_cells:
            raise ValueError(f"duplicate raw-matrix cell_id: {cell_id}")
        seen_cells.add(cell_id)

        gene_reads = as_nonnegative_int(source, "gene_assigned_reads")
        exact_umis = as_nonnegative_int(source, "exact_umis")
        genes = as_nonnegative_int(source, "genes_detected")
        human_reads = as_nonnegative_int(source, "human_assigned_reads")
        pig_reads = as_nonnegative_int(source, "pig_assigned_reads")
        human_umis = as_nonnegative_int(source, "human_exact_umis")
        pig_umis = as_nonnegative_int(source, "pig_exact_umis")
        if human_reads + pig_reads != gene_reads:
            raise ValueError(f"assigned-read species counts do not reconcile: {cell_id}")
        if human_umis + pig_umis != exact_umis:
            raise ValueError(f"exact-UMI species counts do not reconcile: {cell_id}")

        read_label, read_major = classify_species(
            human_reads, pig_reads, args.minimum_major_fraction
        )
        umi_label, umi_major = classify_species(
            human_umis, pig_umis, args.minimum_major_fraction
        )
        passes = int(exact_umis >= args.min_umis and genes >= args.min_genes)
        annotated.append(
            {
                "cell_id": cell_id,
                "sample": source["sample"],
                "well_index": well_index,
                "well_design": well_design(well_index, args.human_only_wells),
                "CB": source["CB"],
                "exact_umis": exact_umis,
                "genes_detected": genes,
                "gene_assigned_reads": gene_reads,
                "human_assigned_reads": human_reads,
                "pig_assigned_reads": pig_reads,
                "human_exact_umis": human_umis,
                "pig_exact_umis": pig_umis,
                "read_human_fraction": fraction(human_reads, human_reads + pig_reads),
                "umi_human_fraction": fraction(human_umis, human_umis + pig_umis),
                "read_species_label": read_label,
                "read_major_fraction": read_major,
                "umi_species_label": umi_label,
                "umi_major_fraction": umi_major,
                "read_umi_status": comparison_status(read_label, umi_label),
                "human_reads_per_human_UMI": reads_per_umi(human_reads, human_umis),
                "pig_reads_per_pig_UMI": reads_per_umi(pig_reads, pig_umis),
                "passes_expression_QC": passes,
                "expression_QC_rule": f"UMI>={args.min_umis};gene>={args.min_genes}",
            }
        )

    if seen_wells != set(range(1, args.expected_wells + 1)):
        missing = sorted(set(range(1, args.expected_wells + 1)) - seen_wells)
        raise ValueError(f"RNA cell-QC does not cover all 96 wells; missing={missing}")

    fields = list(annotated[0].keys()) if annotated else []
    if not fields:
        raise ValueError("no cells in RNA raw matrix")
    annotated_path = output_root / "rna_read_umi_cellline_labels.tsv.gz"
    atomic_tsv(annotated_path, fields, annotated)

    retained = [row for row in annotated if row["passes_expression_QC"] == 1]
    groups = {
        "all_wells": retained,
        "GM12878_only_wells": [row for row in retained if row["well_design"] == "GM12878_only_wells"],
        "GM12878_PK15_mix_wells": [row for row in retained if row["well_design"] == "GM12878_PK15_mix_wells"],
    }
    group_metrics = {group: summarize_group(rows) for group, rows in groups.items()}
    per_group_rows: list[dict[str, object]] = []
    for group in GROUP_ORDER:
        per_group_rows.append({"well_group": group, **group_metrics[group]})
    per_group_path = qc_root / "rna_read_umi_cellline_audit_per_group.tsv"
    atomic_tsv(per_group_path, list(per_group_rows[0].keys()), per_group_rows)

    per_well_rows: list[dict[str, object]] = []
    for well_index in range(1, args.expected_wells + 1):
        subset = [row for row in retained if row["well_index"] == well_index]
        metrics = summarize_group(subset)
        per_well_rows.append(
            {
                "well_index": well_index,
                "well_design": well_design(well_index, args.human_only_wells),
                **metrics,
            }
        )
    per_well_path = qc_root / "rna_read_umi_cellline_audit_per_well.tsv"
    atomic_tsv(per_well_path, list(per_well_rows[0].keys()), per_well_rows)

    gm_metrics = group_metrics["GM12878_only_wells"]
    recommendation = recommendation_from_gm_group(gm_metrics)
    summary_rows = [
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE, "interpretation": "validated implementation signature"},
        {"metric": "status", "value": "validated", "interpretation": "input counts reconcile; biological labels still require review"},
        {"metric": "raw_matrix_cells", "value": len(annotated), "interpretation": "all cells in the unfiltered RNA matrix"},
        {"metric": "expression_QC_rule", "value": f"UMI>={args.min_umis};gene>={args.min_genes}", "interpretation": "diagnostic RNA cell threshold"},
        {"metric": "retained_cells", "value": len(retained), "interpretation": "cells evaluated for RNA label choice and UMAP"},
        {"metric": "minimum_major_fraction", "value": args.minimum_major_fraction, "interpretation": "human/pig label confidence threshold"},
        {"metric": "GM_only_read_human_cell_fraction", "value": gm_metrics["read_GM12878_fraction"], "interpretation": "cell-level read-layer expectation check in wells 1-48"},
        {"metric": "GM_only_UMI_human_cell_fraction", "value": gm_metrics["umi_GM12878_fraction"], "interpretation": "cell-level exact-UMI expectation check in wells 1-48"},
        {"metric": "GM_only_confident_discordant_fraction", "value": gm_metrics["confident_discordant_fraction"], "interpretation": "read and exact-UMI high-confidence labels that disagree"},
        {"metric": "recommended_RNA_species_label_layer", "value": recommendation, "interpretation": "conservative gate before RNA labels validate DNA UMAP"},
        {"metric": "RNA_to_DNA_label_transfer_ready", "value": int(recommendation == "read_and_UMI_species_labels_supported"), "interpretation": "0 means do not transfer RNA species labels yet"},
    ]
    summary_path = qc_root / "rna_read_umi_cellline_audit_summary.tsv"
    atomic_tsv(summary_path, ["metric", "value", "interpretation"], summary_rows)

    plot_path = output_root / "rna_read_umi_cellline_audit.png"
    build_plot(annotated, plot_path)
    print(f"[OK] raw cells={len(annotated):,} retained={len(retained):,}")
    print(f"[OK] recommendation={recommendation}")
    print(f"[OK] summary: {summary_path}")
    print(f"[OK] annotated cells: {annotated_path}")
    print(f"[OK] plot: {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
