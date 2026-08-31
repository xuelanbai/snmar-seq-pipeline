#!/usr/bin/env python3
"""Step 83: build complete per-cell WCG/GCH coverage QC.

This step combines the cell-resolved matrix summaries with the reference
context denominators calculated for the BigWig tracks. No cell filtering or
biological species/doublet threshold is applied here.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


CONTEXTS = ("WCG", "GCH")
SPECIES = ("human", "pig")
QUANTILES = (
    ("minimum", 0.0),
    ("p25", 0.25),
    ("median", 0.5),
    ("p75", 0.75),
    ("p90", 0.9),
    ("p95", 0.95),
    ("p99", 0.99),
    ("maximum", 1.0),
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


def atomic_metrics(path: Path, rows: list[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as handle:
            handle.write("metric\tvalue\n")
            for metric, value in rows:
                handle.write(f"{metric}\t{value}\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_nonnegative_integer(row: dict[str, str], field: str) -> int:
    value = int(row[field])
    if value < 0:
        raise ValueError(f"negative {field}: {row}")
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--input-per-cell", type=Path)
    parser.add_argument("--bigwig-summary", type=Path)
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    input_per_cell = args.input_per_cell or (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete/"
        "dna_single_cytosine_complete_per_cell.tsv.gz"
    )
    bigwig_summary = args.bigwig_summary or (
        project
        / "08_qc_summary/24_dna_single_cytosine_bigwig/"
        "dna_single_base_bigwig_summary.tsv"
    )
    output_root = args.output_root or (
        project / "08_qc_summary/25_dna_single_cell_qc"
    )
    output_root.mkdir(parents=True, exist_ok=True)

    track_metrics = read_metrics(bigwig_summary)
    if track_metrics.get("status") != "validated":
        raise ValueError("single-base BigWig result is not validated")
    denominators: dict[tuple[str, str], int] = {}
    for species in SPECIES:
        for context in CONTEXTS:
            key = f"{species}_{context}_reference_assayable_sites"
            value = int(track_metrics[key])
            if value < 1:
                raise ValueError(f"invalid reference denominator: {key}")
            denominators[(species, context)] = value

    if not input_per_cell.is_file() or input_per_cell.stat().st_size == 0:
        raise FileNotFoundError(input_per_cell)
    output = output_root / "dna_single_cell_coverage_qc.tsv.gz"
    descriptor, temporary = tempfile.mkstemp(
        prefix=output.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    distributions: defaultdict[str, list[float]] = defaultdict(list)
    category_counts: Counter[str] = Counter()
    threshold_counts: Counter[tuple[str, int]] = Counter()
    cell_ids: set[str] = set()
    rows_written = 0
    try:
        with (
            gzip.open(
                input_per_cell, "rt", encoding="utf-8", newline=""
            ) as source,
            gzip.open(
                temporary,
                "wt",
                encoding="utf-8",
                newline="",
                compresslevel=1,
            ) as target,
        ):
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"missing fields: {input_per_cell}")
            required = {
                "cell_id",
                "WCG_observed_sites",
                "WCG_covered_calls",
                "WCG_methylation_fraction",
                "GCH_observed_sites",
                "GCH_covered_calls",
                "GCH_accessibility_fraction",
                "human_signal_covered_calls",
                "pig_signal_covered_calls",
                "dominant_reference_signal",
            }
            required.update(
                f"{species}_{context}_observed_sites"
                for species in SPECIES
                for context in CONTEXTS
            )
            required.update(
                f"{species}_{context}_covered_calls"
                for species in SPECIES
                for context in CONTEXTS
            )
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(
                    "per-cell matrix QC lacks current fields; rerun steps "
                    "68-75: " + ", ".join(sorted(missing))
                )
            new_fields = [
                "total_signal_covered_calls",
                "major_reference_signal_fraction",
                *[
                    f"{species}_{context}_reference_coverage_fraction"
                    for species in SPECIES
                    for context in CONTEXTS
                ],
            ]
            writer = csv.DictWriter(
                target,
                fieldnames=reader.fieldnames + new_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for row in reader:
                cell_id = row["cell_id"]
                if not cell_id or cell_id in cell_ids:
                    raise ValueError(f"duplicate/empty cell_id: {cell_id}")
                cell_ids.add(cell_id)
                context_sites = {
                    context: parse_nonnegative_integer(
                        row, f"{context}_observed_sites"
                    )
                    for context in CONTEXTS
                }
                context_calls = {
                    context: parse_nonnegative_integer(
                        row, f"{context}_covered_calls"
                    )
                    for context in CONTEXTS
                }
                for context in CONTEXTS:
                    separated = sum(
                        parse_nonnegative_integer(
                            row, f"{species}_{context}_observed_sites"
                        )
                        for species in SPECIES
                    )
                    if separated != context_sites[context]:
                        raise ValueError(
                            f"species {context} sites do not reconcile: "
                            f"{cell_id}"
                        )
                    separated_calls = sum(
                        parse_nonnegative_integer(
                            row, f"{species}_{context}_covered_calls"
                        )
                        for species in SPECIES
                    )
                    if separated_calls != context_calls[context]:
                        raise ValueError(
                            f"species {context} calls do not reconcile: "
                            f"{cell_id}"
                        )
                    if context_calls[context] < context_sites[context]:
                        raise ValueError(
                            f"{context} calls below observed sites: {cell_id}"
                        )
                    distributions[f"{context}_observed_sites"].append(
                        context_sites[context]
                    )
                    distributions[f"{context}_covered_calls"].append(
                        context_calls[context]
                    )
                    for threshold in (1_000, 10_000, 100_000):
                        if context_sites[context] >= threshold:
                            threshold_counts[(context, threshold)] += 1
                for field in (
                    "WCG_methylation_fraction",
                    "GCH_accessibility_fraction",
                ):
                    if row[field] != "NA":
                        value = float(row[field])
                        if not 0.0 <= value <= 1.0:
                            raise ValueError(f"invalid {field}: {cell_id}")
                        distributions[field].append(value)

                human_calls = parse_nonnegative_integer(
                    row, "human_signal_covered_calls"
                )
                pig_calls = parse_nonnegative_integer(
                    row, "pig_signal_covered_calls"
                )
                expected_human_calls = sum(
                    parse_nonnegative_integer(
                        row, f"human_{context}_covered_calls"
                    )
                    for context in CONTEXTS
                )
                expected_pig_calls = sum(
                    parse_nonnegative_integer(
                        row, f"pig_{context}_covered_calls"
                    )
                    for context in CONTEXTS
                )
                if (
                    human_calls != expected_human_calls
                    or pig_calls != expected_pig_calls
                ):
                    raise ValueError(
                        f"species signal calls do not reconcile: {cell_id}"
                    )
                total_calls = human_calls + pig_calls
                major_fraction: float | str = (
                    max(human_calls, pig_calls) / total_calls
                    if total_calls
                    else "NA"
                )
                row["total_signal_covered_calls"] = total_calls
                row["major_reference_signal_fraction"] = major_fraction
                if total_calls:
                    distributions["major_reference_signal_fraction"].append(
                        float(major_fraction)
                    )
                category_counts[row["dominant_reference_signal"]] += 1

                for species in SPECIES:
                    for context in CONTEXTS:
                        observed = parse_nonnegative_integer(
                            row, f"{species}_{context}_observed_sites"
                        )
                        coverage_fraction = (
                            observed / denominators[(species, context)]
                        )
                        if not 0.0 <= coverage_fraction <= 1.0:
                            raise ValueError(
                                "reference coverage outside 0..1: "
                                f"{cell_id} {species} {context}"
                            )
                        field = (
                            f"{species}_{context}_"
                            "reference_coverage_fraction"
                        )
                        row[field] = coverage_fraction
                        distributions[field].append(coverage_fraction)
                writer.writerow(row)
                rows_written += 1
        if rows_written == 0:
            raise ValueError("no per-cell rows")
        os.replace(temporary, output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("modality", "DNA_WCG_methylation_and_GCH_accessibility"),
        ("cells", rows_written),
        ("input_per_cell", input_per_cell),
        ("per_cell_qc", output),
        ("coverage_denominator", "reference_assayable_sites_both_strands"),
        ("cell_filtering", "none"),
        ("species_doublet_threshold", "none"),
    ]
    for species in SPECIES:
        for context in CONTEXTS:
            summary_rows.append(
                (
                    f"{species}_{context}_reference_assayable_sites",
                    denominators[(species, context)],
                )
            )
    for context in CONTEXTS:
        summary_rows.append(
            (
                f"cells_with_{context}_observed_sites",
                sum(
                    value > 0
                    for value in distributions[
                        f"{context}_observed_sites"
                    ]
                ),
            )
        )
        for threshold in (1_000, 10_000, 100_000):
            summary_rows.append(
                (
                    f"cells_{context}_observed_sites_ge_{threshold}",
                    threshold_counts[(context, threshold)],
                )
            )
    for category in ("human", "pig", "tie", "no_signal"):
        summary_rows.append(
            (
                f"cells_dominant_reference_{category}",
                category_counts[category],
            )
        )
    for metric in sorted(distributions):
        values = np.asarray(distributions[metric], dtype=np.float64)
        summary_rows.append((f"{metric}_nonmissing_cells", len(values)))
        if len(values):
            for label, probability in QUANTILES:
                summary_rows.append(
                    (
                        f"{metric}_{label}",
                        float(np.quantile(values, probability)),
                    )
                )
    atomic_metrics(
        output_root / "dna_single_cell_coverage_qc_summary.tsv",
        summary_rows,
    )
    print(
        "[OK] DNA single-cell coverage QC: "
        f"cells={rows_written} output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
