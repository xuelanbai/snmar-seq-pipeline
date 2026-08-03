#!/usr/bin/env python3
"""Step 71: merge 96 site summaries into WCG and GCH track shards."""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import heapq
import os
import tempfile
from collections import Counter
from pathlib import Path


EXPECTED_SAMPLES = 96
VALIDATION_SIGNATURE = "snmar_single_cytosine_validation_20260801"
CONTEXT_ORDER = {"WCG": 1, "GCH": 2}
SIGNAL_TYPES = {
    "WCG": "endogenous_methylation",
    "GCH": "GpC_accessibility",
}
FIELDS = [
    "contig_index",
    "species",
    "chrom",
    "start",
    "end",
    "strand",
    "context",
    "signal_type",
    "covered_cells",
    "covered_calls",
    "converted_calls",
    "signal_fraction",
]


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


def atomic_summary(path: Path, rows: list[dict[str, object]]) -> None:
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
                fieldnames=["metric", "value"],
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
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.shard < 1 or args.shard > 8:
        raise ValueError("shard must be 1..8")
    project = args.project
    input_root = project / "08_qc_summary/20_single_cytosine_inputs"
    tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
    contig_index_path = input_root / "dna_single_cytosine_contig_index.tsv"
    for path in (tasks_path, contig_index_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    tasks = [
        row
        for row in read_tsv(tasks_path)
        if int(row["shard"]) == args.shard
    ]
    tasks.sort(key=lambda row: int(row["sample_index"]))
    if (
        len(tasks) != EXPECTED_SAMPLES
        or [int(row["sample_index"]) for row in tasks]
        != list(range(1, EXPECTED_SAMPLES + 1))
    ):
        raise ValueError("track shard does not have 96 sample inputs")
    contigs = {
        int(row["contig_index"]): row
        for row in read_tsv(contig_index_path)
        if int(row["shard"]) == args.shard
    }
    if not contigs:
        raise ValueError("track shard has no indexed contigs")
    inputs = []
    for task in tasks:
        sample = task["sample"]
        shard = int(task["shard"])
        prefix = f"{sample}.shard_{shard:02d}"
        site_path = (
            Path(task["output_dir"]) / f"{prefix}.site_summary.tsv.gz"
        )
        validation = (
            project
            / "08_qc_summary/21_dna_single_cytosine_validation"
            / sample
            / f"{sample}.single_cytosine.validation_summary.tsv"
        )
        for path in (site_path, validation):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
        metrics = read_metrics(validation)
        if (
            metrics.get("status") != "validated"
            or metrics.get("pipeline_signature") != VALIDATION_SIGNATURE
            or metrics.get("sample") != sample
        ):
            raise ValueError(f"sample validation failed: {validation}")
        inputs.append(site_path)

    output_root = (
        args.output_root
        or project
        / "07_postmapping/complete_matrices/dna_single_cytosine/tracks"
    )
    qc_root = (
        args.qc_root
        or project / "08_qc_summary/22_dna_single_cytosine_tracks"
    )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    detailed_path = (
        output_root
        / f"dna_single_cytosine.pseudobulk.shard_{args.shard:02d}.tsv.gz"
    )
    bedgraph_paths = {
        "WCG": output_root
        / f"dna_WCG_methylation.shard_{args.shard:02d}.bedGraph.gz",
        "GCH": output_root
        / f"dna_GCH_accessibility.shard_{args.shard:02d}.bedGraph.gz",
    }
    temporary_paths: dict[str, str] = {}
    descriptor, temporary_paths["detailed"] = tempfile.mkstemp(
        prefix=detailed_path.name + ".", suffix=".tmp", dir=output_root
    )
    os.close(descriptor)
    for context, path in bedgraph_paths.items():
        descriptor, temporary_paths[context] = tempfile.mkstemp(
            prefix=path.name + ".", suffix=".tmp", dir=output_root
        )
        os.close(descriptor)
    site_counts: Counter[str] = Counter()
    covered_cells: Counter[str] = Counter()
    covered: Counter[str] = Counter()
    converted: Counter[str] = Counter()
    try:
        with contextlib.ExitStack() as stack:
            readers = []
            for path in inputs:
                handle = stack.enter_context(
                    gzip.open(path, "rt", encoding="utf-8", newline="")
                )
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != FIELDS:
                    raise ValueError(f"unexpected site-summary fields: {path}")
                readers.append(reader)
            detailed = stack.enter_context(
                gzip.open(
                    temporary_paths["detailed"],
                    "wt",
                    encoding="utf-8",
                    newline="",
                    compresslevel=1,
                )
            )
            bedgraphs = {
                context: stack.enter_context(
                    gzip.open(
                        temporary_paths[context],
                        "wt",
                        encoding="utf-8",
                        newline="",
                        compresslevel=1,
                    )
                )
                for context in ("WCG", "GCH")
            }
            writer = csv.DictWriter(
                detailed,
                fieldnames=FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            heap: list[
                tuple[int, int, int, int, int, dict[str, str]]
            ] = []
            previous: list[tuple[int, int, int, int] | None] = [
                None
            ] * len(readers)

            def push_next(source_index: int) -> None:
                try:
                    row = next(readers[source_index])
                except StopIteration:
                    return
                contig_index = int(row["contig_index"])
                start = int(row["start"])
                strand = 1 if row["strand"] == "+" else -1
                context = row["context"]
                context_order = CONTEXT_ORDER.get(context)
                key = (contig_index, start, strand, context_order or 0)
                if (
                    context_order is None
                    or row["signal_type"] != SIGNAL_TYPES[context]
                    or (
                        previous[source_index] is not None
                        and key <= previous[source_index]
                    )
                ):
                    raise ValueError("sample site summary context/order invalid")
                previous[source_index] = key
                metadata = contigs.get(contig_index)
                if (
                    metadata is None
                    or row["chrom"] != metadata["contig"]
                    or row["species"] != metadata["species"]
                    or int(row["end"]) != start + 1
                    or start < 0
                    or start >= int(metadata["length"])
                ):
                    raise ValueError("site row conflicts with contig index")
                heapq.heappush(
                    heap,
                    (
                        contig_index,
                        start,
                        strand,
                        context_order,
                        source_index,
                        row,
                    ),
                )

            for index in range(len(readers)):
                push_next(index)
            while heap:
                (
                    contig_index,
                    start,
                    strand,
                    context_order,
                    source_index,
                    row,
                ) = heapq.heappop(heap)
                key = (contig_index, start, strand, context_order)
                aggregate_cells = int(row["covered_cells"])
                aggregate_covered = int(row["covered_calls"])
                aggregate_converted = int(row["converted_calls"])
                context = row["context"]
                push_next(source_index)
                while heap and heap[0][:4] == key:
                    *_key, source_index, row = heapq.heappop(heap)
                    aggregate_cells += int(row["covered_cells"])
                    aggregate_covered += int(row["covered_calls"])
                    aggregate_converted += int(row["converted_calls"])
                    push_next(source_index)
                if (
                    aggregate_cells < 1
                    or aggregate_covered < aggregate_cells
                    or aggregate_converted < 0
                    or aggregate_converted > aggregate_covered
                ):
                    raise ValueError("invalid aggregate cytosine counts")
                metadata = contigs[contig_index]
                fraction = aggregate_converted / aggregate_covered
                writer.writerow(
                    {
                        "contig_index": contig_index,
                        "species": metadata["species"],
                        "chrom": metadata["contig"],
                        "start": start,
                        "end": start + 1,
                        "strand": "+" if strand == 1 else "-",
                        "context": context,
                        "signal_type": SIGNAL_TYPES[context],
                        "covered_cells": aggregate_cells,
                        "covered_calls": aggregate_covered,
                        "converted_calls": aggregate_converted,
                        "signal_fraction": fraction,
                    }
                )
                bedgraphs[context].write(
                    f"{metadata['contig']}\t{start}\t{start + 1}\t"
                    f"{fraction:.8g}\n"
                )
                site_counts[context] += 1
                covered_cells[context] += aggregate_cells
                covered[context] += aggregate_covered
                converted[context] += aggregate_converted
        os.replace(temporary_paths["detailed"], detailed_path)
        for context, path in bedgraph_paths.items():
            os.replace(temporary_paths[context], path)
    except BaseException:
        for path in temporary_paths.values():
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
        raise
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "shard", "value": args.shard},
        {"metric": "input_samples", "value": EXPECTED_SAMPLES},
        {"metric": "contigs", "value": len(contigs)},
        {"metric": "detailed_track", "value": detailed_path},
        {"metric": "WCG_bedgraph", "value": bedgraph_paths["WCG"]},
        {"metric": "GCH_bedgraph", "value": bedgraph_paths["GCH"]},
        {
            "metric": "coordinate_system",
            "value": (
                "BED_0-based_reference_C_plus_reference_G_minus;"
                "strand_in_detailed_track"
            ),
        },
        {
            "metric": "context_model",
            "value": "WCG=methylation;GCH=accessibility;GCG=excluded",
        },
    ]
    for context in ("WCG", "GCH"):
        summary_rows.extend(
            [
                {
                    "metric": f"{context}_observed_sites",
                    "value": site_counts[context],
                },
                {
                    "metric": f"{context}_sum_covered_cells_per_site",
                    "value": covered_cells[context],
                },
                {
                    "metric": f"{context}_covered_calls",
                    "value": covered[context],
                },
                {
                    "metric": f"{context}_converted_calls",
                    "value": converted[context],
                },
                {
                    "metric": f"{context}_signal_fraction",
                    "value": converted[context] / covered[context]
                    if covered[context]
                    else "NA",
                },
            ]
        )
    atomic_summary(
        qc_root
        / f"dna_single_cytosine.track_shard_{args.shard:02d}.summary.tsv",
        summary_rows,
    )
    print(
        f"[OK] track shard={args.shard}: "
        f"WCG={covered['WCG']} GCH={covered['GCH']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
