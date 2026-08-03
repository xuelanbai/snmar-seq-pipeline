#!/usr/bin/env python3
"""Step 81: build species-separated per-base tables and BigWig tracks.

Inputs are the eight validated pseudo-bulk single-cytosine track shards made
from CB-aware, PCR-deduplicated DNA templates.  BigWig is necessarily
pseudo-bulk: the complete HDF5 matrix retains the single-cell dimension.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import gzip
import heapq
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np


EXPECTED_SHARDS = 8
DEDUP_SIGNATURE = "snmar_dna_cb_markdup_corrected_20260731"
REFERENCE_CHUNK_BASES = 4 * 1024 * 1024
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
CONTEXT_SIGNAL = {
    "WCG": "endogenous_methylation",
    "GCH": "GpC_accessibility",
}
TRACK_LABEL = {
    "WCG": "methylation",
    "GCH": "accessibility",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def read_metrics(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    metrics = {row["metric"]: row["value"] for row in rows}
    if len(metrics) != len(rows):
        raise ValueError(f"duplicate metric in {path}")
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


def count_contexts_in_block(sequence: bytes) -> Counter[str]:
    """Count strand-specific WCG/GCH/GCG cytosines in a sequence block."""
    counts: Counter[str] = Counter()
    if len(sequence) < 3:
        return counts
    bases = np.frombuffer(sequence, dtype=np.uint8)
    left = bases[:-2]
    center = bases[1:-1]
    right = bases[2:]
    a, c, g, t = map(ord, "ACGT")

    # Forward-strand cytosines use the reference C coordinate.
    counts["WCG"] += int(
        np.count_nonzero(
            (center == c) & (right == g) & ((left == a) | (left == t))
        )
    )
    counts["GCH"] += int(
        np.count_nonzero(
            (center == c)
            & (left == g)
            & ((right == a) | (right == c) | (right == t))
        )
    )
    counts["GCG"] += int(
        np.count_nonzero((center == c) & (left == g) & (right == g))
    )

    # Reverse-strand cytosines use the complementary reference G coordinate.
    counts["WCG"] += int(
        np.count_nonzero(
            (center == g) & (left == c) & ((right == a) | (right == t))
        )
    )
    counts["GCH"] += int(
        np.count_nonzero(
            (center == g)
            & (right == c)
            & ((left == a) | (left == g) | (left == t))
        )
    )
    counts["GCG"] += int(
        np.count_nonzero((center == g) & (left == c) & (right == c))
    )
    return counts


def count_reference_contexts(
    reference_fasta: Path,
    contigs: dict[int, dict[str, str]],
) -> Counter[tuple[str, str]]:
    """Count assayable reference contexts without loading a genome at once."""
    if not reference_fasta.is_file() or reference_fasta.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty reference FASTA: {reference_fasta}")
    wanted = {row["contig"]: row for row in contigs.values()}
    seen: set[str] = set()
    counts: Counter[tuple[str, str]] = Counter()
    current_name: str | None = None
    current_row: dict[str, str] | None = None
    current_length = 0
    buffer = bytearray()

    def consume(final: bool = False) -> None:
        nonlocal buffer
        if current_row is None:
            buffer.clear()
            return
        while len(buffer) >= REFERENCE_CHUNK_BASES + 2:
            block = bytes(buffer[: REFERENCE_CHUNK_BASES + 2])
            local = count_contexts_in_block(block)
            species = current_row["species"]
            for context, value in local.items():
                counts[(species, context)] += value
            del buffer[:REFERENCE_CHUNK_BASES]
        if final:
            local = count_contexts_in_block(bytes(buffer))
            species = current_row["species"]
            for context, value in local.items():
                counts[(species, context)] += value
            buffer.clear()

    def finish_contig() -> None:
        nonlocal current_name, current_row, current_length
        consume(final=True)
        if current_row is not None:
            expected_length = int(current_row["length"])
            if current_length != expected_length:
                raise ValueError(
                    "FASTA/contig-index length mismatch for "
                    f"{current_name}: {current_length} != {expected_length}"
                )
            if current_name in seen:
                raise ValueError(f"duplicate FASTA contig: {current_name}")
            seen.add(current_name)
        current_name = None
        current_row = None
        current_length = 0

    opener = gzip.open if reference_fasta.suffix == ".gz" else open
    with opener(reference_fasta, "rb") as handle:
        for raw_line in handle:
            if raw_line.startswith(b">"):
                finish_contig()
                current_name = raw_line[1:].split(maxsplit=1)[0].decode(
                    "utf-8"
                )
                current_row = wanted.get(current_name)
                continue
            if current_name is None:
                if raw_line.strip():
                    raise ValueError("sequence encountered before FASTA header")
                continue
            if current_row is None:
                continue
            sequence = raw_line.strip().upper()
            current_length += len(sequence)
            buffer.extend(sequence)
            consume()
    finish_contig()
    missing = sorted(set(wanted) - seen)
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(
            f"{len(missing)} indexed contigs missing from FASTA: {preview}"
        )
    return counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--track-root", type=Path)
    parser.add_argument("--input-qc-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--reference-fasta", type=Path)
    parser.add_argument(
        "--bedgraph-to-bigwig", default="bedGraphToBigWig"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    track_root = args.track_root or (
        project
        / "07_postmapping/complete_matrices/dna_single_cytosine/tracks"
    )
    input_qc_root = args.input_qc_root or (
        project / "08_qc_summary/22_dna_single_cytosine_tracks"
    )
    output_root = args.output_root or (
        project
        / "07_postmapping/complete_matrices/dna_single_cytosine/bigwig"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/24_dna_single_cytosine_bigwig"
    )
    reference_fasta = args.reference_fasta or (
        project
        / "02_reference/rna_human_pig_joint/human_pig.prefixed.fa"
    )
    dedup_summary = (
        project
        / "08_qc_summary/15_dna_cb_markdup_full/"
        "dna_cb_markdup_full_summary.tsv"
    )
    dedup = read_metrics(dedup_summary)
    if (
        dedup.get("status") != "validated"
        or dedup.get("pipeline_signature") != DEDUP_SIGNATURE
    ):
        raise ValueError(
            "BigWig requires the corrected, validated DNA dedup result"
        )

    contig_index_path = (
        project
        / "08_qc_summary/20_single_cytosine_inputs/"
        "dna_single_cytosine_contig_index.tsv"
    )
    contig_rows = read_tsv(contig_index_path)
    contigs: dict[int, dict[str, str]] = {}
    for row in contig_rows:
        index = int(row["contig_index"])
        if (
            index in contigs
            or row["species"] not in {"human", "pig"}
            or int(row["length"]) < 1
        ):
            raise ValueError(f"invalid contig-index row: {row}")
        contigs[index] = row
    if sorted(contigs) != list(range(len(contigs))):
        raise ValueError("contig indices must be contiguous from zero")
    reference_contexts = count_reference_contexts(reference_fasta, contigs)

    input_paths: list[Path] = []
    for shard in range(1, EXPECTED_SHARDS + 1):
        summary_path = (
            input_qc_root
            / f"dna_single_cytosine.track_shard_{shard:02d}.summary.tsv"
        )
        metrics = read_metrics(summary_path)
        if (
            metrics.get("status") != "validated"
            or int(metrics.get("shard", "0")) != shard
        ):
            raise ValueError(f"invalid track shard summary: {summary_path}")
        # Prevent old track shards from being mixed with a corrected dedup run.
        if summary_path.stat().st_mtime <= dedup_summary.stat().st_mtime:
            raise ValueError(
                f"track shard predates corrected dedup; rerun steps 53-75: "
                f"{summary_path}"
            )
        detailed = Path(metrics["detailed_track"])
        expected = (
            track_root
            / f"dna_single_cytosine.pseudobulk.shard_{shard:02d}.tsv.gz"
        )
        if detailed != expected:
            raise ValueError(
                f"unexpected detailed-track path: {detailed} != {expected}"
            )
        if not detailed.is_file() or detailed.stat().st_size == 0:
            raise FileNotFoundError(detailed)
        input_paths.append(detailed)

    converter = shutil.which(args.bedgraph_to_bigwig)
    if converter is None:
        raise FileNotFoundError(
            "bedGraphToBigWig is unavailable. Install "
            "bioconda::ucsc-bedgraphtobigwig in snTAPS_align."
        )
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)

    chrom_sizes: dict[str, Path] = {}
    for species in ("human", "pig"):
        path = output_root / f"{species}.chrom.sizes"
        descriptor, temporary = tempfile.mkstemp(
            prefix=path.name + ".", suffix=".tmp", dir=output_root
        )
        try:
            with os.fdopen(
                descriptor, "w", encoding="utf-8", newline=""
            ) as handle:
                for index in sorted(contigs):
                    row = contigs[index]
                    if row["species"] == species:
                        handle.write(f"{row['contig']}\t{row['length']}\n")
            os.replace(temporary, path)
        except BaseException:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise
        chrom_sizes[species] = path

    keys = [
        (species, context)
        for species in ("human", "pig")
        for context in ("WCG", "GCH")
    ]
    per_base_paths = {
        key: output_root
        / (
            f"{key[0]}.{key[1]}_{TRACK_LABEL[key[1]]}."
            "per_base.tsv.gz"
        )
        for key in keys
    }
    ratio_paths = {
        key: output_root
        / f"{key[0]}.{key[1]}_{TRACK_LABEL[key[1]]}_fraction.bw"
        for key in keys
    }
    coverage_paths = {
        key: output_root
        / f"{key[0]}.{key[1]}_{TRACK_LABEL[key[1]]}_coverage.bw"
        for key in keys
    }
    covered_cells_paths = {
        key: output_root
        / f"{key[0]}.{key[1]}_{TRACK_LABEL[key[1]]}_covered_cells.bw"
        for key in keys
    }
    temporary: dict[tuple[str, str, str], str] = {}
    for species, context in keys:
        for kind, final in (
            ("per_base", per_base_paths[(species, context)]),
            ("ratio_bedgraph", ratio_paths[(species, context)]),
            ("coverage_bedgraph", coverage_paths[(species, context)]),
            (
                "covered_cells_bedgraph",
                covered_cells_paths[(species, context)],
            ),
        ):
            descriptor, path = tempfile.mkstemp(
                prefix=final.name + ".", suffix=".tmp", dir=output_root
            )
            os.close(descriptor)
            temporary[(species, context, kind)] = path

    site_counts: Counter[tuple[str, str]] = Counter()
    covered_cell_entries: Counter[tuple[str, str]] = Counter()
    covered_calls: Counter[tuple[str, str]] = Counter()
    converted_calls: Counter[tuple[str, str]] = Counter()
    temporary_bigwigs: list[Path] = []
    try:
        with contextlib.ExitStack() as stack:
            readers: list[csv.DictReader] = []
            for path in input_paths:
                handle = stack.enter_context(
                    gzip.open(path, "rt", encoding="utf-8", newline="")
                )
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != FIELDS:
                    raise ValueError(f"unexpected fields in {path}")
                readers.append(reader)

            per_base_handles = {}
            ratio_handles = {}
            coverage_handles = {}
            covered_cells_handles = {}
            for key in keys:
                per_base = stack.enter_context(
                    gzip.open(
                        temporary[(*key, "per_base")],
                        "wt",
                        encoding="utf-8",
                        newline="",
                        compresslevel=1,
                    )
                )
                per_base.write(
                    "chrom\tstart\tend\tstrand\tcovered_cells\t"
                    "covered_calls\tconverted_calls\tsignal_fraction\n"
                )
                per_base_handles[key] = per_base
                ratio_handles[key] = stack.enter_context(
                    open(
                        temporary[(*key, "ratio_bedgraph")],
                        "wt",
                        encoding="utf-8",
                        newline="",
                    )
                )
                coverage_handles[key] = stack.enter_context(
                    open(
                        temporary[(*key, "coverage_bedgraph")],
                        "wt",
                        encoding="utf-8",
                        newline="",
                    )
                )
                covered_cells_handles[key] = stack.enter_context(
                    open(
                        temporary[(*key, "covered_cells_bedgraph")],
                        "wt",
                        encoding="utf-8",
                        newline="",
                    )
                )

            heap: list[
                tuple[int, int, int, int, int, dict[str, str]]
            ] = []
            previous: list[tuple[int, int, int, int] | None] = [
                None
            ] * len(readers)

            def push_next(source: int) -> None:
                try:
                    row = next(readers[source])
                except StopIteration:
                    return
                context = row["context"]
                context_order = {"WCG": 1, "GCH": 2}.get(context)
                strand = 1 if row["strand"] == "+" else -1
                key = (
                    int(row["contig_index"]),
                    int(row["start"]),
                    strand,
                    context_order or 0,
                )
                if (
                    context_order is None
                    or previous[source] is not None
                    and key <= previous[source]
                ):
                    raise ValueError("input track is not strictly sorted")
                previous[source] = key
                heapq.heappush(heap, (*key, source, row))

            for source in range(len(readers)):
                push_next(source)
            last_output_key: tuple[int, int] | None = None
            while heap:
                contig_index, start, strand, _order, source, row = (
                    heapq.heappop(heap)
                )
                push_next(source)
                metadata = contigs.get(contig_index)
                context = row["context"]
                species = row["species"]
                end = int(row["end"])
                cells = int(row["covered_cells"])
                covered = int(row["covered_calls"])
                converted = int(row["converted_calls"])
                fraction = float(row["signal_fraction"])
                output_key = (contig_index, start)
                if (
                    metadata is None
                    or species != metadata["species"]
                    or row["chrom"] != metadata["contig"]
                    or row["signal_type"] != CONTEXT_SIGNAL[context]
                    or end != start + 1
                    or start < 0
                    or end > int(metadata["length"])
                    or cells < 1
                    or cells > covered
                    or covered < 1
                    or converted < 0
                    or converted > covered
                    or not 0.0 <= fraction <= 1.0
                    or last_output_key is not None
                    and output_key <= last_output_key
                ):
                    raise ValueError(f"invalid/overlapping per-base row: {row}")
                last_output_key = output_key
                key = (species, context)
                per_base_handles[key].write(
                    f"{row['chrom']}\t{start}\t{end}\t{row['strand']}\t"
                    f"{row['covered_cells']}\t{covered}\t{converted}\t"
                    f"{fraction:.8g}\n"
                )
                ratio_handles[key].write(
                    f"{row['chrom']}\t{start}\t{end}\t{fraction:.8g}\n"
                )
                coverage_handles[key].write(
                    f"{row['chrom']}\t{start}\t{end}\t{covered}\n"
                )
                covered_cells_handles[key].write(
                    f"{row['chrom']}\t{start}\t{end}\t{cells}\n"
                )
                site_counts[key] += 1
                covered_cell_entries[key] += cells
                covered_calls[key] += covered
                converted_calls[key] += converted

        for key in keys:
            species, context = key
            if site_counts[key] == 0:
                raise ValueError(f"no per-base observations for {key}")
            for kind, final in (
                ("ratio_bedgraph", ratio_paths[key]),
                ("coverage_bedgraph", coverage_paths[key]),
                ("covered_cells_bedgraph", covered_cells_paths[key]),
            ):
                temporary_bw = final.with_name(
                    final.name + f".tmp.{os.getpid()}"
                )
                temporary_bigwigs.append(temporary_bw)
                subprocess.run(
                    [
                        converter,
                        temporary[(*key, kind)],
                        str(chrom_sizes[species]),
                        str(temporary_bw),
                    ],
                    check=True,
                )
                if (
                    not temporary_bw.is_file()
                    or temporary_bw.stat().st_size == 0
                ):
                    raise RuntimeError(f"empty BigWig: {temporary_bw}")
                os.replace(temporary_bw, final)
            os.replace(
                temporary[(*key, "per_base")], per_base_paths[key]
            )

        summary_rows: list[tuple[str, object]] = [
            ("status", "validated"),
            ("pipeline_signature", "snmar_single_base_bigwig_20260731"),
            ("source_dedup_pipeline_signature", DEDUP_SIGNATURE),
            ("aggregation", "pseudobulk_all_CBs_by_reference_species"),
            ("input_track_shards", EXPECTED_SHARDS),
            ("coordinate_system", "zero_based_half_open"),
            ("fraction_range", "0_to_1"),
            ("WCG_semantics", "endogenous_methylation"),
            ("GCH_semantics", "GpC_accessibility"),
            ("GCG_policy", "excluded_as_ambiguous"),
            (
                "reference_site_definition",
                "both_strands_WCG_or_GCH_cytosines_in_reference",
            ),
            ("reference_fasta", reference_fasta),
            ("bedGraphToBigWig", converter),
        ]
        for species in ("human", "pig"):
            summary_rows.append(
                (
                    f"{species}_GCG_ambiguous_reference_sites",
                    reference_contexts[(species, "GCG")],
                )
            )
        for species, context in keys:
            key = (species, context)
            prefix = f"{species}_{context}"
            reference_sites = reference_contexts[key]
            if reference_sites < 1:
                raise ValueError(f"no assayable reference sites for {key}")
            summary_rows.extend(
                [
                    (f"{prefix}_reference_assayable_sites", reference_sites),
                    (f"{prefix}_observed_sites", site_counts[key]),
                    (
                        f"{prefix}_observed_site_coverage_fraction",
                        site_counts[key] / reference_sites,
                    ),
                    (
                        f"{prefix}_covered_cell_site_entries",
                        covered_cell_entries[key],
                    ),
                    (f"{prefix}_covered_calls", covered_calls[key]),
                    (f"{prefix}_converted_calls", converted_calls[key]),
                    (
                        f"{prefix}_global_fraction",
                        converted_calls[key] / covered_calls[key],
                    ),
                    (f"{prefix}_per_base_table", per_base_paths[key]),
                    (f"{prefix}_fraction_bigwig", ratio_paths[key]),
                    (f"{prefix}_coverage_bigwig", coverage_paths[key]),
                    (
                        f"{prefix}_covered_cells_bigwig",
                        covered_cells_paths[key],
                    ),
                ]
            )
        atomic_metrics(
            qc_root / "dna_single_base_bigwig_summary.tsv",
            summary_rows,
        )
    finally:
        for path in temporary.values():
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
        for path in temporary_bigwigs:
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    print(
        "[OK] species-separated WCG/GCH per-base tables and BigWigs: "
        f"{output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
