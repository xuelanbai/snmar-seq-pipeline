#!/usr/bin/env python3
"""Step 85: build IGV BigWigs for separate human and pig references.

Step 81 deliberately preserves the ``human_``/``pig_`` chromosome prefixes
used by the joint mapping reference.  Those BigWigs are correct for the joint
reference, but appear empty when IGV is displaying an ordinary single-species
reference.  This display-only step rebuilds the same three signals from the
validated per-base tables while removing exactly one species prefix.

No calls, coordinates, coverage values, or fractions are recalculated.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


SPECIES = ("human", "pig")
CONTEXT_LABEL = {
    "WCG": "methylation",
    "GCH": "accessibility",
}
PER_BASE_FIELDS = [
    "chrom",
    "start",
    "end",
    "strand",
    "covered_cells",
    "covered_calls",
    "converted_calls",
    "signal_fraction",
]
SIGNALS = (
    ("fraction", "signal_fraction"),
    ("coverage", "covered_calls"),
    ("covered_cells", "covered_cells"),
)
SOURCE_SIGNATURE = "snmar_single_base_bigwig_20260731"
PIPELINE_SIGNATURE = "snmar_igv_species_bigwig_20260804"
DEFAULT_REFERENCE = {
    "human": (
        "/share/org/YZWL/yzwl_wujx/Projects/FOODIE_project/reference/"
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    ),
    "pig": (
        "/share/org/YZWL/yzwl_wujx/Projects/FOODIE_pig/reference/"
        "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    ),
}


def read_metrics(path: Path) -> dict[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty summary: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    metrics = {row["metric"]: row["value"] for row in rows}
    if not rows or len(metrics) != len(rows):
        raise ValueError(f"empty or duplicate metrics in {path}")
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


def load_chromosomes(
    path: Path, species: str
) -> tuple[list[tuple[str, str, int]], dict[str, tuple[int, str, int]]]:
    """Return source/IGV names in reference order and a lookup table."""
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty chromosome sizes: {path}")
    prefix = species + "_"
    chromosomes: list[tuple[str, str, int]] = []
    lookup: dict[str, tuple[int, str, int]] = {}
    target_names: set[str] = set()
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for order, line in enumerate(handle):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"invalid chromosome-size row: {line!r}")
            source, raw_length = fields
            if not source.startswith(prefix):
                raise ValueError(
                    f"{source!r} does not start with {prefix!r} in {path}"
                )
            target = source[len(prefix) :]
            length = int(raw_length)
            if not target or length < 1:
                raise ValueError(f"invalid chromosome in {path}: {line!r}")
            if source in lookup or target in target_names:
                raise ValueError(
                    f"chromosome name collision after prefix removal: {target}"
                )
            chromosomes.append((source, target, length))
            lookup[source] = (order, target, length)
            target_names.add(target)
    if not chromosomes:
        raise ValueError(f"no chromosomes in {path}")
    return chromosomes, lookup


def atomic_chrom_sizes(
    path: Path, chromosomes: list[tuple[str, str, int]]
) -> None:
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as handle:
            for _source, target, length in chromosomes:
                handle.write(f"{target}\t{length}\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def validate_species_reference(
    reference: Path,
    chromosomes: list[tuple[str, str, int]],
) -> Path:
    """Require exact target names/order/lengths in the single-species FAI."""
    if not reference.is_file() or reference.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty reference FASTA: {reference}")
    fai = Path(str(reference) + ".fai")
    if not fai.is_file() or fai.stat().st_size == 0:
        raise FileNotFoundError(
            f"missing/empty reference FAI: {fai}; run samtools faidx"
        )
    observed: list[tuple[str, int]] = []
    seen: set[str] = set()
    with fai.open("rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"invalid FAI row {line_number}: {line!r}")
            name = fields[0]
            length = int(fields[1])
            if not name or length < 1 or name in seen:
                raise ValueError(
                    f"invalid/duplicate FAI contig at row {line_number}: "
                    f"{name!r}"
                )
            seen.add(name)
            observed.append((name, length))
    expected = [(target, length) for _source, target, length in chromosomes]
    if observed != expected:
        mismatch = next(
            (
                index
                for index, (left, right) in enumerate(
                    zip(observed, expected), 1
                )
                if left != right
            ),
            min(len(observed), len(expected)) + 1,
        )
        observed_value = (
            observed[mismatch - 1] if mismatch <= len(observed) else "EOF"
        )
        expected_value = (
            expected[mismatch - 1] if mismatch <= len(expected) else "EOF"
        )
        raise ValueError(
            "prefix-free chromosome sizes do not match the target FASTA FAI "
            f"at row {mismatch}: observed={observed_value!r}, "
            f"expected={expected_value!r}"
        )
    return fai


def write_igv_bedgraphs(
    input_path: Path,
    species: str,
    chromosome_lookup: dict[str, tuple[int, str, int]],
    bedgraph_paths: dict[str, Path],
) -> int:
    """Strip a species prefix and write the three unchanged signals."""
    handles = {
        name: path.open("wt", encoding="utf-8", newline="")
        for name, path in bedgraph_paths.items()
    }
    rows = 0
    previous: tuple[int, int, int] | None = None
    try:
        with gzip.open(
            input_path, "rt", encoding="utf-8", newline=""
        ) as source:
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames != PER_BASE_FIELDS:
                raise ValueError(
                    f"unexpected per-base fields in {input_path}: "
                    f"{reader.fieldnames}"
                )
            for row in reader:
                metadata = chromosome_lookup.get(row["chrom"])
                if metadata is None:
                    raise ValueError(
                        f"unexpected {species} chromosome: {row['chrom']}"
                    )
                order, target, length = metadata
                start = int(row["start"])
                end = int(row["end"])
                covered_cells = int(row["covered_cells"])
                covered_calls = int(row["covered_calls"])
                converted_calls = int(row["converted_calls"])
                fraction = float(row["signal_fraction"])
                current = (order, start, end)
                if (
                    start < 0
                    or end != start + 1
                    or end > length
                    or covered_cells < 1
                    or covered_cells > covered_calls
                    or converted_calls < 0
                    or converted_calls > covered_calls
                    or not 0.0 <= fraction <= 1.0
                    or previous is not None
                    and current <= previous
                ):
                    raise ValueError(f"invalid or unsorted row: {row}")
                previous = current
                values = {
                    "fraction": f"{fraction:.8g}",
                    "coverage": str(covered_calls),
                    "covered_cells": str(covered_cells),
                }
                for signal, value in values.items():
                    handles[signal].write(
                        f"{target}\t{start}\t{end}\t{value}\n"
                    )
                rows += 1
    finally:
        for handle in handles.values():
            handle.close()
    if rows < 1:
        raise ValueError(f"no observations in {input_path}")
    return rows


def write_igv_readme(path: Path) -> None:
    text = """IGV species-reference BigWigs

These files contain the same values as the joint-reference BigWigs from step
81.  Only chromosome names changed: human_* -> * and pig_* -> *.

Load the matching single-species FASTA/genome in IGV before loading a track:
- human: Homo_sapiens.GRCh38.dna.primary_assembly.fa (Ensembl names)
- pig: Sus_scrofa.Sscrofa11.1.dna.toplevel.fa (Ensembl names)

Do not use these prefix-free tracks with human_pig.prefixed.fa.  Use the
original step-81 BigWigs for that joint reference.  Single-base tracks are
sparse, so inspect a region of about 100 kb to 1 Mb before zooming further.

Signal meanings:
- fraction: converted_calls / covered_calls, range 0-1
- coverage: PCR-deduplicated covered template calls
- covered_cells: cells with at least one covered call at that position
- WCG: endogenous DNA methylation
- GCH: GpC-labelled chromatin accessibility
"""
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
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
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--source-summary", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument(
        "--human-reference",
        type=Path,
        default=Path(DEFAULT_REFERENCE["human"]),
    )
    parser.add_argument(
        "--pig-reference",
        type=Path,
        default=Path(DEFAULT_REFERENCE["pig"]),
    )
    parser.add_argument(
        "--species",
        choices=("all", *SPECIES),
        default="all",
        help="build both species by default, or rebuild one species only",
    )
    parser.add_argument(
        "--bedgraph-to-bigwig", default="bedGraphToBigWig"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    input_root = args.input_root or (
        project
        / "07_postmapping/complete_matrices/dna_single_cytosine/bigwig"
    )
    source_summary = args.source_summary or (
        project
        / "08_qc_summary/24_dna_single_cytosine_bigwig/"
        "dna_single_base_bigwig_summary.tsv"
    )
    output_root = args.output_root or input_root / "igv_species_reference"
    qc_root = args.qc_root or (
        project / "08_qc_summary/26_dna_igv_species_bigwig"
    )

    source_metrics = read_metrics(source_summary)
    if (
        source_metrics.get("status") != "validated"
        or source_metrics.get("pipeline_signature") != SOURCE_SIGNATURE
    ):
        raise ValueError(
            "step 85 requires the validated step-81 BigWig summary"
        )
    converter = shutil.which(args.bedgraph_to_bigwig)
    if converter is None:
        raise FileNotFoundError(
            "bedGraphToBigWig is unavailable; use foodie_env or set "
            "--bedgraph-to-bigwig"
        )

    selected = SPECIES if args.species == "all" else (args.species,)
    references = {
        "human": args.human_reference,
        "pig": args.pig_reference,
    }
    output_root.mkdir(parents=True, exist_ok=True)
    qc_root.mkdir(parents=True, exist_ok=True)
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("source_pipeline_signature", SOURCE_SIGNATURE),
        ("source_summary", source_summary),
        ("source_bigwig_root", input_root),
        ("output_root", output_root),
        ("selected_species", ",".join(selected)),
        ("coordinate_system", "zero_based_half_open"),
        ("coordinate_transform", "remove_exact_species_prefix_only"),
        ("signal_values_recalculated", 0),
        ("bedGraphToBigWig", converter),
    ]
    temporary_paths: list[Path] = []
    try:
        for species in selected:
            source_sizes = input_root / f"{species}.chrom.sizes"
            chromosomes, lookup = load_chromosomes(source_sizes, species)
            target_sizes = (
                output_root / f"{species}.species_reference.chrom.sizes"
            )
            reference = references[species]
            reference_fai = validate_species_reference(
                reference, chromosomes
            )
            atomic_chrom_sizes(target_sizes, chromosomes)
            summary_rows.extend(
                [
                    (f"{species}_reference_fasta", reference),
                    (f"{species}_reference_fai", reference_fai),
                    (
                        f"{species}_reference_compatibility",
                        "exact_contig_order_names_and_lengths",
                    ),
                    (f"{species}_chromosomes", len(chromosomes)),
                    (f"{species}_chrom_sizes", target_sizes),
                ]
            )

            for context, label in CONTEXT_LABEL.items():
                input_path = (
                    input_root
                    / f"{species}.{context}_{label}.per_base.tsv.gz"
                )
                if not input_path.is_file() or input_path.stat().st_size == 0:
                    raise FileNotFoundError(input_path)
                bedgraphs: dict[str, Path] = {}
                for signal, _field in SIGNALS:
                    descriptor, temporary = tempfile.mkstemp(
                        prefix=f"{species}.{context}.{signal}.",
                        suffix=".bedGraph.tmp",
                        dir=output_root,
                    )
                    os.close(descriptor)
                    path = Path(temporary)
                    temporary_paths.append(path)
                    bedgraphs[signal] = path
                observations = write_igv_bedgraphs(
                    input_path, species, lookup, bedgraphs
                )
                summary_rows.append(
                    (f"{species}_{context}_observed_sites", observations)
                )
                for signal, _field in SIGNALS:
                    final = (
                        output_root
                        / f"{species}.{context}_{label}_{signal}.igv.bw"
                    )
                    temporary_bw = final.with_name(
                        final.name + f".tmp.{os.getpid()}"
                    )
                    temporary_paths.append(temporary_bw)
                    subprocess.run(
                        [
                            converter,
                            str(bedgraphs[signal]),
                            str(target_sizes),
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
                    summary_rows.append(
                        (f"{species}_{context}_{signal}_igv_bigwig", final)
                    )
                    bedgraphs[signal].unlink()

        readme = output_root / "README_IGV.txt"
        write_igv_readme(readme)
        summary_rows.append(("igv_readme", readme))
        atomic_metrics(
            qc_root / "dna_igv_species_bigwig_summary.tsv",
            summary_rows,
        )
    finally:
        for path in temporary_paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    print(f"[OK] IGV species-reference BigWigs: {output_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
