#!/usr/bin/env python3
"""Step 65: call one sample/shard WCG-methylation and GCH-accessibility matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import heapq
import os
import re
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Callable, Iterable, Iterator


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
FILTER_FLAGS = 4 + 256 + 512 + 1024 + 2048
CONTEXT_GCG_AMBIGUOUS = 0
CONTEXT_WCG = 1
CONTEXT_GCH = 2
CONTEXT_NAMES = {
    CONTEXT_WCG: "WCG",
    CONTEXT_GCH: "GCH",
}
SIGNAL_TYPES = {
    CONTEXT_WCG: "endogenous_methylation",
    CONTEXT_GCH: "GpC_accessibility",
}
MATRIX_FIELDS = (
    "cell_index",
    "contig_index",
    "cytosine_position0",
    "strand",
    "context_code",
    "covered",
    "converted",
    "signal_fraction",
)


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not operations or "".join(f"{n}{op}" for n, op in operations) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar}")
    return operations


def aligned_bases(
    pos1: int, cigar: str, sequence: str, quality: str
) -> Iterator[tuple[int, int, str, int]]:
    """Yield reference position, query position, base and Phred quality."""
    query = 0
    reference = pos1 - 1
    for length, op in parse_cigar(cigar):
        if op in {"M", "=", "X"}:
            for offset in range(length):
                qpos = query + offset
                if qpos >= len(sequence) or qpos >= len(quality):
                    raise ValueError("CIGAR consumes beyond SEQ/QUAL")
                yield (
                    reference + offset,
                    qpos,
                    sequence[qpos].upper(),
                    ord(quality[qpos]) - 33,
                )
            query += length
            reference += length
        elif op in {"I", "S"}:
            query += length
        elif op in {"D", "N"}:
            reference += length
        elif op in {"H", "P"}:
            continue
    if query != len(sequence) or query != len(quality):
        raise ValueError("CIGAR does not consume complete SEQ/QUAL")


def context_call(
    is_reverse: bool,
    reference: str,
    reference_pos0: int,
    read_base: str,
) -> tuple[int, int, int, int] | None:
    """Return reference position, strand, context and converted state.

    A reverse-strand cytosine is represented at the coordinate of its
    complementary reference G and retains strand=-1. GCG is returned with
    context code 0 so it can be audited and excluded from both main signals.
    """
    if reference_pos0 <= 0 or reference_pos0 + 1 >= len(reference):
        return None
    if not is_reverse:
        if reference[reference_pos0] != "C" or read_base not in {"C", "T"}:
            return None
        previous = reference[reference_pos0 - 1]
        following = reference[reference_pos0 + 1]
        converted = int(read_base == "T")
        if previous in {"A", "T"} and following == "G":
            return reference_pos0, 1, CONTEXT_WCG, converted
        if previous == "G" and following in {"A", "C", "T"}:
            return reference_pos0, 1, CONTEXT_GCH, converted
        if previous == "G" and following == "G":
            return reference_pos0, 1, CONTEXT_GCG_AMBIGUOUS, converted
        return None
    if reference[reference_pos0] != "G" or read_base not in {"G", "A"}:
        return None
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    previous = complement.get(reference[reference_pos0 + 1], "N")
    following = complement.get(reference[reference_pos0 - 1], "N")
    converted = int(read_base == "A")
    if previous in {"A", "T"} and following == "G":
        return reference_pos0, -1, CONTEXT_WCG, converted
    if previous == "G" and following in {"A", "C", "T"}:
        return reference_pos0, -1, CONTEXT_GCH, converted
    if previous == "G" and following == "G":
        return reference_pos0, -1, CONTEXT_GCG_AMBIGUOUS, converted
    return None


def update_template_call(
    pending: dict[
        tuple[int, str, int, int, int], tuple[int, int | None]
    ],
    key: tuple[int, str, int, int, int],
    base_quality: int,
    converted: int,
) -> str:
    previous = pending.get(key)
    if previous is None:
        pending[key] = (base_quality, converted)
        return "new"
    old_quality, old_converted = previous
    if base_quality > old_quality:
        pending[key] = (base_quality, converted)
        return "replaced"
    if base_quality < old_quality:
        return "ignored_lower_quality"
    if old_converted == converted:
        return "same"
    pending[key] = (base_quality, None)
    return "conflict"


def flush_pending_before(
    limit: int | None,
    heap: list[
        tuple[int, int, tuple[int, str, int, int, int]]
    ],
    pending: dict[
        tuple[int, str, int, int, int], tuple[int, int | None]
    ],
    consume: Callable[[int, int, int, int, int], None],
) -> int:
    conflicts = 0
    while heap and (limit is None or heap[0][0] < limit):
        site, _serial, key = heapq.heappop(heap)
        value = pending.pop(key, None)
        if value is None:
            continue
        _quality, converted = value
        if converted is None:
            conflicts += 1
        else:
            consume(key[0], site, key[3], key[4], converted)
    return conflicts


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def run_text(command: list[str]) -> str:
    result = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stderr}"
        )
    return result.stdout


def fetch_contig(samtools: str, fasta: Path, contig: str) -> str:
    text = run_text([samtools, "faidx", str(fasta), contig])
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"reference fetch failed: {contig}")
    sequence = "".join(lines[1:]).upper()
    if not sequence:
        raise ValueError(f"empty reference sequence: {contig}")
    return sequence


def read_contigs(path: Path) -> list[str]:
    contigs = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not contigs or len(contigs) != len(set(contigs)):
        raise ValueError(f"empty/duplicate contig list: {path}")
    return contigs


def read_fai(path: Path) -> dict[str, int]:
    lengths = {}
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def read_cell_indices(path: Path, sample: str) -> dict[str, int]:
    indices: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample"] != sample:
                continue
            cb = row["CB"]
            if cb in indices:
                raise ValueError(f"duplicate CB in {sample}: {cb}")
            indices[cb] = int(row["cell_index"])
    if not indices:
        raise ValueError(f"no indexed cells for sample: {sample}")
    return indices


def read_contig_indices(path: Path) -> dict[str, dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {row["contig"]: row for row in rows}
    if not rows or len(result) != len(rows):
        raise ValueError(f"invalid contig index: {path}")
    return result


def atomic_tsv(
    path: Path, fields: list[str], rows: Iterable[dict[str, object]]
) -> None:
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


def atomic_gzip_tsv(
    path: Path, fields: list[str], rows: Iterable[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as handle:
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


def write_h5(
    path: Path,
    records: list[tuple[int, int, int, int, int, int, int]],
    attributes: dict[str, object],
) -> None:
    try:
        import h5py
        import numpy as np
    except ImportError as error:
        raise RuntimeError(
            "Step 65 requires numpy and h5py in snTAPS_align"
        ) from error
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        arrays = {
            "cell_index": np.fromiter(
                (row[4] for row in records), dtype=np.uint32
            ),
            "contig_index": np.fromiter(
                (row[0] for row in records), dtype=np.uint32
            ),
            "cytosine_position0": np.fromiter(
                (row[1] for row in records), dtype=np.uint32
            ),
            "strand": np.fromiter(
                (row[2] for row in records), dtype=np.int8
            ),
            "context_code": np.fromiter(
                (row[3] for row in records), dtype=np.uint8
            ),
            "covered": np.fromiter(
                (row[5] for row in records), dtype=np.uint32
            ),
            "converted": np.fromiter(
                (row[6] for row in records), dtype=np.uint32
            ),
        }
        covered = arrays["covered"]
        converted = arrays["converted"]
        arrays["signal_fraction"] = np.divide(
            converted,
            covered,
            out=np.zeros(len(records), dtype=np.float32),
            where=covered > 0,
        )
        with h5py.File(temporary, "w") as output:
            group = output.create_group("matrix")
            for field in MATRIX_FIELDS:
                values = arrays[field]
                options = {}
                if len(values):
                    options = {
                        "compression": "gzip",
                        "compression_opts": 1,
                        "shuffle": True,
                        "chunks": True,
                    }
                group.create_dataset(field, data=values, **options)
            for key, value in attributes.items():
                output.attrs[key] = value
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--contig-file", type=Path, required=True)
    parser.add_argument("--cell-index", type=Path, required=True)
    parser.add_argument("--contig-index", type=Path, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--exclude-read-ends", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.shard < 1
        or args.threads < 1
        or args.min_mapq < 0
        or args.min_baseq < 0
        or args.exclude_read_ends < 0
    ):
        raise ValueError("invalid single-cytosine setting")
    fai = Path(str(args.fasta) + ".fai")
    prerequisites = (
        args.bam,
        args.fasta,
        fai,
        args.contig_file,
        args.cell_index,
        args.contig_index,
    )
    for path in prerequisites:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    contigs = read_contigs(args.contig_file)
    lengths = read_fai(fai)
    index_rows = read_contig_indices(args.contig_index)
    cell_indices = read_cell_indices(args.cell_index, args.sample)
    if any(contig not in lengths or contig not in index_rows for contig in contigs):
        raise ValueError("contig absent from FAI or global index")
    matrix_cov: Counter[tuple[int, int, int, int, int]] = Counter()
    matrix_converted: Counter[
        tuple[int, int, int, int, int]
    ] = Counter()
    context_cov: Counter[int] = Counter()
    context_converted: Counter[int] = Counter()
    species_context_cov: Counter[tuple[str, int]] = Counter()
    species_context_converted: Counter[tuple[str, int]] = Counter()
    counters: Counter[str] = Counter()

    for contig in contigs:
        reference = fetch_contig(args.samtools, args.fasta, contig)
        if len(reference) != lengths[contig]:
            raise ValueError(f"FAI/reference length mismatch: {contig}")
        metadata = index_rows[contig]
        contig_index = int(metadata["contig_index"])
        species = metadata["species"]
        pending: dict[
            tuple[int, str, int, int, int], tuple[int, int | None]
        ] = {}
        pending_heap: list[
            tuple[int, int, tuple[int, str, int, int, int]]
        ] = []
        serial = 0
        last_pos0 = -1

        def consume(
            cell_index: int,
            site: int,
            strand: int,
            context_code: int,
            converted: int,
        ) -> None:
            if context_code == CONTEXT_GCG_AMBIGUOUS:
                counters["ambiguous_GCG_template_calls_excluded"] += 1
                counters[
                    "ambiguous_GCG_converted_template_calls_excluded"
                ] += converted
                return
            key = (contig_index, site, strand, context_code, cell_index)
            matrix_cov[key] += 1
            matrix_converted[key] += converted
            context_cov[context_code] += 1
            context_converted[context_code] += converted
            species_context_cov[(species, context_code)] += 1
            species_context_converted[(species, context_code)] += converted

        command = [
            args.samtools,
            "view",
            "-@",
            str(args.threads),
            "-q",
            str(args.min_mapq),
            "-F",
            str(FILTER_FLAGS),
            str(args.bam),
            contig,
        ]
        process = subprocess.Popen(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert process.stdout is not None
        for line in process.stdout:
            counters["alignment_records_after_flag_mapq_filter"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("malformed SAM record")
            qname, flag_text, rname, pos_text, _mapq, cigar = fields[:6]
            sequence, quality = fields[9], fields[10]
            if rname != contig:
                raise ValueError(f"unexpected contig: {rname}")
            pos0 = int(pos_text) - 1
            if pos0 < last_pos0:
                raise ValueError(f"non-coordinate-sorted BAM: {contig}")
            last_pos0 = pos0
            counters["conflicting_template_cytosine_calls_discarded"] += (
                flush_pending_before(pos0, pending_heap, pending, consume)
            )
            if sequence == "*" or quality == "*":
                counters["missing_sequence_or_quality"] += 1
                continue
            tags = parse_tags(fields[11:])
            cb = tags.get("CB", "").split()[0]
            xd = tags.get("XD")
            if not cb or xd is None:
                counters["missing_required_CB_or_XD"] += 1
                continue
            if xd != "0":
                counters["nonzero_XD"] += 1
                continue
            cell_index = cell_indices.get(cb)
            if cell_index is None:
                raise ValueError(f"strict BAM CB absent from cell index: {cb}")
            counters["tag_eligible_alignment_records"] += 1
            is_reverse = bool(int(flag_text) & 16)
            for reference_pos0, qpos, read_base, base_quality in aligned_bases(
                pos0 + 1, cigar, sequence, quality
            ):
                if (
                    qpos < args.exclude_read_ends
                    or qpos >= len(sequence) - args.exclude_read_ends
                ):
                    counters["read_end_aligned_bases_excluded"] += 1
                    continue
                if base_quality < args.min_baseq:
                    counters["low_baseq_aligned_bases"] += 1
                    continue
                call = context_call(
                    is_reverse, reference, reference_pos0, read_base
                )
                if call is None:
                    continue
                site, strand, context_code, converted = call
                context_name = (
                    "GCG_ambiguous"
                    if context_code == CONTEXT_GCG_AMBIGUOUS
                    else CONTEXT_NAMES[context_code]
                )
                counters[
                    f"raw_directional_{context_name}_observations"
                ] += 1
                key = (cell_index, qname, site, strand, context_code)
                outcome = update_template_call(
                    pending, key, base_quality, converted
                )
                if outcome == "new":
                    serial += 1
                    heapq.heappush(pending_heap, (site, serial, key))
                else:
                    counters["overlap_observations_collapsed"] += 1
                if outcome == "conflict":
                    counters["equal_quality_conflicts"] += 1
        assert process.stderr is not None
        stderr = process.stderr.read()
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(
                f"samtools view failed ({return_code}): {stderr}"
            )
        counters["conflicting_template_cytosine_calls_discarded"] += (
            flush_pending_before(None, pending_heap, pending, consume)
        )
        if pending:
            raise RuntimeError("pending-call cache did not flush")
        counters["completed_contigs"] += 1
        counters["reference_bases"] += len(reference)

    if counters["missing_required_CB_or_XD"] or counters["nonzero_XD"]:
        raise ValueError("caller encountered missing/non-strict tags")
    records = [
        (
            key[0],
            key[1],
            key[2],
            key[3],
            key[4],
            covered,
            matrix_converted[key],
        )
        for key, covered in matrix_cov.items()
    ]
    records.sort(key=lambda row: (row[0], row[1], row[2], row[3], row[4]))
    total_cov = sum(row[5] for row in records)
    total_converted = sum(row[6] for row in records)
    if total_converted > total_cov:
        raise ValueError("converted calls exceed covered calls")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{args.sample}.shard_{args.shard:02d}"
    matrix_path = args.out_dir / f"{prefix}.single_cytosine.h5"
    write_h5(
        matrix_path,
        records,
        {
            "status": "validated",
            "sample": args.sample,
            "shard": args.shard,
            "coordinate_system": (
                "0-based_reference_C_for_plus_reference_G_for_minus"
            ),
            "context_code_1": "WCG_endogenous_methylation",
            "context_code_2": "GCH_GpC_accessibility",
            "ambiguous_context_policy": "GCG_excluded_from_both_signals",
            "converted_semantics": (
                "forward_T_or_reverse_A_after_TAPS"
            ),
            "missing_semantics": "absent_coordinate_is_NA_not_zero",
            "minimum_mapq": args.min_mapq,
            "minimum_baseq": args.min_baseq,
            "excluded_bases_each_read_end": args.exclude_read_ends,
        },
    )

    def site_rows() -> Iterator[dict[str, object]]:
        current: tuple[int, int, int, int] | None = None
        covered_cells = covered_calls = converted_calls = 0
        for (
            contig_index,
            position0,
            strand,
            context_code,
            _cell,
            covered,
            converted,
        ) in records:
            key = (contig_index, position0, strand, context_code)
            if current is not None and key != current:
                row = index_rows_by_number[current[0]]
                yield {
                    "contig_index": current[0],
                    "species": row["species"],
                    "chrom": row["contig"],
                    "start": current[1],
                    "end": current[1] + 1,
                    "strand": "+" if current[2] == 1 else "-",
                    "context": CONTEXT_NAMES[current[3]],
                    "signal_type": SIGNAL_TYPES[current[3]],
                    "covered_cells": covered_cells,
                    "covered_calls": covered_calls,
                    "converted_calls": converted_calls,
                    "signal_fraction": converted_calls / covered_calls,
                }
                covered_cells = covered_calls = converted_calls = 0
            current = key
            covered_cells += 1
            covered_calls += covered
            converted_calls += converted
        if current is not None:
            row = index_rows_by_number[current[0]]
            yield {
                "contig_index": current[0],
                "species": row["species"],
                "chrom": row["contig"],
                "start": current[1],
                "end": current[1] + 1,
                "strand": "+" if current[2] == 1 else "-",
                "context": CONTEXT_NAMES[current[3]],
                "signal_type": SIGNAL_TYPES[current[3]],
                "covered_cells": covered_cells,
                "covered_calls": covered_calls,
                "converted_calls": converted_calls,
                "signal_fraction": converted_calls / covered_calls,
            }

    index_rows_by_number = {
        int(row["contig_index"]): row for row in index_rows.values()
    }
    site_path = args.out_dir / f"{prefix}.site_summary.tsv.gz"
    site_fields = [
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
    atomic_gzip_tsv(site_path, site_fields, site_rows())
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "sample", "value": args.sample},
        {"metric": "shard", "value": args.shard},
        {"metric": "input_bam", "value": args.bam},
        {"metric": "matrix_h5", "value": matrix_path},
        {"metric": "site_summary", "value": site_path},
        {"metric": "expected_contigs", "value": len(contigs)},
        {"metric": "observed_cell_cytosine_entries", "value": len(records)},
        {"metric": "covered_template_cytosine_calls", "value": total_cov},
        {"metric": "converted_template_cytosine_calls", "value": total_converted},
        {
            "metric": "combined_signal_fraction_for_QC_only",
            "value": total_converted / total_cov if total_cov else "NA",
        },
        {"metric": "WCG_covered_calls", "value": context_cov[CONTEXT_WCG]},
        {
            "metric": "WCG_converted_calls",
            "value": context_converted[CONTEXT_WCG],
        },
        {
            "metric": "WCG_methylation_fraction",
            "value": context_converted[CONTEXT_WCG]
            / context_cov[CONTEXT_WCG]
            if context_cov[CONTEXT_WCG]
            else "NA",
        },
        {"metric": "GCH_covered_calls", "value": context_cov[CONTEXT_GCH]},
        {
            "metric": "GCH_converted_calls",
            "value": context_converted[CONTEXT_GCH],
        },
        {
            "metric": "GCH_accessibility_fraction",
            "value": context_converted[CONTEXT_GCH]
            / context_cov[CONTEXT_GCH]
            if context_cov[CONTEXT_GCH]
            else "NA",
        },
        {"metric": "minimum_mapq", "value": args.min_mapq},
        {"metric": "minimum_baseq", "value": args.min_baseq},
        {
            "metric": "excluded_bases_each_read_end",
            "value": args.exclude_read_ends,
        },
        {
            "metric": "context_model",
            "value": (
                "WCG=endogenous_methylation;"
                "GCH=GpC_accessibility;GCG=excluded"
            ),
        },
        {
            "metric": "converted_semantics",
            "value": "forward_T_or_reverse_A_after_TAPS",
        },
        {
            "metric": "missing_cell_cytosine_semantics",
            "value": "NA_not_zero_signal",
        },
    ]
    for species in ("human", "pig"):
        for context_code in (CONTEXT_WCG, CONTEXT_GCH):
            context_name = CONTEXT_NAMES[context_code]
            summary_rows.extend(
                [
                    {
                        "metric": f"{species}_{context_name}_covered_calls",
                        "value": species_context_cov[
                            (species, context_code)
                        ],
                    },
                    {
                        "metric": f"{species}_{context_name}_converted_calls",
                        "value": species_context_converted[
                            (species, context_code)
                        ],
                    },
                ]
            )
    summary_rows.extend(
        {"metric": key, "value": value}
        for key, value in sorted(counters.items())
    )
    atomic_tsv(
        args.out_dir / f"{prefix}.summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] {args.sample} shard={args.shard}: "
        f"entries={len(records)} calls={total_cov}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
