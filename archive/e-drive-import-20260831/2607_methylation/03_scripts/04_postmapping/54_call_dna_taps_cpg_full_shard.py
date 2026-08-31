#!/usr/bin/env python3
"""Step 54: call CB-aware TAPS CpGs for one sample/contig shard."""

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


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not operations or "".join(f"{n}{op}" for n, op in operations) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar}")
    return operations


def aligned_bases(
    pos1: int, cigar: str, sequence: str, quality: str
) -> Iterator[tuple[int, str, int]]:
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


def taps_call(
    is_reverse: bool,
    reference: str,
    reference_pos0: int,
    read_base: str,
) -> tuple[int, int] | None:
    if reference_pos0 < 0 or reference_pos0 >= len(reference):
        return None
    if not is_reverse:
        if (
            reference[reference_pos0] == "C"
            and reference_pos0 + 1 < len(reference)
            and reference[reference_pos0 + 1] == "G"
        ):
            if read_base == "T":
                return reference_pos0, 1
            if read_base == "C":
                return reference_pos0, 0
    elif (
        reference[reference_pos0] == "G"
        and reference_pos0 > 0
        and reference[reference_pos0 - 1] == "C"
    ):
        if read_base == "A":
            return reference_pos0 - 1, 1
        if read_base == "G":
            return reference_pos0 - 1, 0
    return None


def update_template_call(
    pending: dict[tuple[str, str, int], tuple[int, int | None]],
    key: tuple[str, str, int],
    base_quality: int,
    modified: int,
) -> str:
    previous = pending.get(key)
    if previous is None:
        pending[key] = (base_quality, modified)
        return "new"
    old_quality, old_modified = previous
    if base_quality > old_quality:
        pending[key] = (base_quality, modified)
        return "replaced"
    if base_quality < old_quality:
        return "ignored_lower_quality"
    if old_modified == modified:
        return "same"
    pending[key] = (base_quality, None)
    return "conflict"


def flush_pending_before(
    limit: int | None,
    heap: list[tuple[int, int, tuple[str, str, int]]],
    pending: dict[tuple[str, str, int], tuple[int, int | None]],
    consume: Callable[[str, int, int], None],
) -> int:
    """Flush sites that no future coordinate-sorted read can cover."""
    conflicts = 0
    while heap and (limit is None or heap[0][0] < limit):
        site, _serial, key = heapq.heappop(heap)
        value = pending.pop(key, None)
        if value is None:
            continue
        _quality, modified = value
        if modified is None:
            conflicts += 1
        else:
            consume(key[0], site, modified)
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


def species_of(contig: str) -> str:
    if contig.startswith("human_"):
        return "human"
    if contig.startswith("pig_"):
        return "pig"
    raise ValueError(f"unprefixed contig: {contig}")


def atomic_tsv(
    path: Path, fields: list[str], rows: list[dict[str, object]]
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
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
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
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--contig-file", type=Path, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--bin-size", type=int, default=100_000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.shard < 1
        or args.threads < 1
        or args.min_mapq < 0
        or args.min_baseq < 0
        or args.bin_size < 1
    ):
        raise ValueError("invalid full-calling setting")
    fai = Path(str(args.fasta) + ".fai")
    for path in (args.bam, args.fasta, fai, args.contig_file):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    contigs = read_contigs(args.contig_file)
    lengths = read_fai(fai)
    if any(contig not in lengths for contig in contigs):
        raise ValueError("contig list contains a name absent from FAI")

    cell_reads: Counter[str] = Counter()
    cell_cov: Counter[str] = Counter()
    cell_mod: Counter[str] = Counter()
    bin_cov: Counter[tuple[str, str, str, int]] = Counter()
    bin_mod: Counter[tuple[str, str, str, int]] = Counter()
    species_cov: Counter[str] = Counter()
    species_mod: Counter[str] = Counter()
    counters: Counter[str] = Counter()

    for contig in contigs:
        reference = fetch_contig(args.samtools, args.fasta, contig)
        if len(reference) != lengths[contig]:
            raise ValueError(f"FAI/reference length mismatch: {contig}")
        species = species_of(contig)
        pending: dict[tuple[str, str, int], tuple[int, int | None]] = {}
        pending_heap: list[tuple[int, int, tuple[str, str, int]]] = []
        serial = 0
        last_pos0 = -1

        def consume(cell_id: str, site: int, modified: int) -> None:
            bin_start = (site // args.bin_size) * args.bin_size
            cell_cov[cell_id] += 1
            cell_mod[cell_id] += modified
            key = (cell_id, species, contig, bin_start)
            bin_cov[key] += 1
            bin_mod[key] += modified
            species_cov[species] += 1
            species_mod[species] += modified

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
                raise ValueError(f"unexpected contig from samtools: {rname}")
            pos0 = int(pos_text) - 1
            if pos0 < last_pos0:
                raise ValueError(f"non-coordinate-sorted BAM at {contig}")
            last_pos0 = pos0
            # A reverse read starting at position p can call CpG C at p-1.
            counters["conflicting_template_cpg_calls_discarded"] += (
                flush_pending_before(
                    pos0 - 1, pending_heap, pending, consume
                )
            )
            if sequence == "*" or quality == "*":
                counters["missing_sequence_or_quality"] += 1
                continue
            tags = parse_tags(fields[11:])
            cb = tags.get("CB")
            xd = tags.get("XD")
            if not cb or xd is None:
                counters["missing_required_CB_or_XD"] += 1
                continue
            if xd != "0":
                counters["nonzero_XD"] += 1
                continue
            cell_id = f"{args.sample}|{cb.split()[0]}"
            cell_reads[cell_id] += 1
            counters["tag_eligible_alignment_records"] += 1
            is_reverse = bool(int(flag_text) & 16)
            for reference_pos0, read_base, base_quality in aligned_bases(
                pos0 + 1, cigar, sequence, quality
            ):
                if base_quality < args.min_baseq:
                    counters["low_baseq_aligned_bases"] += 1
                    continue
                call = taps_call(
                    is_reverse, reference, reference_pos0, read_base
                )
                if call is None:
                    continue
                site, modified = call
                counters["raw_directional_cpg_observations"] += 1
                key = (cell_id, qname, site)
                outcome = update_template_call(
                    pending, key, base_quality, modified
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
        counters["conflicting_template_cpg_calls_discarded"] += (
            flush_pending_before(None, pending_heap, pending, consume)
        )
        if pending:
            raise RuntimeError("pending-call cache did not flush")
        counters["completed_contigs"] += 1
        counters["reference_bases"] += len(reference)

    if counters["missing_required_CB_or_XD"] or counters["nonzero_XD"]:
        raise ValueError("full caller encountered missing/non-strict tags")
    total_cov = sum(cell_cov.values())
    total_mod = sum(cell_mod.values())
    if total_cov == 0:
        raise ValueError("full shard produced zero CpG calls")

    bins_by_cell: Counter[str] = Counter()
    for cell_id, _species, _contig, _start in bin_cov:
        bins_by_cell[cell_id] += 1
    per_cell_rows = []
    for cell_id in sorted(cell_cov, key=lambda value: (-cell_cov[value], value)):
        covered = cell_cov[cell_id]
        modified = cell_mod[cell_id]
        per_cell_rows.append(
            {
                "sample": args.sample,
                "cell_id": cell_id,
                "CB": cell_id.split("|", 1)[1],
                "eligible_alignment_records": cell_reads[cell_id],
                "covered_template_cpg_calls": covered,
                "modified_template_cpg_calls": modified,
                "modified_rate": modified / covered,
                "bins_with_calls": bins_by_cell[cell_id],
            }
        )
    def iter_bin_rows() -> Iterator[dict[str, object]]:
        for key, covered in bin_cov.items():
            cell_id, species, contig, start = key
            modified = bin_mod[key]
            yield {
                "sample": args.sample,
                "cell_id": cell_id,
                "CB": cell_id.split("|", 1)[1],
                "species": species,
                "chrom": contig,
                "start": start,
                "end": min(start + args.bin_size, lengths[contig]),
                "modified_template_cpg_calls": modified,
                "covered_template_cpg_calls": covered,
                "modified_rate": modified / covered,
            }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{args.sample}.shard_{args.shard:02d}"
    atomic_tsv(
        args.out_dir / f"{prefix}.per_cell.tsv",
        list(per_cell_rows[0]),
        per_cell_rows,
    )
    atomic_gzip_tsv(
        args.out_dir / f"{prefix}.bin_counts.tsv.gz",
        [
            "sample",
            "cell_id",
            "CB",
            "species",
            "chrom",
            "start",
            "end",
            "modified_template_cpg_calls",
            "covered_template_cpg_calls",
            "modified_rate",
        ],
        iter_bin_rows(),
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "sample", "value": args.sample},
        {"metric": "shard", "value": args.shard},
        {"metric": "input_bam", "value": args.bam},
        {"metric": "input_bam_bytes", "value": args.bam.stat().st_size},
        {"metric": "input_bam_mtime", "value": int(args.bam.stat().st_mtime)},
        {"metric": "contig_file", "value": args.contig_file},
        {"metric": "expected_contigs", "value": len(contigs)},
        {"metric": "min_mapq", "value": args.min_mapq},
        {"metric": "min_baseq", "value": args.min_baseq},
        {"metric": "bin_size", "value": args.bin_size},
        {
            "metric": "call_model",
            "value": "forward_ref_C_CT_reverse_ref_G_GA",
        },
        {
            "metric": "overlap_policy",
            "value": "one_call_per_sample_CB_QNAME_CpG_highest_baseQ",
        },
        {"metric": "cells_with_cpg_calls", "value": len(cell_cov)},
        {"metric": "covered_template_cpg_calls", "value": total_cov},
        {"metric": "modified_template_cpg_calls", "value": total_mod},
        {"metric": "modified_rate", "value": total_mod / total_cov},
        {"metric": "human_covered_calls", "value": species_cov["human"]},
        {"metric": "human_modified_calls", "value": species_mod["human"]},
        {"metric": "pig_covered_calls", "value": species_cov["pig"]},
        {"metric": "pig_modified_calls", "value": species_mod["pig"]},
    ]
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
        f"contigs={len(contigs)} calls={total_cov}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
