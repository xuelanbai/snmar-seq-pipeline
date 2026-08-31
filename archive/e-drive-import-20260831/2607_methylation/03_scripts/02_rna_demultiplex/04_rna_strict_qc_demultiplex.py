#!/usr/bin/env python3
"""Step 04: probe or strictly demultiplex the current snMAR RNA structure."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


DNA_RE = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
GTCA_INSERT_MOTIF = "CTCCGACGC"
GTCA_INSERT = "GTCA"


@dataclass(frozen=True)
class T7Linker:
    name: str
    sequence: str
    linker_length: int
    linker_sequence: str
    dna_residual_linker_sequence: str


def load_barcodes(path: Path) -> dict[str, str]:
    seq_to_name: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_no}: expected NAME SEQUENCE")
            name, sequence = fields[0], fields[1].upper()
            if not DNA_RE.fullmatch(sequence):
                raise ValueError(f"{path}:{line_no}: invalid sequence {sequence!r}")
            if sequence in seq_to_name and seq_to_name[sequence] != name:
                raise ValueError(f"{path}:{line_no}: duplicate sequence {sequence}")
            seq_to_name[sequence] = name
    if not seq_to_name:
        raise ValueError(f"empty barcode table: {path}")
    if len({len(sequence) for sequence in seq_to_name}) != 1:
        raise ValueError(f"mixed barcode lengths in {path}")
    return seq_to_name


def make_unique_prefix_map(
    seq_to_name: dict[str, str],
    prefix_length: int,
) -> dict[str, tuple[str, str]]:
    """Map a whitelist-unique prefix for diagnostics, not formal assignment."""
    prefix_map: dict[str, tuple[str, str]] = {}
    for sequence, name in seq_to_name.items():
        prefix = sequence[:prefix_length]
        if prefix in prefix_map:
            other_name, other_sequence = prefix_map[prefix]
            raise ValueError(
                f"non-unique {prefix_length}-bp prefix {prefix}: "
                f"{other_name}/{other_sequence} and {name}/{sequence}"
            )
        prefix_map[prefix] = (name, sequence)
    return prefix_map


def make_dna_residual_linker(linker: str, path: Path, line_no: int) -> str:
    if linker.count(GTCA_INSERT_MOTIF) != 1:
        raise ValueError(
            f"{path}:{line_no}: linker must contain one "
            f"{GTCA_INSERT_MOTIF} motif"
        )
    return linker.replace(
        GTCA_INSERT_MOTIF,
        GTCA_INSERT_MOTIF + GTCA_INSERT,
        1,
    )


def load_t7_linkers(path: Path) -> dict[str, T7Linker]:
    result: dict[str, T7Linker] = {}
    names: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"name", "sequence", "linker_length", "linker_sequence"}
        if set(reader.fieldnames or ()) != required:
            raise ValueError(
                f"{path}: expected columns name, sequence, linker_length, linker_sequence"
            )
        for line_no, row in enumerate(reader, 2):
            name = row["name"]
            sequence = row["sequence"].upper()
            linker_sequence = row["linker_sequence"].upper()
            linker_length = int(row["linker_length"])
            if name in names:
                raise ValueError(f"{path}:{line_no}: duplicate T7 name {name}")
            if sequence in result:
                raise ValueError(f"{path}:{line_no}: duplicate T7 sequence {sequence}")
            if (
                len(sequence) != 8
                or not DNA_RE.fullmatch(sequence)
                or linker_length not in (22, 23, 24, 25)
                or len(linker_sequence) != linker_length
                or not DNA_RE.fullmatch(linker_sequence)
            ):
                raise ValueError(f"{path}:{line_no}: invalid T7/linker row")
            names.add(name)
            result[sequence] = T7Linker(
                name=name,
                sequence=sequence,
                linker_length=linker_length,
                linker_sequence=linker_sequence,
                dna_residual_linker_sequence=make_dna_residual_linker(
                    linker_sequence,
                    path,
                    line_no,
                ),
            )
    if not result:
        raise ValueError(f"empty T7-linker table: {path}")
    return result


def fastq_pairs(
    r1_path: Path,
    r2_path: Path,
) -> Iterator[tuple[str, str, str, str, str, str, str, str]]:
    with gzip.open(r1_path, "rt") as r1, gzip.open(r2_path, "rt") as r2:
        record_no = 0
        while True:
            h1 = r1.readline().rstrip()
            h2 = r2.readline().rstrip()
            if not h1 and not h2:
                break
            record_no += 1
            if not h1 or not h2:
                raise ValueError(f"R1/R2 record count differs near pair {record_no}")
            s1 = r1.readline().rstrip().upper()
            p1 = r1.readline().rstrip()
            q1 = r1.readline().rstrip()
            s2 = r2.readline().rstrip().upper()
            p2 = r2.readline().rstrip()
            q2 = r2.readline().rstrip()
            if not all((s1, p1, q1, s2, p2, q2)):
                raise ValueError(f"truncated FASTQ record near pair {record_no}")
            if len(s1) != len(q1) or len(s2) != len(q2):
                raise ValueError(
                    f"sequence/quality length mismatch near pair {record_no}"
                )
            read_id_1 = h1.split(maxsplit=1)[0].removeprefix("@").removesuffix("/1")
            read_id_2 = h2.split(maxsplit=1)[0].removeprefix("@").removesuffix("/2")
            if read_id_1 != read_id_2:
                raise ValueError(f"R1/R2 read IDs differ near pair {record_no}")
            yield h1, s1, p1, q1, h2, s2, p2, q2


def parse_header_indices(header: str) -> tuple[str, str] | None:
    for token in reversed(header.split()):
        candidate = token.rsplit(":", 1)[-1].upper()
        match = re.fullmatch(r"([ACGTN]+)\+([ACGTN]+)", candidate)
        if match:
            return match.group(1), match.group(2)
    match = re.search(r"([ACGTN]+)\+([ACGTN]+)\s*$", header.upper())
    if match:
        return match.group(1), match.group(2)
    return None


def trailing_run(sequence: str, base: str) -> int:
    count = 0
    for nucleotide in reversed(sequence):
        if nucleotide != base:
            break
        count += 1
    return count


def write_record(
    handle: TextIO,
    header: str,
    tags: str,
    sequence: str,
    plus: str,
    quality: str,
) -> None:
    read_id = header.split(maxsplit=1)[0]
    handle.write(f"{read_id} {tags}\n{sequence}\n{plus}\n{quality}\n")


def process(args: argparse.Namespace) -> None:
    t7_linkers = load_t7_linkers(args.t7_linker)
    t5 = load_barcodes(args.t5)
    i7 = load_barcodes(args.i7)
    i5 = load_barcodes(args.i5)
    i7_prefix8 = make_unique_prefix_map(i7, 8)
    t5_length = len(next(iter(t5)))

    stats: Counter[str] = Counter()
    barcode_counts: Counter[tuple[str, str]] = Counter()
    unknown_t7: Counter[str] = Counter()
    unknown_t5: Counter[str] = Counter()
    unknown_header: Counter[str] = Counter()
    prefix8_i7_candidate: Counter[str] = Counter()
    unknown_linker: Counter[str] = Counter()

    out1 = out2 = None
    if args.write:
        if not args.out_r1 or not args.out_r2:
            raise SystemExit("--write requires --out-r1 and --out-r2")
        args.out_r1.parent.mkdir(parents=True, exist_ok=True)
        args.out_r2.parent.mkdir(parents=True, exist_ok=True)
        out1 = gzip.open(args.out_r1, "wt", compresslevel=args.gzip_level)
        out2 = gzip.open(args.out_r2, "wt", compresslevel=args.gzip_level)

    try:
        for record in fastq_pairs(args.r1, args.r2):
            if args.max_reads and stats["total_pairs"] >= args.max_reads:
                break
            h1, s1, p1, q1, h2, s2, p2, q2 = record
            stats["total_pairs"] += 1

            r1_poly_a = trailing_run(s1, "A")
            if r1_poly_a >= 10:
                stats["R1_suffix_polyA_ge10"] += 1
            if r1_poly_a >= 20:
                stats["R1_suffix_polyA_ge20"] += 1

            if len(s2) < 8:
                stats["too_short_T7"] += 1
                stats["unassigned"] += 1
                continue

            observed_t7 = s2[:8]
            t7_record = t7_linkers.get(observed_t7)
            if t7_record is None:
                stats["T7_unmatched"] += 1
                unknown_t7[observed_t7] += 1
                stats["unassigned"] += 1
                continue
            stats["T7_exact"] += 1

            linker_start = 8
            linker_end = linker_start + t7_record.linker_length
            t5_end = linker_end + t5_length
            umi_end = t5_end + args.umi_length
            poly_t_end = umi_end + args.polyt_length
            if len(s2) < umi_end:
                stats["too_short_barcode_umi"] += 1
                stats["unassigned"] += 1
                continue

            observed_linker = s2[linker_start:linker_end]
            linker_exact = observed_linker == t7_record.linker_sequence
            if linker_exact:
                stats["linker_exact_given_T7"] += 1
            else:
                stats["linker_mismatch_given_T7"] += 1
                if s2.startswith(
                    t7_record.dna_residual_linker_sequence,
                    linker_start,
                ):
                    stats["DNA_residual_GTCA_linker_given_T7"] += 1
                else:
                    stats["other_linker_mismatch_given_T7"] += 1
                unknown_linker[
                    f"{t7_record.name}:{observed_linker}"
                ] += 1

            observed_t5 = s2[linker_end:t5_end]
            t5_name = t5.get(observed_t5)
            if t5_name is None:
                stats["T5_unmatched_at_expected_position"] += 1
                unknown_t5[observed_t5] += 1
            else:
                stats["T5_exact_at_expected_position"] += 1

            umi = s2[t5_end:umi_end]
            if DNA_RE.fullmatch(umi):
                stats["UMI_evaluable"] += 1
            else:
                stats["UMI_invalid"] += 1

            if len(s2) >= poly_t_end:
                stats["R2_polyT31_evaluable"] += 1
                poly_t = s2[umi_end:poly_t_end]
                t_fraction = poly_t.count("T") / args.polyt_length
                a_fraction = poly_t.count("A") / args.polyt_length
                if t_fraction == 1:
                    stats["R2_polyT31_all_T"] += 1
                if t_fraction >= 0.90:
                    stats["R2_polyT31_ge90pct_T"] += 1
                if a_fraction >= 0.90:
                    stats["R2_polyT31_ge90pct_A"] += 1
            else:
                stats["R2_polyT31_not_evaluable"] += 1

            indices_r1 = parse_header_indices(h1)
            indices_r2 = parse_header_indices(h2)
            if indices_r1 is None or indices_r2 is None:
                stats["bad_header_index"] += 1
                unknown_header["UNPARSED"] += 1
                stats["unassigned"] += 1
                continue
            if indices_r1 != indices_r2:
                stats["header_index_R1_R2_mismatch"] += 1
                unknown_header[
                    f"{indices_r1[0]}+{indices_r1[1]}|"
                    f"{indices_r2[0]}+{indices_r2[1]}"
                ] += 1
                stats["unassigned"] += 1
                continue
            idx1, idx2 = indices_r1
            stats["header_index_parsed"] += 1
            i7_name = None
            canonical_i7 = None
            i7_mode = None
            prefix8_candidate_name = None
            if len(idx1) == 10:
                i7_name = i7.get(idx1)
                if i7_name is not None:
                    canonical_i7 = idx1
                    i7_mode = "full10"
                    stats["I7_full10_exact"] += 1
            elif len(idx1) == 8:
                prefix_match = i7_prefix8.get(idx1)
                if prefix_match is not None:
                    prefix8_candidate_name, _ = prefix_match
                    stats["I7_prefix8_unique_candidate_rejected"] += 1
                    prefix8_i7_candidate[
                        f"{idx1}+{idx2}:{prefix8_candidate_name}"
                    ] += 1
            i5_name = i5.get(idx2)
            if i7_name is None:
                stats["I7_rejected"] += 1
                if prefix8_candidate_name is None:
                    stats["I7_unmatched"] += 1
            else:
                stats["I7_exact"] += 1
            if i5_name is None:
                stats["I5_unmatched"] += 1
            else:
                stats["I5_exact"] += 1
            if (
                (i7_name is None and prefix8_candidate_name is None)
                or i5_name is None
            ):
                unknown_header[f"{idx1}+{idx2}"] += 1

            structure_exact = linker_exact and t5_name is not None
            if structure_exact:
                stats["R2_T7_linker_T5_exact"] += 1
                if prefix8_candidate_name is not None and i5_name is not None:
                    stats["potential_prefix8_rescue_structure_I5_exact"] += 1

            if not (structure_exact and i7_name and i5_name):
                stats["unassigned"] += 1
                continue

            raw_barcode = observed_t7 + observed_t5 + idx1 + idx2
            if canonical_i7 is None or i7_mode is None:
                raise AssertionError("accepted I7 lacks canonical sequence/mode")
            canonical_barcode = (
                observed_t7 + observed_t5 + canonical_i7 + idx2
            )
            cb_name = f"{t7_record.name}_{t5_name}_{i7_name}_{i5_name}"
            barcode_counts[(cb_name, canonical_barcode)] += 1
            stats["strict_full_barcode"] += 1
            stats["strict_assigned"] += 1
            stats["strict_assigned_I7_full10"] += 1
            stats["assigned_XD0"] += 1
            stats["assigned"] += 1

            if out1 is not None and out2 is not None:
                tags = (
                    f"CR:Z:{raw_barcode}"
                    f"\tCB:Z:{canonical_barcode}"
                    f"\tCN:Z:{cb_name}"
                    f"\tUR:Z:{umi}"
                    f"\tUB:Z:{umi}"
                    f"\tXD:i:0"
                    f"\tXE:Z:0,0,0,0"
                    f"\tXL:i:{t7_record.linker_length}"
                    f"\tXI:Z:{i7_mode}"
                )
                write_record(out1, h1, tags, s1, p1, q1)
                write_record(out2, h2, tags, s2, p2, q2)
    finally:
        if out1 is not None:
            out1.close()
        if out2 is not None:
            out2.close()

    total = stats["total_pairs"] or 1
    ordered_metrics = (
        "total_pairs",
        "T7_exact",
        "T7_unmatched",
        "linker_exact_given_T7",
        "linker_mismatch_given_T7",
        "DNA_residual_GTCA_linker_given_T7",
        "other_linker_mismatch_given_T7",
        "T5_exact_at_expected_position",
        "T5_unmatched_at_expected_position",
        "R2_T7_linker_T5_exact",
        "UMI_evaluable",
        "UMI_invalid",
        "header_index_parsed",
        "bad_header_index",
        "header_index_R1_R2_mismatch",
        "I7_exact",
        "I7_full10_exact",
        "I7_prefix8_unique_candidate_rejected",
        "I7_rejected",
        "I7_unmatched",
        "I5_exact",
        "I5_unmatched",
        "potential_prefix8_rescue_structure_I5_exact",
        "strict_full_barcode",
        "strict_assigned",
        "strict_assigned_I7_full10",
        "assigned_XD0",
        "assigned",
        "unassigned",
        "too_short_T7",
        "too_short_barcode_umi",
        "R1_suffix_polyA_ge10",
        "R1_suffix_polyA_ge20",
        "R2_polyT31_evaluable",
        "R2_polyT31_all_T",
        "R2_polyT31_ge90pct_T",
        "R2_polyT31_ge90pct_A",
        "R2_polyT31_not_evaluable",
    )
    args.stat.parent.mkdir(parents=True, exist_ok=True)
    with args.stat.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pair_id", "sample", "metric", "count", "percent"))
        for metric in ordered_metrics:
            count = stats[metric]
            writer.writerow(
                (
                    args.pair_id,
                    args.sample,
                    metric,
                    count,
                    f"{count / total * 100:.6f}",
                )
            )

    args.bc_count.parent.mkdir(parents=True, exist_ok=True)
    with args.bc_count.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("pair_id", "sample", "barcode_name", "barcode_seq", "count")
        )
        for (name, sequence), count in barcode_counts.most_common():
            writer.writerow((args.pair_id, args.sample, name, sequence, count))

    unknown_payload = {
        "pair_id": args.pair_id,
        "sample": args.sample,
        "R2_first8": unknown_t7.most_common(args.top_n),
        "T5_at_expected_position": unknown_t5.most_common(args.top_n),
        "header_index": unknown_header.most_common(args.top_n),
        "I7_prefix8_candidate_rejected": prefix8_i7_candidate.most_common(
            args.top_n
        ),
        "linker_by_T7": unknown_linker.most_common(args.top_n),
    }
    args.unknown.parent.mkdir(parents=True, exist_ok=True)
    args.unknown.write_text(
        json.dumps(unknown_payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--t7-linker", required=True, type=Path)
    parser.add_argument("--t5", required=True, type=Path)
    parser.add_argument("--i7", required=True, type=Path)
    parser.add_argument("--i5", required=True, type=Path)
    parser.add_argument("--stat", required=True, type=Path)
    parser.add_argument("--bc-count", required=True, type=Path)
    parser.add_argument("--unknown", required=True, type=Path)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--umi-length", type=int, default=8)
    parser.add_argument("--polyt-length", type=int, default=31)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--out-r1", type=Path)
    parser.add_argument("--out-r2", type=Path)
    parser.add_argument("--gzip-level", type=int, choices=range(1, 10), default=4)
    parser.add_argument("--top-n", type=int, default=100)
    return parser.parse_args()


if __name__ == "__main__":
    process(parse_args())
