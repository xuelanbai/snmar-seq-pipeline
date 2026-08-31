#!/usr/bin/env python3
"""Step 04: probe or demultiplex snMAR DNA FASTQ by four barcode segments."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


DNA_RE = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
SAMPLE_NUMBER_RE = re.compile(r"-(\d+)$")


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def load_barcodes(path: Path) -> dict[str, str]:
    seq_to_name: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_no}: expected NAME SEQUENCE")
            name, seq = fields[0], fields[1].upper()
            if not DNA_RE.fullmatch(seq):
                raise ValueError(f"{path}:{line_no}: invalid sequence {seq!r}")
            if seq in seq_to_name and seq_to_name[seq] != name:
                raise ValueError(f"{path}:{line_no}: duplicate sequence {seq}")
            seq_to_name[seq] = name
    if not seq_to_name:
        raise ValueError(f"empty barcode file: {path}")
    lengths = {len(seq) for seq in seq_to_name}
    if len(lengths) != 1:
        raise ValueError(f"mixed barcode lengths in {path}: {sorted(lengths)}")
    return seq_to_name


@dataclass(frozen=True)
class Match:
    status: str
    name: str | None
    corrected: str | None
    distance: int | None


class BarcodeMatcher:
    def __init__(self, seq_to_name: dict[str, str]):
        self.seq_to_name = seq_to_name
        self.sequences = tuple(seq_to_name)
        self.length = len(self.sequences[0])
        neighbours: dict[str, list[str]] = {}
        alphabet = "ACGTN"
        for candidate in self.sequences:
            for position, original in enumerate(candidate):
                for base in alphabet:
                    if base == original:
                        continue
                    observed = candidate[:position] + base + candidate[position + 1 :]
                    neighbours.setdefault(observed, []).append(candidate)
        self.one_mismatch_neighbours = neighbours

    def match(self, observed: str, max_mismatch: int) -> Match:
        observed = observed.upper()
        exact_name = self.seq_to_name.get(observed)
        if exact_name is not None:
            return Match("exact", exact_name, observed, 0)
        if max_mismatch <= 0 or len(observed) != self.length:
            return Match("unmatched", None, None, None)
        best = self.one_mismatch_neighbours.get(observed, [])
        if not best:
            return Match("unmatched", None, None, None)
        unique_best = tuple(dict.fromkeys(best))
        if len(unique_best) != 1:
            return Match("ambiguous", None, None, 1)
        corrected = unique_best[0]
        return Match("rescued", self.seq_to_name[corrected], corrected, 1)


def fastq_pairs(r1_path: Path, r2_path: Path) -> Iterator[tuple[str, str, str, str, str, str, str, str]]:
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
            s1, p1, q1 = r1.readline().rstrip(), r1.readline().rstrip(), r1.readline().rstrip()
            s2, p2, q2 = r2.readline().rstrip(), r2.readline().rstrip(), r2.readline().rstrip()
            if not all((s1, p1, q1, s2, p2, q2)):
                raise ValueError(f"truncated FASTQ record near pair {record_no}")
            if len(s1) != len(q1) or len(s2) != len(q2):
                raise ValueError(f"sequence/quality length mismatch near pair {record_no}")
            yield h1, s1.upper(), p1, q1, h2, s2.upper(), p2, q2


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


def write_record(
    handle: TextIO,
    header: str,
    tags: str,
    seq: str,
    plus: str,
    qual: str,
) -> None:
    read_id = header.split(maxsplit=1)[0]
    handle.write(f"{read_id} {tags}\n{seq}\n{plus}\n{qual}\n")


def sample_expected_t7(sample: str, t7_by_name: dict[str, str]) -> str | None:
    match = SAMPLE_NUMBER_RE.search(sample)
    if not match:
        return None
    return t7_by_name.get(f"T7-{int(match.group(1))}")


def orientation_metrics(
    stats: Counter[str],
    idx1: str,
    idx2: str,
    i7: dict[str, str],
    i5: dict[str, str],
) -> None:
    variants = {
        "header_I7_I5_raw": (idx1, idx2, i7, i5),
        "header_I7rc_I5_raw": (reverse_complement(idx1), idx2, i7, i5),
        "header_I7_I5rc_raw": (idx1, reverse_complement(idx2), i7, i5),
        "header_I7rc_I5rc_raw": (
            reverse_complement(idx1),
            reverse_complement(idx2),
            i7,
            i5,
        ),
        "header_I5_I7_swapped": (idx1, idx2, i5, i7),
        "header_I5rc_I7_swapped": (reverse_complement(idx1), idx2, i5, i7),
        "header_I5_I7rc_swapped": (idx1, reverse_complement(idx2), i5, i7),
        "header_I5rc_I7rc_swapped": (
            reverse_complement(idx1),
            reverse_complement(idx2),
            i5,
            i7,
        ),
    }
    for metric, (left, right, left_map, right_map) in variants.items():
        if left in left_map and right in right_map:
            stats[metric] += 1


def process(args: argparse.Namespace) -> None:
    t5 = load_barcodes(args.t5)
    t7 = load_barcodes(args.t7)
    i5 = load_barcodes(args.i5)
    i7 = load_barcodes(args.i7)
    matchers = {
        "T5": BarcodeMatcher(t5),
        "T7": BarcodeMatcher(t7),
        "I5": BarcodeMatcher(i5),
        "I7": BarcodeMatcher(i7),
    }
    t7_by_name = {name: seq for seq, name in t7.items()}
    expected_t7 = sample_expected_t7(args.sample, t7_by_name)

    stats: Counter[str] = Counter()
    barcode_counts: Counter[tuple[str, str]] = Counter()
    unknown_r1: Counter[str] = Counter()
    unknown_r2: Counter[str] = Counter()
    unknown_header: Counter[str] = Counter()

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
            h1, s1, p1, q1, h2, s2, p2, q2 = record
            if args.max_reads and stats["total_pairs"] >= args.max_reads:
                break
            stats["total_pairs"] += 1

            if len(s1) < matchers["T5"].length or len(s2) < matchers["T7"].length:
                stats["too_short"] += 1
                stats["unassigned"] += 1
                continue

            obs_t5 = s1[: matchers["T5"].length]
            obs_t7 = s2[: matchers["T7"].length]
            if obs_t5 in t5:
                stats["R1_T5_exact"] += 1
            else:
                unknown_r1[obs_t5] += 1
            if obs_t7 in t7:
                stats["R2_T7_exact"] += 1
            else:
                unknown_r2[obs_t7] += 1
            if obs_t5 in t5 and obs_t7 in t7:
                stats["R1_T5_and_R2_T7"] += 1
            if s1[:8] in t7 and s2[:8] in t5:
                stats["R1_T7_and_R2_T5_swapped"] += 1
            if expected_t7 is not None:
                stats["R2_expected_T7_evaluable"] += 1
                if obs_t7 == expected_t7:
                    stats["R2_expected_T7_from_sample_name"] += 1

            indices = parse_header_indices(h1)
            if indices is None:
                stats["bad_header_index"] += 1
                unknown_header["UNPARSED"] += 1
                stats["unassigned"] += 1
                continue
            idx1, idx2 = indices
            stats["header_index_parsed"] += 1
            orientation_metrics(stats, idx1, idx2, i7, i5)
            if not (idx1 in i7 and idx2 in i5):
                unknown_header[f"{idx1}+{idx2}"] += 1

            exact_matches = (
                matchers["T5"].match(obs_t5, 0),
                matchers["T7"].match(obs_t7, 0),
                matchers["I7"].match(idx1, 0),
                matchers["I5"].match(idx2, 0),
            )
            if all(match.status == "exact" for match in exact_matches):
                stats["strict_full_barcode"] += 1

            matches = (
                matchers["T5"].match(obs_t5, args.max_segment_mismatch),
                matchers["T7"].match(obs_t7, args.max_segment_mismatch),
                matchers["I7"].match(idx1, args.max_segment_mismatch),
                matchers["I5"].match(idx2, args.max_segment_mismatch),
            )
            labels = ("T5", "T7", "I7", "I5")
            for label, match in zip(labels, matches):
                stats[f"{label}_{match.status}"] += 1

            if any(match.status == "ambiguous" for match in matches):
                stats["ambiguous_barcode"] += 1
                stats["unassigned"] += 1
                continue
            if any(match.name is None for match in matches):
                stats["unmatched_barcode"] += 1
                stats["unassigned"] += 1
                continue

            total_distance = sum(int(match.distance or 0) for match in matches)
            if total_distance > args.max_total_mismatch:
                stats["over_total_mismatch_limit"] += 1
                stats["unassigned"] += 1
                continue

            t5_match, t7_match, i7_match, i5_match = matches
            assert all(match.corrected and match.name for match in matches)
            raw_barcode = obs_t7 + obs_t5 + idx1 + idx2
            cb_seq = (
                str(t7_match.corrected)
                + str(t5_match.corrected)
                + str(i7_match.corrected)
                + str(i5_match.corrected)
            )
            cb_name = (
                f"{t7_match.name}_{t5_match.name}_{i7_match.name}_{i5_match.name}"
            )
            segment_distances = (
                int(t7_match.distance or 0),
                int(t5_match.distance or 0),
                int(i7_match.distance or 0),
                int(i5_match.distance or 0),
            )
            distance_vector = ",".join(str(distance) for distance in segment_distances)
            barcode_counts[(cb_name, cb_seq)] += 1
            stats["assigned"] += 1
            stats[f"assigned_XD{total_distance}"] += 1
            if total_distance == 0:
                stats["strict_assigned"] += 1
            else:
                stats["rescued_assigned"] += 1
                stats["rescue_unique_total_le2"] += 1

            if out1 is not None and out2 is not None:
                tags = (
                    f"CR:Z:{raw_barcode}"
                    f"\tCB:Z:{cb_seq}"
                    f"\tCN:Z:{cb_name}"
                    f"\tXD:i:{total_distance}"
                    f"\tXE:Z:{distance_vector}"
                )
                write_record(out1, h1, tags, s1, p1, q1)
                write_record(out2, h2, tags, s2, p2, q2)
    finally:
        if out1 is not None:
            out1.close()
        if out2 is not None:
            out2.close()

    total = stats["total_pairs"] or 1
    args.stat.parent.mkdir(parents=True, exist_ok=True)
    with args.stat.open("w", encoding="utf-8") as handle:
        handle.write("pair_id\tsample\tmetric\tcount\tpercent\n")
        ordered = (
            "total_pairs",
            "R1_T5_exact",
            "R2_T7_exact",
            "R1_T5_and_R2_T7",
            "R1_T7_and_R2_T5_swapped",
            "R2_expected_T7_evaluable",
            "R2_expected_T7_from_sample_name",
            "header_index_parsed",
            "header_I7_I5_raw",
            "header_I7rc_I5_raw",
            "header_I7_I5rc_raw",
            "header_I7rc_I5rc_raw",
            "header_I5_I7_swapped",
            "header_I5rc_I7_swapped",
            "header_I5_I7rc_swapped",
            "header_I5rc_I7rc_swapped",
            "strict_full_barcode",
            "strict_assigned",
            "rescued_assigned",
            "rescue_unique_total_le2",
            "assigned_XD0",
            "assigned_XD1",
            "assigned_XD2",
            "assigned",
            "unassigned",
            "ambiguous_barcode",
            "unmatched_barcode",
            "over_total_mismatch_limit",
            "too_short",
            "bad_header_index",
            "T5_exact",
            "T5_rescued",
            "T5_ambiguous",
            "T5_unmatched",
            "T7_exact",
            "T7_rescued",
            "T7_ambiguous",
            "T7_unmatched",
            "I7_exact",
            "I7_rescued",
            "I7_ambiguous",
            "I7_unmatched",
            "I5_exact",
            "I5_rescued",
            "I5_ambiguous",
            "I5_unmatched",
        )
        for metric in ordered:
            count = stats[metric]
            handle.write(
                f"{args.pair_id}\t{args.sample}\t{metric}\t{count}\t"
                f"{count / total * 100:.6f}\n"
            )

    args.bc_count.parent.mkdir(parents=True, exist_ok=True)
    with args.bc_count.open("w", encoding="utf-8") as handle:
        handle.write("pair_id\tsample\tbarcode_name\tbarcode_seq\tcount\n")
        for (name, seq), count in barcode_counts.most_common():
            handle.write(f"{args.pair_id}\t{args.sample}\t{name}\t{seq}\t{count}\n")

    unknown_payload = {
        "pair_id": args.pair_id,
        "sample": args.sample,
        "R1_first8": unknown_r1.most_common(args.top_n),
        "R2_first8": unknown_r2.most_common(args.top_n),
        "header_index": unknown_header.most_common(args.top_n),
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
    parser.add_argument("--t5", required=True, type=Path)
    parser.add_argument("--t7", required=True, type=Path)
    parser.add_argument("--i5", required=True, type=Path)
    parser.add_argument("--i7", required=True, type=Path)
    parser.add_argument("--stat", required=True, type=Path)
    parser.add_argument("--bc-count", required=True, type=Path)
    parser.add_argument("--unknown", required=True, type=Path)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--max-segment-mismatch", type=int, choices=(0, 1), default=0)
    parser.add_argument("--max-total-mismatch", type=int, choices=(0, 1, 2), default=0)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--out-r1", type=Path)
    parser.add_argument("--out-r2", type=Path)
    parser.add_argument("--gzip-level", type=int, choices=range(1, 10), default=1)
    parser.add_argument("--top-n", type=int, default=30)
    return parser.parse_args()


if __name__ == "__main__":
    process(parse_args())
