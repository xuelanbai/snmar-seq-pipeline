#!/usr/bin/env python3
"""Step 11: diagnose RNA linker variants, T5 offsets, and index lengths.

This is a read-only probe. It does not write demultiplexed FASTQ files and does
not change the strict matching rules in 04_rna_strict_qc_demultiplex.py.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


DNA_RE = re.compile(r"^[ACGTN]+$", re.IGNORECASE)
GTCA_INSERT_MOTIF = "CTCCGACGC"
GTCA_INSERT = "GTCA"


@dataclass(frozen=True)
class T7Linker:
    name: str
    sequence: str
    expected_linker: str
    gtca_linker: str


@dataclass(frozen=True)
class PairRow:
    pair_id: str
    sample: str
    r2: str


@dataclass
class PairResult:
    pair_id: str
    sample: str
    stats: Counter[str]
    offsets: Counter[tuple[str, str, str, int]]
    top_sequences: Counter[tuple[str, str]]
    index_lengths: Counter[tuple[int, int]]


ORDERED_METRICS = (
    "total_reads",
    "T7_exact",
    "T7_unmatched",
    "too_short_T7",
    "linker_expected_exact",
    "linker_GTCA_exact",
    "linker_neither_exact",
    "T5_exact_after_expected_linker",
    "T5_exact_after_GTCA_linker",
    "T5_exact_scanned_any_forward",
    "T5_exact_scanned_any_reverse_complement",
    "T5_no_exact_hit_in_scan",
    "polyT31_evaluable_after_expected_structure",
    "polyT31_ge90pct_after_expected_structure",
    "polyT31_evaluable_after_GTCA_structure",
    "polyT31_ge90pct_after_GTCA_structure",
    "header_index_parsed",
    "bad_header_index",
    "I7_full10_exact",
    "I7_prefix8_exact_unique",
    "I7_unmatched",
    "I5_exact",
    "I5_unmatched",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_two_column_index(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
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
            if sequence in result:
                raise ValueError(f"{path}:{line_no}: duplicate sequence {sequence}")
            result[sequence] = name
    if not result:
        raise ValueError(f"empty index table: {path}")
    return result


def make_gtca_linker(linker: str, source: Path, line_no: int) -> str:
    if linker.count(GTCA_INSERT_MOTIF) != 1:
        raise ValueError(
            f"{source}:{line_no}: linker must contain one "
            f"{GTCA_INSERT_MOTIF} motif"
        )
    return linker.replace(
        GTCA_INSERT_MOTIF,
        GTCA_INSERT_MOTIF + GTCA_INSERT,
        1,
    )


def read_t7_linkers(path: Path) -> dict[str, T7Linker]:
    result: dict[str, T7Linker] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"name", "sequence", "linker_length", "linker_sequence"}
        if set(reader.fieldnames or ()) != required:
            raise ValueError(
                f"{path}: expected name, sequence, linker_length, linker_sequence"
            )
        for line_no, row in enumerate(reader, 2):
            sequence = row["sequence"].upper()
            linker = row["linker_sequence"].upper()
            if (
                len(sequence) != 8
                or not DNA_RE.fullmatch(sequence)
                or int(row["linker_length"]) != len(linker)
                or not DNA_RE.fullmatch(linker)
            ):
                raise ValueError(f"{path}:{line_no}: invalid T7/linker row")
            if sequence in result:
                raise ValueError(f"{path}:{line_no}: duplicate T7 {sequence}")
            result[sequence] = T7Linker(
                name=row["name"],
                sequence=sequence,
                expected_linker=linker,
                gtca_linker=make_gtca_linker(linker, path, line_no),
            )
    if not result:
        raise ValueError(f"empty T7/linker table: {path}")
    return result


def read_pair_table(path: Path) -> list[PairRow]:
    rows: list[PairRow] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"pair_id", "sample", "R1", "R2"}
        if set(reader.fieldnames or ()) != required:
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in reader:
            rows.append(PairRow(row["pair_id"], row["sample"], row["R2"]))
    if not rows:
        raise ValueError(f"empty pair table: {path}")
    if len({row.pair_id for row in rows}) != len(rows):
        raise ValueError(f"{path}: duplicate pair_id")
    return rows


def fastq_records(path: Path) -> Iterator[tuple[str, str]]:
    with gzip.open(path, "rt") as handle:
        record_no = 0
        while True:
            header = handle.readline().rstrip()
            if not header:
                break
            record_no += 1
            sequence = handle.readline().rstrip().upper()
            plus = handle.readline().rstrip()
            quality = handle.readline().rstrip()
            if not sequence or not plus or not quality:
                raise ValueError(f"{path}: truncated FASTQ at record {record_no}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"{path}: invalid FASTQ at record {record_no}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: sequence/quality length mismatch at record {record_no}"
                )
            yield header, sequence


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


def prefix_map(
    index: dict[str, str],
    length: int,
) -> tuple[dict[str, str], set[str]]:
    grouped: dict[str, list[str]] = {}
    for sequence, name in index.items():
        grouped.setdefault(sequence[:length], []).append(name)
    unique = {
        prefix: names[0]
        for prefix, names in grouped.items()
        if len(names) == 1
    }
    ambiguous = {prefix for prefix, names in grouped.items() if len(names) > 1}
    return unique, ambiguous


def polyt_ge90(
    sequence: str,
    t5_start: int,
    t5_length: int,
    umi_length: int,
    polyt_length: int,
) -> bool | None:
    start = t5_start + t5_length + umi_length
    end = start + polyt_length
    if len(sequence) < end:
        return None
    return sequence[start:end].count("T") / polyt_length >= 0.90


def process_pair(
    row: PairRow,
    t7_linkers: dict[str, T7Linker],
    t5_forward: dict[str, str],
    t5_reverse: dict[str, str],
    i7: dict[str, str],
    i7_prefix8: dict[str, str],
    i5: dict[str, str],
    max_reads: int,
    scan_left: int,
    scan_right: int,
    umi_length: int,
    polyt_length: int,
) -> PairResult:
    stats: Counter[str] = Counter()
    offsets: Counter[tuple[str, str, str, int]] = Counter()
    top_sequences: Counter[tuple[str, str]] = Counter()
    index_lengths: Counter[tuple[int, int]] = Counter()
    t5_length = len(next(iter(t5_forward)))

    for header, sequence in fastq_records(Path(row.r2)):
        if max_reads and stats["total_reads"] >= max_reads:
            break
        stats["total_reads"] += 1

        indices = parse_header_indices(header)
        if indices is None:
            stats["bad_header_index"] += 1
        else:
            idx1, idx2 = indices
            stats["header_index_parsed"] += 1
            index_lengths[(len(idx1), len(idx2))] += 1
            if idx1 in i7:
                stats["I7_full10_exact"] += 1
            elif len(idx1) == 8 and idx1 in i7_prefix8:
                stats["I7_prefix8_exact_unique"] += 1
            else:
                stats["I7_unmatched"] += 1
            if idx2 in i5:
                stats["I5_exact"] += 1
            else:
                stats["I5_unmatched"] += 1

        if len(sequence) < 8:
            stats["too_short_T7"] += 1
            continue
        t7_record = t7_linkers.get(sequence[:8])
        if t7_record is None:
            stats["T7_unmatched"] += 1
            continue
        stats["T7_exact"] += 1

        linker_start = 8
        expected_start = linker_start + len(t7_record.expected_linker)
        gtca_start = linker_start + len(t7_record.gtca_linker)
        expected_exact = sequence.startswith(
            t7_record.expected_linker,
            linker_start,
        )
        gtca_exact = sequence.startswith(
            t7_record.gtca_linker,
            linker_start,
        )
        if expected_exact:
            stats["linker_expected_exact"] += 1
            linker_class = "expected"
        elif gtca_exact:
            stats["linker_GTCA_exact"] += 1
            linker_class = "GTCA"
        else:
            stats["linker_neither_exact"] += 1
            linker_class = "neither"

        for label, start in (
            ("after_expected_first8", expected_start),
            ("after_expected_first16", expected_start),
            ("after_GTCA_first8", gtca_start),
            ("after_GTCA_first16", gtca_start),
        ):
            length = 16 if label.endswith("16") else 8
            observed = sequence[start : start + length]
            if len(observed) == length:
                top_sequences[(label, observed)] += 1

        expected_t5 = sequence[expected_start : expected_start + t5_length]
        gtca_t5 = sequence[gtca_start : gtca_start + t5_length]
        expected_t5_exact = expected_t5 in t5_forward
        gtca_t5_exact = gtca_t5 in t5_forward
        if expected_exact and expected_t5_exact:
            stats["T5_exact_after_expected_linker"] += 1
        if gtca_exact and gtca_t5_exact:
            stats["T5_exact_after_GTCA_linker"] += 1

        if expected_exact and expected_t5_exact:
            polyt = polyt_ge90(
                sequence,
                expected_start,
                t5_length,
                umi_length,
                polyt_length,
            )
            if polyt is not None:
                stats["polyT31_evaluable_after_expected_structure"] += 1
                if polyt:
                    stats["polyT31_ge90pct_after_expected_structure"] += 1
        if gtca_exact and gtca_t5_exact:
            polyt = polyt_ge90(
                sequence,
                gtca_start,
                t5_length,
                umi_length,
                polyt_length,
            )
            if polyt is not None:
                stats["polyT31_evaluable_after_GTCA_structure"] += 1
                if polyt:
                    stats["polyT31_ge90pct_after_GTCA_structure"] += 1

        hits: list[tuple[int, int, str]] = []
        forward_seen = False
        reverse_seen = False
        for offset in range(-scan_left, scan_right + 1):
            start = expected_start + offset
            if start < 8 or start + t5_length > len(sequence):
                continue
            observed = sequence[start : start + t5_length]
            if observed in t5_forward:
                offsets[(linker_class, "all", "forward", offset)] += 1
                hits.append((abs(offset), offset, "forward"))
                forward_seen = True
            if observed in t5_reverse:
                offsets[
                    (linker_class, "all", "reverse_complement", offset)
                ] += 1
                hits.append((abs(offset), offset, "reverse_complement"))
                reverse_seen = True
        if forward_seen:
            stats["T5_exact_scanned_any_forward"] += 1
        if reverse_seen:
            stats["T5_exact_scanned_any_reverse_complement"] += 1
        if not hits:
            stats["T5_no_exact_hit_in_scan"] += 1
        else:
            _, nearest_offset, orientation = min(
                hits,
                key=lambda hit: (
                    hit[0],
                    hit[1],
                    0 if hit[2] == "forward" else 1,
                ),
            )
            offsets[
                (linker_class, "nearest", orientation, nearest_offset)
            ] += 1

    return PairResult(
        pair_id=row.pair_id,
        sample=row.sample,
        stats=stats,
        offsets=offsets,
        top_sequences=top_sequences,
        index_lengths=index_lengths,
    )


def write_metric_tables(results: list[PairResult], out_dir: Path) -> None:
    total_stats: Counter[str] = Counter()
    for result in results:
        total_stats.update(result.stats)
    total = total_stats["total_reads"] or 1

    with (out_dir / "structure_probe_summary.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "count", "percent_total"))
        for metric in ORDERED_METRICS:
            count = total_stats[metric]
            writer.writerow((metric, count, f"{count / total * 100:.6f}"))

    with (out_dir / "structure_probe_per_pair.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pair_id", "sample", "metric", "count", "percent"))
        for result in sorted(results, key=lambda item: (item.sample, item.pair_id)):
            pair_total = result.stats["total_reads"] or 1
            for metric in ORDERED_METRICS:
                count = result.stats[metric]
                writer.writerow(
                    (
                        result.pair_id,
                        result.sample,
                        metric,
                        count,
                        f"{count / pair_total * 100:.6f}",
                    )
                )


def write_counter_tables(
    results: list[PairResult],
    out_dir: Path,
    top_n: int,
) -> None:
    total = sum(result.stats["total_reads"] for result in results) or 1
    offsets: Counter[tuple[str, str, str, int]] = Counter()
    top_sequences: Counter[tuple[str, str]] = Counter()
    index_lengths: Counter[tuple[int, int]] = Counter()
    for result in results:
        offsets.update(result.offsets)
        top_sequences.update(result.top_sequences)
        index_lengths.update(result.index_lengths)

    with (out_dir / "structure_probe_t5_offsets.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "linker_class",
                "hit_scope",
                "orientation",
                "offset_from_expected",
                "count",
                "percent_total",
            )
        )
        for key, count in sorted(
            offsets.items(),
            key=lambda item: (
                item[0][0],
                item[0][1],
                item[0][2],
                item[0][3],
            ),
        ):
            writer.writerow((*key, count, f"{count / total * 100:.6f}"))

    by_category: dict[str, Counter[str]] = {}
    for (category, sequence), count in top_sequences.items():
        by_category.setdefault(category, Counter())[sequence] += count
    with (out_dir / "structure_probe_top_sequences.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("category", "sequence", "count", "percent_total"))
        for category in sorted(by_category):
            for sequence, count in by_category[category].most_common(top_n):
                writer.writerow(
                    (category, sequence, count, f"{count / total * 100:.6f}")
                )

    with (out_dir / "structure_probe_index_lengths.tsv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("I7_length", "I5_length", "count", "percent_total"))
        for (i7_length, i5_length), count in index_lengths.most_common():
            writer.writerow(
                (
                    i7_length,
                    i5_length,
                    count,
                    f"{count / total * 100:.6f}",
                )
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--t7-linker", required=True, type=Path)
    parser.add_argument("--t5", required=True, type=Path)
    parser.add_argument("--i7", required=True, type=Path)
    parser.add_argument("--i5", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=32)
    parser.add_argument("--max-reads", type=int, default=200000)
    parser.add_argument("--scan-left", type=int, default=8)
    parser.add_argument("--scan-right", type=int, default=64)
    parser.add_argument("--umi-length", type=int, default=8)
    parser.add_argument("--polyt-length", type=int, default=31)
    parser.add_argument("--top-n", type=int, default=100)
    args = parser.parse_args()

    if args.jobs < 1 or args.max_reads < 0:
        raise SystemExit("--jobs must be >=1 and --max-reads must be >=0")
    if min(args.scan_left, args.scan_right, args.umi_length, args.polyt_length) < 0:
        raise SystemExit("scan/UMI/polyT lengths must be non-negative")

    pair_rows = read_pair_table(args.pair_table)
    t7_linkers = read_t7_linkers(args.t7_linker)
    t5_forward = read_two_column_index(args.t5)
    i7 = read_two_column_index(args.i7)
    i5 = read_two_column_index(args.i5)
    if {len(sequence) for sequence in t5_forward} != {8}:
        raise SystemExit("T5 table must contain only 8-bp sequences")
    if {len(sequence) for sequence in i7} != {10}:
        raise SystemExit("I7 table must contain only 10-bp sequences")
    if {len(sequence) for sequence in i5} != {10}:
        raise SystemExit("I5 table must contain only 10-bp sequences")
    i7_prefix8, ambiguous_prefixes = prefix_map(i7, 8)
    if ambiguous_prefixes:
        preview = ",".join(sorted(ambiguous_prefixes)[:10])
        raise SystemExit(f"I7 first-8 prefixes are not unique: {preview}")
    t5_reverse = {
        reverse_complement(sequence): name
        for sequence, name in t5_forward.items()
    }
    if len(t5_reverse) != len(t5_forward):
        raise SystemExit("reverse-complement T5 sequences are not unique")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results: list[PairResult] = []
    workers = min(args.jobs, len(pair_rows))
    print(
        f"[INFO] pairs={len(pair_rows)} jobs={workers} "
        f"max_reads_per_pair={args.max_reads}"
    )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_row = {
            executor.submit(
                process_pair,
                row,
                t7_linkers,
                t5_forward,
                t5_reverse,
                i7,
                i7_prefix8,
                i5,
                args.max_reads,
                args.scan_left,
                args.scan_right,
                args.umi_length,
                args.polyt_length,
            ): row
            for row in pair_rows
        }
        for future in as_completed(future_to_row):
            row = future_to_row[future]
            result = future.result()
            results.append(result)
            print(
                f"[DONE] {row.sample} {row.pair_id}: "
                f"{result.stats['total_reads']} reads"
            )

    write_metric_tables(results, args.out_dir)
    write_counter_tables(results, args.out_dir, args.top_n)
    print(f"[OK] diagnostic tables: {args.out_dir}")
    print("[NOTE] strict demultiplexing code and whitelist tables were not changed.")


if __name__ == "__main__":
    main()
