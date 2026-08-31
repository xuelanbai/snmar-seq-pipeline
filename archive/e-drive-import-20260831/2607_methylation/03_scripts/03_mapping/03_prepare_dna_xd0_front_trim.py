#!/usr/bin/env python3
"""Step 03: keep paired XD=0 DNA reads and remove inline library structure."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path
from typing import Iterator, TextIO


TAG_RE = re.compile(
    r"(?:^|[ \t])"
    r"(?P<tag>(?:CR|CB|CN|UR|UB|XD|XE|XL|XI):[AifZHB]:[^ \t\r\n]+)"
)
TAG_ORDER = ("CR", "CB", "CN", "UR", "UB", "XD", "XE", "XL", "XI")


def paired_records(
    r1_path: Path,
    r2_path: Path,
) -> Iterator[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]]:
    with gzip.open(r1_path, "rt") as r1, gzip.open(r2_path, "rt") as r2:
        record_no = 0
        while True:
            h1 = r1.readline().rstrip()
            h2 = r2.readline().rstrip()
            if not h1 and not h2:
                break
            record_no += 1
            if not h1 or not h2:
                raise ValueError(
                    f"R1/R2 count differs near pair {record_no}: "
                    f"{r1_path}, {r2_path}"
                )
            rec1 = (
                h1,
                r1.readline().rstrip(),
                r1.readline().rstrip(),
                r1.readline().rstrip(),
            )
            rec2 = (
                h2,
                r2.readline().rstrip(),
                r2.readline().rstrip(),
                r2.readline().rstrip(),
            )
            for path, (header, sequence, plus, quality) in (
                (r1_path, rec1),
                (r2_path, rec2),
            ):
                if not all((sequence, plus, quality)):
                    raise ValueError(f"{path}: truncated FASTQ near {record_no}")
                if len(sequence) != len(quality):
                    raise ValueError(
                        f"{path}: sequence/quality mismatch near {record_no}"
                    )
            id1 = rec1[0].split(maxsplit=1)[0].removeprefix("@").removesuffix("/1")
            id2 = rec2[0].split(maxsplit=1)[0].removeprefix("@").removesuffix("/2")
            if id1 != id2:
                raise ValueError(f"R1/R2 read IDs differ near pair {record_no}")
            yield rec1, rec2


def normalize_header(header: str) -> tuple[str, str]:
    tags: dict[str, str] = {}
    for match in TAG_RE.finditer(header):
        tag = match.group("tag")
        tags[tag[:2]] = tag
    if "XD" not in tags:
        raise ValueError(f"FASTQ header lacks XD tag: {header[:200]}")
    read_id = header.split(maxsplit=1)[0]
    ordered = [tags[name] for name in TAG_ORDER if name in tags]
    return f"{read_id} " + "\t".join(ordered), tags["XD"]


def trim_record(
    record: tuple[str, str, str, str],
    normalized_header: str,
    trim_front: int,
) -> tuple[str, str, str, str]:
    _, sequence, plus, quality = record
    return (
        normalized_header,
        sequence[trim_front:],
        plus,
        quality[trim_front:],
    )


def write_record(handle: TextIO, record: tuple[str, str, str, str]) -> None:
    handle.write("\n".join(record) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r2", required=True, type=Path)
    parser.add_argument("--stat", required=True, type=Path)
    parser.add_argument("--trim-front-r1", type=int, default=30)
    parser.add_argument("--trim-front-r2", type=int, default=43)
    parser.add_argument("--max-input-pairs", type=int, default=0)
    parser.add_argument("--gzip-level", type=int, choices=range(1, 10), default=4)
    args = parser.parse_args()

    for value, label in (
        (args.trim_front_r1, "--trim-front-r1"),
        (args.trim_front_r2, "--trim-front-r2"),
    ):
        if value < 0:
            raise SystemExit(f"{label} must be non-negative")

    args.out_r1.parent.mkdir(parents=True, exist_ok=True)
    args.out_r2.parent.mkdir(parents=True, exist_ok=True)
    counts = {
        "total_pairs": 0,
        "XD0_pairs": 0,
        "non_XD0_pairs": 0,
        "too_short_after_front_trim": 0,
        "output_pairs": 0,
    }

    with gzip.open(
        args.out_r1, "wt", compresslevel=args.gzip_level
    ) as out1, gzip.open(
        args.out_r2, "wt", compresslevel=args.gzip_level
    ) as out2:
        for rec1, rec2 in paired_records(args.r1, args.r2):
            if (
                args.max_input_pairs
                and counts["total_pairs"] >= args.max_input_pairs
            ):
                break
            counts["total_pairs"] += 1
            header1, xd1 = normalize_header(rec1[0])
            header2, xd2 = normalize_header(rec2[0])
            if xd1 != xd2:
                raise ValueError(
                    f"R1/R2 XD tags differ near input pair "
                    f"{counts['total_pairs']}: {xd1} vs {xd2}"
                )
            if xd1 != "XD:i:0":
                counts["non_XD0_pairs"] += 1
                continue
            counts["XD0_pairs"] += 1
            trimmed1 = trim_record(
                rec1,
                header1,
                args.trim_front_r1,
            )
            trimmed2 = trim_record(
                rec2,
                header2,
                args.trim_front_r2,
            )
            if not trimmed1[1] or not trimmed2[1]:
                counts["too_short_after_front_trim"] += 1
                continue
            write_record(out1, trimmed1)
            write_record(out2, trimmed2)
            counts["output_pairs"] += 1

    args.stat.parent.mkdir(parents=True, exist_ok=True)
    with args.stat.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "pair_id",
                "sample",
                "total_pairs",
                "XD0_pairs",
                "non_XD0_pairs",
                "too_short_after_front_trim",
                "output_pairs",
                "trim_front_R1",
                "trim_front_R2",
            )
        )
        writer.writerow(
            (
                args.pair_id,
                args.sample,
                counts["total_pairs"],
                counts["XD0_pairs"],
                counts["non_XD0_pairs"],
                counts["too_short_after_front_trim"],
                counts["output_pairs"],
                args.trim_front_r1,
                args.trim_front_r2,
            )
        )
    print(
        f"[OK] {args.pair_id}: input={counts['total_pairs']} "
        f"XD0={counts['XD0_pairs']} output={counts['output_pairs']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
