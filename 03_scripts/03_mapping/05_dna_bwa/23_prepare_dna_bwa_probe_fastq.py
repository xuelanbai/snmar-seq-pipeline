#!/usr/bin/env python3
"""Write a bounded paired DNA FASTQ probe with SAM-safe comment tags.

The demultiplex mismatch vector is called XE in the source FASTQ. BWA owns
several X* auxiliary tags and may emit XE itself, so this script deliberately
renames the source mismatch vector to ZE before ``bwa mem -C`` imports the
comment fields into BAM.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import tempfile
from pathlib import Path
from typing import TextIO

REQUIRED = ("CR", "CB", "CN", "XD", "XE")
TAG_PATTERNS = {
    "CR": re.compile(r"^[ACGTN]+$"),
    "CB": re.compile(r"^[ACGTN]+$"),
    "CN": re.compile(r"^\S+$"),
    "XD": re.compile(r"^[0-9]+$"),
    "XE": re.compile(r"^[0-9]+(?:,[0-9]+)*$"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--out-r1", required=True, type=Path)
    parser.add_argument("--out-r2", required=True, type=Path)
    parser.add_argument("--stats-out", required=True, type=Path)
    parser.add_argument("--max-pairs", required=True, type=int)
    parser.add_argument("--expected-source-pairs", required=True, type=int)
    return parser.parse_args()


def open_input(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="ascii", newline="")
    return path.open("rt", encoding="ascii", newline="")


def read_record(handle: TextIO, label: str, number: int) -> tuple[str, str, str, str] | None:
    lines = tuple(handle.readline() for _ in range(4))
    if all(line == "" for line in lines):
        return None
    if any(line == "" for line in lines):
        raise ValueError(f"{label}: truncated FASTQ record at pair {number:,}")
    header, sequence, plus, quality = (line.rstrip("\r\n") for line in lines)
    if not header.startswith("@"):
        raise ValueError(f"{label}: invalid header at pair {number:,}: {header[:80]!r}")
    if not plus.startswith("+"):
        raise ValueError(f"{label}: invalid plus line at pair {number:,}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"{label}: sequence/quality length mismatch at pair {number:,}: "
            f"{len(sequence)} != {len(quality)}"
        )
    return header, sequence, plus, quality


def canonical_id(header: str) -> str:
    token = header[1:].split(None, 1)[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def parse_tags(header: str, label: str, number: int) -> dict[str, str]:
    fields = header.split()
    tags: dict[str, str] = {}
    for field in fields[1:]:
        parts = field.split(":", 2)
        if len(parts) != 3:
            raise ValueError(f"{label}: malformed FASTQ comment field at pair {number}: {field}")
        key, tag_type, value = parts
        if key in tags:
            raise ValueError(f"{label}: duplicate FASTQ tag {key} at pair {number}")
        if key in REQUIRED:
            expected_type = "i" if key == "XD" else "Z"
            if tag_type != expected_type:
                raise ValueError(
                    f"{label}: {key} type is {tag_type}, expected {expected_type}, "
                    f"at pair {number}"
                )
            if not TAG_PATTERNS[key].fullmatch(value):
                raise ValueError(f"{label}: invalid {key} value at pair {number}: {value}")
            tags[key] = value
    missing = [key for key in REQUIRED if key not in tags]
    if missing:
        raise ValueError(f"{label}: missing tags {missing} at pair {number}")
    return tags


def format_header(header: str, tags: dict[str, str]) -> str:
    read_id = header[1:].split(None, 1)[0]
    comment = "\t".join(
        (
            f"CR:Z:{tags['CR']}",
            f"CB:Z:{tags['CB']}",
            f"CN:Z:{tags['CN']}",
            f"XD:i:{tags['XD']}",
            f"ZE:Z:{tags['XE']}",
        )
    )
    return f"@{read_id} {comment}"


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_name, path)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    args = parse_args()
    if args.max_pairs < 1 or args.expected_source_pairs < args.max_pairs:
        raise ValueError(
            "max-pairs must be positive and no larger than expected-source-pairs"
        )
    for path in (args.r1, args.r2):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty input FASTQ: {path}")
    for path in (args.out_r1, args.out_r2):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            raise FileExistsError(f"refusing to overwrite probe FASTQ: {path}")

    count = 0
    tag_lengths: dict[str, set[int]] = {key: set() for key in REQUIRED}
    try:
        with (
            open_input(args.r1) as r1_handle,
            open_input(args.r2) as r2_handle,
            args.out_r1.open("wt", encoding="ascii", newline="\n") as out1,
            args.out_r2.open("wt", encoding="ascii", newline="\n") as out2,
        ):
            while count < args.max_pairs:
                number = count + 1
                record1 = read_record(r1_handle, "R1", number)
                record2 = read_record(r2_handle, "R2", number)
                if record1 is None or record2 is None:
                    raise ValueError(
                        f"input ended before requested probe size at pair {number:,}"
                    )
                if canonical_id(record1[0]) != canonical_id(record2[0]):
                    raise ValueError(
                        f"R1/R2 read-ID mismatch at pair {number:,}: "
                        f"{canonical_id(record1[0])!r} != {canonical_id(record2[0])!r}"
                    )
                tags1 = parse_tags(record1[0], "R1", number)
                tags2 = parse_tags(record2[0], "R2", number)
                if tags1 != tags2:
                    raise ValueError(f"R1/R2 tag mismatch at pair {number:,}")
                for key, value in tags1.items():
                    tag_lengths[key].add(len(value))

                out1.write(
                    f"{format_header(record1[0], tags1)}\n"
                    f"{record1[1]}\n+\n{record1[3]}\n"
                )
                out2.write(
                    f"{format_header(record2[0], tags2)}\n"
                    f"{record2[1]}\n+\n{record2[3]}\n"
                )
                count += 1
    except BaseException:
        for path in (args.out_r1, args.out_r2):
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        raise

    payload = {
        "input_r1": str(args.r1),
        "input_r2": str(args.r2),
        "expected_source_pairs": args.expected_source_pairs,
        "written_probe_pairs": count,
        "source_to_bam_tag_rename": "XE->ZE",
        "observed_source_tag_lengths": {
            key: sorted(values) for key, values in tag_lengths.items()
        },
    }
    atomic_json(args.stats_out, payload)
    print(f"[OK] wrote {count:,} synchronized DNA probe read pairs")
    print("[OK] source FASTQ XE mismatch-vector tag will enter BAM as ZE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
