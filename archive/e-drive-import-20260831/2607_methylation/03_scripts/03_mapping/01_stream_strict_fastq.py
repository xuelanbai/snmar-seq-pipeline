#!/usr/bin/env python3
"""Stream XD=0 FASTQ records and normalize comments to SAM-compatible tags."""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from typing import TextIO


TAG_RE = re.compile(
    r"(?:^|[ \t])"
    r"(?P<tag>(?:CR|CB|CN|UR|UB|XD|XE|XL|XI):[AifZHB]:[^ \t\r\n]+)"
)
TAG_ORDER = ("CR", "CB", "CN", "UR", "UB", "XD", "XE", "XL", "XI")


def open_fastq(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8")


def normalize_header(header: str) -> tuple[str, bool]:
    tags_by_name: dict[str, str] = {}
    for match in TAG_RE.finditer(header):
        tag = match.group("tag")
        tags_by_name[tag[:2]] = tag
    strict = tags_by_name.get("XD") == "XD:i:0"
    read_id = header.split(maxsplit=1)[0]
    ordered = [tags_by_name[name] for name in TAG_ORDER if name in tags_by_name]
    if not ordered:
        return read_id, False
    return f"{read_id} " + "\t".join(ordered), strict


def main() -> int:
    if len(sys.argv) != 2:
        print(
            "Usage: python3 01_stream_strict_fastq.py INPUT.fastq[.gz]",
            file=sys.stderr,
        )
        return 2
    path = Path(sys.argv[1])
    total = tagged = written = 0
    with open_fastq(path) as handle:
        while True:
            header = handle.readline().rstrip()
            if not header:
                break
            sequence = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            quality = handle.readline().rstrip()
            total += 1
            if not all((sequence, plus, quality)):
                raise ValueError(f"{path}: truncated FASTQ record near {total}")
            if len(sequence) != len(quality):
                raise ValueError(f"{path}: sequence/quality mismatch near {total}")
            normalized, strict = normalize_header(header)
            if "XD:i:" in header:
                tagged += 1
            if not strict:
                continue
            sys.stdout.write(f"{normalized}\n{sequence}\n{plus}\n{quality}\n")
            written += 1

    if total and tagged != total:
        raise ValueError(
            f"{path}: XD tag missing from {total - tagged}/{total} record(s)"
        )
    print(
        f"[STRICT_STREAM] file={path} total={total} written_XD0={written}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
