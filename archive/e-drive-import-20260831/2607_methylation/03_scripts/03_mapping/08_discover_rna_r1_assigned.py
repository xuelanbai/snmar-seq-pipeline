#!/usr/bin/env python3
"""Step 08: discover strict10 RNA R1 FASTQ files for R1-only processing."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


R1_SUFFIX = ".assigned_R1.fq.gz"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--expected-files", type=int, default=96)
    args = parser.parse_args()

    if args.expected_files < 0:
        parser.error("--expected-files must be non-negative")
    if not args.split_root.is_dir():
        raise NotADirectoryError(args.split_root)

    rows: list[tuple[str, str, Path]] = []
    seen_ids: set[str] = set()
    for r1 in sorted(args.split_root.rglob(f"*{R1_SUFFIX}")):
        if not r1.is_file() or r1.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty strict RNA R1 FASTQ: {r1}")
        pair_id = r1.name.removesuffix(R1_SUFFIX)
        if pair_id in seen_ids:
            raise ValueError(
                f"duplicate pair_id across RNA strict10 tree: {pair_id}"
            )
        seen_ids.add(pair_id)
        rows.append((pair_id, r1.parent.name, r1.resolve()))

    if not rows:
        raise FileNotFoundError(
            f"no *{R1_SUFFIX} files found below {args.split_root}"
        )
    if args.expected_files and len(rows) != args.expected_files:
        raise ValueError(
            f"found {len(rows)} strict RNA R1 files, "
            f"expected {args.expected_files}"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pair_id", "sample", "R1"))
        writer.writerows(rows)

    print(f"[OK] strict RNA R1 files={len(rows)}; table={args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
