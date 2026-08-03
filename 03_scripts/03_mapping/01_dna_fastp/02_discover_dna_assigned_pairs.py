#!/usr/bin/env python3
"""Step 02: discover DNA assigned FASTQ pairs and choose an even probe subset."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


R1_SUFFIX = ".assigned_R1.fq.gz"
R2_SUFFIX = ".assigned_R2.fq.gz"


def write_table(path: Path, rows: list[tuple[str, str, Path, Path]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pair_id", "sample", "R1", "R2"))
        for pair_id, sample, r1, r2 in rows:
            writer.writerow((pair_id, sample, r1, r2))


def evenly_spaced(
    rows: list[tuple[str, str, Path, Path]],
    count: int,
) -> list[tuple[str, str, Path, Path]]:
    if count <= 0:
        raise ValueError("--probe-count must be positive")
    if count >= len(rows):
        return rows
    if count == 1:
        return [rows[len(rows) // 2]]
    positions = {
        round(index * (len(rows) - 1) / (count - 1))
        for index in range(count)
    }
    if len(positions) != count:
        raise AssertionError("even probe selection produced duplicate positions")
    return [rows[position] for position in sorted(positions)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--probe-out", type=Path)
    parser.add_argument("--probe-count", type=int, default=16)
    parser.add_argument("--expected-pairs", type=int, default=768)
    args = parser.parse_args()

    rows: list[tuple[str, str, Path, Path]] = []
    seen_ids: set[str] = set()
    for r1 in sorted(args.split_root.rglob(f"*{R1_SUFFIX}")):
        pair_id = r1.name.removesuffix(R1_SUFFIX)
        r2 = r1.with_name(pair_id + R2_SUFFIX)
        if not r2.is_file():
            raise FileNotFoundError(f"R2 missing for {r1}: expected {r2}")
        if pair_id in seen_ids:
            raise ValueError(f"duplicate pair_id across split tree: {pair_id}")
        seen_ids.add(pair_id)
        rows.append((pair_id, r1.parent.name, r1.resolve(), r2.resolve()))

    if not rows:
        raise FileNotFoundError(
            f"no *{R1_SUFFIX} files found below {args.split_root}"
        )
    if args.expected_pairs and len(rows) != args.expected_pairs:
        raise ValueError(
            f"found {len(rows)} DNA pairs, expected {args.expected_pairs}"
        )

    write_table(args.out, rows)
    if args.probe_out is not None:
        probe_rows = evenly_spaced(rows, args.probe_count)
        write_table(args.probe_out, probe_rows)
        print(
            f"[OK] DNA pairs={len(rows)}; "
            f"evenly_spaced_probe_pairs={len(probe_rows)}"
        )
    else:
        print(f"[OK] DNA pairs={len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
