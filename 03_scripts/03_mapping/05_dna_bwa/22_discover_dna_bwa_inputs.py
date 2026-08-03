#!/usr/bin/env python3
"""Validate DNA clean FASTQs and write deterministic BWA input manifests."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastp-per-pair", required=True, type=Path)
    parser.add_argument("--clean-root", required=True, type=Path)
    parser.add_argument("--pair-out", required=True, type=Path)
    parser.add_argument("--sample-out", required=True, type=Path)
    parser.add_argument("--probe-out", required=True, type=Path)
    parser.add_argument("--expected-pairs-files", type=int, default=768)
    parser.add_argument("--expected-samples", type=int, default=96)
    parser.add_argument("--expected-lanes-per-sample", type=int, default=8)
    parser.add_argument("--expected-total-clean-pairs", type=int, default=16_947_758_776)
    parser.add_argument("--probe-files-per-group", type=int, default=4)
    parser.add_argument("--probe-read-pairs", type=int, default=1_000_000)
    return parser.parse_args()


def atomic_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, path)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def as_int(row: dict[str, str], key: str) -> int:
    try:
        return int(row[key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"invalid integer in column {key!r}: {row.get(key)!r}") from exc


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"invalid float in column {key!r}: {row.get(key)!r}") from exc


def group_for_sample(sample: str) -> str:
    if "-GM12878-" in sample:
        return "GM12878"
    if "-cellMix-" in sample:
        return "cellMix"
    raise ValueError(f"cannot determine biological group from sample: {sample}")


def choose_probe(group_rows: list[dict[str, object]], count: int) -> list[tuple[dict[str, object], str]]:
    """Choose size/retention extremes, then fill evenly if any choices collide."""
    if len(group_rows) < count:
        raise ValueError(f"group has only {len(group_rows)} rows; cannot choose {count}")
    candidates = [
        (min(group_rows, key=lambda x: int(x["expected_clean_pairs"])), "smallest_clean"),
        (max(group_rows, key=lambda x: int(x["expected_clean_pairs"])), "largest_clean"),
        (min(group_rows, key=lambda x: float(x["fastp_retention_percent"])), "lowest_retention"),
        (max(group_rows, key=lambda x: float(x["fastp_retention_percent"])), "highest_retention"),
    ]
    chosen: list[tuple[dict[str, object], str]] = []
    used: set[str] = set()
    for row, reason in candidates:
        pair_id = str(row["pair_id"])
        if pair_id not in used:
            chosen.append((row, reason))
            used.add(pair_id)
        if len(chosen) == count:
            return chosen

    ordered = sorted(group_rows, key=lambda x: int(x["expected_clean_pairs"]))
    if count == 1:
        positions = [len(ordered) // 2]
    else:
        positions = [round(i * (len(ordered) - 1) / (count - 1)) for i in range(count)]
    for position in positions + list(range(len(ordered))):
        row = ordered[position]
        pair_id = str(row["pair_id"])
        if pair_id not in used:
            chosen.append((row, f"size_rank_{position + 1}"))
            used.add(pair_id)
        if len(chosen) == count:
            return chosen
    raise AssertionError("probe selection failed")


def main() -> int:
    args = parse_args()
    if min(
        args.expected_pairs_files,
        args.expected_samples,
        args.expected_lanes_per_sample,
        args.expected_total_clean_pairs,
        args.probe_files_per_group,
        args.probe_read_pairs,
    ) < 1:
        raise ValueError("all count arguments must be positive")

    with args.fastp_per_pair.open(newline="", encoding="utf-8-sig") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(source_rows) != args.expected_pairs_files:
        raise ValueError(
            f"expected {args.expected_pairs_files} fastp rows, observed {len(source_rows)}"
        )

    required = {
        "pair_id",
        "sample",
        "fastp_output_pairs",
        "fastp_read_retention_percent",
        "checked_strict_paired_headers",
    }
    missing = required.difference(source_rows[0] if source_rows else {})
    if missing:
        raise ValueError(f"fastp table is missing columns: {sorted(missing)}")

    pair_rows: list[dict[str, object]] = []
    seen_pair_ids: set[str] = set()
    for row in source_rows:
        pair_id = row["pair_id"].strip()
        sample = row["sample"].strip()
        if not pair_id or not sample:
            raise ValueError("empty pair_id or sample")
        if pair_id in seen_pair_ids:
            raise ValueError(f"duplicate pair_id: {pair_id}")
        seen_pair_ids.add(pair_id)
        if as_int(row, "checked_strict_paired_headers") < 1:
            raise ValueError(f"strict paired-header check failed: {pair_id}")

        r1 = args.clean_root / sample / f"{pair_id}.clean_R1.fq.gz"
        r2 = args.clean_root / sample / f"{pair_id}.clean_R2.fq.gz"
        for path in (r1, r2):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"missing/empty clean FASTQ: {path}")
        if "," in str(r1) or "," in str(r2):
            raise ValueError(f"commas are not supported in FASTQ paths: {pair_id}")

        pair_rows.append(
            {
                "pair_index": 0,
                "pair_id": pair_id,
                "sample": sample,
                "group": group_for_sample(sample),
                "r1": str(r1),
                "r2": str(r2),
                "expected_clean_pairs": as_int(row, "fastp_output_pairs"),
                "fastp_retention_percent": (
                    f"{as_float(row, 'fastp_read_retention_percent'):.6f}"
                ),
                "r1_bytes": r1.stat().st_size,
                "r2_bytes": r2.stat().st_size,
            }
        )

    pair_rows.sort(key=lambda x: (str(x["sample"]), str(x["pair_id"])))
    for index, row in enumerate(pair_rows, start=1):
        row["pair_index"] = index

    by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        by_sample[str(row["sample"])].append(row)
    if len(by_sample) != args.expected_samples:
        raise ValueError(f"expected {args.expected_samples} samples, observed {len(by_sample)}")

    sample_rows: list[dict[str, object]] = []
    for sample_index, sample in enumerate(sorted(by_sample), start=1):
        rows = by_sample[sample]
        if len(rows) != args.expected_lanes_per_sample:
            raise ValueError(
                f"{sample}: expected {args.expected_lanes_per_sample} lane pairs, "
                f"observed {len(rows)}"
            )
        groups = {str(x["group"]) for x in rows}
        if len(groups) != 1:
            raise ValueError(f"{sample}: inconsistent groups: {sorted(groups)}")
        sample_rows.append(
            {
                "sample_index": sample_index,
                "sample": sample,
                "group": next(iter(groups)),
                "lane_pairs": len(rows),
                "expected_clean_pairs": sum(int(x["expected_clean_pairs"]) for x in rows),
                "pair_ids": ",".join(str(x["pair_id"]) for x in rows),
                "r1_files": ",".join(str(x["r1"]) for x in rows),
                "r2_files": ",".join(str(x["r2"]) for x in rows),
            }
        )

    total = sum(int(x["expected_clean_pairs"]) for x in pair_rows)
    if total != args.expected_total_clean_pairs:
        raise ValueError(
            f"expected {args.expected_total_clean_pairs:,} total clean pairs, observed {total:,}"
        )

    probe_rows: list[dict[str, object]] = []
    probe_index = 0
    for group in ("GM12878", "cellMix"):
        group_rows = [x for x in pair_rows if x["group"] == group]
        for row, reason in choose_probe(group_rows, args.probe_files_per_group):
            probe_index += 1
            probe_rows.append(
                {
                    "probe_index": probe_index,
                    "run_id": f"probe_{probe_index:02d}_{row['pair_id']}",
                    "pair_id": row["pair_id"],
                    "sample": row["sample"],
                    "group": row["group"],
                    "r1": row["r1"],
                    "r2": row["r2"],
                    "source_clean_pairs": row["expected_clean_pairs"],
                    "expected_pairs": min(
                        args.probe_read_pairs, int(row["expected_clean_pairs"])
                    ),
                    "selection_reason": reason,
                }
            )

    atomic_tsv(args.pair_out, list(pair_rows[0]), pair_rows)
    atomic_tsv(args.sample_out, list(sample_rows[0]), sample_rows)
    atomic_tsv(args.probe_out, list(probe_rows[0]), probe_rows)
    group_counts = defaultdict(int)
    for row in sample_rows:
        group_counts[str(row["group"])] += 1
    print(
        f"[OK] DNA BWA manifests: {len(pair_rows)} lane pairs, "
        f"{len(sample_rows)} samples, {total:,} clean read pairs"
    )
    print(
        "[OK] sample groups: "
        + ", ".join(f"{key}={group_counts[key]}" for key in sorted(group_counts))
    )
    print(f"[OK] probe lane pairs: {len(probe_rows)} x <= {args.probe_read_pairs:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
