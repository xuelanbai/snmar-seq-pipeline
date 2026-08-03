#!/usr/bin/env python3
"""Step 05: merge per-pair long-format barcode QC statistics."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def read_stat(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stat-dir", required=True, type=Path)
    parser.add_argument("--unknown-dir", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--per-pair", required=True, type=Path)
    parser.add_argument("--unknown-summary", required=True, type=Path)
    parser.add_argument("--pair-table", required=True, type=Path)
    args = parser.parse_args()

    with args.pair_table.open(encoding="utf-8", newline="") as handle:
        expected_pair_ids = {
            row["pair_id"] for row in csv.DictReader(handle, delimiter="\t")
        }
    stat_files = sorted(
        path
        for path in args.stat_dir.glob("*.stat.tsv")
        if path.name.removesuffix(".stat.tsv") in expected_pair_ids
    )
    if not stat_files:
        raise SystemExit(f"[ERROR] no .stat.tsv files in {args.stat_dir}")

    all_rows: list[dict[str, str]] = []
    totals: Counter[str] = Counter()
    pair_totals: dict[str, int] = {}
    for path in stat_files:
        rows = read_stat(path)
        all_rows.extend(rows)
        for row in rows:
            count = int(row["count"])
            totals[row["metric"]] += count
            if row["metric"] == "total_pairs":
                pair_totals[row["pair_id"]] = count

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    total_pairs = totals["total_pairs"] or 1
    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "count", "percent"))
        for metric, count in totals.items():
            writer.writerow((metric, count, f"{count / total_pairs * 100:.6f}"))

    with args.per_pair.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("pair_id", "sample", "metric", "count", "percent"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(all_rows)

    unknown_totals = {
        "R1_first8": Counter(),
        "R2_first8": Counter(),
        "header_index": Counter(),
    }
    for path in sorted(
        path
        for path in args.unknown_dir.glob("*.unknown.json")
        if path.name.removesuffix(".unknown.json") in expected_pair_ids
    ):
        payload = json.loads(path.read_text(encoding="utf-8"))
        for category in unknown_totals:
            unknown_totals[category].update(
                {sequence: int(count) for sequence, count in payload.get(category, [])}
            )

    with args.unknown_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("category", "sequence", "count"))
        for category, counts in unknown_totals.items():
            for sequence, count in counts.most_common(100):
                writer.writerow((category, sequence, count))

    strict = totals["strict_assigned"]
    rescued = totals["rescued_assigned"]
    print(f"[RESULT] total_pairs={totals['total_pairs']}")
    print(f"[RESULT] strict_assigned={strict} ({strict / total_pairs * 100:.4f}%)")
    print(f"[RESULT] rescued_assigned={rescued} ({rescued / total_pairs * 100:.4f}%)")
    print(
        "[RESULT] final_assigned="
        f"{strict + rescued} ({(strict + rescued) / total_pairs * 100:.4f}%)"
    )
    print(f"[INFO] summary: {args.summary}")


if __name__ == "__main__":
    main()
