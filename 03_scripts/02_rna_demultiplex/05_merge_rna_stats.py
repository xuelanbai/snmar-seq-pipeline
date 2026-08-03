#!/usr/bin/env python3
"""Step 05: merge RNA per-pair long-format statistics after completeness checks."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


UNKNOWN_CATEGORIES = (
    "R2_first8",
    "T5_at_expected_position",
    "header_index",
    "I7_prefix8_candidate_rejected",
    "linker_by_T7",
)


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
        pair_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_pair_ids = {row["pair_id"] for row in pair_rows}
    stat_by_id = {
        path.name.removesuffix(".stat.tsv"): path
        for path in args.stat_dir.glob("*.stat.tsv")
        if path.name.removesuffix(".stat.tsv") in expected_pair_ids
    }
    missing_stats = sorted(expected_pair_ids - set(stat_by_id))
    if missing_stats:
        raise SystemExit(
            f"[ERROR] missing {len(missing_stats)} per-pair stat file(s); "
            f"first={missing_stats[:10]}"
        )

    all_rows: list[dict[str, str]] = []
    totals: Counter[str] = Counter()
    metric_order: list[str] = []
    for pair_id in sorted(expected_pair_ids):
        with stat_by_id[pair_id].open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        all_rows.extend(rows)
        for row in rows:
            metric = row["metric"]
            if metric not in metric_order:
                metric_order.append(metric)
            totals[metric] += int(row["count"])

    total_pairs = totals["total_pairs"] or 1
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "count", "percent"))
        for metric in metric_order:
            count = totals[metric]
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

    unknown_totals = {category: Counter() for category in UNKNOWN_CATEGORIES}
    missing_unknown: list[str] = []
    for pair_id in sorted(expected_pair_ids):
        path = args.unknown_dir / f"{pair_id}.unknown.json"
        if not path.is_file():
            missing_unknown.append(pair_id)
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        for category in UNKNOWN_CATEGORIES:
            unknown_totals[category].update(
                {
                    sequence: int(count)
                    for sequence, count in payload.get(category, [])
                }
            )
    if missing_unknown:
        raise SystemExit(
            f"[ERROR] missing {len(missing_unknown)} unknown JSON file(s); "
            f"first={missing_unknown[:10]}"
        )

    with args.unknown_summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("category", "sequence", "count"))
        for category, counter in unknown_totals.items():
            for sequence, count in counter.most_common(100):
                writer.writerow((category, sequence, count))

    strict = totals["strict_assigned"]
    print(f"[RESULT] total_pairs={totals['total_pairs']}")
    print(f"[RESULT] strict_assigned={strict} ({strict / total_pairs * 100:.4f}%)")
    print(f"[INFO] completed_pairs={len(expected_pair_ids)}")
    print(f"[INFO] summary={args.summary}")


if __name__ == "__main__":
    main()
