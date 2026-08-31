#!/usr/bin/env python3
"""Step 14: validate RNA clean R1 files and select a mapping probe."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


OUTPUT_FIELDS = (
    "pair_id",
    "sample",
    "clean_R1",
    "expected_reads",
    "fastp_R1_retention_percent",
    "selection_reason",
)


def validate_token(label: str, value: str) -> None:
    if not value or any(character.isspace() for character in value):
        raise ValueError(f"{label} contains whitespace or is empty: {value!r}")


def atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with tmp.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=OUTPUT_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def select_probe(
    rows: list[dict[str, object]],
    requested: int,
) -> list[dict[str, object]]:
    by_retention = sorted(
        rows,
        key=lambda row: (
            float(row["fastp_R1_retention_percent"]),
            str(row["pair_id"]),
        ),
    )
    reasons: dict[str, list[str]] = {}
    selected: dict[str, dict[str, object]] = {}

    def add(row: dict[str, object], reason: str) -> None:
        pair_id = str(row["pair_id"])
        selected[pair_id] = row
        reasons.setdefault(pair_id, []).append(reason)

    retention_positions = (
        (0.00, "lowest_retention"),
        (0.25, "retention_q25"),
        (0.50, "median_retention"),
        (0.75, "retention_q75"),
        (1.00, "highest_retention"),
    )
    for fraction, reason in retention_positions:
        index = round((len(by_retention) - 1) * fraction)
        add(by_retention[index], reason)

    by_reads = sorted(
        rows,
        key=lambda row: (int(row["expected_reads"]), str(row["pair_id"])),
    )
    add(by_reads[0], "smallest_clean_file")
    add(by_reads[-1], "largest_clean_file")
    by_input = sorted(
        rows,
        key=lambda row: (
            int(row["fastp_input_R1_reads"]),
            str(row["pair_id"]),
        ),
    )
    add(by_input[-1], "largest_fastp_input")

    if len(selected) < requested:
        for index in range(requested):
            rank = round((len(by_retention) - 1) * index / max(requested - 1, 1))
            add(by_retention[rank], f"even_retention_rank_{index + 1}")
            if len(selected) >= requested:
                break
    if len(selected) < requested:
        for row in sorted(rows, key=lambda item: str(item["pair_id"])):
            add(row, "deterministic_fill")
            if len(selected) >= requested:
                break

    priority = {str(row["pair_id"]): index for index, row in enumerate(by_retention)}
    chosen_ids = sorted(
        selected,
        key=lambda pair_id: (
            priority[pair_id],
            pair_id,
        ),
    )
    if len(chosen_ids) > requested:
        mandatory = {
            str(by_retention[0]["pair_id"]),
            str(by_retention[-1]["pair_id"]),
            str(by_reads[-1]["pair_id"]),
            str(by_input[-1]["pair_id"]),
        }
        kept = [pair_id for pair_id in chosen_ids if pair_id in mandatory]
        for pair_id in chosen_ids:
            if pair_id not in kept:
                kept.append(pair_id)
            if len(kept) == requested:
                break
        chosen_ids = kept

    output: list[dict[str, object]] = []
    for pair_id in chosen_ids:
        row = dict(selected[pair_id])
        row["selection_reason"] = ",".join(dict.fromkeys(reasons[pair_id]))
        output.append(row)
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastp-per-file", type=Path, required=True)
    parser.add_argument("--clean-root", type=Path, required=True)
    parser.add_argument("--all-out", type=Path, required=True)
    parser.add_argument("--probe-out", type=Path, required=True)
    parser.add_argument("--expected-files", type=int, default=96)
    parser.add_argument("--probe-files", type=int, default=8)
    parser.add_argument("--expected-total-reads", type=int)
    args = parser.parse_args()

    if args.expected_files < 1:
        parser.error("--expected-files must be positive")
    if not 1 <= args.probe_files <= args.expected_files:
        parser.error("--probe-files must be within expected file count")
    if args.expected_total_reads is not None and args.expected_total_reads < 1:
        parser.error("--expected-total-reads must be positive")
    if not args.fastp_per_file.is_file():
        raise FileNotFoundError(args.fastp_per_file)

    with args.fastp_per_file.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(source_rows) != args.expected_files:
        raise ValueError(
            f"{args.fastp_per_file}: rows={len(source_rows)}, "
            f"expected={args.expected_files}"
        )

    seen: set[str] = set()
    rows: list[dict[str, object]] = []
    required = {
        "pair_id",
        "sample",
        "fastp_input_R1_reads",
        "fastp_output_R1_reads",
        "fastp_R1_retention_percent",
        "checked_strict_R1_headers",
    }
    missing_columns = required - set(source_rows[0])
    if missing_columns:
        raise ValueError(
            f"{args.fastp_per_file}: missing columns {sorted(missing_columns)}"
        )

    for source in source_rows:
        pair_id = source["pair_id"]
        sample = source["sample"]
        validate_token("pair_id", pair_id)
        validate_token("sample", sample)
        if pair_id in seen:
            raise ValueError(f"duplicate pair_id: {pair_id}")
        seen.add(pair_id)

        expected_reads = int(source["fastp_output_R1_reads"])
        input_reads = int(source["fastp_input_R1_reads"])
        checked_headers = int(source["checked_strict_R1_headers"])
        retention = float(source["fastp_R1_retention_percent"])
        if expected_reads < 1 or input_reads < expected_reads:
            raise ValueError(f"{pair_id}: invalid fastp input/output counts")
        if checked_headers < min(expected_reads, 1000):
            raise ValueError(
                f"{pair_id}: strict header checks={checked_headers}, "
                f"expected at least {min(expected_reads, 1000)}"
            )
        clean_r1 = (
            args.clean_root / sample / f"{pair_id}.clean_R1.fq.gz"
        )
        validate_token("clean_R1", str(clean_r1))
        if not clean_r1.is_file() or clean_r1.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty clean R1: {clean_r1}")

        rows.append(
            {
                "pair_id": pair_id,
                "sample": sample,
                "clean_R1": str(clean_r1),
                "expected_reads": expected_reads,
                "fastp_R1_retention_percent": f"{retention:.6f}",
                "selection_reason": "all",
                "fastp_input_R1_reads": input_reads,
            }
        )

    rows.sort(key=lambda row: (str(row["sample"]), str(row["pair_id"])))
    total_reads = sum(int(row["expected_reads"]) for row in rows)
    if (
        args.expected_total_reads is not None
        and total_reads != args.expected_total_reads
    ):
        raise ValueError(
            f"clean R1 total={total_reads}, "
            f"expected={args.expected_total_reads}"
        )
    probe = select_probe(rows, args.probe_files)
    all_rows = [
        {field: row[field] for field in OUTPUT_FIELDS} for row in rows
    ]
    probe_rows = [
        {field: row[field] for field in OUTPUT_FIELDS} for row in probe
    ]
    atomic_write(args.all_out, all_rows)
    atomic_write(args.probe_out, probe_rows)
    print(
        f"[OK] RNA STAR inputs={len(all_rows)}; "
        f"probe={len(probe_rows)}; clean_reads="
        f"{total_reads}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
