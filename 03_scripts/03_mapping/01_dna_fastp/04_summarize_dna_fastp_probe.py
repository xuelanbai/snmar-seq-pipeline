#!/usr/bin/env python3
"""Step 04: merge DNA preparation/fastp metrics and verify strict tags."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly one data row")
    return rows[0]


def verify_headers(path: Path, max_records: int) -> int:
    checked = 0
    with gzip.open(path, "rt") as handle:
        while checked < max_records:
            header = handle.readline().rstrip()
            if not header:
                break
            sequence = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            quality = handle.readline().rstrip()
            if not all((sequence, plus, quality)):
                raise ValueError(f"{path}: truncated FASTQ near record {checked + 1}")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: sequence/quality mismatch near record {checked + 1}"
                )
            for token in ("CR:Z:", "CB:Z:", "CN:Z:", "XD:i:0", "XE:Z:"):
                if token not in header:
                    raise ValueError(
                        f"{path}: header missing {token}: {header[:200]}"
                    )
            checked += 1
    return checked


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe-table", required=True, type=Path)
    parser.add_argument("--stat-root", required=True, type=Path)
    parser.add_argument("--clean-root", required=True, type=Path)
    parser.add_argument("--json-root", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--per-pair", required=True, type=Path)
    parser.add_argument("--check-headers", type=int, default=1000)
    args = parser.parse_args()

    with args.probe_table.open(encoding="utf-8", newline="") as handle:
        probe_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not probe_rows:
        raise ValueError(f"empty probe table: {args.probe_table}")

    output_rows: list[dict[str, object]] = []
    totals = {
        "input_pairs_to_prepare": 0,
        "XD0_pairs": 0,
        "front_trim_output_pairs": 0,
        "fastp_before_reads": 0,
        "fastp_after_reads": 0,
        "fastp_after_q20_bases": 0,
        "fastp_after_q30_bases": 0,
        "fastp_after_total_bases": 0,
    }

    for pair in probe_rows:
        pair_id = pair["pair_id"]
        sample = pair["sample"]
        stat_path = args.stat_root / sample / f"{pair_id}.prepare.tsv"
        json_path = args.json_root / sample / f"{pair_id}.fastp.json"
        clean_r1 = args.clean_root / sample / f"{pair_id}.clean_R1.fq.gz"
        clean_r2 = args.clean_root / sample / f"{pair_id}.clean_R2.fq.gz"
        for path in (stat_path, json_path, clean_r1, clean_r2):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"missing/empty probe output: {path}")

        prepared = read_single_row(stat_path)
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        before = payload["summary"]["before_filtering"]
        after = payload["summary"]["after_filtering"]
        checked_r1 = verify_headers(clean_r1, args.check_headers)
        checked_r2 = verify_headers(clean_r2, args.check_headers)
        if checked_r1 != checked_r2:
            raise ValueError(
                f"{pair_id}: checked R1/R2 record counts differ "
                f"({checked_r1} vs {checked_r2})"
            )

        row = {
            "pair_id": pair_id,
            "sample": sample,
            "input_pairs_to_prepare": int(prepared["total_pairs"]),
            "XD0_pairs": int(prepared["XD0_pairs"]),
            "front_trim_output_pairs": int(prepared["output_pairs"]),
            "trim_front_R1": int(prepared["trim_front_R1"]),
            "trim_front_R2": int(prepared["trim_front_R2"]),
            "fastp_before_reads": int(before["total_reads"]),
            "fastp_after_reads": int(after["total_reads"]),
            "fastp_read_retention_percent": (
                int(after["total_reads"]) / max(int(before["total_reads"]), 1) * 100
            ),
            "fastp_after_q20_rate": float(after["q20_rate"]),
            "fastp_after_q30_rate": float(after["q30_rate"]),
            "checked_strict_headers_per_read": checked_r1,
        }
        output_rows.append(row)
        for key in (
            "input_pairs_to_prepare",
            "XD0_pairs",
            "front_trim_output_pairs",
            "fastp_before_reads",
            "fastp_after_reads",
        ):
            totals[key] += int(row[key])
        totals["fastp_after_q20_bases"] += int(after["q20_bases"])
        totals["fastp_after_q30_bases"] += int(after["q30_bases"])
        totals["fastp_after_total_bases"] += int(after["total_bases"])

    args.per_pair.parent.mkdir(parents=True, exist_ok=True)
    with args.per_pair.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(output_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    total_bases = max(totals["fastp_after_total_bases"], 1)
    summary_rows = (
        ("probe_pairs", len(output_rows), ""),
        (
            "XD0_percent_in_probe_input",
            totals["XD0_pairs"] / max(totals["input_pairs_to_prepare"], 1) * 100,
            "%",
        ),
        (
            "fastp_read_retention_percent",
            totals["fastp_after_reads"] / max(totals["fastp_before_reads"], 1) * 100,
            "%",
        ),
        (
            "fastp_after_Q20_percent",
            totals["fastp_after_q20_bases"] / total_bases * 100,
            "%",
        ),
        (
            "fastp_after_Q30_percent",
            totals["fastp_after_q30_bases"] / total_bases * 100,
            "%",
        ),
        (
            "strict_header_check_records_per_read_max",
            args.check_headers,
            "records",
        ),
    )
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value", "unit"))
        for metric, value, unit in summary_rows:
            if isinstance(value, float):
                value = f"{value:.6f}"
            writer.writerow((metric, value, unit))
    print(
        f"[OK] DNA fastp probe pairs={len(output_rows)}; "
        f"read retention="
        f"{totals['fastp_after_reads'] / max(totals['fastp_before_reads'], 1) * 100:.3f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
