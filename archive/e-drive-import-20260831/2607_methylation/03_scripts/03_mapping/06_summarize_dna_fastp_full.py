#!/usr/bin/env python3
"""Reconcile DNA XD0 demultiplex, paired pre-trimming, and fastp results."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path


REQUIRED_TAGS = ("CR:Z:", "CB:Z:", "CN:Z:", "XD:i:0", "XE:Z:")


def load_expected_xd0(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("metric") != "assigned_XD0":
                continue
            pair_id = row["pair_id"]
            if pair_id in result:
                raise ValueError(f"{path}: duplicate assigned_XD0 for {pair_id}")
            result[pair_id] = int(row["count"])
    if not result:
        raise ValueError(f"{path}: no assigned_XD0 rows")
    return result


def load_pretrim(path: Path, pair_id: str, sample: str) -> dict[str, int]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty pretrim stat: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected one row, found {len(rows)}")
    row = rows[0]
    if row.get("pair_id") != pair_id or row.get("sample") != sample:
        raise ValueError(f"{path}: pair_id/sample mismatch")
    fields = (
        "total_pairs",
        "XD0_pairs",
        "non_XD0_pairs",
        "too_short_after_front_trim",
        "output_pairs",
        "trim_front_R1",
        "trim_front_R2",
    )
    result: dict[str, int] = {}
    for field in fields:
        try:
            result[field] = int(row[field])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path}: invalid {field}={row.get(field)!r}") from exc
    return result


def canonical_id(header: str) -> str:
    return (
        header.split(maxsplit=1)[0]
        .removeprefix("@")
        .removesuffix("/1")
        .removesuffix("/2")
    )


def read_record(handle, path: Path, record_no: int):
    header = handle.readline().rstrip()
    if not header:
        return None
    sequence = handle.readline().rstrip()
    plus = handle.readline().rstrip()
    quality = handle.readline().rstrip()
    if not all((sequence, plus, quality)):
        raise ValueError(f"{path}: truncated FASTQ near record {record_no}")
    if not plus.startswith("+"):
        raise ValueError(f"{path}: invalid plus line near record {record_no}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"{path}: sequence/quality mismatch near record {record_no}"
        )
    for tag in REQUIRED_TAGS:
        if tag not in header:
            raise ValueError(f"{path}: header missing {tag}: {header[:200]}")
    return header, sequence, plus, quality


def verify_paired_headers(r1_path: Path, r2_path: Path, limit: int) -> int:
    checked = 0
    with gzip.open(r1_path, "rt") as r1, gzip.open(r2_path, "rt") as r2:
        while checked < limit:
            rec1 = read_record(r1, r1_path, checked + 1)
            rec2 = read_record(r2, r2_path, checked + 1)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise ValueError(f"R1/R2 counts differ near record {checked + 1}")
            if canonical_id(rec1[0]) != canonical_id(rec2[0]):
                raise ValueError(f"R1/R2 IDs differ near record {checked + 1}")
            checked += 1
    return checked


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--demux-per-pair", required=True, type=Path)
    parser.add_argument("--pretrim-stat-root", required=True, type=Path)
    parser.add_argument("--clean-root", required=True, type=Path)
    parser.add_argument("--json-root", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--per-pair", required=True, type=Path)
    parser.add_argument("--check-headers", type=int, default=1000)
    args = parser.parse_args()

    if args.check_headers < 0:
        raise ValueError("--check-headers must be non-negative")

    with args.pair_table.open(encoding="utf-8", newline="") as handle:
        pairs = list(csv.DictReader(handle, delimiter="\t"))
    if not pairs:
        raise ValueError(f"empty pair table: {args.pair_table}")
    pair_ids = [row["pair_id"] for row in pairs]
    if len(pair_ids) != len(set(pair_ids)):
        raise ValueError(f"{args.pair_table}: duplicate pair_id")

    expected_xd0 = load_expected_xd0(args.demux_per_pair)
    missing = sorted(set(pair_ids) - set(expected_xd0))
    if missing:
        raise ValueError(
            f"missing assigned_XD0 for {len(missing)} pair(s): {missing[:10]}"
        )

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    totals: Counter[str] = Counter()
    weighted_r1_length = 0.0
    weighted_r2_length = 0.0

    for pair in pairs:
        pair_id = pair["pair_id"]
        sample = pair["sample"]
        clean_r1 = args.clean_root / sample / f"{pair_id}.clean_R1.fq.gz"
        clean_r2 = args.clean_root / sample / f"{pair_id}.clean_R2.fq.gz"
        json_path = args.json_root / sample / f"{pair_id}.fastp.json"
        pretrim_path = (
            args.pretrim_stat_root / sample / f"{pair_id}.pretrim.tsv"
        )
        try:
            for path in (clean_r1, clean_r2, json_path):
                if not path.is_file() or path.stat().st_size == 0:
                    raise FileNotFoundError(f"missing/empty output: {path}")

            expected = expected_xd0[pair_id]
            pretrim = load_pretrim(pretrim_path, pair_id, sample)
            if pretrim["XD0_pairs"] != expected:
                raise ValueError(
                    f"pretrim XD0={pretrim['XD0_pairs']}, expected={expected}"
                )
            if pretrim["too_short_after_front_trim"] != 0:
                raise ValueError(
                    "pretrim too_short="
                    f"{pretrim['too_short_after_front_trim']}"
                )
            if pretrim["output_pairs"] != expected:
                raise ValueError(
                    f"pretrim output={pretrim['output_pairs']}, "
                    f"expected={expected}"
                )

            payload = json.loads(json_path.read_text(encoding="utf-8"))
            before = payload["summary"]["before_filtering"]
            after = payload["summary"]["after_filtering"]
            filtering = payload.get("filtering_result", {})
            adapter = payload.get("adapter_cutting", {})
            before_reads = int(before["total_reads"])
            after_reads = int(after["total_reads"])
            if before_reads != pretrim["output_pairs"] * 2:
                raise ValueError(
                    f"fastp before_reads={before_reads}, "
                    f"expected={pretrim['output_pairs'] * 2}"
                )
            if before_reads % 2 or after_reads % 2:
                raise ValueError("odd read count in paired fastp JSON")

            checked = verify_paired_headers(
                clean_r1, clean_r2, args.check_headers
            )
            output_pairs = after_reads // 2
            r1_mean = float(after.get("read1_mean_length", 0.0))
            r2_mean = float(after.get("read2_mean_length", 0.0))
            row: dict[str, object] = {
                "pair_id": pair_id,
                "sample": sample,
                "demux_XD0_pairs": expected,
                "pretrim_XD0_pairs": pretrim["XD0_pairs"],
                "pretrim_output_pairs": pretrim["output_pairs"],
                "pretrim_R1_front_trim": pretrim["trim_front_R1"],
                "pretrim_R2_front_trim": pretrim["trim_front_R2"],
                "fastp_input_pairs": before_reads // 2,
                "fastp_output_pairs": output_pairs,
                "fastp_read_retention_percent": (
                    after_reads / max(before_reads, 1) * 100
                ),
                "fastp_after_q20_rate": float(after["q20_rate"]),
                "fastp_after_q30_rate": float(after["q30_rate"]),
                "fastp_after_R1_mean_length": r1_mean,
                "fastp_after_R2_mean_length": r2_mean,
                "low_quality_reads": int(
                    filtering.get("low_quality_reads", 0)
                ),
                "too_many_N_reads": int(
                    filtering.get("too_many_N_reads", 0)
                ),
                "too_short_reads": int(filtering.get("too_short_reads", 0)),
                "adapter_trimmed_reads": int(
                    adapter.get("adapter_trimmed_reads", 0)
                ),
                "checked_strict_paired_headers": checked,
            }
            rows.append(row)

            totals["input_XD0_pairs"] += expected
            totals["pretrim_output_pairs"] += pretrim["output_pairs"]
            totals["before_reads"] += before_reads
            totals["after_reads"] += after_reads
            totals["fastp_output_pairs"] += output_pairs
            totals["after_total_bases"] += int(after["total_bases"])
            totals["after_q20_bases"] += int(after["q20_bases"])
            totals["after_q30_bases"] += int(after["q30_bases"])
            totals["low_quality_reads"] += int(row["low_quality_reads"])
            totals["too_many_N_reads"] += int(row["too_many_N_reads"])
            totals["too_short_reads"] += int(row["too_short_reads"])
            totals["adapter_trimmed_reads"] += int(
                row["adapter_trimmed_reads"]
            )
            weighted_r1_length += r1_mean * output_pairs
            weighted_r2_length += r2_mean * output_pairs
        except Exception as exc:
            errors.append(f"{pair_id}: {exc}")

    if errors:
        shown = "\n".join(f"  - {error}" for error in errors[:30])
        extra = len(errors) - min(len(errors), 30)
        suffix = f"\n  ... and {extra} more" if extra else ""
        raise ValueError(
            f"validation failed for {len(errors)} pair(s):\n{shown}{suffix}"
        )

    args.per_pair.parent.mkdir(parents=True, exist_ok=True)
    with args.per_pair.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    before_reads = max(totals["before_reads"], 1)
    after_bases = max(totals["after_total_bases"], 1)
    output_pairs = max(totals["fastp_output_pairs"], 1)
    summary_rows = (
        ("completed_pairs", len(rows), "FASTQ pairs"),
        ("input_XD0_pairs", totals["input_XD0_pairs"], "read pairs"),
        ("pretrim_output_pairs", totals["pretrim_output_pairs"], "read pairs"),
        ("fastp_input_pairs", totals["before_reads"] // 2, "read pairs"),
        ("fastp_output_pairs", totals["fastp_output_pairs"], "read pairs"),
        (
            "fastp_read_retention_percent",
            totals["after_reads"] / before_reads * 100,
            "%",
        ),
        (
            "fastp_after_Q20_percent",
            totals["after_q20_bases"] / after_bases * 100,
            "%",
        ),
        (
            "fastp_after_Q30_percent",
            totals["after_q30_bases"] / after_bases * 100,
            "%",
        ),
        (
            "fastp_after_R1_weighted_mean_length",
            weighted_r1_length / output_pairs,
            "bp",
        ),
        (
            "fastp_after_R2_weighted_mean_length",
            weighted_r2_length / output_pairs,
            "bp",
        ),
        (
            "low_quality_reads_percent",
            totals["low_quality_reads"] / before_reads * 100,
            "%",
        ),
        (
            "too_many_N_reads_percent",
            totals["too_many_N_reads"] / before_reads * 100,
            "%",
        ),
        (
            "too_short_reads_percent",
            totals["too_short_reads"] / before_reads * 100,
            "%",
        ),
        (
            "adapter_trimmed_reads_percent",
            totals["adapter_trimmed_reads"] / before_reads * 100,
            "%",
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
        f"[OK] full DNA fastp pairs={len(rows)}; "
        f"retention={totals['after_reads'] / before_reads * 100:.3f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
