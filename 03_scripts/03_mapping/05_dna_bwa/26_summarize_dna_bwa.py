#!/usr/bin/env python3
"""Combine validated per-run DNA BWA JSON files into QC tables."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("probe", "full"), required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--validation-root", required=True, type=Path)
    parser.add_argument("--summary-out", required=True, type=Path)
    parser.add_argument("--per-run-out", required=True, type=Path)
    parser.add_argument("--expected-runs", required=True, type=int)
    parser.add_argument("--expected-total-pairs", type=int)
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


def percent(numerator: int, denominator: int) -> str:
    return f"{100.0 * numerator / denominator:.6f}" if denominator else "0.000000"


def main() -> int:
    args = parse_args()
    if args.expected_runs < 1:
        raise ValueError("expected-runs must be positive")
    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if len(manifest) != args.expected_runs:
        raise ValueError(
            f"expected {args.expected_runs} manifest rows, observed {len(manifest)}"
        )

    if args.mode == "probe":
        needed = {"run_id", "sample", "group", "expected_pairs", "selection_reason"}
        run_specs = [
            {
                "run_id": row["run_id"],
                "sample": row["sample"],
                "group": row["group"],
                "expected_pairs": int(row["expected_pairs"]),
                "selection_reason": row["selection_reason"],
            }
            for row in manifest
        ]
    else:
        needed = {"sample", "group", "expected_clean_pairs", "lane_pairs"}
        run_specs = [
            {
                "run_id": row["sample"],
                "sample": row["sample"],
                "group": row["group"],
                "expected_pairs": int(row["expected_clean_pairs"]),
                "selection_reason": "full_sample_8_lanes",
            }
            for row in manifest
        ]
    missing = needed.difference(manifest[0] if manifest else {})
    if missing:
        raise ValueError(f"manifest is missing columns: {sorted(missing)}")

    per_run: list[dict[str, object]] = []
    missing_json: list[str] = []
    totals = {
        "expected_read_pairs": 0,
        "primary_records": 0,
        "primary_mapped": 0,
        "properly_paired": 0,
        "primary_duplicates": 0,
        "human_primary_mapped": 0,
        "pig_primary_mapped": 0,
        "bam_bytes": 0,
    }
    for spec in run_specs:
        run_id = str(spec["run_id"])
        json_path = args.validation_root / f"{run_id}.validation.json"
        if not json_path.is_file() or json_path.stat().st_size == 0:
            missing_json.append(str(json_path))
            continue
        with json_path.open(encoding="utf-8") as handle:
            data = json.load(handle)
        if data.get("status") != "validated":
            raise ValueError(f"{run_id}: validation status is not 'validated'")
        for key in ("run_id", "sample"):
            if str(data.get(key)) != str(spec[key]):
                raise ValueError(
                    f"{run_id}: {key} mismatch: {data.get(key)!r} != {spec[key]!r}"
                )
        if int(data["expected_read_pairs"]) != int(spec["expected_pairs"]):
            raise ValueError(f"{run_id}: expected-pair count differs from manifest")

        row = {
            "mode": args.mode,
            "run_id": run_id,
            "sample": spec["sample"],
            "group": spec["group"],
            "selection_reason": spec["selection_reason"],
            "expected_read_pairs": int(data["expected_read_pairs"]),
            "primary_records": int(data["primary_records"]),
            "primary_mapped": int(data["primary_mapped"]),
            "primary_mapping_percent": f"{float(data['primary_mapping_percent']):.6f}",
            "properly_paired": int(data["properly_paired"]),
            "properly_paired_percent": f"{float(data['properly_paired_percent']):.6f}",
            "primary_duplicates": int(data["primary_duplicates"]),
            "primary_duplicate_percent": f"{float(data['primary_duplicate_percent']):.6f}",
            "human_primary_mapped": int(data["human_primary_mapped"]),
            "pig_primary_mapped": int(data["pig_primary_mapped"]),
            "human_percent_of_primary_mapped": (
                f"{float(data['human_percent_of_primary_mapped']):.6f}"
            ),
            "pig_percent_of_primary_mapped": (
                f"{float(data['pig_percent_of_primary_mapped']):.6f}"
            ),
            "bam_bytes": int(data["bam_bytes"]),
            "bam": data["bam"],
            "validation_json": str(json_path),
            "status": data["status"],
        }
        per_run.append(row)
        for key in totals:
            totals[key] += int(row[key])

    if missing_json:
        preview = "\n".join(f"  {path}" for path in missing_json[:20])
        raise FileNotFoundError(
            f"{len(missing_json)} validation JSON files are missing; mapping is not complete:\n"
            f"{preview}"
        )
    if len(per_run) != args.expected_runs:
        raise AssertionError("validated run count does not match expected-runs")
    if (
        args.expected_total_pairs is not None
        and totals["expected_read_pairs"] != args.expected_total_pairs
    ):
        raise ValueError(
            f"expected {args.expected_total_pairs:,} total input pairs, "
            f"observed {totals['expected_read_pairs']:,}"
        )

    summary_metrics: list[tuple[str, object]] = [
        ("mode", args.mode),
        ("status", "validated"),
        ("completed_runs", len(per_run)),
        ("expected_read_pairs", totals["expected_read_pairs"]),
        ("primary_records", totals["primary_records"]),
        ("primary_mapped", totals["primary_mapped"]),
        (
            "primary_mapping_percent",
            percent(totals["primary_mapped"], totals["primary_records"]),
        ),
        ("properly_paired", totals["properly_paired"]),
        (
            "properly_paired_percent",
            percent(totals["properly_paired"], totals["primary_records"]),
        ),
        ("primary_duplicates", totals["primary_duplicates"]),
        (
            "primary_duplicate_percent",
            percent(totals["primary_duplicates"], totals["primary_records"]),
        ),
        ("human_primary_mapped", totals["human_primary_mapped"]),
        ("pig_primary_mapped", totals["pig_primary_mapped"]),
        (
            "human_percent_of_primary_mapped",
            percent(totals["human_primary_mapped"], totals["primary_mapped"]),
        ),
        (
            "pig_percent_of_primary_mapped",
            percent(totals["pig_primary_mapped"], totals["primary_mapped"]),
        ),
        ("total_bam_bytes", totals["bam_bytes"]),
        ("source_fastq_tag_rename", "XE_to_ZE"),
        ("global_duplicate_marking", "disabled_multi_cell_CB_aware_required"),
    ]
    summary_rows = [{"metric": key, "value": value} for key, value in summary_metrics]
    atomic_tsv(args.per_run_out, list(per_run[0]), per_run)
    atomic_tsv(args.summary_out, ["metric", "value"], summary_rows)
    print(
        f"[OK] {args.mode} DNA BWA summary: {len(per_run)} runs, "
        f"{totals['expected_read_pairs']:,} input pairs"
    )
    print(
        f"[OK] primary mapping={percent(totals['primary_mapped'], totals['primary_records'])}% "
        f"properly paired={percent(totals['properly_paired'], totals['primary_records'])}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
