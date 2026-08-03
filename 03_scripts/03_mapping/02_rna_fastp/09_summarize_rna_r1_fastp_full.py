#!/usr/bin/env python3
"""Step 09: verify and summarize full RNA R1-only fastp outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


PIPELINE_SIGNATURE = "rna_r1_se_fastp_v9"
REQUIRED_TAGS = (
    "CR:Z:",
    "CB:Z:",
    "CN:Z:",
    "UR:Z:",
    "UB:Z:",
    "XD:i:0",
    "XE:Z:",
    "XL:i:",
    "XI:Z:full10",
)


def load_expected_assigned(path: Path) -> dict[str, int]:
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


def read_record(handle, path: Path, record_no: int):
    raw_header = handle.readline()
    if raw_header == "":
        return None
    header = raw_header.rstrip("\r\n")
    if not header.startswith("@"):
        raise ValueError(f"{path}: invalid header near record {record_no}")
    sequence = handle.readline().rstrip("\r\n")
    plus = handle.readline().rstrip("\r\n")
    quality = handle.readline().rstrip("\r\n")
    if not sequence or not plus or not quality:
        raise ValueError(f"{path}: truncated FASTQ near record {record_no}")
    if not plus.startswith("+"):
        raise ValueError(f"{path}: invalid plus line near record {record_no}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"{path}: sequence/quality length mismatch near record {record_no}"
        )
    for tag in REQUIRED_TAGS:
        if tag not in header:
            raise ValueError(f"{path}: header missing {tag}: {header[:240]}")
    return header, sequence, plus, quality


def verify_clean_r1(path: Path, limit: int | None) -> int:
    checked = 0
    with gzip.open(path, "rt") as handle:
        while limit is None or checked < limit:
            record = read_record(handle, path, checked + 1)
            if record is None:
                break
            checked += 1
    return checked


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(4 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def load_done_manifest(
    path: Path,
    pair_id: str,
    sample: str,
    raw_r1: Path,
    clean_r1: Path,
    json_path: Path,
    html_path: Path,
) -> dict[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty done manifest: {path}")

    result: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"{path}:{line_no}: expected key<TAB>value")
            key, value = fields
            if key in result:
                raise ValueError(f"{path}:{line_no}: duplicate key {key}")
            result[key] = value

    expected_text = {
        "status": "complete",
        "pipeline_signature": PIPELINE_SIGNATURE,
        "pair_id": pair_id,
        "sample": sample,
        "raw_R1_path": str(raw_r1),
        "raw_R1_bytes": str(raw_r1.stat().st_size),
        "raw_R1_mtime": str(int(raw_r1.stat().st_mtime)),
    }
    for key, expected in expected_text.items():
        if result.get(key) != expected:
            raise ValueError(
                f"{path}: {key}={result.get(key)!r}, expected {expected!r}"
            )

    size_paths = {
        "clean_R1_bytes": clean_r1,
        "json_bytes": json_path,
        "html_bytes": html_path,
    }
    for key, file_path in size_paths.items():
        try:
            recorded = int(result[key])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"{path}: invalid/missing {key}") from exc
        if recorded != file_path.stat().st_size:
            raise ValueError(
                f"{path}: {key}={recorded}, "
                f"actual={file_path.stat().st_size}"
            )

    checksum = sha256_file(clean_r1)
    if result.get("clean_R1_sha256", "").lower() != checksum:
        raise ValueError(f"{path}: clean R1 SHA256 mismatch")

    for key in (
        "expected_input_R1_reads",
        "fastp_before_reads",
        "fastp_after_reads",
        "verified_output_R1_reads",
    ):
        try:
            int(result[key])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"{path}: invalid/missing {key}") from exc
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1-table", type=Path)
    parser.add_argument("--demux-per-pair", type=Path)
    parser.add_argument("--clean-root", type=Path)
    parser.add_argument("--json-root", type=Path)
    parser.add_argument("--done-root", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--per-file", type=Path)
    parser.add_argument("--check-headers", type=int, default=1000)
    parser.add_argument("--verify-clean-r1", type=Path)
    parser.add_argument("--expected-output-reads", type=int)
    args = parser.parse_args()

    if args.check_headers < 0:
        parser.error("--check-headers must be non-negative")

    verify_requested = (
        args.verify_clean_r1 is not None
        or args.expected_output_reads is not None
    )
    if verify_requested:
        if (
            args.verify_clean_r1 is None
            or args.expected_output_reads is None
        ):
            parser.error(
                "--verify-clean-r1 and --expected-output-reads "
                "must be used together"
            )
        if args.expected_output_reads < 1:
            parser.error("--expected-output-reads must be positive")
        checked = verify_clean_r1(args.verify_clean_r1, None)
        if checked != args.expected_output_reads:
            raise ValueError(
                f"clean R1 reads={checked}, "
                f"expected={args.expected_output_reads}"
            )
        print(checked)
        return 0

    required = (
        "r1_table",
        "demux_per_pair",
        "clean_root",
        "json_root",
        "done_root",
        "summary",
        "per_file",
    )
    missing = [
        f"--{name.replace('_', '-')}"
        for name in required
        if getattr(args, name) is None
    ]
    if missing:
        parser.error("missing required argument(s): " + ", ".join(missing))

    assert args.r1_table is not None
    assert args.demux_per_pair is not None
    assert args.clean_root is not None
    assert args.json_root is not None
    assert args.done_root is not None
    assert args.summary is not None
    assert args.per_file is not None

    with args.r1_table.open(encoding="utf-8", newline="") as handle:
        input_rows = list(csv.DictReader(handle, delimiter="\t"))
    if not input_rows:
        raise ValueError(f"empty R1 table: {args.r1_table}")

    pair_ids = [row["pair_id"] for row in input_rows]
    if len(pair_ids) != len(set(pair_ids)):
        raise ValueError(f"{args.r1_table}: duplicate pair_id")

    expected = load_expected_assigned(args.demux_per_pair)
    missing_expected = sorted(set(pair_ids) - set(expected))
    extra_expected = sorted(set(expected) - set(pair_ids))
    if missing_expected or extra_expected:
        raise ValueError(
            "pair IDs differ between discovery and demultiplex tables: "
            f"missing_expected={missing_expected[:10]} "
            f"extra_expected={extra_expected[:10]}"
        )

    rows: list[dict[str, object]] = []
    totals: Counter[str] = Counter()

    for input_row in input_rows:
        pair_id = input_row["pair_id"]
        sample = input_row["sample"]
        raw_r1 = Path(input_row["R1"])
        clean_r1 = (
            args.clean_root / sample / f"{pair_id}.clean_R1.fq.gz"
        )
        json_path = args.json_root / sample / f"{pair_id}.fastp.json"
        html_path = args.json_root / sample / f"{pair_id}.fastp.html"
        done_path = args.done_root / f"{pair_id}.done"
        for file_path in (raw_r1, clean_r1, json_path, html_path):
            if not file_path.is_file() or file_path.stat().st_size == 0:
                raise FileNotFoundError(f"missing/empty RNA R1 file: {file_path}")

        done = load_done_manifest(
            done_path,
            pair_id,
            sample,
            raw_r1,
            clean_r1,
            json_path,
            html_path,
        )
        expected_reads = expected[pair_id]

        payload = json.loads(json_path.read_text(encoding="utf-8"))
        before = payload["summary"]["before_filtering"]
        after = payload["summary"]["after_filtering"]
        filtering = payload.get("filtering_result", {})
        adapter = payload.get("adapter_cutting", {})
        before_reads = int(before["total_reads"])
        after_reads = int(after["total_reads"])

        if before_reads != expected_reads:
            raise ValueError(
                f"{pair_id}: fastp input R1 reads={before_reads}, "
                f"demultiplex assigned_XD0={expected_reads}"
            )
        if after_reads <= 0 or after_reads > before_reads:
            raise ValueError(
                f"{pair_id}: invalid fastp output reads={after_reads}"
            )
        if int(done["expected_input_R1_reads"]) != expected_reads:
            raise ValueError(f"{pair_id}: done/demultiplex input mismatch")
        if int(done["fastp_before_reads"]) != before_reads:
            raise ValueError(f"{pair_id}: done/JSON before mismatch")
        if int(done["fastp_after_reads"]) != after_reads:
            raise ValueError(f"{pair_id}: done/JSON after mismatch")
        if int(done["verified_output_R1_reads"]) != after_reads:
            raise ValueError(f"{pair_id}: verified/JSON output mismatch")

        checked = verify_clean_r1(clean_r1, args.check_headers)
        after_bases = int(after["total_bases"])
        row: dict[str, object] = {
            "pair_id": pair_id,
            "sample": sample,
            "demultiplex_assigned_XD0_R1_reads": expected_reads,
            "fastp_input_R1_reads": before_reads,
            "fastp_output_R1_reads": after_reads,
            "fastp_R1_retention_percent": (
                after_reads / max(before_reads, 1) * 100
            ),
            "fastp_after_q20_rate": float(after["q20_rate"]),
            "fastp_after_q30_rate": float(after["q30_rate"]),
            "fastp_after_R1_mean_length": float(
                after.get("read1_mean_length", 0.0)
            ),
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
            "checked_strict_R1_headers": checked,
        }
        rows.append(row)

        totals["expected_reads"] += expected_reads
        totals["before_reads"] += before_reads
        totals["after_reads"] += after_reads
        totals["after_bases"] += after_bases
        totals["after_q20_bases"] += int(after["q20_bases"])
        totals["after_q30_bases"] += int(after["q30_bases"])
        totals["low_quality_reads"] += int(row["low_quality_reads"])
        totals["too_many_N_reads"] += int(row["too_many_N_reads"])
        totals["too_short_reads"] += int(row["too_short_reads"])
        totals["adapter_trimmed_reads"] += int(
            row["adapter_trimmed_reads"]
        )

    if totals["expected_reads"] != totals["before_reads"]:
        raise ValueError("global demultiplex/fastp input R1 count mismatch")

    args.per_file.parent.mkdir(parents=True, exist_ok=True)
    with args.per_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    before_reads = max(totals["before_reads"], 1)
    after_reads = max(totals["after_reads"], 1)
    after_bases = max(totals["after_bases"], 1)
    summary_rows: list[tuple[str, object, str]] = [
        ("completed_R1_FASTQ_files", len(rows), "R1 FASTQ files"),
        (
            "demultiplex_assigned_XD0_R1_reads",
            totals["expected_reads"],
            "R1 reads",
        ),
        ("fastp_input_R1_reads", totals["before_reads"], "R1 reads"),
        ("fastp_output_R1_reads", totals["after_reads"], "R1 reads"),
        (
            "fastp_R1_retention_percent",
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
            totals["after_bases"] / after_reads,
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
    ]

    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with args.summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value", "unit"))
        for metric, value, unit in summary_rows:
            if isinstance(value, float):
                value = f"{value:.6f}"
            writer.writerow((metric, value, unit))

    print(
        f"[OK] RNA R1-only fastp files={len(rows)}; "
        f"input={totals['before_reads']}; "
        f"output={totals['after_reads']}; "
        f"retention={totals['after_reads'] / before_reads * 100:.3f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
