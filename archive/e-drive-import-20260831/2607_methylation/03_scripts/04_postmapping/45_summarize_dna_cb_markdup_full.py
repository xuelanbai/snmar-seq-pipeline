#!/usr/bin/env python3
"""Step 45: validate and summarize all 96 DNA CB-aware dedup BAMs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from pathlib import Path


PROJECT_DEFAULT = Path(
    "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
)
INTEGER_RE = re.compile(r"^-?[0-9]+$")
PIPELINE_SIGNATURE = "snmar_dna_cb_markdup_corrected_20260731"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows in TSV: {path}")
    return rows


def read_markdup_stats(path: Path) -> dict[str, int | str]:
    stats: dict[str, int | str] = {}
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or ":" not in line:
                continue
            key, value = line.split(":", 1)
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()
            stats[key] = int(value) if INTEGER_RE.fullmatch(value) else value
    return stats


def sample_tags(
    bam_path: Path, samtools: str, limit: int
) -> dict[str, int]:
    command = [
        samtools,
        "view",
        "-@",
        "2",
        "-F",
        "2308",
        str(bam_path),
    ]
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    assert proc.stdout is not None
    checked = 0
    missing = 0
    wrong_xd = 0
    duplicate_flags = 0
    unique_cb: set[str] = set()
    try:
        for raw in proc.stdout:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"malformed SAM record: {bam_path}")
            checked += 1
            flag = int(fields[1])
            if flag & 0x400:
                duplicate_flags += 1
            tags: dict[str, str] = {}
            for field in fields[11:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            required = ("CR", "CB", "CN", "XD", "ZE")
            if any(tag not in tags for tag in required):
                missing += 1
            else:
                unique_cb.add(tags["CB"].split()[0])
                if int(tags["XD"]) != 0:
                    wrong_xd += 1
            if checked >= limit:
                break
    finally:
        proc.stdout.close()
        if proc.poll() is None:
            proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
    return {
        "checked_primary_mapped_records": checked,
        "missing_required_tag_records": missing,
        "nonzero_XD_records": wrong_xd,
        "duplicate_flag_records": duplicate_flags,
        "sampled_distinct_CB": len(unique_cb),
    }


def atomic_table(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with tmp.open("wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=PROJECT_DEFAULT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--sample-tags", type=int, default=1000)
    args = parser.parse_args()

    project = args.project
    manifest = args.manifest or (
        project
        / "08_qc_summary/14_postmapping_full_inputs/"
        "dna_cb_markdup_full_manifest.tsv"
    )
    output_root = args.output_root or (
        project / "07_postmapping/dna_cb_dedup_full"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/15_dna_cb_markdup_full"
    )
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools unavailable: {args.samtools}")
    if args.sample_tags < 1:
        parser.error("--sample-tags must be positive")

    manifest_rows = read_tsv(manifest)
    if len(manifest_rows) != 96:
        raise ValueError(
            f"expected 96 manifest rows, found {len(manifest_rows)}"
        )
    per_sample: list[dict[str, object]] = []
    for row in manifest_rows:
        run_id = row["run_id"]
        bam = output_root / run_id / f"{run_id}.cb_dedup.bam"
        bai = Path(str(bam) + ".bai")
        stats_path = qc_root / "by_run" / f"{run_id}.markdup.stats.txt"
        done = qc_root / "by_run" / f"{run_id}.done.tsv"
        for path in (bam, bai, stats_path, done):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"missing/empty output: {path}")
        done_metrics = {
            item["metric"]: item["value"] for item in read_tsv(done)
        }
        if (
            done_metrics.get("status") != "completed"
            or done_metrics.get("pipeline_signature")
            != PIPELINE_SIGNATURE
            or done_metrics.get("duplicate_key")
            != "CB+reference+unclipped_template_endpoints+orientation"
            or done_metrics.get("marked_bam_retained") != "0"
        ):
            raise ValueError(f"obsolete/invalid dedup completion file: {done}")
        stats = read_markdup_stats(stats_path)
        examined = int(stats.get("examined", -1))
        duplicate_total = int(stats.get("duplicate_total", -1))
        duplicate_primary = int(
            stats.get("duplicate_primary_total", duplicate_total)
        )
        expected_examined = int(row["input_primary_mapped"])
        if examined != expected_examined:
            raise ValueError(
                f"{run_id}: markdup examined={examined}, "
                f"expected primary mapped={expected_examined}"
            )
        if duplicate_primary < 0 or duplicate_primary > examined:
            raise ValueError(f"{run_id}: invalid primary duplicate count")
        tags = sample_tags(bam, samtools, args.sample_tags)
        if tags["checked_primary_mapped_records"] == 0:
            raise ValueError(f"{run_id}: no primary mapped output records")
        if (
            tags["missing_required_tag_records"]
            or tags["nonzero_XD_records"]
            or tags["duplicate_flag_records"]
        ):
            raise ValueError(f"{run_id}: BAM tag/flag validation failed")
        retained = examined - duplicate_primary
        per_sample.append(
            {
                "sample_index": row["sample_index"],
                "run_id": run_id,
                "sample": row["sample"],
                "group": row["group"],
                "input_primary_mapped": expected_examined,
                "duplicate_primary_records_removed": duplicate_primary,
                "duplicate_total_records_flagged": duplicate_total,
                "retained_primary_mapped": retained,
                "duplicate_fraction": duplicate_primary / examined,
                "retained_fraction": retained / examined,
                "input_bam_bytes": int(row["input_bam_bytes"]),
                "output_bam_bytes": bam.stat().st_size,
                "output_bam": str(bam),
                **tags,
                "status": "validated",
            }
        )

    per_path = qc_root / "dna_cb_markdup_full_per_sample.tsv"
    atomic_table(per_path, list(per_sample[0]), per_sample)

    total_input = sum(
        int(row["input_primary_mapped"]) for row in per_sample
    )
    total_duplicate = sum(
        int(row["duplicate_primary_records_removed"]) for row in per_sample
    )
    total_retained = sum(
        int(row["retained_primary_mapped"]) for row in per_sample
    )
    total_input_bytes = sum(
        int(row["input_bam_bytes"]) for row in per_sample
    )
    total_output_bytes = sum(
        int(row["output_bam_bytes"]) for row in per_sample
    )
    summary_rows = [
        {"metric": "status", "value": "validated"},
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE},
        {"metric": "completed_samples", "value": 96},
        {
            "metric": "duplicate_scope",
            "value": "within_CB_and_reference_species",
        },
        {
            "metric": "duplicate_key",
            "value": (
                "CB+reference+unclipped_template_endpoints+orientation"
            ),
        },
        {"metric": "markdup_mode", "value": "template"},
        {
            "metric": "markdup_expected_max_read_length",
            "value": "default_300",
        },
        {
            "metric": "reported_unit",
            "value": "primary_alignment_records",
        },
        {"metric": "duplicates_removed", "value": 1},
        {"metric": "input_primary_mapped", "value": total_input},
        {
            "metric": "duplicate_primary_records_removed",
            "value": total_duplicate,
        },
        {"metric": "retained_primary_mapped", "value": total_retained},
        {
            "metric": "duplicate_fraction",
            "value": total_duplicate / total_input,
        },
        {
            "metric": "retained_fraction",
            "value": total_retained / total_input,
        },
        {"metric": "input_bam_bytes", "value": total_input_bytes},
        {"metric": "output_bam_bytes", "value": total_output_bytes},
        {
            "metric": "output_to_input_byte_fraction",
            "value": total_output_bytes / total_input_bytes,
        },
    ]
    atomic_table(
        qc_root / "dna_cb_markdup_full_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] DNA full CB-aware dedup validated: "
        f"retained={total_retained}/{total_input}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
