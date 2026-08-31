#!/usr/bin/env python3
"""Step 38: summarize and validate DNA CB-aware markdup Probe outputs."""

from __future__ import annotations

import argparse
import csv
import json
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
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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
    checked = 0
    missing = 0
    wrong_xd = 0
    duplicate_flags = 0
    unique_cb: set[str] = set()
    command = [
        samtools,
        "view",
        "-@",
        "4",
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
    try:
        for raw in proc.stdout:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"malformed SAM record from {bam_path}")
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
        if proc.stdout is not None:
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


def atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        tmp.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=PROJECT_DEFAULT)
    parser.add_argument("--sample-tags", type=int, default=10000)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    args = parser.parse_args()

    project = args.project
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools unavailable: {args.samtools}")
    manifest = args.manifest or (
        project
        / "08_qc_summary/13_postmapping_probe/"
        "dna_cb_markdup_probe_manifest.tsv"
    )
    output_root = args.output_root or (
        project / "07_postmapping/dna_cb_dedup_probe"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/13_postmapping_probe/dna_cb_markdup"
    )
    rows = read_tsv(manifest)
    if len(rows) != 4:
        raise ValueError(f"expected four manifest rows, found {len(rows)}")

    per_sample: list[dict[str, object]] = []
    for row in rows:
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
        ):
            raise ValueError(f"obsolete/invalid dedup completion file: {done}")
        marked_bam = output_root / run_id / f"{run_id}.cb_marked.bam"
        marked_bai = Path(str(marked_bam) + ".bai")
        if done_metrics.get("marked_bam_retained") != "1":
            raise ValueError(f"Probe must retain its marked BAM: {done}")
        for path in (marked_bam, marked_bai):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"missing/empty Probe output: {path}")
        stats = read_markdup_stats(stats_path)
        tags = sample_tags(bam, samtools, args.sample_tags)
        if tags["checked_primary_mapped_records"] == 0:
            raise ValueError(f"no primary mapped records sampled: {bam}")
        if tags["missing_required_tag_records"] != 0:
            raise ValueError(f"required BAM tags missing: {bam}")
        if tags["nonzero_XD_records"] != 0:
            raise ValueError(f"non-strict XD observed: {bam}")
        if tags["duplicate_flag_records"] != 0:
            raise ValueError(f"duplicate flag remains after filtering: {bam}")
        marked_tags = sample_tags(
            marked_bam, samtools, args.sample_tags
        )
        if (
            marked_tags["checked_primary_mapped_records"] == 0
            or marked_tags["missing_required_tag_records"] != 0
            or marked_tags["nonzero_XD_records"] != 0
        ):
            raise ValueError(f"marked BAM validation failed: {marked_bam}")
        examined = int(stats.get("examined", 0))
        duplicate_total = int(stats.get("duplicate_total", 0))
        duplicate_primary = int(
            stats.get("duplicate_primary_total", duplicate_total)
        )
        rate = duplicate_primary / examined if examined else 0.0
        result: dict[str, object] = {
            "probe_index": row["probe_index"],
            "run_id": run_id,
            "sample": row["sample"],
            "group": row["group"],
            "selection_reason": row["selection_reason"],
            "input_bam_bytes": int(row["input_bam_bytes"]),
            "output_bam_bytes": bam.stat().st_size,
            "markdup_examined": examined,
            "markdup_duplicate_total": duplicate_total,
            "markdup_duplicate_primary_total": duplicate_primary,
            "markdup_duplicate_fraction": rate,
            "output_bam": str(bam),
            "marked_bam": str(marked_bam),
            "marked_sample_duplicate_flag_records": marked_tags[
                "duplicate_flag_records"
            ],
            **tags,
            "status": "validated",
        }
        per_sample.append(result)
        atomic_write_json(
            qc_root / "by_run" / f"{run_id}.validation.json", result
        )

    qc_root.mkdir(parents=True, exist_ok=True)
    per_path = qc_root / "dna_cb_markdup_probe_per_sample.tsv"
    tmp_per = per_path.with_name(per_path.name + f".tmp.{os.getpid()}")
    with tmp_per.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(per_sample[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(per_sample)
    os.replace(tmp_per, per_path)

    total_examined = sum(int(row["markdup_examined"]) for row in per_sample)
    total_duplicate = sum(
        int(row["markdup_duplicate_primary_total"]) for row in per_sample
    )
    summary_rows = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("completed_samples", len(per_sample)),
        ("duplicate_scope", "within_CB_and_reference_species"),
        (
            "duplicate_key",
            "CB+reference+unclipped_template_endpoints+orientation",
        ),
        ("markdup_mode", "template"),
        ("markdup_expected_max_read_length", "default_300"),
        ("duplicates_removed", 1),
        ("markdup_examined", total_examined),
        ("markdup_duplicate_primary_total", total_duplicate),
        # Backward-compatible name used by step 42 and the completed v12 Probe.
        ("markdup_duplicate_total", total_duplicate),
        (
            "markdup_duplicate_fraction",
            total_duplicate / total_examined if total_examined else 0.0,
        ),
    ]
    summary_path = qc_root / "dna_cb_markdup_probe_summary.tsv"
    tmp_summary = summary_path.with_name(
        summary_path.name + f".tmp.{os.getpid()}"
    )
    with tmp_summary.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("metric\tvalue\n")
        for key, value in summary_rows:
            handle.write(f"{key}\t{value}\n")
    os.replace(tmp_summary, summary_path)
    print(f"[OK] validated {len(per_sample)} DNA CB-aware dedup samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
