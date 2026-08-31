#!/usr/bin/env python3
"""Step 55: merge validated contig shards for one full DNA sample."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import tempfile
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def read_metrics(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    metrics = {row["metric"]: row["value"] for row in rows}
    if len(metrics) != len(rows):
        raise ValueError(f"duplicate metric: {path}")
    return metrics


def atomic_tsv(
    path: Path, fields: list[str], rows: list[dict[str, object]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(
            descriptor, "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--shard-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--expected-shards", type=int, default=8)
    parser.add_argument("--cleanup-shards", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.expected_shards < 1:
        raise ValueError("expected-shards must be positive")
    if not args.bam.is_file() or args.bam.stat().st_size == 0:
        raise FileNotFoundError(args.bam)
    shard_summaries: list[tuple[int, dict[str, str], Path]] = []
    aggregate: Counter[str] = Counter()
    total_cov = 0
    total_mod = 0
    human_cov = 0
    human_mod = 0
    pig_cov = 0
    pig_mod = 0
    expected_bin_size: str | None = None
    expected_call_model: str | None = None
    expected_overlap_policy: str | None = None

    for shard in range(1, args.expected_shards + 1):
        prefix = f"{args.sample}.shard_{shard:02d}"
        shard_dir = args.shard_root / f"shard_{shard:02d}"
        summary_path = shard_dir / f"{prefix}.summary.tsv"
        summary = read_metrics(summary_path)
        if (
            summary.get("status") != "validated"
            or summary.get("sample") != args.sample
            or int(summary.get("shard", "0")) != shard
            or summary.get("input_bam") != str(args.bam)
        ):
            raise ValueError(f"invalid shard summary: {summary_path}")
        if int(summary.get("expected_contigs", "0")) != int(
            summary.get("completed_contigs", "-1")
        ):
            raise ValueError(f"incomplete contigs: {summary_path}")
        if int(summary.get("missing_required_CB_or_XD", "0")) != 0:
            raise ValueError(f"missing tags in shard {shard}")
        if int(summary.get("nonzero_XD", "0")) != 0:
            raise ValueError(f"non-strict XD in shard {shard}")
        for key, expected in (
            ("bin_size", expected_bin_size),
            ("call_model", expected_call_model),
            ("overlap_policy", expected_overlap_policy),
        ):
            if expected is not None and summary[key] != expected:
                raise ValueError(f"inconsistent {key} across shards")
        expected_bin_size = expected_bin_size or summary["bin_size"]
        expected_call_model = expected_call_model or summary["call_model"]
        expected_overlap_policy = (
            expected_overlap_policy or summary["overlap_policy"]
        )
        covered = int(summary["covered_template_cpg_calls"])
        modified = int(summary["modified_template_cpg_calls"])
        if not (covered > 0 and 0 <= modified <= covered):
            raise ValueError(f"invalid CpG totals in shard {shard}")
        total_cov += covered
        total_mod += modified
        human_cov += int(summary["human_covered_calls"])
        human_mod += int(summary["human_modified_calls"])
        pig_cov += int(summary["pig_covered_calls"])
        pig_mod += int(summary["pig_modified_calls"])
        for key in (
            "alignment_records_after_flag_mapq_filter",
            "tag_eligible_alignment_records",
            "raw_directional_cpg_observations",
            "overlap_observations_collapsed",
            "equal_quality_conflicts",
            "conflicting_template_cpg_calls_discarded",
            "completed_contigs",
            "reference_bases",
        ):
            aggregate[key] += int(summary.get(key, "0"))
        shard_summaries.append((shard, summary, shard_dir))

    cell_reads: Counter[str] = Counter()
    cell_cov: Counter[str] = Counter()
    cell_mod: Counter[str] = Counter()
    cell_bins: Counter[str] = Counter()
    for shard, _summary, shard_dir in shard_summaries:
        path = (
            shard_dir
            / f"{args.sample}.shard_{shard:02d}.per_cell.tsv"
        )
        for row in read_tsv(path):
            if row["sample"] != args.sample:
                raise ValueError(f"sample mismatch in: {path}")
            cell_id = row["cell_id"]
            if not cell_id.startswith(args.sample + "|"):
                raise ValueError(f"cell namespace mismatch: {cell_id}")
            cell_reads[cell_id] += int(row["eligible_alignment_records"])
            cell_cov[cell_id] += int(row["covered_template_cpg_calls"])
            cell_mod[cell_id] += int(row["modified_template_cpg_calls"])
            cell_bins[cell_id] += int(row["bins_with_calls"])
    if sum(cell_cov.values()) != total_cov or sum(cell_mod.values()) != total_mod:
        raise ValueError("per-cell totals do not equal shard summaries")

    sample_dir = args.output_root / args.sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    per_cell_rows = []
    for cell_id in sorted(cell_cov, key=lambda value: (-cell_cov[value], value)):
        covered = cell_cov[cell_id]
        modified = cell_mod[cell_id]
        per_cell_rows.append(
            {
                "sample": args.sample,
                "cell_id": cell_id,
                "CB": cell_id.split("|", 1)[1],
                "eligible_alignment_records": cell_reads[cell_id],
                "covered_template_cpg_calls": covered,
                "modified_template_cpg_calls": modified,
                "modified_rate": modified / covered,
                "bins_with_calls": cell_bins[cell_id],
            }
        )
    atomic_tsv(
        sample_dir / f"{args.sample}.taps_cpg_full.per_cell.tsv",
        list(per_cell_rows[0]),
        per_cell_rows,
    )

    bin_fields = [
        "sample",
        "cell_id",
        "CB",
        "species",
        "chrom",
        "start",
        "end",
        "modified_template_cpg_calls",
        "covered_template_cpg_calls",
        "modified_rate",
    ]
    bin_output = sample_dir / f"{args.sample}.taps_cpg_full.bin_counts.tsv.gz"
    descriptor, temporary = tempfile.mkstemp(
        prefix=bin_output.name + ".", suffix=".tmp", dir=sample_dir
    )
    os.close(descriptor)
    merged_bin_cov = 0
    merged_bin_mod = 0
    try:
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as out:
            writer = csv.DictWriter(
                out,
                fieldnames=bin_fields,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for shard, _summary, shard_dir in shard_summaries:
                path = (
                    shard_dir
                    / f"{args.sample}.shard_{shard:02d}.bin_counts.tsv.gz"
                )
                with gzip.open(
                    path, "rt", encoding="utf-8", newline=""
                ) as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    if reader.fieldnames != bin_fields:
                        raise ValueError(f"unexpected bin header: {path}")
                    for row in reader:
                        if row["sample"] != args.sample:
                            raise ValueError(f"sample mismatch in: {path}")
                        covered = int(row["covered_template_cpg_calls"])
                        modified = int(row["modified_template_cpg_calls"])
                        if not (covered > 0 and 0 <= modified <= covered):
                            raise ValueError(f"invalid bin counts in: {path}")
                        merged_bin_cov += covered
                        merged_bin_mod += modified
                        writer.writerow(row)
        if merged_bin_cov != total_cov or merged_bin_mod != total_mod:
            raise ValueError("merged bin totals do not equal shard summaries")
        os.replace(temporary, bin_output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "sample", "value": args.sample},
        {"metric": "input_bam", "value": args.bam},
        {"metric": "input_bam_bytes", "value": args.bam.stat().st_size},
        {"metric": "input_bam_mtime", "value": int(args.bam.stat().st_mtime)},
        {"metric": "completed_shards", "value": len(shard_summaries)},
        {"metric": "completed_contigs", "value": aggregate["completed_contigs"]},
        {"metric": "reference_bases", "value": aggregate["reference_bases"]},
        {"metric": "bin_size", "value": expected_bin_size},
        {"metric": "call_model", "value": expected_call_model},
        {"metric": "overlap_policy", "value": expected_overlap_policy},
        {"metric": "cells_with_cpg_calls", "value": len(cell_cov)},
        {"metric": "covered_template_cpg_calls", "value": total_cov},
        {"metric": "modified_template_cpg_calls", "value": total_mod},
        {"metric": "modified_rate", "value": total_mod / total_cov},
        {"metric": "human_covered_calls", "value": human_cov},
        {"metric": "human_modified_calls", "value": human_mod},
        {
            "metric": "human_modified_rate",
            "value": human_mod / human_cov if human_cov else 0.0,
        },
        {"metric": "pig_covered_calls", "value": pig_cov},
        {"metric": "pig_modified_calls", "value": pig_mod},
        {
            "metric": "pig_modified_rate",
            "value": pig_mod / pig_cov if pig_cov else 0.0,
        },
    ]
    summary_rows.extend(
        {"metric": key, "value": aggregate[key]}
        for key in (
            "alignment_records_after_flag_mapq_filter",
            "tag_eligible_alignment_records",
            "raw_directional_cpg_observations",
            "overlap_observations_collapsed",
            "equal_quality_conflicts",
            "conflicting_template_cpg_calls_discarded",
        )
    )
    summary_path = sample_dir / f"{args.sample}.taps_cpg_full.summary.tsv"
    atomic_tsv(summary_path, ["metric", "value"], summary_rows)

    if args.cleanup_shards:
        resolved = args.shard_root.resolve()
        if (
            not resolved.is_dir()
            or resolved.name != args.sample
            or resolved == resolved.parent
        ):
            raise ValueError(f"unsafe shard cleanup target: {resolved}")
        shutil.rmtree(resolved)
    print(
        f"[OK] {args.sample}: full CpG calls={total_cov}, "
        f"cells={len(cell_cov)}, shards={len(shard_summaries)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
