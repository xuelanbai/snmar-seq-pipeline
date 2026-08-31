#!/usr/bin/env python3
"""Step 58: validate and summarize all 96 full DNA TAPS outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


EXPECTED = 96


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
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    manifest = (
        args.manifest
        or project
        / "08_qc_summary/17_downstream_full_inputs/dna_taps_cpg_full_manifest.tsv"
    )
    output_root = (
        args.output_root or project / "07_postmapping/dna_taps_cpg_full"
    )
    qc_root = (
        args.qc_root or project / "08_qc_summary/18_dna_taps_cpg_full"
    )
    manifest_rows = read_tsv(manifest)
    if len(manifest_rows) != EXPECTED:
        raise ValueError(f"expected 96 DNA samples, found {len(manifest_rows)}")

    per_sample = []
    total_cov = 0
    total_mod = 0
    human_cov = 0
    human_mod = 0
    pig_cov = 0
    pig_mod = 0
    total_cells = 0
    total_contigs: int | None = None
    total_reference_bases: int | None = None
    expected_call_model: str | None = None
    expected_overlap_policy: str | None = None
    per_cell_sources: list[tuple[str, Path]] = []
    for row in manifest_rows:
        sample = row["sample"]
        sample_dir = output_root / sample
        summary_path = sample_dir / f"{sample}.taps_cpg_full.summary.tsv"
        summary = read_metrics(summary_path)
        if (
            summary.get("status") != "validated"
            or summary.get("sample") != sample
            or summary.get("input_bam") != row["input_bam"]
        ):
            raise ValueError(f"invalid full DNA summary: {sample}")
        for key, expected in (
            ("completed_contigs", total_contigs),
            ("reference_bases", total_reference_bases),
            ("call_model", expected_call_model),
            ("overlap_policy", expected_overlap_policy),
        ):
            observed: str | int
            if key in {"completed_contigs", "reference_bases"}:
                observed = int(summary[key])
            else:
                observed = summary[key]
            if expected is not None and observed != expected:
                raise ValueError(f"inconsistent {key}: {sample}")
        total_contigs = total_contigs or int(summary["completed_contigs"])
        total_reference_bases = total_reference_bases or int(
            summary["reference_bases"]
        )
        expected_call_model = expected_call_model or summary["call_model"]
        expected_overlap_policy = (
            expected_overlap_policy or summary["overlap_policy"]
        )
        covered = int(summary["covered_template_cpg_calls"])
        modified = int(summary["modified_template_cpg_calls"])
        cells = int(summary["cells_with_cpg_calls"])
        if not (covered > 0 and 0 <= modified <= covered and cells > 0):
            raise ValueError(f"invalid full DNA counters: {sample}")
        sample_human_cov = int(summary["human_covered_calls"])
        sample_human_mod = int(summary["human_modified_calls"])
        sample_pig_cov = int(summary["pig_covered_calls"])
        sample_pig_mod = int(summary["pig_modified_calls"])
        per_sample.append(
            {
                "sample_index": row["sample_index"],
                "sample": sample,
                "group": row["group"],
                "cells_with_cpg_calls": cells,
                "covered_template_cpg_calls": covered,
                "modified_template_cpg_calls": modified,
                "modified_rate": modified / covered,
                "human_covered_calls": sample_human_cov,
                "human_modified_rate": (
                    sample_human_mod / sample_human_cov
                    if sample_human_cov
                    else 0.0
                ),
                "pig_covered_calls": sample_pig_cov,
                "pig_modified_rate": (
                    sample_pig_mod / sample_pig_cov
                    if sample_pig_cov
                    else 0.0
                ),
                "overlap_observations_collapsed": int(
                    summary.get("overlap_observations_collapsed", "0")
                ),
                "status": "validated",
            }
        )
        total_cov += covered
        total_mod += modified
        total_cells += cells
        human_cov += sample_human_cov
        human_mod += sample_human_mod
        pig_cov += sample_pig_cov
        pig_mod += sample_pig_mod
        per_cell_sources.append(
            (
                sample,
                sample_dir / f"{sample}.taps_cpg_full.per_cell.tsv",
            )
        )
    if human_cov == 0 or pig_cov == 0:
        raise ValueError("full DNA calling lacks human or pig CpG calls")

    qc_root.mkdir(parents=True, exist_ok=True)
    per_sample.sort(key=lambda value: int(value["sample_index"]))
    atomic_tsv(
        qc_root / "dna_taps_cpg_full_per_sample.tsv",
        list(per_sample[0]),
        per_sample,
    )

    per_cell_output = qc_root / "dna_taps_cpg_full_per_cell.tsv.gz"
    descriptor, temporary = tempfile.mkstemp(
        prefix=per_cell_output.name + ".",
        suffix=".tmp",
        dir=qc_root,
    )
    os.close(descriptor)
    per_cell_fields: list[str] | None = None
    merged_cells = 0
    merged_cell_cov = 0
    merged_cell_mod = 0
    try:
        with gzip.open(
            temporary, "wt", encoding="utf-8", newline="", compresslevel=1
        ) as out:
            writer: csv.DictWriter | None = None
            for sample, source in per_cell_sources:
                with source.open(
                    "rt", encoding="utf-8", newline=""
                ) as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    if not reader.fieldnames:
                        raise ValueError(f"missing header: {source}")
                    if per_cell_fields is None:
                        per_cell_fields = reader.fieldnames
                        writer = csv.DictWriter(
                            out,
                            fieldnames=per_cell_fields,
                            delimiter="\t",
                            lineterminator="\n",
                        )
                        writer.writeheader()
                    elif reader.fieldnames != per_cell_fields:
                        raise ValueError(f"inconsistent header: {source}")
                    assert writer is not None
                    for cell_row in reader:
                        if cell_row["sample"] != sample:
                            raise ValueError(f"sample mismatch: {source}")
                        merged_cells += 1
                        merged_cell_cov += int(
                            cell_row["covered_template_cpg_calls"]
                        )
                        merged_cell_mod += int(
                            cell_row["modified_template_cpg_calls"]
                        )
                        writer.writerow(cell_row)
        if (
            merged_cells != total_cells
            or merged_cell_cov != total_cov
            or merged_cell_mod != total_mod
        ):
            raise ValueError("full per-cell totals do not reconcile")
        os.replace(temporary, per_cell_output)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "completed_samples", "value": len(per_sample)},
        {"metric": "contigs_per_sample", "value": total_contigs},
        {"metric": "reference_bases_per_sample", "value": total_reference_bases},
        {"metric": "call_model", "value": expected_call_model},
        {"metric": "overlap_policy", "value": expected_overlap_policy},
        {"metric": "cells_with_cpg_calls_sum", "value": total_cells},
        {"metric": "covered_template_cpg_calls", "value": total_cov},
        {"metric": "modified_template_cpg_calls", "value": total_mod},
        {"metric": "modified_rate", "value": total_mod / total_cov},
        {"metric": "human_covered_calls", "value": human_cov},
        {"metric": "human_modified_rate", "value": human_mod / human_cov},
        {"metric": "pig_covered_calls", "value": pig_cov},
        {"metric": "pig_modified_rate", "value": pig_mod / pig_cov},
    ]
    atomic_tsv(
        qc_root / "dna_taps_cpg_full_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] full DNA TAPS: samples={len(per_sample)} "
        f"calls={total_cov} cells={total_cells}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
