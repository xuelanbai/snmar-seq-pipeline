#!/usr/bin/env python3
"""Step 40: validate and summarize RNA GeneFull/UMI Probe results."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path


PROJECT_DEFAULT = Path(
    "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
)
MODES = ("forward", "reverse", "unstranded")


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=PROJECT_DEFAULT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--result-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    args = parser.parse_args()

    project = args.project
    manifest = args.manifest or (
        project
        / "08_qc_summary/13_postmapping_probe/"
        "rna_genefull_probe_manifest.tsv"
    )
    result_root = args.result_root or (
        project / "07_postmapping/rna_genefull_probe"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/13_postmapping_probe/rna_genefull"
    )
    rows = read_manifest(manifest)
    if len(rows) != 8:
        raise ValueError(f"expected eight RNA Probe rows, found {len(rows)}")

    per_file: list[dict[str, object]] = []
    aggregate = {
        mode: {
            "assigned_reads": 0,
            "exact_umis": 0,
            "eligible": 0,
        }
        for mode in MODES
    }
    for row in rows:
        pair_id = row["pair_id"]
        path = result_root / pair_id / f"{pair_id}.summary.json"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty RNA Probe summary: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") != "completed":
            raise ValueError(f"incomplete RNA Probe result: {path}")
        counters = payload["counters"]
        eligible = int(counters.get("tag_eligible_unique_primary", 0))
        if eligible <= 0:
            raise ValueError(f"no tag-eligible unique reads: {pair_id}")
        if int(counters.get("missing_required_tags", 0)) != 0:
            raise ValueError(f"missing required tags: {pair_id}")
        if int(counters.get("nonzero_XD", 0)) != 0:
            raise ValueError(f"non-strict XD observed: {pair_id}")
        for mode in MODES:
            metrics = payload["modes"][mode]
            assigned = int(metrics["assigned_reads"])
            umis = int(metrics["exact_umis"])
            aggregate[mode]["assigned_reads"] += assigned
            aggregate[mode]["exact_umis"] += umis
            aggregate[mode]["eligible"] += eligible
            per_file.append(
                {
                    "probe_index": row["probe_index"],
                    "pair_id": pair_id,
                    "sample": row["sample"],
                    "selection_reason": row["selection_reason"],
                    "mode": mode,
                    "tag_eligible_unique_primary": eligible,
                    **metrics,
                    "status": "validated",
                }
            )

    qc_root.mkdir(parents=True, exist_ok=True)
    per_path = qc_root / "rna_genefull_probe_per_file_mode.tsv"
    tmp_per = per_path.with_name(per_path.name + f".tmp.{os.getpid()}")
    with tmp_per.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(per_file[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(per_file)
    os.replace(tmp_per, per_path)

    forward = aggregate["forward"]["assigned_reads"]
    reverse = aggregate["reverse"]["assigned_reads"]
    if reverse == 0:
        ratio: float | str = "inf" if forward else "undefined"
    else:
        ratio = forward / reverse
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("completed_files", 8),
        ("counting_feature", "GeneFull_gene_body"),
        ("umi_deduplication", "exact_within_sample_CB_gene"),
        ("strand_decision", "pending_probe_review"),
        ("forward_to_reverse_assigned_read_ratio", ratio),
    ]
    for mode in MODES:
        values = aggregate[mode]
        summary_rows.extend(
            (
                (f"{mode}_tag_eligible_unique_primary", values["eligible"]),
                (f"{mode}_assigned_reads", values["assigned_reads"]),
                (
                    f"{mode}_assigned_fraction",
                    values["assigned_reads"] / values["eligible"],
                ),
                (f"{mode}_exact_umis", values["exact_umis"]),
            )
        )
    summary_path = qc_root / "rna_genefull_probe_summary.tsv"
    tmp_summary = summary_path.with_name(
        summary_path.name + f".tmp.{os.getpid()}"
    )
    with tmp_summary.open("wt", encoding="utf-8", newline="") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in summary_rows:
            handle.write(f"{metric}\t{value}\n")
    os.replace(tmp_summary, summary_path)
    print("[OK] validated and summarized eight RNA GeneFull Probe files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
