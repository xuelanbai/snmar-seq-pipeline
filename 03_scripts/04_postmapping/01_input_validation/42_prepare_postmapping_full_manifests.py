#!/usr/bin/env python3
"""Step 42: validate both Probes and prepare full DNA/RNA manifests."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path


PROJECT_DEFAULT = Path(
    "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
)
TRAILING_NUMBER_RE = re.compile(r"([0-9]+)$")


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows in TSV: {path}")
    return rows


def read_metrics(path: Path) -> dict[str, str]:
    rows = read_table(path)
    if set(rows[0]) != {"metric", "value"}:
        raise ValueError(f"expected metric/value TSV: {path}")
    metrics = {row["metric"]: row["value"] for row in rows}
    if len(metrics) != len(rows):
        raise ValueError(f"duplicate metric in {path}")
    return metrics


def trailing_number(value: str) -> int:
    match = TRAILING_NUMBER_RE.search(value)
    if not match:
        raise ValueError(f"identifier lacks trailing number: {value}")
    return int(match.group(1))


def atomic_tsv(
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
    parser.add_argument("--dna-full-table", type=Path)
    parser.add_argument("--rna-full-table", type=Path)
    parser.add_argument("--dna-probe-summary", type=Path)
    parser.add_argument("--rna-probe-summary", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--rna-shards", type=int, default=4)
    args = parser.parse_args()

    project = args.project
    dna_full_table = args.dna_full_table or (
        project
        / "08_qc_summary/12_dna_bwa_full/dna_bwa_full_per_sample.tsv"
    )
    rna_full_table = args.rna_full_table or (
        project
        / "08_qc_summary/09_rna_star_full/rna_star_full_per_file.tsv"
    )
    dna_probe_summary = args.dna_probe_summary or (
        project
        / "08_qc_summary/13_postmapping_probe/dna_cb_markdup/"
        "dna_cb_markdup_probe_summary.tsv"
    )
    rna_probe_summary = args.rna_probe_summary or (
        project
        / "08_qc_summary/13_postmapping_probe/rna_genefull/"
        "rna_genefull_probe_summary.tsv"
    )
    out_dir = args.out_dir or (
        project / "08_qc_summary/14_postmapping_full_inputs"
    )
    if args.rna_shards < 1 or args.rna_shards > 16:
        parser.error("--rna-shards must be between 1 and 16")

    dna_probe = read_metrics(dna_probe_summary)
    if (
        dna_probe.get("status") != "validated"
        or int(dna_probe.get("completed_samples", "0")) != 4
        or dna_probe.get("pipeline_signature")
        != "snmar_dna_cb_markdup_corrected_20260731"
        or dna_probe.get("duplicate_scope")
        != "within_CB_and_reference_species"
        or dna_probe.get("duplicate_key")
        != "CB+reference+unclipped_template_endpoints+orientation"
        or int(dna_probe.get("duplicates_removed", "0")) != 1
    ):
        raise ValueError("DNA CB-aware markdup Probe has not passed")
    dna_duplicate_fraction = float(
        dna_probe["markdup_duplicate_fraction"]
    )
    if not 0.0 < dna_duplicate_fraction < 1.0:
        raise ValueError("invalid DNA Probe duplicate fraction")

    rna_probe = read_metrics(rna_probe_summary)
    if (
        rna_probe.get("status") != "validated"
        or int(rna_probe.get("completed_files", "0")) != 8
        or rna_probe.get("counting_feature") != "GeneFull_gene_body"
        or rna_probe.get("umi_deduplication")
        != "exact_within_sample_CB_gene"
    ):
        raise ValueError("RNA GeneFull/UMI Probe has not passed")
    forward_reads = int(rna_probe["forward_assigned_reads"])
    reverse_reads = int(rna_probe["reverse_assigned_reads"])
    forward_umis = int(rna_probe["forward_exact_umis"])
    reverse_umis = int(rna_probe["reverse_exact_umis"])
    if reverse_reads <= forward_reads or reverse_umis <= forward_umis:
        raise ValueError(
            "RNA Probe does not support the required reverse-strand decision"
        )

    dna_rows = read_table(dna_full_table)
    rna_rows = read_table(rna_full_table)
    if len(dna_rows) != 96 or len(rna_rows) != 96:
        raise ValueError(
            f"expected 96 DNA and 96 RNA rows; "
            f"found {len(dna_rows)} and {len(rna_rows)}"
        )
    if any(row.get("status") != "validated" for row in dna_rows):
        raise ValueError("all DNA full mapping rows must be validated")
    if len({row["run_id"] for row in dna_rows}) != 96:
        raise ValueError("DNA run_id values are not unique")
    if len({row["pair_id"] for row in rna_rows}) != 96:
        raise ValueError("RNA pair_id values are not unique")

    dna_sorted = sorted(
        dna_rows, key=lambda row: trailing_number(row["sample"])
    )
    dna_manifest: list[dict[str, object]] = []
    for index, row in enumerate(dna_sorted, 1):
        if trailing_number(row["sample"]) != index:
            raise ValueError(
                f"DNA sample numbering is not complete at index {index}"
            )
        dna_manifest.append(
            {
                "sample_index": index,
                "run_id": row["run_id"],
                "sample": row["sample"],
                "group": row["group"],
                "input_primary_mapped": row["primary_mapped"],
                "input_bam_bytes": row["bam_bytes"],
                "input_bam": row["bam"],
            }
        )

    # Greedy size balancing keeps the four RNA jobs close in expected work.
    shard_totals = [0 for _ in range(args.rna_shards)]
    shard_members: list[list[dict[str, str]]] = [
        [] for _ in range(args.rna_shards)
    ]
    for row in sorted(
        rna_rows, key=lambda item: int(item["star_input_reads"]), reverse=True
    ):
        shard_index = min(
            range(args.rna_shards), key=lambda index: shard_totals[index]
        )
        shard_members[shard_index].append(row)
        shard_totals[shard_index] += int(row["star_input_reads"])

    rna_manifest: list[dict[str, object]] = []
    global_index = 0
    for shard_index, members in enumerate(shard_members, 1):
        members.sort(key=lambda row: trailing_number(row["sample"]))
        for shard_order, row in enumerate(members, 1):
            global_index += 1
            rna_manifest.append(
                {
                    "manifest_index": global_index,
                    "shard": shard_index,
                    "shard_order": shard_order,
                    "pair_id": row["pair_id"],
                    "sample": row["sample"],
                    "star_input_reads": row["star_input_reads"],
                    "primary_mapped_reads": row["primary_mapped_reads"],
                    "primary_mapping_percent": row[
                        "primary_mapping_percent"
                    ],
                    "input_bam": row["bam_path"],
                }
            )
    if sorted(len(members) for members in shard_members) != [
        24,
        24,
        24,
        24,
    ]:
        raise ValueError(
            "four RNA shards must contain exactly 24 files each"
        )

    atomic_tsv(
        out_dir / "dna_cb_markdup_full_manifest.tsv",
        list(dna_manifest[0]),
        dna_manifest,
    )
    atomic_tsv(
        out_dir / "rna_reverse_genefull_full_manifest.tsv",
        list(rna_manifest[0]),
        rna_manifest,
    )
    decisions = [
        {
            "modality": "DNA",
            "decision": "duplicate_scope",
            "value": "CB_tag",
            "evidence": f"probe_duplicate_fraction={dna_duplicate_fraction}",
        },
        {
            "modality": "RNA",
            "decision": "strand",
            "value": "reverse",
            "evidence": (
                f"assigned_reads={reverse_reads}>{forward_reads};"
                f"exact_umis={reverse_umis}>{forward_umis}"
            ),
        },
        {
            "modality": "RNA",
            "decision": "feature",
            "value": "GeneFull_gene_body",
            "evidence": "single_nucleus_RNA",
        },
        {
            "modality": "RNA",
            "decision": "umi_deduplication",
            "value": "exact_within_sample_CB_gene",
            "evidence": "validated_probe",
        },
    ]
    atomic_tsv(
        out_dir / "postmapping_full_decisions.tsv",
        ["modality", "decision", "value", "evidence"],
        decisions,
    )
    shard_text = ", ".join(
        f"{index + 1}:{total}" for index, total in enumerate(shard_totals)
    )
    print(
        "[OK] prepared 96 DNA samples and 96 RNA files; "
        f"RNA shard input reads={shard_text}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
