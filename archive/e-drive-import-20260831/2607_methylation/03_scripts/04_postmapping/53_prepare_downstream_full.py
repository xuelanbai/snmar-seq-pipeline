#!/usr/bin/env python3
"""Step 53: validate full dedup/count outputs and prepare full jobs."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from pathlib import Path


EXPECTED = 96
DNA_DEDUP_PIPELINE_SIGNATURE = "snmar_dna_cb_markdup_corrected_20260731"


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
        raise ValueError(f"duplicate metric in: {path}")
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


def atomic_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as out:
            for line in lines:
                out.write(line.rstrip("\n") + "\n")
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def validate_dna(
    summary_path: Path, per_sample_path: Path
) -> list[dict[str, str]]:
    summary = read_metrics(summary_path)
    rows = read_tsv(per_sample_path)
    if (
        summary.get("status") != "validated"
        or summary.get("pipeline_signature")
        != DNA_DEDUP_PIPELINE_SIGNATURE
        or int(summary.get("completed_samples", "0")) != EXPECTED
        or len(rows) != EXPECTED
        or any(row["status"] != "validated" for row in rows)
    ):
        raise ValueError("DNA full CB-aware dedup validation failed")
    input_primary = sum(int(row["input_primary_mapped"]) for row in rows)
    removed = sum(
        int(row["duplicate_primary_records_removed"]) for row in rows
    )
    retained = sum(int(row["retained_primary_mapped"]) for row in rows)
    if input_primary - removed != retained:
        raise ValueError("DNA input - duplicate != retained")
    if (
        input_primary != int(summary["input_primary_mapped"])
        or removed != int(summary["duplicate_primary_records_removed"])
        or retained != int(summary["retained_primary_mapped"])
    ):
        raise ValueError("DNA full totals do not reconcile")
    for row in rows:
        if any(
            int(row[key]) != 0
            for key in (
                "missing_required_tag_records",
                "nonzero_XD_records",
                "duplicate_flag_records",
            )
        ):
            raise ValueError(f"DNA tag/flag validation failed: {row['sample']}")
        if not row["output_bam"].endswith(".cb_dedup.bam"):
            raise ValueError(f"unexpected DNA BAM: {row['output_bam']}")
    return sorted(rows, key=lambda row: int(row["sample_index"]))


def validate_rna(
    summary_path: Path,
    per_file_path: Path,
    cell_qc_path: Path,
) -> dict[str, str]:
    summary = read_metrics(summary_path)
    rows = read_tsv(per_file_path)
    expected = {
        "status": "validated",
        "completed_files": str(EXPECTED),
        "strand": "reverse",
        "counting_feature": "GeneFull_gene_body",
        "umi_deduplication": "exact_within_sample_CB_gene",
        "cell_filtering": "not_applied_raw_matrix",
    }
    if any(summary.get(key) != value for key, value in expected.items()):
        raise ValueError("RNA full raw-matrix validation failed")
    if len(rows) != EXPECTED or any(
        row["status"] != "validated" for row in rows
    ):
        raise ValueError("RNA per-file validation failed")
    exact_umis = sum(int(row["exact_umis"]) for row in rows)
    if exact_umis != int(summary["raw_matrix_umi_sum"]):
        raise ValueError("RNA per-file and raw-matrix UMI totals differ")
    if sum(
        int(row["human_exact_umis"]) + int(row["pig_exact_umis"])
        for row in rows
    ) != exact_umis:
        raise ValueError("RNA species UMI totals differ")

    cell_rows = 0
    raw_cells = 0
    cell_umis = 0
    with gzip.open(cell_qc_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cell_rows += 1
            raw_cells += int(row["in_raw_matrix"])
            cell_umis += int(row["exact_umis"])
    if cell_rows != sum(int(row["cells_tag_eligible"]) for row in rows):
        raise ValueError("RNA cell-QC row total differs")
    if raw_cells != int(summary["raw_matrix_cells"]):
        raise ValueError("RNA raw-cell total differs")
    if cell_umis != exact_umis:
        raise ValueError("RNA cell-QC UMI total differs")
    return summary


def read_fai(path: Path) -> list[dict[str, object]]:
    contigs = []
    with path.open("rt", encoding="utf-8") as handle:
        for order, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"malformed FAI row {order}: {path}")
            name = fields[0]
            length = int(fields[1])
            if name.startswith("human_"):
                species = "human"
            elif name.startswith("pig_"):
                species = "pig"
            else:
                raise ValueError(f"unexpected unprefixed contig: {name}")
            contigs.append(
                {
                    "reference_order": order,
                    "contig": name,
                    "length": length,
                    "species": species,
                }
            )
    if not contigs:
        raise ValueError(f"empty FAI: {path}")
    return contigs


def balance_contigs(
    contigs: list[dict[str, object]], shards: int
) -> list[list[dict[str, object]]]:
    if shards < 1 or len(contigs) < shards:
        raise ValueError("invalid contig-shard count")
    assignments: list[list[dict[str, object]]] = [[] for _ in range(shards)]
    loads = [0] * shards
    for row in sorted(
        contigs,
        key=lambda value: (
            -int(value["length"]),
            int(value["reference_order"]),
        ),
    ):
        shard = min(range(shards), key=lambda index: (loads[index], index))
        assignments[shard].append(row)
        loads[shard] += int(row["length"])
    for rows in assignments:
        rows.sort(key=lambda value: int(value["reference_order"]))
    return assignments


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
        ),
    )
    parser.add_argument("--joint-fasta", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--contig-shards", type=int, default=8)
    parser.add_argument("--rna-min-umis", type=int, default=10)
    parser.add_argument("--rna-min-genes", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    project = args.project
    if (
        args.contig_shards < 1
        or args.rna_min_umis < 1
        or args.rna_min_genes < 1
    ):
        raise ValueError("invalid full-processing setting")
    dna_root = project / "08_qc_summary/15_dna_cb_markdup_full"
    rna_root = project / "08_qc_summary/16_rna_reverse_genefull_full"
    raw_matrix = (
        project
        / "07_postmapping/rna_reverse_genefull_full/raw_feature_bc_matrix"
    )
    fasta = (
        args.joint_fasta
        or project
        / "02_reference/rna_human_pig_joint/human_pig.prefixed.fa"
    )
    out_dir = (
        args.out_dir or project / "08_qc_summary/17_downstream_full_inputs"
    )
    prerequisites = (
        dna_root / "dna_cb_markdup_full_summary.tsv",
        dna_root / "dna_cb_markdup_full_per_sample.tsv",
        rna_root / "rna_reverse_genefull_full_summary.tsv",
        rna_root / "rna_reverse_genefull_full_per_file.tsv",
        rna_root / "rna_reverse_cell_qc.tsv.gz",
        fasta,
        Path(str(fasta) + ".fai"),
        raw_matrix / "matrix.mtx.gz",
        raw_matrix / "barcodes.tsv.gz",
        raw_matrix / "features.tsv.gz",
        raw_matrix / "feature_metadata.tsv.gz",
    )
    for path in prerequisites:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty prerequisite: {path}")

    dna_rows = validate_dna(prerequisites[0], prerequisites[1])
    rna_summary = validate_rna(
        prerequisites[2], prerequisites[3], prerequisites[4]
    )
    contigs = read_fai(prerequisites[6])
    assignments = balance_contigs(contigs, args.contig_shards)

    manifest_rows = [
        {
            "sample_index": row["sample_index"],
            "run_id": row["run_id"],
            "sample": row["sample"],
            "group": row["group"],
            "retained_primary_mapped": row["retained_primary_mapped"],
            "input_bam": row["output_bam"],
        }
        for row in dna_rows
    ]
    atomic_tsv(
        out_dir / "dna_taps_cpg_full_manifest.tsv",
        list(manifest_rows[0]),
        manifest_rows,
    )

    shard_rows = []
    shard_dir = out_dir / "dna_contig_shards"
    for shard_index, rows in enumerate(assignments, 1):
        shard_file = shard_dir / f"shard_{shard_index:02d}.contigs.txt"
        atomic_lines(shard_file, [str(row["contig"]) for row in rows])
        for row in rows:
            shard_rows.append(
                {
                    "shard": shard_index,
                    "reference_order": row["reference_order"],
                    "contig": row["contig"],
                    "species": row["species"],
                    "length": row["length"],
                    "contig_file": shard_file,
                }
            )
    atomic_tsv(
        out_dir / "dna_contig_shards.tsv",
        [
            "shard",
            "reference_order",
            "contig",
            "species",
            "length",
            "contig_file",
        ],
        shard_rows,
    )
    if len({str(row["contig"]) for row in shard_rows}) != len(contigs):
        raise ValueError("contig assignment is not one-to-one")

    decisions: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "mode", "value": "direct_full_no_data_probe"},
        {"metric": "dna_samples", "value": len(manifest_rows)},
        {"metric": "dna_contigs", "value": len(contigs)},
        {"metric": "dna_contig_shards_per_sample", "value": args.contig_shards},
        {
            "metric": "dna_call_model",
            "value": "forward_ref_C_CT_reverse_ref_G_GA",
        },
        {
            "metric": "dna_overlap_policy",
            "value": "one_call_per_sample_CB_QNAME_CpG_highest_baseQ",
        },
        {
            "metric": "rna_raw_matrix_cells",
            "value": rna_summary["raw_matrix_cells"],
        },
        {
            "metric": "rna_raw_matrix_umi_sum",
            "value": rna_summary["raw_matrix_umi_sum"],
        },
        {"metric": "rna_min_exact_umis", "value": args.rna_min_umis},
        {"metric": "rna_min_genes", "value": args.rna_min_genes},
        {
            "metric": "rna_filter_rule",
            "value": "exact_umis>=minimum AND genes_detected>=minimum",
        },
    ]
    atomic_tsv(
        out_dir / "downstream_full_decisions.tsv",
        ["metric", "value"],
        decisions,
    )
    print(
        f"[OK] full inputs: DNA={len(manifest_rows)} samples, "
        f"contigs={len(contigs)} x {args.contig_shards} shards, "
        f"RNA cells={rna_summary['raw_matrix_cells']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
