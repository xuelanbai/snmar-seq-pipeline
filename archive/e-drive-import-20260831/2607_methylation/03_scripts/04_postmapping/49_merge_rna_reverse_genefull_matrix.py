#!/usr/bin/env python3
"""Step 49: validate 96 RNA counts and create a raw sparse UMI matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
from collections import Counter
from pathlib import Path


PROJECT_DEFAULT = Path(
    "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
)
TRAILING_NUMBER_RE = re.compile(r"([0-9]+)$")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"no rows in TSV: {path}")
    return rows


def read_gzip_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty gzip TSV: {path}")
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_number(sample: str) -> int:
    match = TRAILING_NUMBER_RE.search(sample)
    if not match:
        raise ValueError(f"sample lacks trailing number: {sample}")
    return int(match.group(1))


def atomic_gzip_lines(path: Path, lines) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with gzip.open(tmp, "wt", encoding="utf-8", compresslevel=1) as out:
            for line in lines:
                out.write(line)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


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
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--by-file-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument(
        "--expected-files",
        type=int,
        default=96,
        help="Keep the production default at 96; smaller values are for tests.",
    )
    args = parser.parse_args()

    project = args.project
    manifest = args.manifest or (
        project
        / "08_qc_summary/14_postmapping_full_inputs/"
        "rna_reverse_genefull_full_manifest.tsv"
    )
    by_file_root = args.by_file_root or (
        project / "07_postmapping/rna_reverse_genefull_full/by_file"
    )
    output_root = args.output_root or (
        project
        / "07_postmapping/rna_reverse_genefull_full/raw_feature_bc_matrix"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/16_rna_reverse_genefull_full"
    )

    if args.expected_files < 1:
        parser.error("--expected-files must be positive")
    manifest_rows = read_tsv(manifest)
    if len(manifest_rows) != args.expected_files:
        raise ValueError(
            f"expected {args.expected_files} manifest rows, "
            f"found {len(manifest_rows)}"
        )
    features: dict[str, tuple[str, str]] = {}
    cells: dict[str, dict[str, str]] = {}
    matrix_counts: dict[tuple[str, str], int] = {}
    per_file: list[dict[str, object]] = []
    aggregate: Counter[str] = Counter()

    for row in manifest_rows:
        pair_id = row["pair_id"]
        file_root = by_file_root / pair_id
        summary_path = file_root / f"{pair_id}.summary.json"
        if not summary_path.is_file() or summary_path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty summary: {summary_path}")
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        if (
            payload.get("status") != "validated"
            or payload.get("strand") != "reverse"
            or payload.get("counting_feature") != "GeneFull_gene_body"
            or payload.get("umi_deduplication")
            != "exact_within_sample_CB_gene"
        ):
            raise ValueError(f"invalid full RNA result: {summary_path}")
        counters = payload["counters"]
        if int(counters["primary_records"]) != int(row["star_input_reads"]):
            raise ValueError(f"{pair_id}: primary-record reconciliation failed")
        if int(counters["primary_mapped"]) != int(
            row["primary_mapped_reads"]
        ):
            raise ValueError(f"{pair_id}: primary-mapped reconciliation failed")
        if int(counters.get("missing_required_tags", 0)) != 0:
            raise ValueError(f"{pair_id}: required tags are missing")
        if int(counters.get("nonzero_XD", 0)) != 0:
            raise ValueError(f"{pair_id}: non-strict XD observed")

        counts_path = Path(payload["cell_gene_counts"])
        cell_qc_path = Path(payload["cell_qc"])
        count_rows = read_gzip_tsv(counts_path)
        qc_rows = read_gzip_tsv(cell_qc_path)
        per_file_umis = 0
        per_file_species_umis: Counter[str] = Counter()
        for count_row in count_rows:
            if count_row["pair_id"] != pair_id:
                raise ValueError(f"{pair_id}: count-row pair_id mismatch")
            cell_id = count_row["cell_id"]
            gene_id = count_row["gene_id"]
            umi_count = int(count_row["exact_umis"])
            if umi_count <= 0:
                raise ValueError(f"{pair_id}: non-positive matrix count")
            token = (cell_id, gene_id)
            if token in matrix_counts:
                raise ValueError(
                    f"duplicate cell/gene matrix entry: {cell_id}, {gene_id}"
                )
            matrix_counts[token] = umi_count
            per_file_umis += umi_count
            per_file_species_umis[count_row["species"]] += umi_count
            feature_value = (
                count_row["gene_name"],
                count_row["species"],
            )
            previous = features.setdefault(gene_id, feature_value)
            if previous != feature_value:
                raise ValueError(f"inconsistent gene metadata: {gene_id}")

        qc_umi_sum = 0
        for qc_row in qc_rows:
            if qc_row["pair_id"] != pair_id:
                raise ValueError(f"{pair_id}: cell-QC pair_id mismatch")
            cell_id = qc_row["cell_id"]
            if cell_id in cells:
                raise ValueError(f"duplicate cell_id across files: {cell_id}")
            cells[cell_id] = qc_row
            qc_umi_sum += int(qc_row["exact_umis"])
        expected_umis = int(payload["exact_umis"])
        if per_file_umis != expected_umis or qc_umi_sum != expected_umis:
            raise ValueError(f"{pair_id}: exact-UMI reconciliation failed")

        assigned_reads = int(counters.get("gene_assigned", 0))
        eligible = int(counters.get("tag_eligible_unique_primary", 0))
        file_result = {
            "manifest_index": row["manifest_index"],
            "shard": row["shard"],
            "pair_id": pair_id,
            "sample": row["sample"],
            "star_input_reads": int(row["star_input_reads"]),
            "primary_mapped_reads": int(row["primary_mapped_reads"]),
            "unique_primary_reads": int(counters.get("unique_primary", 0)),
            "tag_eligible_unique_primary": eligible,
            "gene_assigned_reads": assigned_reads,
            "gene_assignment_fraction": (
                assigned_reads / eligible if eligible else 0.0
            ),
            "gene_ambiguous_reads": int(
                counters.get("gene_ambiguous", 0)
            ),
            "gene_no_feature_reads": int(
                counters.get("gene_no_feature", 0)
            ),
            "exact_umis": expected_umis,
            "human_exact_umis": per_file_species_umis["human"],
            "pig_exact_umis": per_file_species_umis["pig"],
            "cells_tag_eligible": int(payload["cells_tag_eligible"]),
            "cells_with_exact_umi": int(payload["cells_with_exact_umi"]),
            "detected_cell_gene_pairs": int(
                payload["detected_cell_gene_pairs"]
            ),
            "status": "validated",
        }
        per_file.append(file_result)
        for key in (
            "star_input_reads",
            "primary_mapped_reads",
            "unique_primary_reads",
            "tag_eligible_unique_primary",
            "gene_assigned_reads",
            "gene_ambiguous_reads",
            "gene_no_feature_reads",
            "exact_umis",
            "human_exact_umis",
            "pig_exact_umis",
        ):
            aggregate[key] += int(file_result[key])

    matrix_cells = {
        cell_id for cell_id, _gene_id in matrix_counts
    }
    cell_order = sorted(
        matrix_cells,
        key=lambda cell_id: (
            sample_number(cells[cell_id]["sample"]),
            cells[cell_id]["CB"],
        ),
    )
    feature_order = sorted(
        features,
        key=lambda gene_id: (features[gene_id][1], gene_id),
    )
    cell_index = {
        cell_id: index for index, cell_id in enumerate(cell_order, 1)
    }
    feature_index = {
        gene_id: index for index, gene_id in enumerate(feature_order, 1)
    }
    matrix_total = sum(matrix_counts.values())
    if matrix_total != aggregate["exact_umis"]:
        raise ValueError("matrix total does not equal aggregated exact UMIs")

    def matrix_lines():
        yield "%%MatrixMarket matrix coordinate integer general\n"
        yield "% snMAR reverse GeneFull exact-UMI raw matrix\n"
        yield (
            f"{len(feature_order)} {len(cell_order)} "
            f"{len(matrix_counts)}\n"
        )
        for (cell_id, gene_id), count in sorted(
            matrix_counts.items(),
            key=lambda item: (
                cell_index[item[0][0]],
                feature_index[item[0][1]],
            ),
        ):
            yield (
                f"{feature_index[gene_id]} {cell_index[cell_id]} {count}\n"
            )

    atomic_gzip_lines(output_root / "matrix.mtx.gz", matrix_lines())
    atomic_gzip_lines(
        output_root / "barcodes.tsv.gz",
        (f"{cell_id}\n" for cell_id in cell_order),
    )
    atomic_gzip_lines(
        output_root / "features.tsv.gz",
        (
            f"{gene_id}\t{features[gene_id][0]}\tGene Expression\n"
            for gene_id in feature_order
        ),
    )
    atomic_gzip_lines(
        output_root / "feature_metadata.tsv.gz",
        (
            "gene_id\tgene_name\tspecies\n",
            *(
                f"{gene_id}\t{features[gene_id][0]}\t"
                f"{features[gene_id][1]}\n"
                for gene_id in feature_order
            ),
        ),
    )

    cell_qc_fields = [
        "pair_id",
        "sample",
        "cell_id",
        "CB",
        "CN",
        "tag_eligible_unique_primary",
        "gene_assigned_reads",
        "gene_assignment_fraction",
        "exact_umis",
        "genes_detected",
        "human_assigned_reads",
        "pig_assigned_reads",
        "human_exact_umis",
        "pig_exact_umis",
        "in_raw_matrix",
    ]

    def cell_qc_lines():
        yield "\t".join(cell_qc_fields) + "\n"
        for cell_id in sorted(
            cells,
            key=lambda value: (
                sample_number(cells[value]["sample"]),
                cells[value]["CB"],
            ),
        ):
            row = dict(cells[cell_id])
            row["in_raw_matrix"] = "1" if cell_id in matrix_cells else "0"
            yield "\t".join(row[field] for field in cell_qc_fields) + "\n"

    atomic_gzip_lines(qc_root / "rna_reverse_cell_qc.tsv.gz", cell_qc_lines())
    per_file.sort(key=lambda row: int(row["manifest_index"]))
    atomic_tsv(
        qc_root / "rna_reverse_genefull_full_per_file.tsv",
        list(per_file[0]),
        per_file,
    )
    summary_rows = [
        {"metric": "status", "value": "validated"},
        {"metric": "completed_files", "value": args.expected_files},
        {"metric": "strand", "value": "reverse"},
        {"metric": "counting_feature", "value": "GeneFull_gene_body"},
        {
            "metric": "umi_deduplication",
            "value": "exact_within_sample_CB_gene",
        },
        {"metric": "cell_filtering", "value": "not_applied_raw_matrix"},
        *(
            {"metric": key, "value": value}
            for key, value in aggregate.items()
        ),
        {
            "metric": "gene_assignment_fraction",
            "value": (
                aggregate["gene_assigned_reads"]
                / aggregate["tag_eligible_unique_primary"]
            ),
        },
        {"metric": "raw_matrix_cells", "value": len(cell_order)},
        {"metric": "raw_matrix_features", "value": len(feature_order)},
        {"metric": "raw_matrix_nonzero_entries", "value": len(matrix_counts)},
        {"metric": "raw_matrix_umi_sum", "value": matrix_total},
    ]
    atomic_tsv(
        qc_root / "rna_reverse_genefull_full_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] raw RNA matrix: genes={len(feature_order)} "
        f"cells={len(cell_order)} nnz={len(matrix_counts)} "
        f"umis={matrix_total}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
