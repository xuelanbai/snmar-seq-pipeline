#!/usr/bin/env python3
"""Step 118: build a complete exon-supported RNA salvage matrix.

This creates a new, explicitly low-depth matrix from the exon_reverse counts
emitted by step 113.  It never overwrites the original GeneFull matrix.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_rna_exon_salvage_matrix_20260807"
RECOUNT_SIGNATURE = "snmar_rna_expression_evidence_recount_20260806"
MODE = "exon_reverse"


def opener_for(path: Path):
    return gzip.open if path.suffix == ".gz" else open


def iter_tsv(path: Path) -> Iterator[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with opener_for(path)(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"TSV has no header: {path}")
        yield from reader


def read_tsv(path: Path) -> list[dict[str, str]]:
    return list(iter_tsv(path))


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"sample lacks trailing physical-well number: {value}")
    return int(match.group(1))


def dominant_species(human: int, pig: int) -> tuple[str, float | None]:
    total = human + pig
    if total == 0:
        return "no_signal", None
    if human == pig:
        return "tie", 0.5
    return ("human" if human > pig else "pig"), max(human, pig) / total


def atomic_tsv(
    path: Path, fields: list[str], rows: Iterable[dict[str, object]], gzip_output: bool = False
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        opener = gzip.open if gzip_output else open
        with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def atomic_gzip_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=1) as handle:
            for line in lines:
                handle.write(line)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--by-file-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--expected-wells", type=int, default=96)
    parser.add_argument("--min-umis", type=int, default=5)
    parser.add_argument("--min-genes", type=int, default=3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.expected_wells != 96 or args.min_umis < 1 or args.min_genes < 1:
        raise ValueError("invalid expected-well count or filter")
    project = args.project
    manifest_path = args.manifest or (
        project / "08_qc_summary/14_postmapping_full_inputs/rna_reverse_genefull_full_manifest.tsv"
    )
    by_file_root = args.by_file_root or (
        project / "10_downstream_analysis/05_rna_expression_evidence_audit/by_file"
    )
    output_root = args.output_root or (
        project / "10_downstream_analysis/06_rna_exon_salvage_matrix"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/36_rna_exon_salvage_matrix"
    )
    manifests = read_tsv(manifest_path)
    if len(manifests) != args.expected_wells:
        raise ValueError(f"expected 96 manifest rows, found {len(manifests)}")

    features: dict[str, tuple[str, str]] = {}
    cells: dict[str, dict[str, object]] = {}
    counts: dict[tuple[str, str], int] = {}
    expected_total = 0
    seen_wells: set[int] = set()
    for manifest in manifests:
        pair_id = manifest["pair_id"]
        sample = manifest["sample"]
        well_index = trailing_number(sample)
        if well_index in seen_wells:
            raise ValueError(f"duplicate well: {well_index}")
        seen_wells.add(well_index)
        root = by_file_root / pair_id
        summary_path = root / f"{pair_id}.evidence_mode_summary.json"
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
        if (
            payload.get("status") != "validated"
            or payload.get("pipeline_signature") != RECOUNT_SIGNATURE
            or payload.get("pair_id") != pair_id
            or payload.get("sample") != sample
        ):
            raise ValueError(f"invalid recount summary: {summary_path}")
        expected_file_total = int(payload["mode_totals"][MODE]["exact_umis"])
        observed_file_total = 0
        counts_path = Path(payload["counts"])
        for row in iter_tsv(counts_path):
            if row["mode"] != MODE:
                continue
            if row["pair_id"] != pair_id or row["sample"] != sample:
                raise ValueError(f"count identity mismatch: {pair_id}")
            cell_id = row["cell_id"]
            gene_id = row["gene_id"]
            value = int(row["exact_umis"])
            if value <= 0 or value > int(row["assigned_reads"]):
                raise ValueError(f"invalid exon UMI count: {row}")
            key = (cell_id, gene_id)
            if key in counts:
                raise ValueError(f"duplicate cell/gene count: {key}")
            counts[key] = value
            observed_file_total += value
            feature = (row["gene_name"], row["species"])
            if gene_id in features and features[gene_id] != feature:
                raise ValueError(f"inconsistent feature metadata: {gene_id}")
            features[gene_id] = feature
            if cell_id not in cells:
                cells[cell_id] = {
                    "cell_id": cell_id,
                    "sample": sample,
                    "well_index": well_index,
                    "well_group": "GM12878_only_wells" if well_index <= 48 else "GM12878_PK15_mix_wells",
                    "CB": row["CB"],
                    "CN": row["CN"],
                }
        if observed_file_total != expected_file_total:
            raise ValueError(
                f"{pair_id}: exon UMI total {observed_file_total} != {expected_file_total}"
            )
        expected_total += expected_file_total

    if seen_wells != set(range(1, 97)):
        raise ValueError("recount does not cover physical wells 1-96 exactly once")
    if not counts or sum(counts.values()) != expected_total:
        raise ValueError("complete exon matrix UMI reconciliation failed")

    cell_umis: Counter[str] = Counter()
    cell_genes: Counter[str] = Counter()
    cell_species_umis: Counter[tuple[str, str]] = Counter()
    for (cell_id, gene_id), value in counts.items():
        species = features[gene_id][1]
        if species not in {"human", "pig"}:
            raise ValueError(f"unexpected feature species: {gene_id} {species}")
        cell_umis[cell_id] += value
        cell_genes[cell_id] += 1
        cell_species_umis[(cell_id, species)] += value

    cell_order = sorted(cells, key=lambda cell: (int(cells[cell]["well_index"]), str(cells[cell]["CB"])))
    feature_order = sorted(features, key=lambda gene: (features[gene][1], gene))
    filtered_cells = [
        cell for cell in cell_order if cell_umis[cell] >= args.min_umis and cell_genes[cell] >= args.min_genes
    ]
    raw_cell_index = {cell: index for index, cell in enumerate(cell_order, 1)}
    filtered_cell_index = {cell: index for index, cell in enumerate(filtered_cells, 1)}
    feature_index = {gene: index for index, gene in enumerate(feature_order, 1)}

    def matrix_lines(selected_index: dict[str, int]):
        selected_entries = [
            ((cell, gene), value)
            for (cell, gene), value in counts.items()
            if cell in selected_index
        ]
        yield "%%MatrixMarket matrix coordinate integer general\n"
        yield f"% snMAR exon_reverse salvage matrix; {PIPELINE_SIGNATURE}\n"
        yield f"{len(feature_order)} {len(selected_index)} {len(selected_entries)}\n"
        for (cell, gene), value in sorted(
            selected_entries,
            key=lambda item: (selected_index[item[0][0]], feature_index[item[0][1]]),
        ):
            yield f"{feature_index[gene]} {selected_index[cell]} {value}\n"

    raw_root = output_root / "raw_feature_bc_matrix"
    filtered_root = output_root / f"filtered_feature_bc_matrix_minUMI{args.min_umis}_minGene{args.min_genes}"
    for root, selected, index in (
        (raw_root, cell_order, raw_cell_index),
        (filtered_root, filtered_cells, filtered_cell_index),
    ):
        atomic_gzip_lines(root / "matrix.mtx.gz", matrix_lines(index))
        atomic_gzip_lines(root / "barcodes.tsv.gz", (f"{cell}\n" for cell in selected))
        atomic_gzip_lines(
            root / "features.tsv.gz",
            (
                f"{gene}\t{features[gene][0]}\tGene Expression\n"
                for gene in feature_order
            ),
        )
        atomic_gzip_lines(
            root / "feature_metadata.tsv.gz",
            iter(
                ["gene_id\tgene_name\tspecies\n"]
                + [f"{gene}\t{features[gene][0]}\t{features[gene][1]}\n" for gene in feature_order]
            ),
        )

    metadata_rows: list[dict[str, object]] = []
    for cell in cell_order:
        human = cell_species_umis[(cell, "human")]
        pig = cell_species_umis[(cell, "pig")]
        label, major = dominant_species(human, pig)
        metadata_rows.append(
            {
                **cells[cell],
                "exact_umis": cell_umis[cell],
                "genes_detected": cell_genes[cell],
                "human_exact_umis": human,
                "pig_exact_umis": pig,
                "dominant_species": label,
                "major_species_fraction": major if major is not None else "NA",
                "passes_filter": int(cell in filtered_cell_index),
                "filter_rule": f"UMI>={args.min_umis};gene>={args.min_genes}",
            }
        )
    metadata_fields = list(metadata_rows[0].keys())
    atomic_tsv(output_root / "rna_exon_salvage_cell_metadata.tsv.gz", metadata_fields, metadata_rows, gzip_output=True)

    try:
        import h5py
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("h5py and numpy are required for complete matrix packaging") from exc
    h5_path = output_root / "rna_exon_salvage_complete_sparse.h5"
    temporary_h5 = h5_path.with_name(h5_path.name + f".tmp.{os.getpid()}")
    try:
        with h5py.File(temporary_h5, "w") as handle:
            handle.attrs["pipeline_signature"] = PIPELINE_SIGNATURE
            handle.attrs["mode"] = MODE
            handle.attrs["matrix_orientation"] = "feature_by_cell"
            handle.create_dataset("feature_id", data=np.asarray(feature_order, dtype="S"))
            for group_name, selected, index in (
                ("raw", cell_order, raw_cell_index),
                ("filtered", filtered_cells, filtered_cell_index),
            ):
                group = handle.create_group(group_name)
                selected_entries = [
                    (feature_index[gene] - 1, index[cell] - 1, value)
                    for (cell, gene), value in counts.items()
                    if cell in index
                ]
                group.attrs["shape"] = (len(feature_order), len(selected))
                group.create_dataset("cell_id", data=np.asarray(selected, dtype="S"))
                group.create_dataset("row", data=np.asarray([row for row, _col, _value in selected_entries], dtype=np.int64), compression="gzip")
                group.create_dataset("col", data=np.asarray([col for _row, col, _value in selected_entries], dtype=np.int64), compression="gzip")
                group.create_dataset("data", data=np.asarray([value for _row, _col, value in selected_entries], dtype=np.int64), compression="gzip")
        os.replace(temporary_h5, h5_path)
    finally:
        temporary_h5.unlink(missing_ok=True)

    filtered_set = set(filtered_cells)
    per_well_rows: list[dict[str, object]] = []
    for well in range(1, 97):
        raw_cells = [cell for cell in cell_order if int(cells[cell]["well_index"]) == well]
        kept = [cell for cell in raw_cells if cell in filtered_set]
        per_well_rows.append(
            {
                "well_index": well,
                "well_group": "GM12878_only_wells" if well <= 48 else "GM12878_PK15_mix_wells",
                "raw_exon_cells": len(raw_cells),
                "filtered_exon_cells": len(kept),
                "filtered_exon_umis": sum(cell_umis[cell] for cell in kept),
                "filtered_human_cells": sum(dominant_species(cell_species_umis[(cell, 'human')], cell_species_umis[(cell, 'pig')])[0] == "human" for cell in kept),
                "filtered_pig_cells": sum(dominant_species(cell_species_umis[(cell, 'human')], cell_species_umis[(cell, 'pig')])[0] == "pig" for cell in kept),
            }
        )
    qc_root.mkdir(parents=True, exist_ok=True)
    atomic_tsv(qc_root / "rna_exon_salvage_matrix_per_well.tsv", list(per_well_rows[0].keys()), per_well_rows)
    filtered_umis = sum(cell_umis[cell] for cell in filtered_cells)
    summary_rows = [
        {"metric": "pipeline_signature", "value": PIPELINE_SIGNATURE, "interpretation": "validated implementation signature"},
        {"metric": "status", "value": "validated", "interpretation": "96 recounts and complete sparse matrix reconcile"},
        {"metric": "mode", "value": MODE, "interpretation": "reverse-strand unique exon assignment"},
        {"metric": "raw_cells", "value": len(cell_order), "interpretation": "cells with at least one exon-supported UMI"},
        {"metric": "raw_features", "value": len(feature_order), "interpretation": "genes with at least one exon-supported UMI"},
        {"metric": "raw_exact_umis", "value": expected_total, "interpretation": "all exon-supported exact UMIs"},
        {"metric": "filter_rule", "value": f"UMI>={args.min_umis};gene>={args.min_genes}", "interpretation": "exploratory low-depth salvage filter"},
        {"metric": "filtered_cells", "value": len(filtered_cells), "interpretation": "cells entering exploratory RNA UMAP"},
        {"metric": "filtered_exact_umis", "value": filtered_umis, "interpretation": "UMIs retained after cell filtering"},
        {"metric": "filtered_UMI_retention", "value": filtered_umis / expected_total, "interpretation": "fraction of exon-supported UMIs retained"},
        {"metric": "formal_expression_clustering_ready", "value": 0, "interpretation": "depth remains below conventional expression-clustering QC"},
        {"metric": "intended_use", "value": "exploratory_UMAP_and_cross_modal_species_labels", "interpretation": "not a replacement for a high-quality RNA expression matrix"},
        {"metric": "complete_h5", "value": str(h5_path), "interpretation": "feature-by-cell COO groups raw and filtered"},
    ]
    atomic_tsv(qc_root / "rna_exon_salvage_matrix_summary.tsv", ["metric", "value", "interpretation"], summary_rows)
    print(f"[OK] exon salvage matrix: raw={len(cell_order)} filtered={len(filtered_cells)} UMIs={filtered_umis}/{expected_total}")
    print(f"[OK] complete HDF5: {h5_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
