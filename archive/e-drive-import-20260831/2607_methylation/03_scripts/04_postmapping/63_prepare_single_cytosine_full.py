#!/usr/bin/env python3
"""Step 63: prepare 768 DNA single-cytosine tasks and global indices."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable


EXPECTED_SAMPLES = 96
EXPECTED_SHARDS = 8
EXPECTED_FASTQ_PAIRS = 768


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
    path: Path,
    fields: list[str],
    rows: Iterable[dict[str, object]],
    gzip_output: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    try:
        opener = gzip.open if gzip_output else open
        with opener(
            temporary, "wt", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                delimiter="\t",
                lineterminator="\n",
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
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--min-mapq", type=int, default=30)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--exclude-read-ends", type=int, default=10)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.min_mapq < 0
        or args.min_baseq < 0
        or args.exclude_read_ends < 0
    ):
        raise ValueError("invalid single-cytosine calling setting")
    project = args.project
    downstream = project / "08_qc_summary/17_downstream_full_inputs"
    legacy_qc = project / "08_qc_summary/18_dna_taps_cpg_full"
    manifest_path = downstream / "dna_taps_cpg_full_manifest.tsv"
    contig_shards_path = downstream / "dna_contig_shards.tsv"
    legacy_summary_path = legacy_qc / "dna_taps_cpg_full_summary.tsv"
    legacy_per_sample_path = legacy_qc / "dna_taps_cpg_full_per_sample.tsv"
    legacy_per_cell_path = legacy_qc / "dna_taps_cpg_full_per_cell.tsv.gz"
    barcode_count_root = (
        project
        / "08_qc_summary/01_demultiplex/write_rescue_barcode_counts"
    )
    fasta = (
        args.fasta
        or project
        / "02_reference/rna_human_pig_joint/human_pig.prefixed.fa"
    )
    fai = Path(str(fasta) + ".fai")
    out_dir = (
        args.out_dir
        or project / "08_qc_summary/20_single_cytosine_inputs"
    )
    for path in (
        manifest_path,
        contig_shards_path,
        legacy_summary_path,
        legacy_per_sample_path,
        legacy_per_cell_path,
        barcode_count_root,
        fasta,
        fai,
    ):
        if (
            (path == barcode_count_root and not path.is_dir())
            or (
                path != barcode_count_root
                and (not path.is_file() or path.stat().st_size == 0)
            )
        ):
            raise FileNotFoundError(f"missing/empty prerequisite: {path}")

    legacy_summary = read_metrics(legacy_summary_path)
    if (
        legacy_summary.get("status") != "validated"
        or int(legacy_summary.get("completed_samples", "0"))
        != EXPECTED_SAMPLES
    ):
        raise ValueError("legacy DNA call summary is not validated")
    manifest = read_tsv(manifest_path)
    per_sample = read_tsv(legacy_per_sample_path)
    if len(manifest) != EXPECTED_SAMPLES or len(per_sample) != EXPECTED_SAMPLES:
        raise ValueError("expected 96 DNA samples")
    manifest_by_sample = {row["sample"]: row for row in manifest}
    per_sample_by_sample = {row["sample"]: row for row in per_sample}
    if (
        len(manifest_by_sample) != EXPECTED_SAMPLES
        or set(manifest_by_sample) != set(per_sample_by_sample)
        or sorted(int(row["sample_index"]) for row in manifest)
        != list(range(1, EXPECTED_SAMPLES + 1))
    ):
        raise ValueError("DNA manifest/per-sample identities differ")
    for row in manifest:
        bam = Path(row["input_bam"])
        bai_candidates = (
            Path(str(bam) + ".bai"),
            bam.with_suffix(".bai"),
        )
        if not bam.is_file() or bam.stat().st_size == 0:
            raise FileNotFoundError(bam)
        if not any(path.is_file() and path.stat().st_size > 0 for path in bai_candidates):
            raise FileNotFoundError(f"BAM index absent for {bam}")

    shard_rows = read_tsv(contig_shards_path)
    fai_rows: list[tuple[str, int]] = []
    with fai.open("rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"invalid FAI row: {line.rstrip()}")
            fai_rows.append((fields[0], int(fields[1])))
    contigs: dict[str, dict[str, object]] = {}
    shard_files: dict[int, Path] = {}
    for row in shard_rows:
        shard = int(row["shard"])
        order = int(row["reference_order"])
        contig = row["contig"]
        if shard < 1 or shard > EXPECTED_SHARDS:
            raise ValueError(f"invalid shard: {shard}")
        if contig in contigs:
            raise ValueError(f"duplicate contig: {contig}")
        contigs[contig] = {
            "contig_index": order - 1,
            "reference_order": order,
            "shard": shard,
            "species": row["species"],
            "contig": contig,
            "length": int(row["length"]),
        }
        path = Path(row["contig_file"])
        previous = shard_files.setdefault(shard, path)
        if previous != path:
            raise ValueError(f"shard {shard} has multiple contig files")
    if (
        len(contigs) != int(legacy_summary["contigs_per_sample"])
        or set(shard_files) != set(range(1, EXPECTED_SHARDS + 1))
        or len(fai_rows) != len(contigs)
    ):
        raise ValueError("contig/shard count does not reconcile")
    if sorted(int(row["contig_index"]) for row in contigs.values()) != list(
        range(len(contigs))
    ):
        raise ValueError("contig index is not zero-based/contiguous")
    for index, (contig, length) in enumerate(fai_rows):
        row = contigs.get(contig)
        if (
            row is None
            or int(row["contig_index"]) != index
            or int(row["length"]) != length
        ):
            raise ValueError("contig index differs from joint-reference FAI")
    for path in shard_files.values():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)

    barcode_files = sorted(
        barcode_count_root.glob("*.barcode_count.tsv")
    )
    if len(barcode_files) != EXPECTED_FASTQ_PAIRS:
        raise ValueError(
            f"expected {EXPECTED_FASTQ_PAIRS} barcode-count files, "
            f"found {len(barcode_files)}"
        )
    cells_by_id: dict[str, dict[str, object]] = {}
    for path in barcode_files:
        for row in read_tsv(path):
            sample = row["sample"]
            cb = row["barcode_seq"]
            count = int(row["count"])
            if sample not in manifest_by_sample or not cb or count < 1:
                raise ValueError(f"invalid demultiplex barcode row: {path}")
            cell_id = f"{sample}|{cb}"
            source = manifest_by_sample[sample]
            cell = cells_by_id.setdefault(
                cell_id,
                {
                    "sample_index": int(source["sample_index"]),
                    "sample": sample,
                    "group": source["group"],
                    "cell_id": cell_id,
                    "CB": cb,
                    "barcode_name": row["barcode_name"],
                    "demultiplex_assigned_reads_all_XD": 0,
                    "_pair_ids": set(),
                    "legacy_covered_template_cpg_calls": 0,
                    "legacy_modified_template_cpg_calls": 0,
                    "legacy_bins_with_calls": 0,
                },
            )
            if (
                cell["sample"] != sample
                or cell["CB"] != cb
                or cell["barcode_name"] != row["barcode_name"]
            ):
                raise ValueError(f"inconsistent barcode identity: {cell_id}")
            cell["demultiplex_assigned_reads_all_XD"] = int(
                cell["demultiplex_assigned_reads_all_XD"]
            ) + count
            pair_ids = cell["_pair_ids"]
            assert isinstance(pair_ids, set)
            pair_ids.add(row["pair_id"])

    sample_cell_counts: Counter[str] = Counter()
    sample_covered: Counter[str] = Counter()
    sample_modified: Counter[str] = Counter()
    legacy_seen: set[str] = set()
    with gzip.open(
        legacy_per_cell_path, "rt", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = row["sample"]
            cell_id = row["cell_id"]
            if sample not in manifest_by_sample or cell_id in legacy_seen:
                raise ValueError(f"unknown/duplicate DNA cell: {cell_id}")
            legacy_seen.add(cell_id)
            cell = cells_by_id.get(cell_id)
            if cell is None or cell["CB"] != row["CB"]:
                raise ValueError(
                    f"legacy CpG cell absent from demultiplex CB index: {cell_id}"
                )
            covered = int(row["covered_template_cpg_calls"])
            modified = int(row["modified_template_cpg_calls"])
            cell["legacy_covered_template_cpg_calls"] = covered
            cell["legacy_modified_template_cpg_calls"] = modified
            cell["legacy_bins_with_calls"] = int(row["bins_with_calls"])
            sample_cell_counts[sample] += 1
            sample_covered[sample] += covered
            sample_modified[sample] += modified
    cells = []
    for cell in cells_by_id.values():
        pair_ids = cell.pop("_pair_ids")
        assert isinstance(pair_ids, set)
        cell["demultiplex_pair_files"] = len(pair_ids)
        cells.append(cell)
    cells.sort(
        key=lambda row: (
            int(row["sample_index"]),
            str(row["cell_id"]),
        )
    )
    for cell_index, row in enumerate(cells):
        row["cell_index"] = cell_index
    if len(legacy_seen) != int(legacy_summary["cells_with_cpg_calls_sum"]):
        raise ValueError("legacy DNA cell count does not reconcile")
    for sample, qc in per_sample_by_sample.items():
        if (
            sample_cell_counts[sample] != int(qc["cells_with_cpg_calls"])
            or sample_covered[sample]
            != int(qc["covered_template_cpg_calls"])
            or sample_modified[sample]
            != int(qc["modified_template_cpg_calls"])
        ):
            raise ValueError(f"{sample}: legacy totals differ")

    task_rows = []
    task_index = 0
    matrix_root = (
        project / "07_postmapping/complete_matrices/dna_single_cytosine"
    )
    for sample in sorted(
        manifest_by_sample,
        key=lambda value: int(manifest_by_sample[value]["sample_index"]),
    ):
        source = manifest_by_sample[sample]
        for shard in range(1, EXPECTED_SHARDS + 1):
            task_index += 1
            task_rows.append(
                {
                    "task_index": task_index,
                    "sample_index": int(source["sample_index"]),
                    "sample": sample,
                    "group": source["group"],
                    "shard": shard,
                    "input_bam": source["input_bam"],
                    "reference_fasta": fasta,
                    "contig_file": shard_files[shard],
                    "cell_index_path": out_dir
                    / "dna_single_cytosine_cell_index.tsv.gz",
                    "contig_index_path": out_dir
                    / "dna_single_cytosine_contig_index.tsv",
                    "output_dir": matrix_root
                    / "shards"
                    / sample
                    / f"shard_{shard:02d}",
                    "minimum_mapq": args.min_mapq,
                    "minimum_baseq": args.min_baseq,
                    "excluded_bases_each_read_end": args.exclude_read_ends,
                }
            )
    if len(task_rows) != EXPECTED_SAMPLES * EXPECTED_SHARDS:
        raise ValueError("single-cytosine task count is not 768")

    contig_output = sorted(
        contigs.values(), key=lambda row: int(row["contig_index"])
    )
    atomic_tsv(
        out_dir / "dna_single_cytosine_cell_index.tsv.gz",
        ["cell_index"]
        + [key for key in cells[0] if key != "cell_index"],
        cells,
        gzip_output=True,
    )
    atomic_tsv(
        out_dir / "dna_single_cytosine_contig_index.tsv",
        list(contig_output[0]),
        contig_output,
    )
    atomic_tsv(
        out_dir / "dna_single_cytosine_tasks.tsv",
        list(task_rows[0]),
        task_rows,
    )
    summary_rows: list[dict[str, object]] = [
        {"metric": "status", "value": "validated"},
        {"metric": "mode", "value": "single_cytosine_context_aware"},
        {"metric": "dna_samples", "value": EXPECTED_SAMPLES},
        {"metric": "contig_shards_per_sample", "value": EXPECTED_SHARDS},
        {"metric": "single_cytosine_tasks", "value": len(task_rows)},
        {"metric": "indexed_cells", "value": len(cells)},
        {
            "metric": "legacy_cells_with_CpG_calls",
            "value": len(legacy_seen),
        },
        {
            "metric": "cell_axis_source",
            "value": (
                "768_write_rescue_barcode_count_files_all_XD_candidates;"
                "caller_retains_XD0_only"
            ),
        },
        {"metric": "indexed_contigs", "value": len(contigs)},
        {"metric": "reference_fasta", "value": fasta},
        {"metric": "minimum_mapq", "value": args.min_mapq},
        {"metric": "minimum_baseq", "value": args.min_baseq},
        {"metric": "excluded_bases_each_read_end", "value": args.exclude_read_ends},
        {
            "metric": "context_model",
            "value": "WCG=methylation;GCH=accessibility;GCG=excluded",
        },
        {
            "metric": "strand_model",
            "value": (
                "plus_reference_C_CT;"
                "minus_reference_G_GA_with_strand_retained"
            ),
        },
        {
            "metric": "overlap_policy",
            "value": (
                "one_call_per_sample_CB_QNAME_cytosine_strand_context_"
                "highest_baseQ"
            ),
        },
        {
            "metric": "legacy_covered_calls_without_end_exclusion",
            "value": legacy_summary["covered_template_cpg_calls"],
        },
        {
            "metric": "legacy_modified_calls_without_end_exclusion",
            "value": legacy_summary["modified_template_cpg_calls"],
        },
        {
            "metric": "missing_cell_cytosine_semantics",
            "value": "NA_not_zero_signal",
        },
    ]
    atomic_tsv(
        out_dir / "single_cytosine_prepare_summary.tsv",
        ["metric", "value"],
        summary_rows,
    )
    print(
        f"[OK] single-cytosine inputs: {len(task_rows)} tasks, "
        f"{len(cells)} cells, {len(contigs)} contigs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
