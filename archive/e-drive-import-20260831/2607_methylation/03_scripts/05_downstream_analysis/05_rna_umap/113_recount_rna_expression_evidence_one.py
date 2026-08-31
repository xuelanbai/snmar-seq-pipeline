#!/usr/bin/env python3
"""Step 113: recount one RNA BAM with stricter transcript evidence.

Three modes are emitted from one BAM scan:

* genefull: the existing reverse-strand gene-body definition (reconciliation);
* exon_reverse: reverse-strand reads overlapping exons of exactly one gene;
* junction_genefull: GeneFull-assigned reads whose CIGAR contains an N operation.

This is a salvage audit, not an in-place replacement of the validated raw RNA
matrix.  The per-cell and per-gene outputs allow a later mode to be selected
without scanning the BAM files again.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


PIPELINE_SIGNATURE = "snmar_rna_expression_evidence_recount_20260806"
MODES = ("genefull", "exon_reverse", "junction_genefull")
ATTRIBUTE_RE = re.compile(r'(?:^|;\s*)([A-Za-z0-9_]+) "([^"]*)"')
CIGAR_RE = re.compile(r"([0-9]+)([MIDNSHP=X])")


@dataclass(frozen=True, slots=True)
class Gene:
    gene_id: str
    gene_name: str
    contig: str
    start: int
    end: int
    strand: str
    species: str


class IntervalBins:
    def __init__(self, bin_size: int = 100_000) -> None:
        self.bin_size = bin_size
        self.by_contig: dict[str, dict[int, list[tuple[int, int, str]]]] = defaultdict(
            lambda: defaultdict(list)
        )

    def add(self, contig: str, start: int, end: int, gene_id: str) -> None:
        if start < 0 or end <= start:
            raise ValueError(f"invalid interval: {contig}:{start}-{end} {gene_id}")
        first = start // self.bin_size
        last = (end - 1) // self.bin_size
        interval = (start, end, gene_id)
        for index in range(first, last + 1):
            self.by_contig[contig][index].append(interval)

    def query(self, contig: str, blocks: list[tuple[int, int]]) -> set[str]:
        bins = self.by_contig.get(contig)
        if bins is None:
            return set()
        hits: set[str] = set()
        for block_start, block_end in blocks:
            if block_end <= block_start:
                continue
            tested: set[tuple[int, int, str]] = set()
            first = block_start // self.bin_size
            last = (block_end - 1) // self.bin_size
            for index in range(first, last + 1):
                for interval in bins.get(index, ()):
                    if interval in tested:
                        continue
                    tested.add(interval)
                    start, end, gene_id = interval
                    if start < block_end and block_start < end:
                        hits.add(gene_id)
        return hits


@dataclass(slots=True)
class GeneModels:
    genes: dict[str, Gene]
    bodies: IntervalBins
    exons: IntervalBins
    exon_intervals: int


def parse_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTRIBUTE_RE.findall(text)}


def species_from(contig: str, gene_id: str) -> str:
    if contig.startswith("human_") or gene_id.startswith("HUMAN_"):
        return "human"
    if contig.startswith("pig_") or gene_id.startswith("PIG_"):
        return "pig"
    return "other"


def load_gene_models(gtf: Path) -> GeneModels:
    genes: dict[str, Gene] = {}
    exon_sets: dict[str, set[tuple[int, int]]] = defaultdict(set)
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gtf}:{line_no}: expected 9 GTF fields")
            feature = fields[2]
            if feature not in {"gene", "exon"}:
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                raise ValueError(f"{gtf}:{line_no}: missing gene_id")
            start = int(fields[3]) - 1
            end = int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"{gtf}:{line_no}: invalid interval")
            if feature == "gene":
                if gene_id in genes:
                    raise ValueError(f"duplicate gene_id: {gene_id}")
                genes[gene_id] = Gene(
                    gene_id=gene_id,
                    gene_name=attrs.get("gene_name", gene_id),
                    contig=fields[0],
                    start=start,
                    end=end,
                    strand=fields[6],
                    species=species_from(fields[0], gene_id),
                )
            else:
                exon_sets[gene_id].add((start, end))
    if not genes:
        raise ValueError(f"no gene features in {gtf}")
    bodies = IntervalBins()
    exons = IntervalBins()
    for gene in genes.values():
        bodies.add(gene.contig, gene.start, gene.end, gene.gene_id)
    exon_count = 0
    for gene_id, intervals in exon_sets.items():
        gene = genes.get(gene_id)
        if gene is None:
            raise ValueError(f"exon references absent gene_id: {gene_id}")
        for start, end in sorted(intervals):
            exons.add(gene.contig, start, end, gene_id)
            exon_count += 1
    if exon_count == 0:
        raise ValueError(f"no exon features in {gtf}")
    return GeneModels(genes=genes, bodies=bodies, exons=exons, exon_intervals=exon_count)


def cigar_blocks(start: int, cigar: str) -> tuple[list[tuple[int, int]], bool]:
    if cigar == "*":
        return [], False
    blocks: list[tuple[int, int]] = []
    reference = start
    parsed = 0
    is_spliced = False
    for match in CIGAR_RE.finditer(cigar):
        if match.start() != parsed:
            raise ValueError(f"malformed CIGAR: {cigar}")
        parsed = match.end()
        length = int(match.group(1))
        operation = match.group(2)
        if length <= 0:
            raise ValueError(f"zero-length CIGAR operation: {cigar}")
        if operation in ("M", "=", "X"):
            blocks.append((reference, reference + length))
            reference += length
        elif operation in ("D", "N"):
            if operation == "N":
                is_spliced = True
            reference += length
        elif operation in ("I", "S", "H", "P"):
            continue
    if parsed != len(cigar):
        raise ValueError(f"malformed CIGAR: {cigar}")
    return blocks, is_spliced


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def reverse_unique_gene(
    gene_ids: Iterable[str], read_strand: str, genes: dict[str, Gene]
) -> tuple[str, Gene | None]:
    selected = [genes[gene_id] for gene_id in gene_ids if genes[gene_id].strand != read_strand]
    if not selected:
        return "no_feature", None
    if len(selected) != 1:
        return "ambiguous", None
    return "assigned", selected[0]


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_gzip_rows(path: Path, header: str, rows: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with gzip.open(temporary, "wt", encoding="utf-8", compresslevel=1) as handle:
            handle.write(header)
            for row in rows:
                handle.write(row)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--pair-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--manifest-index", type=int, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--expected-input-reads", type=int, required=True)
    parser.add_argument("--expected-primary-mapped", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for path in (args.bam, Path(str(args.bam) + ".bai"), args.gtf):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty input: {path}")
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools unavailable: {args.samtools}")
    if args.threads < 1 or args.expected_input_reads < 1 or args.expected_primary_mapped < 1:
        raise ValueError("invalid thread or expected-read value")

    models = load_gene_models(args.gtf)
    counters: Counter[str] = Counter()
    cell_names: dict[str, str] = {}
    cell_eligible: Counter[str] = Counter()
    mode_gene_reads: dict[str, Counter[tuple[str, str]]] = {mode: Counter() for mode in MODES}
    mode_gene_umis: dict[str, dict[tuple[str, str], set[str]]] = {
        mode: defaultdict(set) for mode in MODES
    }

    command = [samtools, "view", "-@", str(args.threads), str(args.bam)]
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        bufsize=1024 * 1024,
    )
    assert proc.stdout is not None
    try:
        for raw in proc.stdout:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"malformed SAM record from {args.bam}")
            flag = int(fields[1])
            counters["alignment_records"] += 1
            if flag & (0x100 | 0x800):
                counters["secondary_or_supplementary"] += 1
                continue
            counters["primary_records"] += 1
            if flag & 0x4:
                counters["primary_unmapped"] += 1
                continue
            counters["primary_mapped"] += 1
            if flag & (0x200 | 0x400):
                counters["qcfail_or_duplicate"] += 1
                continue
            tags = parse_tags(fields[11:])
            if "NH" not in tags:
                counters["missing_NH"] += 1
                continue
            if int(tags["NH"]) != 1:
                counters["nonunique_primary"] += 1
                continue
            counters["unique_primary"] += 1
            required = ("CB", "CN", "UB", "XD")
            if any(tag not in tags for tag in required):
                counters["missing_required_tags"] += 1
                continue
            if int(tags["XD"]) != 0:
                counters["nonzero_XD"] += 1
                continue
            cb = tags["CB"].split()[0]
            cn = tags["CN"].split()[0]
            ub = tags["UB"].split()[0]
            if not cb or not ub:
                counters["empty_CB_or_UB"] += 1
                continue
            counters["tag_eligible_unique_primary"] += 1
            cell_names.setdefault(cb, cn)
            cell_eligible[cb] += 1

            contig = fields[2]
            start = int(fields[3]) - 1
            blocks, is_spliced = cigar_blocks(start, fields[5])
            read_strand = "-" if flag & 0x10 else "+"
            body_status, body_gene = reverse_unique_gene(
                models.bodies.query(contig, blocks), read_strand, models.genes
            )
            exon_status, exon_gene = reverse_unique_gene(
                models.exons.query(contig, blocks), read_strand, models.genes
            )
            counters[f"genefull_{body_status}_reads"] += 1
            counters[f"exon_reverse_{exon_status}_reads"] += 1

            assignments: list[tuple[str, Gene]] = []
            if body_gene is not None:
                assignments.append(("genefull", body_gene))
                if is_spliced:
                    assignments.append(("junction_genefull", body_gene))
            if exon_gene is not None:
                assignments.append(("exon_reverse", exon_gene))
            if is_spliced:
                counters["eligible_spliced_reads"] += 1
            for mode, gene in assignments:
                counters[f"{mode}_assigned_reads"] += 1
                counters[f"{mode}_{gene.species}_assigned_reads"] += 1
                mode_gene_reads[mode][(cb, gene.gene_id)] += 1
                mode_gene_umis[mode][(cb, gene.gene_id)].add(ub)
    finally:
        proc.stdout.close()
        stderr = ""
        if proc.stderr is not None:
            stderr = proc.stderr.read()
            proc.stderr.close()
        return_code = proc.wait()
        if return_code != 0:
            raise RuntimeError(f"samtools view failed ({return_code}): {stderr.strip()}")

    if counters["primary_records"] != args.expected_input_reads:
        raise ValueError(
            f"primary records={counters['primary_records']}, expected={args.expected_input_reads}"
        )
    if counters["primary_mapped"] != args.expected_primary_mapped:
        raise ValueError(
            f"primary mapped={counters['primary_mapped']}, expected={args.expected_primary_mapped}"
        )
    if counters["missing_required_tags"] or counters["nonzero_XD"]:
        raise ValueError("required tag or strict-XD validation failed")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    counts_path = args.out_dir / f"{args.pair_id}.evidence_mode_counts.tsv.gz"
    cell_qc_path = args.out_dir / f"{args.pair_id}.evidence_mode_cell_qc.tsv.gz"
    summary_path = args.out_dir / f"{args.pair_id}.evidence_mode_summary.json"

    count_lines: list[str] = []
    mode_cell_umis: dict[str, Counter[str]] = {mode: Counter() for mode in MODES}
    mode_cell_genes: dict[str, Counter[str]] = {mode: Counter() for mode in MODES}
    mode_cell_species_reads: dict[str, Counter[tuple[str, str]]] = {
        mode: Counter() for mode in MODES
    }
    mode_cell_species_umis: dict[str, Counter[tuple[str, str]]] = {
        mode: Counter() for mode in MODES
    }
    mode_totals: dict[str, dict[str, int]] = {}
    for mode in MODES:
        total_umis = 0
        for (cb, gene_id), umi_set in sorted(mode_gene_umis[mode].items()):
            gene = models.genes[gene_id]
            reads = mode_gene_reads[mode][(cb, gene_id)]
            umis = len(umi_set)
            if umis > reads:
                raise ValueError(f"{mode}: UMI count exceeds reads: {cb} {gene_id}")
            total_umis += umis
            mode_cell_umis[mode][cb] += umis
            mode_cell_genes[mode][cb] += 1
            mode_cell_species_reads[mode][(cb, gene.species)] += reads
            mode_cell_species_umis[mode][(cb, gene.species)] += umis
            count_lines.append(
                f"{args.pair_id}\t{args.sample}\t{args.sample}|{cb}\t{cb}\t"
                f"{cell_names.get(cb, '')}\t{mode}\t{gene_id}\t{gene.gene_name}\t"
                f"{gene.species}\t{reads}\t{umis}\n"
            )
        mode_totals[mode] = {
            "assigned_reads": sum(mode_gene_reads[mode].values()),
            "exact_umis": total_umis,
            "cell_gene_entries": len(mode_gene_umis[mode]),
            "cells_with_umi": len(mode_cell_umis[mode]),
            "human_exact_umis": sum(
                value for (cb, species), value in mode_cell_species_umis[mode].items() if species == "human"
            ),
            "pig_exact_umis": sum(
                value for (cb, species), value in mode_cell_species_umis[mode].items() if species == "pig"
            ),
        }
    write_gzip_rows(
        counts_path,
        "pair_id\tsample\tcell_id\tCB\tCN\tmode\tgene_id\tgene_name\tspecies\tassigned_reads\texact_umis\n",
        count_lines,
    )

    cell_lines: list[str] = []
    for cb in sorted(cell_eligible):
        for mode in MODES:
            human_reads = mode_cell_species_reads[mode][(cb, "human")]
            pig_reads = mode_cell_species_reads[mode][(cb, "pig")]
            human_umis = mode_cell_species_umis[mode][(cb, "human")]
            pig_umis = mode_cell_species_umis[mode][(cb, "pig")]
            if human_umis > human_reads or pig_umis > pig_reads:
                raise ValueError(f"{mode}: species UMI count exceeds assigned reads: {cb}")
            cell_lines.append(
                f"{args.pair_id}\t{args.sample}\t{args.sample}|{cb}\t{cb}\t"
                f"{cell_names.get(cb, '')}\t{mode}\t{cell_eligible[cb]}\t"
                f"{human_reads + pig_reads}\t{mode_cell_umis[mode][cb]}\t"
                f"{mode_cell_genes[mode][cb]}\t{human_reads}\t{pig_reads}\t"
                f"{human_umis}\t{pig_umis}\n"
            )
    write_gzip_rows(
        cell_qc_path,
        "pair_id\tsample\tcell_id\tCB\tCN\tmode\ttag_eligible_unique_primary\t"
        "assigned_reads\texact_umis\tgenes_detected\thuman_assigned_reads\t"
        "pig_assigned_reads\thuman_exact_umis\tpig_exact_umis\n",
        cell_lines,
    )

    payload: dict[str, object] = {
        "pipeline_signature": PIPELINE_SIGNATURE,
        "status": "validated",
        "pair_id": args.pair_id,
        "sample": args.sample,
        "manifest_index": args.manifest_index,
        "shard": args.shard,
        "bam": str(args.bam),
        "gtf": str(args.gtf),
        "strand": "reverse",
        "modes": list(MODES),
        "gtf_genes": len(models.genes),
        "gtf_unique_exon_intervals": models.exon_intervals,
        "counters": dict(counters),
        "mode_totals": mode_totals,
        "counts": str(counts_path),
        "cell_qc": str(cell_qc_path),
    }
    atomic_json(summary_path, payload)
    print(
        f"[OK] {args.pair_id}: "
        + " ".join(f"{mode}={mode_totals[mode]['exact_umis']}UMI" for mode in MODES)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
