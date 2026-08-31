#!/usr/bin/env python3
"""Step 46: full reverse-strand GeneFull exact-UMI count for one BAM."""

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


class GeneBins:
    def __init__(self, bin_size: int = 100_000) -> None:
        self.bin_size = bin_size
        self.by_contig: dict[str, dict[int, list[Gene]]] = defaultdict(
            lambda: defaultdict(list)
        )
        self.genes = 0

    def add(self, gene: Gene) -> None:
        first = gene.start // self.bin_size
        last = max(gene.start, gene.end - 1) // self.bin_size
        for index in range(first, last + 1):
            self.by_contig[gene.contig][index].append(gene)
        self.genes += 1

    def query(
        self, contig: str, blocks: list[tuple[int, int]]
    ) -> dict[str, Gene]:
        bins = self.by_contig.get(contig)
        if bins is None:
            return {}
        hits: dict[str, Gene] = {}
        for block_start, block_end in blocks:
            if block_end <= block_start:
                continue
            first = block_start // self.bin_size
            last = (block_end - 1) // self.bin_size
            seen_gene_ids: set[str] = set()
            for index in range(first, last + 1):
                for gene in bins.get(index, ()):
                    if gene.gene_id in seen_gene_ids:
                        continue
                    seen_gene_ids.add(gene.gene_id)
                    if gene.start < block_end and block_start < gene.end:
                        hits[gene.gene_id] = gene
        return hits


def parse_attributes(text: str) -> dict[str, str]:
    return {key: value for key, value in ATTRIBUTE_RE.findall(text)}


def load_genes(gtf: Path) -> GeneBins:
    genes = GeneBins()
    seen_ids: set[str] = set()
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gtf}:{line_no}: expected 9 GTF fields")
            if fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                raise ValueError(f"{gtf}:{line_no}: missing gene_id")
            if gene_id in seen_ids:
                raise ValueError(f"duplicate gene_id: {gene_id}")
            seen_ids.add(gene_id)
            contig = fields[0]
            if contig.startswith("human_") or gene_id.startswith("HUMAN_"):
                species = "human"
            elif contig.startswith("pig_") or gene_id.startswith("PIG_"):
                species = "pig"
            else:
                species = "other"
            start = int(fields[3]) - 1
            end = int(fields[4])
            if start < 0 or end <= start:
                raise ValueError(f"{gtf}:{line_no}: invalid gene interval")
            genes.add(
                Gene(
                    gene_id=gene_id,
                    gene_name=attrs.get("gene_name", gene_id),
                    contig=contig,
                    start=start,
                    end=end,
                    strand=fields[6],
                    species=species,
                )
            )
    if genes.genes == 0:
        raise ValueError(f"no gene features in {gtf}")
    return genes


def cigar_blocks(start: int, cigar: str) -> list[tuple[int, int]]:
    if cigar == "*":
        return []
    blocks: list[tuple[int, int]] = []
    reference = start
    parsed = 0
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
            reference += length
        elif operation in ("I", "S", "H", "P"):
            continue
    if parsed != len(cigar):
        raise ValueError(f"malformed CIGAR: {cigar}")
    return blocks


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def reverse_gene(
    candidates: dict[str, Gene], read_strand: str
) -> tuple[str, Gene | None]:
    selected = {
        gene.gene_id: gene
        for gene in candidates.values()
        if gene.strand != read_strand
    }
    if not selected:
        return "no_feature", None
    if len(selected) != 1:
        return "ambiguous", None
    return "assigned", next(iter(selected.values()))


def atomic_json(path: Path, payload: dict[str, object]) -> None:
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
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pair-id", required=True)
    parser.add_argument("--expected-input-reads", type=int, required=True)
    parser.add_argument("--expected-primary-mapped", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    for path in (args.bam, Path(str(args.bam) + ".bai"), args.gtf):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty input: {path}")
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools unavailable: {args.samtools}")
    if (
        args.threads < 1
        or args.expected_input_reads < 1
        or args.expected_primary_mapped < 1
    ):
        parser.error("invalid thread or expected-read value")

    genes = load_genes(args.gtf)
    counters: Counter[str] = Counter()
    cell_names: dict[str, str] = {}
    cell_unique_reads: Counter[str] = Counter()
    cell_assigned_reads: Counter[str] = Counter()
    cell_gene_reads: Counter[tuple[str, str]] = Counter()
    cell_species_reads: Counter[tuple[str, str]] = Counter()
    cell_gene_umis: dict[tuple[str, str], set[str]] = defaultdict(set)
    gene_lookup: dict[str, Gene] = {}

    command = [
        samtools,
        "view",
        "-@",
        str(args.threads),
        str(args.bam),
    ]
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
            cell_unique_reads[cb] += 1

            contig = fields[2]
            start = int(fields[3]) - 1
            blocks = cigar_blocks(start, fields[5])
            candidates = genes.query(contig, blocks)
            read_strand = "-" if flag & 0x10 else "+"
            status, gene = reverse_gene(candidates, read_strand)
            counters[f"gene_{status}"] += 1
            if gene is None:
                continue
            gene_lookup[gene.gene_id] = gene
            counters[f"{gene.species}_assigned_reads"] += 1
            cell_assigned_reads[cb] += 1
            cell_gene_reads[(cb, gene.gene_id)] += 1
            cell_species_reads[(cb, gene.species)] += 1
            cell_gene_umis[(cb, gene.gene_id)].add(ub)
    finally:
        proc.stdout.close()
        stderr = ""
        if proc.stderr is not None:
            stderr = proc.stderr.read()
            proc.stderr.close()
        return_code = proc.wait()
        if return_code != 0:
            raise RuntimeError(
                f"samtools view failed ({return_code}): {stderr.strip()}"
            )

    if counters["primary_records"] != args.expected_input_reads:
        raise ValueError(
            f"primary records={counters['primary_records']}, "
            f"expected STAR input={args.expected_input_reads}"
        )
    if counters["primary_mapped"] != args.expected_primary_mapped:
        raise ValueError(
            f"primary mapped={counters['primary_mapped']}, "
            f"expected={args.expected_primary_mapped}"
        )
    if counters["missing_required_tags"] or counters["nonzero_XD"]:
        raise ValueError("required tag or strict-XD validation failed")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    counts_path = (
        args.out_dir / f"{args.pair_id}.cell_gene_counts.tsv.gz"
    )
    cell_qc_path = args.out_dir / f"{args.pair_id}.cell_qc.tsv.gz"
    summary_path = args.out_dir / f"{args.pair_id}.summary.json"

    cell_umis: Counter[str] = Counter()
    cell_genes: Counter[str] = Counter()
    cell_species_umis: Counter[tuple[str, str]] = Counter()
    tmp_counts = counts_path.with_name(
        counts_path.name + f".tmp.{os.getpid()}"
    )
    with gzip.open(
        tmp_counts, "wt", encoding="utf-8", compresslevel=1
    ) as out:
        out.write(
            "pair_id\tsample\tcell_id\tCB\tCN\tgene_id\tgene_name"
            "\tspecies\tassigned_reads\texact_umis\n"
        )
        for (cb, gene_id), umi_set in sorted(cell_gene_umis.items()):
            gene = gene_lookup[gene_id]
            umi_count = len(umi_set)
            cell_umis[cb] += umi_count
            cell_genes[cb] += 1
            cell_species_umis[(cb, gene.species)] += umi_count
            out.write(
                f"{args.pair_id}\t{args.sample}\t{args.sample}|{cb}\t"
                f"{cb}\t{cell_names.get(cb, '')}\t{gene_id}\t"
                f"{gene.gene_name}\t{gene.species}\t"
                f"{cell_gene_reads[(cb, gene_id)]}\t{umi_count}\n"
            )
    os.replace(tmp_counts, counts_path)

    tmp_qc = cell_qc_path.with_name(
        cell_qc_path.name + f".tmp.{os.getpid()}"
    )
    with gzip.open(tmp_qc, "wt", encoding="utf-8", compresslevel=1) as out:
        out.write(
            "pair_id\tsample\tcell_id\tCB\tCN"
            "\ttag_eligible_unique_primary\tgene_assigned_reads"
            "\tgene_assignment_fraction\texact_umis\tgenes_detected"
            "\thuman_assigned_reads\tpig_assigned_reads"
            "\thuman_exact_umis\tpig_exact_umis\n"
        )
        for cb in sorted(cell_unique_reads):
            eligible = cell_unique_reads[cb]
            assigned = cell_assigned_reads[cb]
            out.write(
                f"{args.pair_id}\t{args.sample}\t{args.sample}|{cb}\t"
                f"{cb}\t{cell_names.get(cb, '')}\t{eligible}\t{assigned}\t"
                f"{assigned / eligible if eligible else 0.0:.8f}\t"
                f"{cell_umis[cb]}\t{cell_genes[cb]}\t"
                f"{cell_species_reads[(cb, 'human')]}\t"
                f"{cell_species_reads[(cb, 'pig')]}\t"
                f"{cell_species_umis[(cb, 'human')]}\t"
                f"{cell_species_umis[(cb, 'pig')]}\n"
            )
    os.replace(tmp_qc, cell_qc_path)

    total_exact_umis = sum(len(values) for values in cell_gene_umis.values())
    payload: dict[str, object] = {
        "pipeline_signature": "snmar_rna_reverse_genefull_full_v13",
        "status": "validated",
        "pair_id": args.pair_id,
        "sample": args.sample,
        "bam": str(args.bam),
        "gtf": str(args.gtf),
        "strand": "reverse",
        "counting_feature": "GeneFull_gene_body",
        "umi_deduplication": "exact_within_sample_CB_gene",
        "gtf_gene_features": genes.genes,
        "counters": dict(counters),
        "exact_umis": total_exact_umis,
        "cells_tag_eligible": len(cell_unique_reads),
        "cells_with_exact_umi": len(cell_umis),
        "detected_cell_gene_pairs": len(cell_gene_umis),
        "cell_gene_counts": str(counts_path),
        "cell_qc": str(cell_qc_path),
    }
    atomic_json(summary_path, payload)
    print(
        f"[OK] {args.pair_id}: assigned={counters['gene_assigned']} "
        f"exact_umis={total_exact_umis}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
