#!/usr/bin/env python3
"""Step 39: count exact GeneFull UMIs by CB for one RNA Probe BAM."""

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
MODES = ("forward", "reverse", "unstranded")


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
        seen: set[tuple[str, int, int]] = set()
        for block_start, block_end in blocks:
            if block_end <= block_start:
                continue
            first = block_start // self.bin_size
            last = (block_end - 1) // self.bin_size
            for index in range(first, last + 1):
                for gene in bins.get(index, ()):
                    token = (gene.gene_id, block_start, block_end)
                    if token in seen:
                        continue
                    seen.add(token)
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
                raise ValueError(f"{gtf}:{line_no}: gene lacks gene_id")
            if gene_id in seen_ids:
                raise ValueError(f"duplicate gene_id in GTF: {gene_id}")
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
        raise ValueError(f"no gene features found in {gtf}")
    return genes


def assign_mode(
    candidates: dict[str, Gene], read_strand: str, mode: str
) -> tuple[str, Gene | None]:
    if mode == "forward":
        selected = [
            gene for gene in candidates.values() if gene.strand == read_strand
        ]
    elif mode == "reverse":
        selected = [
            gene for gene in candidates.values() if gene.strand != read_strand
        ]
    elif mode == "unstranded":
        selected = list(candidates.values())
    else:
        raise ValueError(mode)
    if not selected:
        return "no_feature", None
    unique = {gene.gene_id: gene for gene in selected}
    if len(unique) != 1:
        return "ambiguous", None
    return "assigned", next(iter(unique.values()))


def cigar_blocks(start: int, cigar: str) -> list[tuple[int, int]]:
    """Return 0-based reference blocks containing aligned query bases."""
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
        else:
            raise ValueError(f"unsupported CIGAR operation: {operation}")
    if parsed != len(cigar):
        raise ValueError(f"malformed CIGAR: {cigar}")
    return blocks


def parse_sam_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument(
        "--max-primary-records",
        type=int,
        default=0,
        help="0 processes the complete BAM; nonzero is only for code tests.",
    )
    args = parser.parse_args()

    for path in (args.bam, Path(str(args.bam) + ".bai"), args.gtf):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing/empty input: {path}")
    if args.threads < 1 or args.max_primary_records < 0:
        parser.error("invalid thread/record limit")
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools unavailable: {args.samtools}")

    gene_bins = load_genes(args.gtf)
    counters: Counter[str] = Counter()
    mode_counts: dict[str, Counter[str]] = {
        mode: Counter() for mode in MODES
    }
    molecules: dict[str, set[tuple[str, str, str]]] = {
        mode: set() for mode in MODES
    }
    cell_gene_pairs: dict[str, set[tuple[str, str]]] = {
        mode: set() for mode in MODES
    }
    cell_assigned_reads: dict[str, Counter[str]] = {
        mode: Counter() for mode in MODES
    }
    cell_names: dict[str, str] = {}
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
    stopped_early = False
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
            if args.max_primary_records and (
                counters["primary_records"] > args.max_primary_records
            ):
                stopped_early = True
                break
            if flag & 0x4:
                counters["primary_unmapped"] += 1
                continue
            if flag & (0x200 | 0x400):
                counters["qcfail_or_duplicate"] += 1
                continue
            counters["primary_mapped"] += 1
            tags = parse_sam_tags(fields[11:])
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

            contig = fields[2]
            start = int(fields[3]) - 1
            blocks = cigar_blocks(start, fields[5])
            candidates = gene_bins.query(contig, blocks)
            read_strand = "-" if flag & 0x10 else "+"
            for mode in MODES:
                status, gene = assign_mode(candidates, read_strand, mode)
                mode_counts[mode][status] += 1
                if gene is None:
                    continue
                gene_lookup[gene.gene_id] = gene
                mode_counts[mode][f"{gene.species}_assigned_reads"] += 1
                cell_assigned_reads[mode][cb] += 1
                molecules[mode].add((cb, gene.gene_id, ub))
                cell_gene_pairs[mode].add((cb, gene.gene_id))
    finally:
        proc.stdout.close()
        if stopped_early:
            proc.terminate()
        stderr = ""
        if proc.stderr is not None:
            stderr = proc.stderr.read()
            proc.stderr.close()
        return_code = proc.wait()
        if not stopped_early and return_code != 0:
            raise RuntimeError(
                f"samtools view failed ({return_code}): {stderr.strip()}"
            )

    per_cell_path = args.out_dir / f"{args.pair_id}.per_cell.tsv.gz"
    per_gene_path = args.out_dir / f"{args.pair_id}.per_gene.tsv.gz"
    summary_path = args.out_dir / f"{args.pair_id}.summary.json"
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cell_umis: dict[str, Counter[str]] = {
        mode: Counter() for mode in MODES
    }
    cell_genes: dict[str, Counter[str]] = {
        mode: Counter() for mode in MODES
    }
    gene_umis: dict[str, Counter[str]] = {
        mode: Counter() for mode in MODES
    }
    for mode in MODES:
        for cb, gene_id, _umi in molecules[mode]:
            cell_umis[mode][cb] += 1
            gene_umis[mode][gene_id] += 1
        for cb, _gene_id in cell_gene_pairs[mode]:
            cell_genes[mode][cb] += 1

    tmp_cell = per_cell_path.with_name(
        per_cell_path.name + f".tmp.{os.getpid()}"
    )
    with gzip.open(tmp_cell, "wt", encoding="utf-8", compresslevel=1) as out:
        out.write(
            "pair_id\tsample\tmode\tcell_id\tCB\tCN\tassigned_reads"
            "\texact_umis\tgenes_detected\n"
        )
        for mode in MODES:
            for cb in sorted(cell_assigned_reads[mode]):
                out.write(
                    f"{args.pair_id}\t{args.sample}\t{mode}\t"
                    f"{args.sample}|{cb}\t{cb}\t{cell_names.get(cb, '')}\t"
                    f"{cell_assigned_reads[mode][cb]}\t"
                    f"{cell_umis[mode][cb]}\t{cell_genes[mode][cb]}\n"
                )
    os.replace(tmp_cell, per_cell_path)

    tmp_gene = per_gene_path.with_name(
        per_gene_path.name + f".tmp.{os.getpid()}"
    )
    with gzip.open(tmp_gene, "wt", encoding="utf-8", compresslevel=1) as out:
        out.write(
            "pair_id\tsample\tmode\tgene_id\tgene_name\tspecies"
            "\texact_umis\n"
        )
        for mode in MODES:
            for gene_id, count in sorted(gene_umis[mode].items()):
                gene = gene_lookup[gene_id]
                out.write(
                    f"{args.pair_id}\t{args.sample}\t{mode}\t{gene_id}\t"
                    f"{gene.gene_name}\t{gene.species}\t{count}\n"
                )
    os.replace(tmp_gene, per_gene_path)

    mode_payload: dict[str, dict[str, int | float]] = {}
    eligible = counters["tag_eligible_unique_primary"]
    for mode in MODES:
        counts = mode_counts[mode]
        human_umis = sum(
            count
            for gene_id, count in gene_umis[mode].items()
            if gene_lookup[gene_id].species == "human"
        )
        pig_umis = sum(
            count
            for gene_id, count in gene_umis[mode].items()
            if gene_lookup[gene_id].species == "pig"
        )
        mode_payload[mode] = {
            "assigned_reads": counts["assigned"],
            "assigned_fraction_of_tag_eligible": (
                counts["assigned"] / eligible if eligible else 0.0
            ),
            "no_feature_reads": counts["no_feature"],
            "ambiguous_reads": counts["ambiguous"],
            "human_assigned_reads": counts["human_assigned_reads"],
            "pig_assigned_reads": counts["pig_assigned_reads"],
            "exact_umis": len(molecules[mode]),
            "human_exact_umis": human_umis,
            "pig_exact_umis": pig_umis,
            "cells_with_assigned_reads": len(cell_assigned_reads[mode]),
            "detected_cell_gene_pairs": len(cell_gene_pairs[mode]),
        }
    payload: dict[str, object] = {
        "pipeline_signature": "snmar_rna_genefull_exact_umi_probe_v12",
        "status": "completed",
        "pair_id": args.pair_id,
        "sample": args.sample,
        "bam": str(args.bam),
        "gtf": str(args.gtf),
        "gtf_gene_features": gene_bins.genes,
        "counting_feature": "GeneFull_gene_body",
        "umi_deduplication": "exact_within_sample_CB_gene",
        "counters": dict(counters),
        "modes": mode_payload,
        "per_cell": str(per_cell_path),
        "per_gene": str(per_gene_path),
    }
    atomic_json(summary_path, payload)
    print(
        f"[OK] {args.pair_id}: eligible={eligible} "
        f"forward={mode_payload['forward']['assigned_reads']} "
        f"reverse={mode_payload['reverse']['assigned_reads']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
