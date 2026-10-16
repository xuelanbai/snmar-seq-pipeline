#!/usr/bin/env python3
"""Step 98: audit one RNA BAM by species, uniqueness, and splice status.

This is a read-only diagnostic.  It does not alter BAMs or matrices.  The
purpose is to distinguish an RNA-specific species-label problem from residual
DNA-like molecules that entered the strict RNA output.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Iterable


PIPELINE_SIGNATURE = "snmar_rna_species_alignment_audit_20260806"
SPECIES = ("human", "pig", "other")


def trailing_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if not match:
        raise ValueError(f"sample lacks trailing well number: {value}")
    return int(match.group(1))


def species_for_contig(contig: str) -> str:
    if contig.startswith("human_"):
        return "human"
    if contig.startswith("pig_"):
        return "pig"
    return "other"


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def audit_sam_lines(lines: Iterable[str]) -> Counter[str]:
    """Return counters matching the filters used by RNA GeneFull counting."""
    counts: Counter[str] = Counter()
    for line_no, raw in enumerate(lines, 1):
        if not raw.strip() or raw.startswith("@"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise ValueError(f"malformed SAM record at line {line_no}")
        flag = int(fields[1])
        if flag & (0x100 | 0x800):
            counts["secondary_or_supplementary"] += 1
            continue
        counts["primary_records"] += 1
        if flag & 0x4 or fields[2] == "*":
            counts["primary_unmapped"] += 1
            continue

        species = species_for_contig(fields[2])
        counts["primary_mapped"] += 1
        counts[f"primary_{species}"] += 1
        if flag & (0x200 | 0x400):
            counts["qcfail_or_duplicate"] += 1
            continue

        tags = parse_tags(fields[11:])
        if "NH" not in tags:
            counts["missing_NH"] += 1
            continue
        nh = int(tags["NH"])
        if nh != 1:
            counts["nonunique_primary"] += 1
            counts[f"nonunique_{species}"] += 1
            continue
        counts["unique_primary"] += 1
        counts[f"unique_{species}"] += 1

        required = ("CB", "CN", "UB", "XD")
        if any(tag not in tags for tag in required):
            counts["missing_required_tags"] += 1
            continue
        if int(tags["XD"]) != 0:
            counts["nonzero_XD"] += 1
            continue
        counts["tag_eligible_unique_primary"] += 1
        counts[f"eligible_{species}"] += 1
        splice_class = "spliced" if "N" in fields[5] else "unspliced"
        counts[f"eligible_{splice_class}_{species}"] += 1
    return counts


def ratio(numerator: int, denominator: int) -> float | str:
    return numerator / denominator if denominator else "NA"


def atomic_tsv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(row),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerow(row)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--pair-id", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--manifest-index", type=int, required=True)
    parser.add_argument("--shard", type=int, required=True)
    parser.add_argument("--expected-primary-mapped", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--human-only-wells", type=int, default=48)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.manifest_index < 1
        or args.shard < 1
        or args.expected_primary_mapped < 1
        or args.threads < 1
        or args.human_only_wells < 1
    ):
        raise ValueError("invalid positive integer argument")
    if not args.bam.is_file() or args.bam.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty BAM: {args.bam}")
    samtools = shutil.which(args.samtools)
    if samtools is None:
        raise FileNotFoundError(f"samtools unavailable: {args.samtools}")

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
        counts = audit_sam_lines(proc.stdout)
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
    if counts["primary_mapped"] != args.expected_primary_mapped:
        raise ValueError(
            f"primary mapped={counts['primary_mapped']}, "
            f"expected={args.expected_primary_mapped}"
        )
    if counts["missing_NH"]:
        raise ValueError(f"BAM records lack NH tags: {counts['missing_NH']}")
    if counts["missing_required_tags"] or counts["nonzero_XD"]:
        raise ValueError("strict RNA tag validation failed")
    if sum(counts[f"primary_{name}"] for name in SPECIES) != counts["primary_mapped"]:
        raise AssertionError("primary species counts do not reconcile")
    if sum(counts[f"unique_{name}"] for name in SPECIES) != counts["unique_primary"]:
        raise AssertionError("unique species counts do not reconcile")
    if (
        sum(counts[f"eligible_{name}"] for name in SPECIES)
        != counts["tag_eligible_unique_primary"]
    ):
        raise AssertionError("eligible species counts do not reconcile")

    well_index = trailing_number(args.sample)
    row: dict[str, object] = {
        "status": "validated",
        "pipeline_signature": PIPELINE_SIGNATURE,
        "manifest_index": args.manifest_index,
        "shard": args.shard,
        "pair_id": args.pair_id,
        "sample": args.sample,
        "well_index": well_index,
        "expected_well_set": (
            "GM12878" if well_index <= args.human_only_wells else "cellMix"
        ),
        "bam": args.bam,
        "primary_records": counts["primary_records"],
        "primary_mapped": counts["primary_mapped"],
        "primary_human": counts["primary_human"],
        "primary_pig": counts["primary_pig"],
        "primary_other": counts["primary_other"],
        "primary_human_fraction": ratio(
            counts["primary_human"], counts["primary_human"] + counts["primary_pig"]
        ),
        "unique_primary": counts["unique_primary"],
        "unique_human": counts["unique_human"],
        "unique_pig": counts["unique_pig"],
        "unique_other": counts["unique_other"],
        "unique_human_fraction": ratio(
            counts["unique_human"], counts["unique_human"] + counts["unique_pig"]
        ),
        "nonunique_primary": counts["nonunique_primary"],
        "tag_eligible_unique_primary": counts["tag_eligible_unique_primary"],
        "eligible_human": counts["eligible_human"],
        "eligible_pig": counts["eligible_pig"],
        "eligible_other": counts["eligible_other"],
        "eligible_spliced_human": counts["eligible_spliced_human"],
        "eligible_unspliced_human": counts["eligible_unspliced_human"],
        "eligible_spliced_pig": counts["eligible_spliced_pig"],
        "eligible_unspliced_pig": counts["eligible_unspliced_pig"],
        "eligible_human_spliced_fraction": ratio(
            counts["eligible_spliced_human"], counts["eligible_human"]
        ),
        "eligible_pig_spliced_fraction": ratio(
            counts["eligible_spliced_pig"], counts["eligible_pig"]
        ),
        "qcfail_or_duplicate": counts["qcfail_or_duplicate"],
        "missing_NH": counts["missing_NH"],
        "missing_required_tags": counts["missing_required_tags"],
        "nonzero_XD": counts["nonzero_XD"],
    }
    atomic_tsv(args.output, row)
    print(
        f"[OK] RNA species alignment audit: {args.sample} "
        f"unique_human={counts['unique_human']} "
        f"unique_pig={counts['unique_pig']} output={args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
