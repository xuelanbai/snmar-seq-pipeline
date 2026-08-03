#!/usr/bin/env python3
"""Validate one coordinate-sorted DNA BWA BAM."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

REQUIRED_TAGS = {"CR", "CB", "CN", "XD", "ZE"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--expected-pairs", required=True, type=int)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--tag-records", type=int, default=1000)
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def run_text(command: list[str], *, check: bool = True) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if check and result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stderr[-4000:]}"
        )
    return result.stdout


def parse_flagstat(text: str) -> dict[str, int]:
    metrics: dict[str, int] = {}
    for line in text.splitlines():
        match = re.match(r"^(\d+) \+ (\d+) (.+)$", line.strip())
        if not match:
            continue
        passed = int(match.group(1))
        failed = int(match.group(2))
        label = match.group(3)
        value = passed + failed
        if label.startswith("in total "):
            metrics["total_records"] = value
        elif label == "primary":
            metrics["primary_records"] = value
        elif label == "secondary":
            metrics["secondary_records"] = value
        elif label == "supplementary":
            metrics["supplementary_records"] = value
        elif label.startswith("primary mapped "):
            metrics["primary_mapped"] = value
        elif label.startswith("mapped ") and "primary_mapped" not in metrics:
            metrics["mapped_all"] = value
        elif label.startswith("primary duplicates"):
            metrics["primary_duplicates"] = value
        elif label.startswith("duplicates") and "primary_duplicates" not in metrics:
            metrics["duplicates_all"] = value
        elif label.startswith("paired in sequencing"):
            metrics["paired_records"] = value
        elif label.startswith("properly paired"):
            metrics["properly_paired"] = value
        elif label.startswith("singletons"):
            metrics["singletons"] = value
    required = {"total_records", "paired_records", "properly_paired"}
    missing = required.difference(metrics)
    if missing:
        raise ValueError(f"could not parse flagstat metrics: {sorted(missing)}")
    if "primary_records" not in metrics:
        metrics["primary_records"] = (
            metrics["total_records"]
            - metrics.get("secondary_records", 0)
            - metrics.get("supplementary_records", 0)
        )
    if "primary_mapped" not in metrics:
        metrics["primary_mapped"] = metrics.get("mapped_all", -1)
    if "primary_duplicates" not in metrics:
        metrics["primary_duplicates"] = metrics.get("duplicates_all", -1)
    return metrics


def parse_header_rg(header: str, sample: str, run_id: str) -> None:
    sq_count = 0
    human = 0
    pig = 0
    matching_rg = False
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            sq_count += 1
            fields = dict(
                field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
            )
            name = fields.get("SN", "")
            human += int(name.startswith("human_"))
            pig += int(name.startswith("pig_"))
        elif line.startswith("@RG\t"):
            fields = dict(
                field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
            )
            if fields.get("SM") == sample and fields.get("ID") == run_id:
                matching_rg = True
    if sq_count == 0 or human == 0 or pig == 0:
        raise ValueError(
            f"BAM header lacks prefixed joint reference: SQ={sq_count}, "
            f"human={human}, pig={pig}"
        )
    if not matching_rg:
        raise ValueError(f"BAM header lacks expected RG ID={run_id!r}, SM={sample!r}")


def contigs_and_species(
    samtools: str, bam: Path, threads: int
) -> tuple[list[str], list[str], list[str], int, int, int]:
    text = run_text([samtools, "idxstats", str(bam)])
    human_contigs: list[str] = []
    pig_contigs: list[str] = []
    other_contigs: list[str] = []
    indexed_human = 0
    indexed_pig = 0
    indexed_other = 0
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 4:
            continue
        name = fields[0]
        mapped = int(fields[2])
        if name == "*":
            continue
        if name.startswith("human_"):
            human_contigs.append(name)
            indexed_human += mapped
        elif name.startswith("pig_"):
            pig_contigs.append(name)
            indexed_pig += mapped
        else:
            other_contigs.append(name)
            indexed_other += mapped
    if not human_contigs or not pig_contigs:
        raise ValueError("BAM index lacks human_ or pig_ contigs")
    if other_contigs:
        raise ValueError(f"unexpected unprefixed BAM contigs: {other_contigs[:5]}")
    return (
        human_contigs,
        pig_contigs,
        other_contigs,
        indexed_human,
        indexed_pig,
        indexed_other,
    )


def primary_species_count(
    samtools: str, bam: Path, threads: int, contigs: list[str]
) -> int:
    command = [
        samtools,
        "view",
        "-@",
        str(threads),
        "-c",
        "-F",
        "2308",
        str(bam),
        *contigs,
    ]
    return int(run_text(command).strip())


def sample_tags(
    samtools: str, bam: Path, threads: int, wanted: int
) -> dict[str, object]:
    process = subprocess.Popen(
        [samtools, "view", "-@", str(threads), "-F", "2308", str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    observed = 0
    missing: dict[str, int] = {key: 0 for key in REQUIRED_TAGS}
    xd_values: set[str] = set()
    cb_lengths: set[int] = set()
    cr_lengths: set[int] = set()
    ze_values: set[str] = set()
    try:
        assert process.stdout is not None
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")
            tags: dict[str, tuple[str, str]] = {}
            for field in fields[11:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = (parts[1], parts[2])
            for key in REQUIRED_TAGS:
                if key not in tags:
                    missing[key] += 1
            if "XD" in tags:
                xd_values.add(tags["XD"][1])
            if "CB" in tags:
                cb_lengths.add(len(tags["CB"][1]))
            if "CR" in tags:
                cr_lengths.add(len(tags["CR"][1]))
            if "ZE" in tags:
                ze_values.add(tags["ZE"][1])
            observed += 1
            if observed >= wanted:
                break
    finally:
        if process.stdout is not None:
            process.stdout.close()
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
    if observed == 0:
        raise ValueError("no primary records available for BAM tag validation")
    failed = {key: value for key, value in missing.items() if value}
    if failed:
        raise ValueError(f"required BAM tags are missing in sampled records: {failed}")
    if xd_values != {"0"}:
        raise ValueError(f"strict DNA BAM contains nonzero XD values: {sorted(xd_values)}")
    return {
        "sampled_primary_records": observed,
        "required_tags": sorted(REQUIRED_TAGS),
        "missing_tag_counts": missing,
        "xd_values": sorted(xd_values),
        "cr_lengths": sorted(cr_lengths),
        "cb_lengths": sorted(cb_lengths),
        "ze_examples": sorted(ze_values)[:10],
    }


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent, text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_name, path)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    args = parse_args()
    if args.expected_pairs < 1 or args.threads < 1 or args.tag_records < 1:
        raise ValueError("counts and threads must be positive")
    if not args.bam.is_file() or args.bam.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty BAM: {args.bam}")
    bai = Path(f"{args.bam}.bai")
    if not bai.is_file() or bai.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty BAM index: {bai}")

    run_text([args.samtools, "quickcheck", "-v", str(args.bam)])
    header = run_text([args.samtools, "view", "-H", str(args.bam)])
    parse_header_rg(header, args.sample, args.run_id)
    flagstat_text = run_text(
        [args.samtools, "flagstat", "-@", str(args.threads), str(args.bam)]
    )
    flags = parse_flagstat(flagstat_text)
    expected_primary = args.expected_pairs * 2
    if flags["primary_records"] != expected_primary:
        raise ValueError(
            f"primary-record count mismatch: expected {expected_primary:,}, "
            f"observed {flags['primary_records']:,}"
        )
    if flags["paired_records"] < expected_primary:
        raise ValueError(
            f"too few records carry the paired flag: expected at least "
            f"{expected_primary:,}, observed {flags['paired_records']:,}"
        )
    if flags["primary_duplicates"] != 0:
        raise ValueError(
            "global duplicate flags are not allowed in the multi-cell mapping BAM: "
            f"observed {flags['primary_duplicates']:,} primary duplicates"
        )

    (
        human_contigs,
        pig_contigs,
        _,
        indexed_human,
        indexed_pig,
        indexed_other,
    ) = contigs_and_species(args.samtools, args.bam, args.threads)
    human_primary = primary_species_count(
        args.samtools, args.bam, args.threads, human_contigs
    )
    pig_primary = primary_species_count(args.samtools, args.bam, args.threads, pig_contigs)
    species_primary = human_primary + pig_primary
    mapped_primary = flags["primary_mapped"]
    if mapped_primary >= 0 and species_primary != mapped_primary:
        raise ValueError(
            "primary species counts do not sum to primary mapped records: "
            f"human={human_primary:,}, pig={pig_primary:,}, mapped={mapped_primary:,}"
        )
    if mapped_primary < 0:
        mapped_primary = species_primary
        flags["primary_mapped"] = mapped_primary

    tags = sample_tags(args.samtools, args.bam, args.threads, args.tag_records)
    denominator = expected_primary
    payload: dict[str, object] = {
        "status": "validated",
        "bam": str(args.bam),
        "bam_bytes": args.bam.stat().st_size,
        "sample": args.sample,
        "run_id": args.run_id,
        "expected_read_pairs": args.expected_pairs,
        "expected_primary_records": expected_primary,
        "global_duplicate_marking": "disabled_multi_cell_CB_aware_required",
        **flags,
        "primary_mapping_percent": (
            100.0 * mapped_primary / denominator if mapped_primary >= 0 else None
        ),
        "properly_paired_percent": 100.0 * flags["properly_paired"] / denominator,
        "primary_duplicate_percent": (
            100.0 * flags["primary_duplicates"] / denominator
            if flags["primary_duplicates"] >= 0
            else None
        ),
        "human_primary_mapped": human_primary,
        "pig_primary_mapped": pig_primary,
        "human_percent_of_primary_mapped": (
            100.0 * human_primary / mapped_primary if mapped_primary > 0 else 0.0
        ),
        "pig_percent_of_primary_mapped": (
            100.0 * pig_primary / mapped_primary if mapped_primary > 0 else 0.0
        ),
        "indexed_human_alignment_records": indexed_human,
        "indexed_pig_alignment_records": indexed_pig,
        "indexed_other_alignment_records": indexed_other,
        "tag_validation": tags,
    }
    atomic_json(args.out, payload)
    print(
        f"[OK] {args.run_id}: primary={expected_primary:,}, "
        f"mapped={mapped_primary:,} ({payload['primary_mapping_percent']:.2f}%), "
        f"proper={flags['properly_paired']:,} "
        f"({payload['properly_paired_percent']:.2f}%)"
    )
    print(
        f"[OK] species primary mapped: human={human_primary:,}, pig={pig_primary:,}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
