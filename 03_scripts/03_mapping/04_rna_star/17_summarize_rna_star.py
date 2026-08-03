#!/usr/bin/env python3
"""Step 17: validate and summarize RNA STAR probe/full mapping outputs."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


PER_FILE_FIELDS = (
    "pair_id",
    "sample",
    "expected_R1_reads",
    "star_input_reads",
    "unique_mapped_reads",
    "unique_mapped_percent",
    "multi_mapped_reads",
    "multi_mapped_percent",
    "primary_mapped_reads",
    "primary_mapping_percent",
    "primary_human_reads",
    "primary_pig_reads",
    "primary_other_contig_reads",
    "human_percent_of_primary_mapped",
    "pig_percent_of_primary_mapped",
    "unmapped_too_many_mismatches_percent",
    "unmapped_too_short_percent",
    "unmapped_other_percent",
    "checked_primary_tag_records",
    "bam_alignment_records",
    "bam_path",
)


def read_two_column(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not fields[0]:
                raise ValueError(f"{path}:{line_no}: expected key<TAB>value")
            key, value = fields
            if key == "metric" and value == "value":
                continue
            if key in result:
                raise ValueError(f"{path}:{line_no}: duplicate key {key}")
            result[key] = value
    return result


def parse_star_log(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = line.split("|", 1)
            result[key.strip()] = value.strip()
    return result


def as_int(mapping: dict[str, str], key: str, path: Path) -> int:
    try:
        return int(mapping[key].replace(",", ""))
    except (KeyError, ValueError) as exc:
        raise ValueError(f"{path}: invalid/missing integer {key!r}") from exc


def as_percent(mapping: dict[str, str], key: str, path: Path) -> float:
    try:
        return float(mapping[key].rstrip("%"))
    except (KeyError, ValueError) as exc:
        raise ValueError(f"{path}: invalid/missing percent {key!r}") from exc


def atomic_tsv(
    path: Path,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
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
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def atomic_summary(path: Path, rows: list[tuple[str, object, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with tmp.open("wt", encoding="utf-8", newline="\n") as handle:
            handle.write("metric\tvalue\tunit\n")
            for metric, value, unit in rows:
                handle.write(f"{metric}\t{value}\t{unit}\n")
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("probe", "full"), required=True)
    parser.add_argument("--input-table", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--done-root", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--per-file", type=Path, required=True)
    parser.add_argument("--expected-files", type=int, required=True)
    args = parser.parse_args()

    if args.expected_files < 1:
        parser.error("--expected-files must be positive")
    with args.input_table.open(encoding="utf-8", newline="") as handle:
        input_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(input_rows) != args.expected_files:
        raise ValueError(
            f"{args.input_table}: rows={len(input_rows)}, "
            f"expected={args.expected_files}"
        )

    rows: list[dict[str, object]] = []
    totals = {
        "input": 0,
        "unique": 0,
        "multi": 0,
        "primary_mapped": 0,
        "human": 0,
        "pig": 0,
        "other_contig": 0,
        "bam_records": 0,
        "checked_tags": 0,
    }
    weighted_unmapped = {
        "too_many_mismatches": 0.0,
        "too_short": 0.0,
        "other": 0.0,
    }
    seen: set[str] = set()

    for source in input_rows:
        pair_id = source["pair_id"]
        sample = source["sample"]
        expected = int(source["expected_reads"])
        if pair_id in seen:
            raise ValueError(f"duplicate pair_id: {pair_id}")
        seen.add(pair_id)

        out_dir = args.output_root / sample / pair_id
        bam = out_dir / f"{pair_id}.rna_joint.sorted.bam"
        bai = Path(str(bam) + ".bai")
        qc_path = out_dir / f"{pair_id}.bam_qc.tsv"
        log_path = out_dir / "star.Log.final.out"
        gene_counts = out_dir / "star.ReadsPerGene.out.tab"
        done_path = args.done_root / f"{pair_id}.done"
        failed_path = args.done_root / f"{pair_id}.failed"
        if failed_path.exists():
            raise RuntimeError(f"failed marker exists: {failed_path}")
        for path in (bam, bai, qc_path, log_path, gene_counts, done_path):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"missing/empty mapping output: {path}")

        done = read_two_column(done_path)
        expected_done = {
            "status": "complete",
            "pipeline_signature": "snmar_rna_star_joint_v10",
            "mode": args.mode,
            "pair_id": pair_id,
            "sample": sample,
            "clean_R1_path": source["clean_R1"],
            "expected_R1_reads": str(expected),
            "bam_path": str(bam),
            "bam_bytes": str(bam.stat().st_size),
        }
        for key, value in expected_done.items():
            if done.get(key) != value:
                raise ValueError(
                    f"{done_path}: {key}={done.get(key)!r}, "
                    f"expected={value!r}"
                )

        qc = read_two_column(qc_path)
        if qc.get("pipeline_signature") != "snmar_rna_bam_validate_v10":
            raise ValueError(f"{qc_path}: wrong validation signature")
        primary = as_int(qc, "primary_records", qc_path)
        primary_mapped = as_int(qc, "primary_mapped_records", qc_path)
        human = as_int(qc, "primary_human_records", qc_path)
        pig = as_int(qc, "primary_pig_records", qc_path)
        other_contig = as_int(
            qc, "primary_other_contig_records", qc_path
        )
        checked_tags = as_int(
            qc, "checked_primary_tag_records", qc_path
        )
        bam_records = as_int(qc, "bam_alignment_records", qc_path)
        if primary != expected:
            raise ValueError(
                f"{qc_path}: primary={primary}, expected={expected}"
            )
        if human + pig + other_contig != primary_mapped:
            raise ValueError(f"{qc_path}: species counts do not close")
        if checked_tags != min(expected, 1000):
            raise ValueError(f"{qc_path}: incomplete tag validation")

        star = parse_star_log(log_path)
        star_input = as_int(star, "Number of input reads", log_path)
        unique = as_int(star, "Uniquely mapped reads number", log_path)
        unique_percent = as_percent(
            star, "Uniquely mapped reads %", log_path
        )
        multi = as_int(
            star, "Number of reads mapped to multiple loci", log_path
        )
        multi_percent = as_percent(
            star, "% of reads mapped to multiple loci", log_path
        )
        mismatch_percent = as_percent(
            star, "% of reads unmapped: too many mismatches", log_path
        )
        too_short_percent = as_percent(
            star, "% of reads unmapped: too short", log_path
        )
        other_percent = as_percent(
            star, "% of reads unmapped: other", log_path
        )
        if star_input != expected:
            raise ValueError(
                f"{log_path}: STAR input={star_input}, expected={expected}"
            )
        if unique + multi != primary_mapped:
            raise ValueError(
                f"{pair_id}: STAR unique+multi={unique + multi}, "
                f"BAM primary mapped={primary_mapped}"
            )

        mapped_denominator = max(primary_mapped, 1)
        row: dict[str, object] = {
            "pair_id": pair_id,
            "sample": sample,
            "expected_R1_reads": expected,
            "star_input_reads": star_input,
            "unique_mapped_reads": unique,
            "unique_mapped_percent": f"{unique_percent:.6f}",
            "multi_mapped_reads": multi,
            "multi_mapped_percent": f"{multi_percent:.6f}",
            "primary_mapped_reads": primary_mapped,
            "primary_mapping_percent": (
                f"{primary_mapped / expected * 100:.6f}"
            ),
            "primary_human_reads": human,
            "primary_pig_reads": pig,
            "primary_other_contig_reads": other_contig,
            "human_percent_of_primary_mapped": (
                f"{human / mapped_denominator * 100:.6f}"
            ),
            "pig_percent_of_primary_mapped": (
                f"{pig / mapped_denominator * 100:.6f}"
            ),
            "unmapped_too_many_mismatches_percent": (
                f"{mismatch_percent:.6f}"
            ),
            "unmapped_too_short_percent": f"{too_short_percent:.6f}",
            "unmapped_other_percent": f"{other_percent:.6f}",
            "checked_primary_tag_records": checked_tags,
            "bam_alignment_records": bam_records,
            "bam_path": str(bam),
        }
        rows.append(row)

        totals["input"] += expected
        totals["unique"] += unique
        totals["multi"] += multi
        totals["primary_mapped"] += primary_mapped
        totals["human"] += human
        totals["pig"] += pig
        totals["other_contig"] += other_contig
        totals["bam_records"] += bam_records
        totals["checked_tags"] += checked_tags
        weighted_unmapped["too_many_mismatches"] += (
            mismatch_percent * expected
        )
        weighted_unmapped["too_short"] += too_short_percent * expected
        weighted_unmapped["other"] += other_percent * expected

    rows.sort(key=lambda row: (str(row["sample"]), str(row["pair_id"])))
    atomic_tsv(args.per_file, PER_FILE_FIELDS, rows)

    input_total = max(totals["input"], 1)
    mapped_total = max(totals["primary_mapped"], 1)
    summary_rows: list[tuple[str, object, str]] = [
        ("mode", args.mode, "mapping mode"),
        ("completed_R1_files", len(rows), "R1 files"),
        ("STAR_input_R1_reads", totals["input"], "R1 reads"),
        ("unique_mapped_R1_reads", totals["unique"], "R1 reads"),
        (
            "unique_mapping_percent",
            f"{totals['unique'] / input_total * 100:.6f}",
            "%",
        ),
        ("multi_mapped_R1_reads", totals["multi"], "R1 reads"),
        (
            "multi_mapping_percent",
            f"{totals['multi'] / input_total * 100:.6f}",
            "%",
        ),
        (
            "primary_mapped_R1_reads",
            totals["primary_mapped"],
            "R1 reads",
        ),
        (
            "primary_mapping_percent",
            f"{totals['primary_mapped'] / input_total * 100:.6f}",
            "%",
        ),
        ("primary_human_R1_reads", totals["human"], "R1 reads"),
        (
            "human_percent_of_primary_mapped",
            f"{totals['human'] / mapped_total * 100:.6f}",
            "%",
        ),
        ("primary_pig_R1_reads", totals["pig"], "R1 reads"),
        (
            "pig_percent_of_primary_mapped",
            f"{totals['pig'] / mapped_total * 100:.6f}",
            "%",
        ),
        (
            "primary_other_contig_R1_reads",
            totals["other_contig"],
            "R1 reads",
        ),
        (
            "unmapped_too_many_mismatches_percent",
            f"{weighted_unmapped['too_many_mismatches'] / input_total:.6f}",
            "%",
        ),
        (
            "unmapped_too_short_percent",
            f"{weighted_unmapped['too_short'] / input_total:.6f}",
            "%",
        ),
        (
            "unmapped_other_percent",
            f"{weighted_unmapped['other'] / input_total:.6f}",
            "%",
        ),
        (
            "checked_primary_tag_records",
            totals["checked_tags"],
            "primary records",
        ),
        (
            "BAM_alignment_records",
            totals["bam_records"],
            "alignment records",
        ),
    ]
    atomic_summary(args.summary, summary_rows)
    print(
        f"[OK] RNA STAR {args.mode}: files={len(rows)} "
        f"input={totals['input']} "
        f"primary_mapping="
        f"{totals['primary_mapped'] / input_total * 100:.3f}% "
        f"human={totals['human'] / mapped_total * 100:.3f}% "
        f"pig={totals['pig'] / mapped_total * 100:.3f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
