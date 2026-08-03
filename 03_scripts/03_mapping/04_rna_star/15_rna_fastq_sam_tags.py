#!/usr/bin/env python3
"""Step 15: protect RNA barcode/UMI tags across FASTQ -> STAR -> BAM."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from pathlib import Path
from typing import TextIO


SIGNATURE = "|||snm10|"
TAG_ORDER = ("CR", "CB", "CN", "UR", "UB", "XD", "XE", "XL", "XI")
TAG_TYPES = {
    "CR": "Z",
    "CB": "Z",
    "CN": "Z",
    "UR": "Z",
    "UB": "Z",
    "XD": "i",
    "XE": "Z",
    "XL": "i",
    "XI": "Z",
}
TAG_RE = re.compile(
    r"(?:^|[ \t])"
    r"(?P<name>CR|CB|CN|UR|UB|XD|XE|XL|XI):"
    r"(?P<type>[AifZHB]):(?P<value>[^ \t\r\n]+)"
)
MAX_QNAME_LENGTH = 254


def open_input(path_text: str) -> TextIO:
    if path_text == "-":
        return sys.stdin
    path = Path(path_text)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def open_output(path_text: str) -> TextIO:
    if path_text == "-":
        return sys.stdout
    path = Path(path_text)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("xt", encoding="utf-8", newline="\n")


def extract_tags(header: str) -> tuple[str, list[str]]:
    if not header.startswith("@"):
        raise ValueError(f"FASTQ header does not start with @: {header[:120]}")
    read_id = header[1:].split(maxsplit=1)[0]
    if not read_id or SIGNATURE in read_id:
        raise ValueError(f"invalid/already encoded read ID: {read_id!r}")

    tags: dict[str, str] = {}
    for match in TAG_RE.finditer(header):
        name = match.group("name")
        tag_type = match.group("type")
        value = match.group("value")
        expected_type = TAG_TYPES[name]
        if tag_type != expected_type:
            raise ValueError(
                f"{read_id}: {name} has type {tag_type}, "
                f"expected {expected_type}"
            )
        if name in tags:
            raise ValueError(f"{read_id}: duplicate tag {name}")
        if "|" in value or any(character.isspace() for character in value):
            raise ValueError(f"{read_id}: unsafe value for {name}: {value!r}")
        tags[name] = f"{name}:{tag_type}:{value}"

    missing = [name for name in TAG_ORDER if name not in tags]
    if missing:
        raise ValueError(f"{read_id}: missing FASTQ tag(s): {missing}")
    if tags["XD"] != "XD:i:0":
        raise ValueError(f"{read_id}: non-strict XD tag: {tags['XD']}")
    if tags["XI"] != "XI:Z:full10":
        raise ValueError(f"{read_id}: unexpected XI tag: {tags['XI']}")
    return read_id, [tags[name] for name in TAG_ORDER]


def encode_qname(read_id: str, tags: list[str]) -> str:
    encoded = read_id + SIGNATURE + "|".join(tags)
    if len(encoded) > MAX_QNAME_LENGTH:
        raise ValueError(
            f"encoded QNAME length={len(encoded)} exceeds "
            f"{MAX_QNAME_LENGTH}: {read_id}"
        )
    return encoded


def decode_qname(qname: str) -> tuple[str, list[str]]:
    if qname.count(SIGNATURE) != 1:
        raise ValueError(f"QNAME lacks unique snMAR v10 signature: {qname[:120]}")
    read_id, payload = qname.split(SIGNATURE, 1)
    encoded_tags = payload.split("|")
    if len(encoded_tags) != len(TAG_ORDER):
        raise ValueError(
            f"{qname[:120]}: encoded tag count={len(encoded_tags)}, "
            f"expected={len(TAG_ORDER)}"
        )
    parsed: dict[str, str] = {}
    for tag in encoded_tags:
        fields = tag.split(":", 2)
        if len(fields) != 3:
            raise ValueError(f"{qname[:120]}: malformed encoded tag {tag!r}")
        name, tag_type, value = fields
        if name not in TAG_TYPES or tag_type != TAG_TYPES[name] or not value:
            raise ValueError(f"{qname[:120]}: invalid encoded tag {tag!r}")
        if name in parsed:
            raise ValueError(f"{qname[:120]}: duplicate encoded tag {name}")
        parsed[name] = tag
    missing = [name for name in TAG_ORDER if name not in parsed]
    if missing:
        raise ValueError(f"{qname[:120]}: missing encoded tags {missing}")
    if parsed["XD"] != "XD:i:0" or parsed["XI"] != "XI:Z:full10":
        raise ValueError(f"{qname[:120]}: non-strict encoded metadata")
    return read_id, [parsed[name] for name in TAG_ORDER]


def read_fastq_record(
    handle: TextIO,
    source: str,
    record_no: int,
) -> tuple[str, str, str] | None:
    raw_header = handle.readline()
    if raw_header == "":
        return None
    header = raw_header.rstrip("\r\n")
    sequence = handle.readline().rstrip("\r\n")
    plus = handle.readline().rstrip("\r\n")
    quality = handle.readline().rstrip("\r\n")
    if not sequence or not plus or not quality:
        raise ValueError(f"{source}: truncated FASTQ near record {record_no}")
    if not plus.startswith("+"):
        raise ValueError(f"{source}: invalid plus line near record {record_no}")
    if len(sequence) != len(quality):
        raise ValueError(
            f"{source}: sequence/quality mismatch near record {record_no}"
        )
    return header, sequence, quality


def write_stats(path_text: str | None, rows: list[tuple[str, object]]) -> None:
    if path_text is None:
        return
    path = Path(path_text)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        with tmp.open("wt", encoding="utf-8", newline="\n") as out:
            out.write("metric\tvalue\n")
            for metric, value in rows:
                out.write(f"{metric}\t{value}\n")
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def encode_fastq(args: argparse.Namespace) -> int:
    count = 0
    max_qname = 0
    source = args.input
    with open_input(source) as source_handle, open_output(args.output) as out:
        while True:
            record = read_fastq_record(source_handle, source, count + 1)
            if record is None:
                break
            header, sequence, quality = record
            read_id, tags = extract_tags(header)
            qname = encode_qname(read_id, tags)
            out.write(f"@{qname}\n{sequence}\n+\n{quality}\n")
            count += 1
            max_qname = max(max_qname, len(qname))
    if count == 0:
        raise ValueError(f"{source}: no FASTQ records")
    write_stats(
        args.stats,
        [
            ("pipeline_signature", "snmar_rna_qname_encode_v10"),
            ("encoded_reads", count),
            ("max_encoded_qname_length", max_qname),
        ],
    )
    print(
        f"[QNAME_ENCODE] reads={count} max_qname_length={max_qname}",
        file=sys.stderr,
    )
    return 0


def restore_sam(args: argparse.Namespace) -> int:
    alignments = 0
    with open_input(args.input) as source, open_output(args.output) as out:
        for line_no, line in enumerate(source, 1):
            if line.startswith("@"):
                out.write(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValueError(
                    f"SAM line {line_no}: expected at least 11 columns"
                )
            read_id, tags = decode_qname(fields[0])
            existing = {
                field[:2]
                for field in fields[11:]
                if len(field) >= 5 and field[2] == ":" and field[4] == ":"
            }
            duplicates = existing.intersection(TAG_ORDER)
            if duplicates:
                raise ValueError(
                    f"SAM line {line_no}: duplicate tag(s) {sorted(duplicates)}"
                )
            fields[0] = read_id
            fields.extend(tags)
            out.write("\t".join(fields) + "\n")
            alignments += 1
    if alignments == 0:
        raise ValueError("SAM stream contained no alignments")
    write_stats(
        args.stats,
        [
            ("pipeline_signature", "snmar_rna_qname_restore_v10"),
            ("restored_alignment_records", alignments),
        ],
    )
    print(
        f"[QNAME_RESTORE] alignment_records={alignments}",
        file=sys.stderr,
    )
    return 0


def parse_aux_tags(fields: list[str], line_no: int) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in fields:
        pieces = field.split(":", 2)
        if len(pieces) != 3:
            continue
        name, tag_type, value = pieces
        if name in TAG_TYPES:
            if name in result:
                raise ValueError(f"SAM line {line_no}: duplicate tag {name}")
            if tag_type != TAG_TYPES[name] or not value:
                raise ValueError(
                    f"SAM line {line_no}: invalid {name} field {field!r}"
                )
            result[name] = field
    return result


def validate_sam(args: argparse.Namespace) -> int:
    metrics = {
        "bam_alignment_records": 0,
        "primary_records": 0,
        "primary_mapped_records": 0,
        "primary_unmapped_records": 0,
        "primary_human_records": 0,
        "primary_pig_records": 0,
        "primary_other_contig_records": 0,
        "secondary_records": 0,
        "supplementary_records": 0,
        "checked_primary_tag_records": 0,
    }
    coordinate_header = False
    with open_input(args.input) as source:
        for line_no, line in enumerate(source, 1):
            if line.startswith("@HD"):
                coordinate_header = "SO:coordinate" in line.rstrip("\r\n").split(
                    "\t"
                )
                continue
            if line.startswith("@"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValueError(
                    f"SAM line {line_no}: expected at least 11 columns"
                )
            if SIGNATURE in fields[0]:
                raise ValueError(
                    f"SAM line {line_no}: encoded QNAME was not restored"
                )
            flag = int(fields[1])
            metrics["bam_alignment_records"] += 1
            if flag & 0x100:
                metrics["secondary_records"] += 1
            if flag & 0x800:
                metrics["supplementary_records"] += 1
            if flag & (0x100 | 0x800):
                continue

            metrics["primary_records"] += 1
            if (
                metrics["checked_primary_tag_records"]
                < args.check_primary_tags
            ):
                tags = parse_aux_tags(fields[11:], line_no)
                missing = [name for name in TAG_ORDER if name not in tags]
                if missing:
                    raise ValueError(
                        f"SAM line {line_no}: missing restored tags {missing}"
                    )
                if tags["XD"] != "XD:i:0" or tags["XI"] != "XI:Z:full10":
                    raise ValueError(
                        f"SAM line {line_no}: non-strict restored tags"
                    )
                metrics["checked_primary_tag_records"] += 1

            rname = fields[2]
            if flag & 0x4 or rname == "*":
                metrics["primary_unmapped_records"] += 1
            else:
                metrics["primary_mapped_records"] += 1
                if rname.startswith("human_"):
                    metrics["primary_human_records"] += 1
                elif rname.startswith("pig_"):
                    metrics["primary_pig_records"] += 1
                else:
                    metrics["primary_other_contig_records"] += 1

    if not coordinate_header:
        raise ValueError("BAM/SAM header is not coordinate sorted")
    if metrics["primary_records"] != args.expected_primary_reads:
        raise ValueError(
            f"primary records={metrics['primary_records']}, "
            f"expected={args.expected_primary_reads}"
        )
    expected_checks = min(args.expected_primary_reads, args.check_primary_tags)
    if metrics["checked_primary_tag_records"] != expected_checks:
        raise ValueError(
            f"checked primary tag records="
            f"{metrics['checked_primary_tag_records']}, "
            f"expected={expected_checks}"
        )
    if (
        metrics["primary_human_records"]
        + metrics["primary_pig_records"]
        + metrics["primary_other_contig_records"]
        != metrics["primary_mapped_records"]
    ):
        raise AssertionError("species primary mapping counts do not close")

    rows: list[tuple[str, object]] = [
        ("pipeline_signature", "snmar_rna_bam_validate_v10")
    ]
    rows.extend(metrics.items())
    write_stats(args.stats, rows)
    print(
        "[BAM_VALIDATE] "
        f"primary={metrics['primary_records']} "
        f"mapped={metrics['primary_mapped_records']} "
        f"human={metrics['primary_human_records']} "
        f"pig={metrics['primary_pig_records']} "
        f"tags_checked={metrics['checked_primary_tag_records']}",
        file=sys.stderr,
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    encode = subparsers.add_parser("encode-fastq")
    encode.add_argument("--input", required=True)
    encode.add_argument("--output", required=True)
    encode.add_argument("--stats")
    encode.set_defaults(function=encode_fastq)

    restore = subparsers.add_parser("restore-sam")
    restore.add_argument("--input", default="-")
    restore.add_argument("--output", default="-")
    restore.add_argument("--stats")
    restore.set_defaults(function=restore_sam)

    validate = subparsers.add_parser("validate-sam")
    validate.add_argument("--input", default="-")
    validate.add_argument("--expected-primary-reads", type=int, required=True)
    validate.add_argument("--check-primary-tags", type=int, default=1000)
    validate.add_argument("--stats", required=True)
    validate.set_defaults(function=validate_sam)

    args = parser.parse_args()
    if getattr(args, "expected_primary_reads", 1) < 1:
        parser.error("--expected-primary-reads must be positive")
    if getattr(args, "check_primary_tags", 0) < 0:
        parser.error("--check-primary-tags must be non-negative")
    return args.function(args)


if __name__ == "__main__":
    raise SystemExit(main())
