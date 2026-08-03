#!/usr/bin/env python3
"""Step 12: build a species-prefixed human+pig FASTA/GTF for RNA mapping."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
from pathlib import Path
from typing import TextIO


ATTRIBUTE_RE = re.compile(
    r'(?P<key>gene_id|transcript_id|exon_id|protein_id) "'
    r'(?P<value>[^"]+)"'
)
PREFIX_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*_$")


def open_text(path: Path) -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def prefixed_attributes(attributes: str, id_prefix: str) -> str:
    def replace(match: re.Match[str]) -> str:
        key = match.group("key")
        value = match.group("value")
        if not value.startswith(id_prefix):
            value = id_prefix + value
        return f'{key} "{value}"'

    return ATTRIBUTE_RE.sub(replace, attributes)


def write_fasta(
    source: Path,
    destination: TextIO,
    contig_prefix: str,
    observed_names: set[str],
    digest: hashlib._Hash,
) -> tuple[int, int]:
    sequences = 0
    bases = 0
    current_name: str | None = None
    with open_text(source) as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith(">"):
                header = line[1:].rstrip("\r\n")
                if not header:
                    raise ValueError(f"{source}:{line_no}: empty FASTA header")
                fields = header.split(maxsplit=1)
                source_name = fields[0]
                current_name = contig_prefix + source_name
                if current_name in observed_names:
                    raise ValueError(
                        f"duplicate combined contig name: {current_name}"
                    )
                observed_names.add(current_name)
                description = f" {fields[1]}" if len(fields) == 2 else ""
                output = f">{current_name}{description}\n"
                destination.write(output)
                digest.update(output.encode())
                sequences += 1
                continue
            if current_name is None:
                raise ValueError(
                    f"{source}:{line_no}: sequence before first FASTA header"
                )
            sequence = line.strip()
            if not sequence:
                continue
            output = sequence + "\n"
            destination.write(output)
            digest.update(output.encode())
            bases += len(sequence)
    if sequences == 0 or bases == 0:
        raise ValueError(f"{source}: no FASTA sequences")
    return sequences, bases


def write_gtf(
    source: Path,
    destination: TextIO,
    contig_prefix: str,
    id_prefix: str,
    fasta_names: set[str],
    digest: hashlib._Hash,
) -> tuple[int, int, set[str]]:
    records = 0
    exons = 0
    missing_contigs: set[str] = set()
    with open_text(source) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{source}:{line_no}: expected 9 GTF columns, "
                    f"observed {len(fields)}"
                )
            fields[0] = contig_prefix + fields[0]
            if fields[0] not in fasta_names:
                missing_contigs.add(fields[0])
            fields[8] = prefixed_attributes(fields[8], id_prefix)
            output = "\t".join(fields) + "\n"
            destination.write(output)
            digest.update(output.encode())
            records += 1
            if fields[2] == "exon":
                exons += 1
    if records == 0 or exons == 0:
        raise ValueError(f"{source}: no usable GTF records/exons")
    return records, exons, missing_contigs


def source_identity(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_size, int(stat.st_mtime)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human-fasta", type=Path, required=True)
    parser.add_argument("--human-gtf", type=Path, required=True)
    parser.add_argument("--pig-fasta", type=Path, required=True)
    parser.add_argument("--pig-gtf", type=Path, required=True)
    parser.add_argument("--out-fasta", type=Path, required=True)
    parser.add_argument("--out-gtf", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--human-contig-prefix", default="human_")
    parser.add_argument("--pig-contig-prefix", default="pig_")
    parser.add_argument("--human-id-prefix", default="HUMAN_")
    parser.add_argument("--pig-id-prefix", default="PIG_")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    for name in (
        "human_contig_prefix",
        "pig_contig_prefix",
        "human_id_prefix",
        "pig_id_prefix",
    ):
        value = getattr(args, name)
        if not PREFIX_RE.fullmatch(value):
            parser.error(
                f"--{name.replace('_', '-')} must end in '_' and contain "
                f"only letters, digits, dot, dash, or underscore: {value!r}"
            )
    if args.human_contig_prefix == args.pig_contig_prefix:
        parser.error("human and pig contig prefixes must differ")
    if args.human_id_prefix == args.pig_id_prefix:
        parser.error("human and pig ID prefixes must differ")

    sources = {
        "human_fasta": args.human_fasta,
        "human_gtf": args.human_gtf,
        "pig_fasta": args.pig_fasta,
        "pig_gtf": args.pig_gtf,
    }
    for label, path in sources.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"{label} missing/empty: {path}")

    outputs = (args.out_fasta, args.out_gtf, args.manifest)
    existing = [path for path in outputs if path.exists()]
    if existing and not args.force:
        raise FileExistsError(
            "output already exists; use --force only after reviewing: "
            + ", ".join(str(path) for path in existing)
        )

    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    fasta_tmp = args.out_fasta.with_name(
        args.out_fasta.name + f".tmp.{os.getpid()}"
    )
    gtf_tmp = args.out_gtf.with_name(
        args.out_gtf.name + f".tmp.{os.getpid()}"
    )
    manifest_tmp = args.manifest.with_name(
        args.manifest.name + f".tmp.{os.getpid()}"
    )

    fasta_digest = hashlib.sha256()
    gtf_digest = hashlib.sha256()
    all_names: set[str] = set()
    metrics: dict[str, int | str] = {}

    try:
        with fasta_tmp.open("wt", encoding="utf-8", newline="\n") as out:
            for species, source, prefix in (
                ("human", args.human_fasta, args.human_contig_prefix),
                ("pig", args.pig_fasta, args.pig_contig_prefix),
            ):
                sequences, bases = write_fasta(
                    source,
                    out,
                    prefix,
                    all_names,
                    fasta_digest,
                )
                metrics[f"{species}_fasta_sequences"] = sequences
                metrics[f"{species}_fasta_bases"] = bases

        with gtf_tmp.open("wt", encoding="utf-8", newline="\n") as out:
            header = (
                "# snMAR joint RNA reference; contigs and stable feature IDs "
                "are species-prefixed\n"
            )
            out.write(header)
            gtf_digest.update(header.encode())
            missing_all: set[str] = set()
            for species, source, contig_prefix, id_prefix in (
                (
                    "human",
                    args.human_gtf,
                    args.human_contig_prefix,
                    args.human_id_prefix,
                ),
                (
                    "pig",
                    args.pig_gtf,
                    args.pig_contig_prefix,
                    args.pig_id_prefix,
                ),
            ):
                records, exons, missing = write_gtf(
                    source,
                    out,
                    contig_prefix,
                    id_prefix,
                    all_names,
                    gtf_digest,
                )
                metrics[f"{species}_gtf_records"] = records
                metrics[f"{species}_gtf_exons"] = exons
                missing_all.update(missing)
            if missing_all:
                preview = ", ".join(sorted(missing_all)[:20])
                raise ValueError(
                    "GTF contig(s) absent from combined FASTA: " + preview
                )

        manifest_rows: list[tuple[str, str]] = [
            ("pipeline_signature", "snmar_joint_rna_reference_v10"),
            ("human_contig_prefix", args.human_contig_prefix),
            ("pig_contig_prefix", args.pig_contig_prefix),
            ("human_id_prefix", args.human_id_prefix),
            ("pig_id_prefix", args.pig_id_prefix),
        ]
        for label, path in sources.items():
            size, mtime = source_identity(path)
            manifest_rows.extend(
                (
                    (f"{label}_path", str(path)),
                    (f"{label}_bytes", str(size)),
                    (f"{label}_mtime", str(mtime)),
                )
            )
        for key, value in metrics.items():
            manifest_rows.append((key, str(value)))
        manifest_rows.extend(
            (
                ("joint_fasta_sha256", fasta_digest.hexdigest()),
                ("joint_gtf_sha256", gtf_digest.hexdigest()),
            )
        )
        with manifest_tmp.open("wt", encoding="utf-8", newline="\n") as out:
            out.write("metric\tvalue\n")
            for key, value in manifest_rows:
                out.write(f"{key}\t{value}\n")

        os.replace(fasta_tmp, args.out_fasta)
        os.replace(gtf_tmp, args.out_gtf)
        os.replace(manifest_tmp, args.manifest)
    except BaseException:
        for path in (fasta_tmp, gtf_tmp, manifest_tmp):
            path.unlink(missing_ok=True)
        raise

    print(
        "[OK] joint RNA reference prepared: "
        f"contigs={len(all_names)} "
        f"bases={int(metrics['human_fasta_bases']) + int(metrics['pig_fasta_bases'])} "
        f"gtf_records={int(metrics['human_gtf_records']) + int(metrics['pig_gtf_records'])}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
