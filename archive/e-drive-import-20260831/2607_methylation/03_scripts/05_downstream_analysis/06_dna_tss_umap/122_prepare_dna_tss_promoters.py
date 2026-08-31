#!/usr/bin/env python3
"""Step 122: build clipped, non-overlapping TSS promoter clusters.

Each gene contributes a [TSS-flank, TSS+flank+1) genomic interval. Windows
are clipped to contig boundaries and overlapping windows on the same contig
are merged so a cytosine is never counted more than once. A conservative
cross-species feature map is also created from single-gene clusters whose
case-normalized gene name occurs exactly once in human and once in pig.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable


PIPELINE_SIGNATURE = "snmar_dna_tss_promoter_index_20260808"
SOURCE_SCOPE = "validated_complete_single_cytosine_96wells_768shards"
ATTRIBUTE_PATTERN = re.compile(r'^\s*([^\s=]+)\s+(?:"([^"]*)"|([^;\s]+))\s*$')


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return path.open(mode + "t", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with open_text(path, "r") as handle:
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
        with opener(temporary, "wt", encoding="utf-8", newline="") as handle:
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


def atomic_metrics(path: Path, rows: list[tuple[str, object]]) -> None:
    atomic_tsv(
        path,
        ["metric", "value"],
        ({"metric": metric, "value": value} for metric, value in rows),
    )


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for raw in text.strip().strip(";").split(";"):
        part = raw.strip()
        if not part:
            continue
        match = ATTRIBUTE_PATTERN.match(part)
        if match:
            key = match.group(1)
            value = match.group(2) if match.group(2) is not None else match.group(3)
        elif "=" in part:
            key, value = part.split("=", 1)
            key, value = key.strip(), value.strip().strip('"')
        else:
            continue
        if key and key not in attributes:
            attributes[key] = value
    return attributes


def load_contigs(path: Path) -> list[dict[str, object]]:
    rows = sorted(read_tsv(path), key=lambda row: int(row["contig_index"]))
    if [int(row["contig_index"]) for row in rows] != list(range(len(rows))):
        raise ValueError("contig indices must be zero-based and contiguous")
    result: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in rows:
        contig = row["contig"]
        species = row["species"]
        length = int(row["length"])
        if contig in seen or species not in {"human", "pig"} or length < 1:
            raise ValueError(f"invalid contig row: {row}")
        seen.add(contig)
        result.append(
            {
                "contig": contig,
                "species": species,
                "length": length,
                "contig_index": int(row["contig_index"]),
            }
        )
    return result


def load_gene_windows(
    gtf: Path,
    contigs: list[dict[str, object]],
    flank: int,
) -> list[dict[str, object]]:
    if flank < 1:
        raise ValueError("TSS flank must be positive")
    by_name = {str(row["contig"]): row for row in contigs}
    genes: list[dict[str, object]] = []
    seen_gene_ids: set[str] = set()
    with open_text(gtf, "r") as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gtf}:{line_no}: expected 9 GTF fields")
            if fields[2] != "gene":
                continue
            chrom, start_text, end_text, strand = (
                fields[0], fields[3], fields[4], fields[6]
            )
            if chrom not in by_name:
                raise ValueError(f"{gtf}:{line_no}: unknown contig {chrom}")
            if strand not in {"+", "-"}:
                raise ValueError(f"{gtf}:{line_no}: invalid gene strand")
            start1, end1 = int(start_text), int(end_text)
            contig = by_name[chrom]
            length = int(contig["length"])
            if start1 < 1 or end1 < start1 or end1 > length:
                raise ValueError(f"{gtf}:{line_no}: gene outside contig")
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("gene_id", "")
            if not gene_id or gene_id in seen_gene_ids:
                raise ValueError(f"{gtf}:{line_no}: missing/duplicate gene_id")
            seen_gene_ids.add(gene_id)
            tss0 = start1 - 1 if strand == "+" else end1 - 1
            start0 = max(0, tss0 - flank)
            end0 = min(length, tss0 + flank + 1)
            genes.append(
                {
                    "species": contig["species"],
                    "contig_index": contig["contig_index"],
                    "chrom": chrom,
                    "contig_length": length,
                    "start0": start0,
                    "end0": end0,
                    "tss0": tss0,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_name": attributes.get("gene_name", ""),
                    "gene_biotype": attributes.get(
                        "gene_biotype", attributes.get("gene_type", "")
                    ),
                    "clipped_left": int(start0 == 0 and tss0 < flank),
                    "clipped_right": int(
                        end0 == length and tss0 + flank + 1 > length
                    ),
                }
            )
    if not genes:
        raise ValueError(f"no gene features found in {gtf}")
    return genes


def build_promoter_clusters(
    windows: list[dict[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in windows:
        grouped[int(row["contig_index"])].append(row)
    clusters: list[dict[str, object]] = []
    for contig_index in sorted(grouped):
        ordered = sorted(
            grouped[contig_index],
            key=lambda row: (int(row["start0"]), int(row["end0"]), str(row["gene_id"])),
        )
        members: list[dict[str, object]] = []
        cluster_start = cluster_end = -1
        for row in ordered:
            start0, end0 = int(row["start0"]), int(row["end0"])
            if not members or start0 < cluster_end:
                if not members:
                    cluster_start = start0
                    cluster_end = end0
                else:
                    cluster_end = max(cluster_end, end0)
                members.append(row)
                continue
            clusters.append(_finalize_cluster(cluster_start, cluster_end, members))
            cluster_start, cluster_end, members = start0, end0, [row]
        if members:
            clusters.append(_finalize_cluster(cluster_start, cluster_end, members))
    clusters.sort(
        key=lambda row: (int(row["contig_index"]), int(row["start0"]), int(row["end0"]))
    )
    for feature_index, row in enumerate(clusters):
        row["feature_index"] = feature_index
        row["feature_id"] = f"TSSP{feature_index + 1:07d}"
        row["ortholog_index"] = -1
        row["ortholog_key"] = ""
    if not clusters:
        raise ValueError("no promoter clusters generated")
    return clusters


def _finalize_cluster(
    start0: int,
    end0: int,
    members: list[dict[str, object]],
) -> dict[str, object]:
    first = members[0]
    gene_ids = sorted({str(row["gene_id"]) for row in members})
    gene_names = sorted(
        {str(row["gene_name"]) for row in members if str(row["gene_name"])}
    )
    biotypes = sorted(
        {str(row["gene_biotype"]) for row in members if str(row["gene_biotype"])}
    )
    return {
        "feature_index": -1,
        "feature_id": "",
        "species": first["species"],
        "contig_index": first["contig_index"],
        "chrom": first["chrom"],
        "start0": start0,
        "end0": end0,
        "width": end0 - start0,
        "source_tss_windows": len(members),
        "gene_count": len(gene_ids),
        "gene_ids": ",".join(gene_ids),
        "gene_names": ",".join(gene_names),
        "gene_biotypes": ",".join(biotypes),
        "tss_positions0": ",".join(
            map(str, sorted({int(row["tss0"]) for row in members}))
        ),
        "strands": ",".join(sorted({str(row["strand"]) for row in members})),
        "clipped_left": max(int(row["clipped_left"]) for row in members),
        "clipped_right": max(int(row["clipped_right"]) for row in members),
        "ortholog_index": -1,
        "ortholog_key": "",
    }


def add_conservative_ortholog_map(
    clusters: list[dict[str, object]],
) -> list[dict[str, object]]:
    candidates: dict[str, dict[str, list[int]]] = defaultdict(
        lambda: {"human": [], "pig": []}
    )
    for row in clusters:
        names = [value for value in str(row["gene_names"]).split(",") if value]
        if int(row["gene_count"]) != 1 or len(names) != 1:
            continue
        key = names[0].strip().upper()
        species = str(row["species"])
        if key and species in {"human", "pig"}:
            candidates[key][species].append(int(row["feature_index"]))
    pairs: list[dict[str, object]] = []
    for key in sorted(candidates):
        human = candidates[key]["human"]
        pig = candidates[key]["pig"]
        if len(human) != 1 or len(pig) != 1:
            continue
        ortholog_index = len(pairs)
        for feature_index in (human[0], pig[0]):
            clusters[feature_index]["ortholog_index"] = ortholog_index
            clusters[feature_index]["ortholog_key"] = key
        pairs.append(
            {
                "ortholog_index": ortholog_index,
                "ortholog_key": key,
                "human_feature_index": human[0],
                "human_feature_id": clusters[human[0]]["feature_id"],
                "human_gene_id": clusters[human[0]]["gene_ids"],
                "pig_feature_index": pig[0],
                "pig_feature_id": clusters[pig[0]]["feature_id"],
                "pig_gene_id": clusters[pig[0]]["gene_ids"],
                "pairing_rule": "single_gene_cluster_unique_case_normalized_gene_name",
            }
        )
    if not pairs:
        raise ValueError("no conservative human/pig promoter pairs found")
    return pairs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--flank", type=int, default=1500)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.flank < 1:
        raise ValueError("TSS flank must be positive")
    project = args.project
    gtf = args.gtf or (
        project / "02_reference/rna_human_pig_joint/human_pig.prefixed.gtf"
    )
    contig_path = (
        project
        / "08_qc_summary/20_single_cytosine_inputs/"
        "dna_single_cytosine_contig_index.tsv"
    )
    source_summary = (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete/"
        "dna_single_cytosine_complete_summary.tsv"
    )
    source = read_metrics(source_summary)
    if (
        source.get("status") != "validated"
        or int(source.get("validated_samples", "0")) != 96
        or int(source.get("source_h5_shards", "0")) != 768
    ):
        raise ValueError("complete single-cytosine source is not current/validated")
    for path in (gtf, contig_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    output_root = args.output_root or (
        project / "10_downstream_analysis/08_dna_tss_promoters"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/38_dna_tss_promoters"
    )
    summary_path = qc_root / "dna_tss_promoter_index_summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and int(existing.get("tss_flank_bp", "0")) == args.flank
        ):
            print(f"[SKIP] validated TSS index exists: {summary_path}")
            return 0
        raise FileExistsError(f"inspect or rerun with --force: {summary_path}")

    contigs = load_contigs(contig_path)
    windows = load_gene_windows(gtf, contigs, args.flank)
    clusters = build_promoter_clusters(windows)
    orthologs = add_conservative_ortholog_map(clusters)
    index_path = output_root / "dna_tss_promoter_index.tsv.gz"
    ortholog_path = output_root / "dna_tss_conservative_ortholog_index.tsv.gz"
    atomic_tsv(index_path, list(clusters[0]), clusters, gzip_output=True)
    atomic_tsv(ortholog_path, list(orthologs[0]), orthologs, gzip_output=True)
    by_species = {
        species: [row for row in clusters if row["species"] == species]
        for species in ("human", "pig")
    }
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("source_single_cytosine_scope", SOURCE_SCOPE),
        ("source_complete_summary", source_summary),
        ("joint_gtf", gtf),
        ("contig_index", contig_path),
        ("tss_flank_bp", args.flank),
        ("window_definition", "[TSS-flank,TSS+flank+1)_clipped_to_contig"),
        ("overlap_policy", "merge_overlapping_windows_within_contig"),
        ("touching_window_policy", "retain_as_separate_nonoverlapping_clusters"),
        ("gene_features", len(windows)),
        ("human_gene_features", sum(row["species"] == "human" for row in windows)),
        ("pig_gene_features", sum(row["species"] == "pig" for row in windows)),
        ("promoter_clusters", len(clusters)),
        ("human_promoter_clusters", len(by_species["human"])),
        ("pig_promoter_clusters", len(by_species["pig"])),
        ("left_clipped_clusters", sum(int(row["clipped_left"]) for row in clusters)),
        ("right_clipped_clusters", sum(int(row["clipped_right"]) for row in clusters)),
        ("conservative_ortholog_features", len(orthologs)),
        ("ortholog_pairing_rule", "single_gene_cluster_unique_case_normalized_gene_name"),
        ("promoter_index", index_path),
        ("ortholog_index", ortholog_path),
    ]
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] TSS promoter index: genes={len(windows)} "
        f"clusters={len(clusters)} orthologs={len(orthologs)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
