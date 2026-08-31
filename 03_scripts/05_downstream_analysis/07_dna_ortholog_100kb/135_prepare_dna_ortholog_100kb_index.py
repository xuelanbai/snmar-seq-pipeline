#!/usr/bin/env python3
"""Step 135: build a label-blind human/pig reciprocal-best 100 kb axis.

Human GRCh38 is divided into non-overlapping 100 kb windows.  Alignment
blocks from the UCSC hg38--susScr11 reciprocal-best chain are projected to
Sscrofa11.1.  Windows are retained only when one chain covers most of the
human window, the pig span has a plausible length, and retained pig spans do
not overlap.  Human and pig intervals are stored as two physical features
that share one ``orthobin_index``; downstream cells are never assigned a
species before matrix construction.

Reference WCG opportunities are counted on both strands.  GCG is excluded by
construction: plus-strand WCG is [AT]CG and minus-strand WCG appears as
CG[AT] on the forward reference.  This matches the single-cytosine caller.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator


PIPELINE_SIGNATURE = "snmar_dna_ortholog_100kb_index_20260825"
CHAIN_SOURCE = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsSusScr11/"
    "reciprocalBest/hg38.susScr11.rbest.chain.gz"
)


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
        ({"metric": key, "value": value} for key, value in rows),
    )


def load_contigs(path: Path) -> list[dict[str, object]]:
    rows = sorted(read_tsv(path), key=lambda row: int(row["contig_index"]))
    if [int(row["contig_index"]) for row in rows] != list(range(len(rows))):
        raise ValueError("contig indices must be zero-based and contiguous")
    result: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in rows:
        name = row["contig"]
        species = row["species"]
        length = int(row["length"])
        if name in seen or species not in {"human", "pig"} or length < 1:
            raise ValueError(f"invalid contig row: {row}")
        seen.add(name)
        result.append(
            {
                "contig": name,
                "species": species,
                "length": length,
                "contig_index": int(row["contig_index"]),
            }
        )
    return result


def ucsc_candidates(joint_name: str, species: str) -> list[str]:
    prefix = species + "_"
    if not joint_name.startswith(prefix):
        return []
    raw = joint_name[len(prefix) :]
    values = [raw]
    if raw in {"MT", "M", "chrMT"}:
        values.extend(["chrM", "MT"])
    elif raw.startswith("chr"):
        values.append(raw)
    else:
        values.append("chr" + raw)
    return list(dict.fromkeys(values))


def build_alias_lookup(
    contigs: list[dict[str, object]], species: str
) -> dict[tuple[str, int], dict[str, object]]:
    lookup: dict[tuple[str, int], dict[str, object]] = {}
    for row in contigs:
        if row["species"] != species:
            continue
        for alias in ucsc_candidates(str(row["contig"]), species):
            key = (alias, int(row["length"]))
            previous = lookup.setdefault(key, row)
            if previous is not row:
                raise ValueError(f"ambiguous UCSC alias: {key}")
    return lookup


def parse_chain_header(line: str) -> dict[str, object]:
    fields = line.split()
    if len(fields) != 13 or fields[0] != "chain":
        raise ValueError(f"invalid chain header: {line!r}")
    return {
        "score": int(fields[1]),
        "t_name": fields[2],
        "t_size": int(fields[3]),
        "t_strand": fields[4],
        "t_start": int(fields[5]),
        "t_end": int(fields[6]),
        "q_name": fields[7],
        "q_size": int(fields[8]),
        "q_strand": fields[9],
        "q_start": int(fields[10]),
        "q_end": int(fields[11]),
        "chain_id": int(fields[12]),
    }


def iter_chain_blocks(path: Path) -> Iterator[dict[str, object]]:
    """Yield aligned blocks in forward-reference coordinates."""
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty reciprocal-best chain: {path}")
    with open_text(path, "r") as handle:
        header: dict[str, object] | None = None
        t_position = q_position = 0
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                header = None
                continue
            if line.startswith("chain "):
                header = parse_chain_header(line)
                if header["t_strand"] != "+" or header["q_strand"] not in {"+", "-"}:
                    raise ValueError(f"unsupported chain strand at line {line_number}")
                t_position = int(header["t_start"])
                q_position = int(header["q_start"])
                continue
            if header is None:
                raise ValueError(f"chain block without header at line {line_number}")
            fields = line.split()
            if len(fields) not in {1, 3}:
                raise ValueError(f"invalid chain block at line {line_number}")
            size = int(fields[0])
            if size < 1:
                raise ValueError(f"non-positive chain block at line {line_number}")
            q_size = int(header["q_size"])
            if header["q_strand"] == "+":
                q_start0, q_end0 = q_position, q_position + size
            else:
                q_start0 = q_size - (q_position + size)
                q_end0 = q_size - q_position
            yield {
                **header,
                "t_block_start0": t_position,
                "t_block_end0": t_position + size,
                "q_block_start0": q_start0,
                "q_block_end0": q_end0,
                "block_size": size,
            }
            if len(fields) == 3:
                t_position += size + int(fields[1])
                q_position += size + int(fields[2])
            else:
                t_position += size
                q_position += size


def project_segment(block: dict[str, object], start0: int, end0: int) -> tuple[int, int]:
    offset = start0 - int(block["t_block_start0"])
    length = end0 - start0
    if block["q_strand"] == "+":
        q_start = int(block["q_block_start0"]) + offset
        return q_start, q_start + length
    q_end = int(block["q_block_end0"]) - offset
    return q_end - length, q_end


def collect_candidates(
    chain: Path,
    human_lookup: dict[tuple[str, int], dict[str, object]],
    pig_lookup: dict[tuple[str, int], dict[str, object]],
    bin_size: int,
) -> tuple[dict[tuple[int, int], list[dict[str, object]]], dict[str, int]]:
    accumulators: dict[
        tuple[int, int, int], dict[str, object]
    ] = {}
    counters = defaultdict(int)
    for block in iter_chain_blocks(chain):
        counters["chain_blocks"] += 1
        human = human_lookup.get((str(block["t_name"]), int(block["t_size"])))
        pig = pig_lookup.get((str(block["q_name"]), int(block["q_size"])))
        if human is None or pig is None:
            counters["blocks_without_joint_contig"] += 1
            continue
        block_start = int(block["t_block_start0"])
        block_end = int(block["t_block_end0"])
        first_bin = block_start // bin_size
        last_bin = (block_end - 1) // bin_size
        for ordinal in range(first_bin, last_bin + 1):
            human_start = ordinal * bin_size
            human_end = min(human_start + bin_size, int(human["length"]))
            overlap_start = max(block_start, human_start)
            overlap_end = min(block_end, human_end)
            if overlap_end <= overlap_start:
                continue
            q_start, q_end = project_segment(block, overlap_start, overlap_end)
            key = (int(human["contig_index"]), human_start, int(block["chain_id"]))
            state = accumulators.setdefault(
                key,
                {
                    "human": human,
                    "pig": pig,
                    "human_start0": human_start,
                    "human_end0": human_end,
                    "chain_id": int(block["chain_id"]),
                    "chain_score": int(block["score"]),
                    "q_strand": block["q_strand"],
                    "aligned_bases": 0,
                    "pig_start0": q_start,
                    "pig_end0": q_end,
                    "alignment_blocks": 0,
                },
            )
            if state["pig"] is not pig or state["q_strand"] != block["q_strand"]:
                raise ValueError("one chain candidate changes pig contig/strand")
            state["aligned_bases"] = int(state["aligned_bases"]) + overlap_end - overlap_start
            state["pig_start0"] = min(int(state["pig_start0"]), q_start)
            state["pig_end0"] = max(int(state["pig_end0"]), q_end)
            state["alignment_blocks"] = int(state["alignment_blocks"]) + 1
    by_bin: dict[tuple[int, int], list[dict[str, object]]] = defaultdict(list)
    for (contig_index, start0, _chain_id), row in accumulators.items():
        by_bin[(contig_index, start0)].append(row)
    counters["candidate_human_bins"] = len(by_bin)
    return by_bin, dict(counters)


def choose_pairs(
    candidates: dict[tuple[int, int], list[dict[str, object]]],
    min_aligned_fraction: float,
    min_span_ratio: float,
    max_span_ratio: float,
) -> tuple[list[dict[str, object]], dict[str, int]]:
    counters = defaultdict(int)
    eligible: list[dict[str, object]] = []
    for values in candidates.values():
        values.sort(
            key=lambda row: (
                -int(row["aligned_bases"]),
                -int(row["chain_score"]),
                int(row["chain_id"]),
            )
        )
        row = values[0]
        human_width = int(row["human_end0"]) - int(row["human_start0"])
        pig_width = int(row["pig_end0"]) - int(row["pig_start0"])
        aligned_fraction = int(row["aligned_bases"]) / human_width
        span_ratio = pig_width / human_width
        row["aligned_fraction"] = aligned_fraction
        row["pig_span_ratio"] = span_ratio
        if human_width <= 0 or pig_width <= 0:
            counters["invalid_span"] += 1
        elif aligned_fraction < min_aligned_fraction:
            counters["low_aligned_fraction"] += 1
        elif not min_span_ratio <= span_ratio <= max_span_ratio:
            counters["implausible_pig_span"] += 1
        else:
            eligible.append(row)
    # Greedy by alignment quality prevents any cytosine entering two axes.
    accepted: list[dict[str, object]] = []
    by_contig: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for row in sorted(
        eligible,
        key=lambda value: (
            -float(value["aligned_fraction"]),
            -int(value["aligned_bases"]),
            -int(value["chain_score"]),
        ),
    ):
        pig = row["pig"]
        contig = int(pig["contig_index"])
        start0, end0 = int(row["pig_start0"]), int(row["pig_end0"])
        if any(start0 < old_end and end0 > old_start for old_start, old_end in by_contig[contig]):
            counters["overlapping_pig_span"] += 1
            continue
        by_contig[contig].append((start0, end0))
        accepted.append(row)
    accepted.sort(
        key=lambda row: (
            int(row["human"]["contig_index"]),
            int(row["human_start0"]),
        )
    )
    counters["accepted_before_WCG_filter"] = len(accepted)
    return accepted, dict(counters)


def count_wcg_in_sequence(sequence: str, start_offset: int, width: int) -> int:
    """Count strand-specific, non-GCG WCG cytosines centred in a region."""
    sequence = sequence.upper()
    count = 0
    for index in range(start_offset, start_offset + width):
        if index <= 0 or index + 1 >= len(sequence):
            continue
        previous, base, following = sequence[index - 1 : index + 2]
        if base == "C" and previous in "AT" and following == "G":
            count += 1
        elif base == "G" and previous == "C" and following in "AT":
            count += 1
    return count


def count_reference_wcg(
    fasta: Path,
    intervals: list[dict[str, object]],
    samtools: str,
) -> list[int]:
    if not fasta.is_file() or fasta.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty FASTA: {fasta}")
    fai = Path(str(fasta) + ".fai")
    if not fai.is_file() or fai.stat().st_size == 0:
        raise FileNotFoundError(f"missing FASTA index: {fai}")
    with tempfile.TemporaryDirectory(prefix="snmar_wcg_") as temporary:
        regions = Path(temporary) / "regions.txt"
        requests: list[tuple[int, int, int]] = []
        with regions.open("wt", encoding="utf-8", newline="") as handle:
            for row in intervals:
                start0, end0 = int(row["start0"]), int(row["end0"])
                length = int(row["contig_length"])
                fetch_start = max(0, start0 - 1)
                fetch_end = min(length, end0 + 1)
                handle.write(f"{row['fasta_contig']}:{fetch_start + 1}-{fetch_end}\n")
                requests.append((start0 - fetch_start, end0 - start0, fetch_end - fetch_start))
        process = subprocess.Popen(
            [samtools, "faidx", "-r", str(regions), str(fasta)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        assert process.stdout is not None
        counts: list[int] = []
        sequence_parts: list[str] = []
        request_index = -1
        for raw in process.stdout:
            if raw.startswith(">"):
                if request_index >= 0:
                    offset, width, expected = requests[request_index]
                    sequence = "".join(sequence_parts)
                    if len(sequence) != expected:
                        raise ValueError("samtools faidx returned an unexpected sequence length")
                    counts.append(count_wcg_in_sequence(sequence, offset, width))
                request_index += 1
                sequence_parts = []
            else:
                sequence_parts.append(raw.strip())
        if request_index >= 0:
            offset, width, expected = requests[request_index]
            sequence = "".join(sequence_parts)
            if len(sequence) != expected:
                raise ValueError("samtools faidx returned an unexpected final sequence length")
            counts.append(count_wcg_in_sequence(sequence, offset, width))
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(f"samtools faidx failed ({return_code}): {stderr}")
    if len(counts) != len(intervals):
        raise ValueError(f"reference sequence count differs: {len(counts)} != {len(intervals)}")
    return counts


def fasta_contig_name(joint_contig: str, species: str) -> str:
    prefix = species + "_"
    if not joint_contig.startswith(prefix):
        raise ValueError(f"unexpected {species} joint contig: {joint_contig}")
    return joint_contig[len(prefix) :]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--chain", type=Path)
    parser.add_argument("--human-fasta", type=Path)
    parser.add_argument("--pig-fasta", type=Path)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--bin-size", type=int, default=100_000)
    parser.add_argument("--min-aligned-fraction", type=float, default=0.80)
    parser.add_argument("--min-span-ratio", type=float, default=0.70)
    parser.add_argument("--max-span-ratio", type=float, default=1.30)
    parser.add_argument("--min-reference-wcg", type=int, default=20)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.bin_size < 1
        or not 0 < args.min_aligned_fraction <= 1
        or not 0 < args.min_span_ratio <= args.max_span_ratio
        or args.min_reference_wcg < 1
    ):
        raise ValueError("invalid ortholog-bin parameter")
    project = args.project
    chain = args.chain or (
        project / "02_reference/human_pig_synteny/hg38.susScr11.rbest.chain.gz"
    )
    human_fasta = args.human_fasta or Path(
        "/share/org/YZWL/yzwl_wujx/Projects/FOODIE_project/reference/"
        "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    )
    pig_fasta = args.pig_fasta or Path(
        "/share/org/YZWL/yzwl_wujx/Projects/FOODIE_pig/reference/"
        "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    )
    contig_path = (
        project
        / "08_qc_summary/20_single_cytosine_inputs/"
        "dna_single_cytosine_contig_index.tsv"
    )
    complete_summary = (
        project
        / "08_qc_summary/23_dna_single_cytosine_complete/"
        "dna_single_cytosine_complete_summary.tsv"
    )
    complete = read_metrics(complete_summary)
    if (
        complete.get("status") != "validated"
        or int(complete.get("validated_samples", "0")) != 96
        or int(complete.get("source_h5_shards", "0")) != 768
    ):
        raise ValueError("complete single-cytosine input is not validated")
    output_root = args.output_root or (
        project / "10_downstream_analysis/10_dna_ortholog_100kb"
    )
    qc_root = args.qc_root or (
        project / "08_qc_summary/40_dna_ortholog_100kb"
    )
    index_path = output_root / "dna_ortholog_100kb_physical_features.tsv.gz"
    pair_path = output_root / "dna_ortholog_100kb_axis.tsv.gz"
    summary_path = qc_root / "dna_ortholog_100kb_index_summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and Path(existing.get("physical_feature_index", "")).is_file()
            and Path(existing.get("orthobin_axis", "")).is_file()
        ):
            print(f"[SKIP] validated ortholog 100 kb index exists: {summary_path}")
            return 0
        raise FileExistsError(f"inspect or rerun with --force: {summary_path}")
    contigs = load_contigs(contig_path)
    candidates, chain_metrics = collect_candidates(
        chain,
        build_alias_lookup(contigs, "human"),
        build_alias_lookup(contigs, "pig"),
        args.bin_size,
    )
    pairs, filter_metrics = choose_pairs(
        candidates,
        args.min_aligned_fraction,
        args.min_span_ratio,
        args.max_span_ratio,
    )
    if not pairs:
        raise ValueError("no reciprocal-best 100 kb pairs passed structural filters")

    human_intervals: list[dict[str, object]] = []
    pig_intervals: list[dict[str, object]] = []
    for row in pairs:
        human, pig = row["human"], row["pig"]
        human_intervals.append(
            {
                "fasta_contig": fasta_contig_name(str(human["contig"]), "human"),
                "start0": int(row["human_start0"]),
                "end0": int(row["human_end0"]),
                "contig_length": int(human["length"]),
            }
        )
        pig_intervals.append(
            {
                "fasta_contig": fasta_contig_name(str(pig["contig"]), "pig"),
                "start0": int(row["pig_start0"]),
                "end0": int(row["pig_end0"]),
                "contig_length": int(pig["length"]),
            }
        )
    human_wcg = count_reference_wcg(human_fasta, human_intervals, args.samtools)
    pig_wcg = count_reference_wcg(pig_fasta, pig_intervals, args.samtools)
    retained = [
        (row, human_count, pig_count)
        for row, human_count, pig_count in zip(pairs, human_wcg, pig_wcg)
        if human_count >= args.min_reference_wcg
        and pig_count >= args.min_reference_wcg
    ]
    if not retained:
        raise ValueError("no ortholog bins pass the bilateral reference-WCG filter")

    physical: list[dict[str, object]] = []
    axes: list[dict[str, object]] = []
    for orthobin_index, (row, human_count, pig_count) in enumerate(retained):
        human, pig = row["human"], row["pig"]
        pair_features: dict[str, int] = {}
        for species, contig, start_key, end_key, reference_count in (
            ("human", human, "human_start0", "human_end0", human_count),
            ("pig", pig, "pig_start0", "pig_end0", pig_count),
        ):
            feature_index = len(physical)
            pair_features[species] = feature_index
            physical.append(
                {
                    "feature_index": feature_index,
                    "feature_id": f"ORTHOBIN{orthobin_index + 1:06d}_{species}",
                    "orthobin_index": orthobin_index,
                    "orthobin_id": f"ORTHOBIN{orthobin_index + 1:06d}",
                    "species": species,
                    "contig_index": int(contig["contig_index"]),
                    "chrom": contig["contig"],
                    "start0": int(row[start_key]),
                    "end0": int(row[end_key]),
                    "width": int(row[end_key]) - int(row[start_key]),
                    "reference_WCG_sites_both_strands": reference_count,
                }
            )
        axes.append(
            {
                "orthobin_index": orthobin_index,
                "orthobin_id": f"ORTHOBIN{orthobin_index + 1:06d}",
                "human_feature_index": pair_features["human"],
                "human_chrom": human["contig"],
                "human_start0": row["human_start0"],
                "human_end0": row["human_end0"],
                "human_reference_WCG_sites": human_count,
                "pig_feature_index": pair_features["pig"],
                "pig_chrom": pig["contig"],
                "pig_start0": row["pig_start0"],
                "pig_end0": row["pig_end0"],
                "pig_reference_WCG_sites": pig_count,
                "aligned_bases": row["aligned_bases"],
                "human_aligned_fraction": row["aligned_fraction"],
                "pig_span_ratio": row["pig_span_ratio"],
                "chain_id": row["chain_id"],
                "chain_score": row["chain_score"],
                "chain_query_strand": row["q_strand"],
                "pairing_rule": "UCSC_hg38_susScr11_reciprocal_best_chain",
            }
        )
    atomic_tsv(index_path, list(physical[0]), physical, gzip_output=True)
    atomic_tsv(pair_path, list(axes[0]), axes, gzip_output=True)
    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("mode", "label_blind_shared_human_pig_orthobin_axis"),
        ("bin_size", args.bin_size),
        ("min_aligned_fraction", args.min_aligned_fraction),
        ("min_span_ratio", args.min_span_ratio),
        ("max_span_ratio", args.max_span_ratio),
        ("min_reference_WCG_sites_each_species", args.min_reference_wcg),
        ("orthobin_axes", len(axes)),
        ("physical_features", len(physical)),
        ("pairs_removed_by_reference_WCG_filter", len(pairs) - len(retained)),
        ("ambiguous_context_policy", "GCG_excluded_from_reference_WCG_count"),
        ("reference_WCG_definition", "plus_[AT]CG_and_minus_CG[AT]"),
        ("chain", chain),
        ("chain_source", CHAIN_SOURCE),
        ("human_fasta", human_fasta),
        ("pig_fasta", pig_fasta),
        ("physical_feature_index", index_path),
        ("orthobin_axis", pair_path),
        ("source_contig_index", contig_path),
        ("source_complete_summary", complete_summary),
    ]
    for key, value in sorted({**chain_metrics, **filter_metrics}.items()):
        summary_rows.append((key, value))
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] ortholog 100 kb index: axes={len(axes)} "
        f"physical_features={len(physical)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
