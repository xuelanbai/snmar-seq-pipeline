#!/usr/bin/env python3
"""Step 33: validate full mapping tables and select downstream Probe inputs."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
from statistics import median


PROJECT_DEFAULT = Path(
    "/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"missing/empty TSV: {path}")
    with path.open("rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"TSV has no data rows: {path}")
    return rows


def as_float(row: dict[str, str], key: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"invalid numeric field {key!r}: {row}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite field {key!r}: {value}")
    return value


def add_choice(
    choices: dict[str, tuple[dict[str, str], str]],
    row: dict[str, str],
    key: str,
    reason: str,
) -> None:
    identifier = row[key]
    if identifier not in choices:
        choices[identifier] = (row, reason)


def closest(
    rows: list[dict[str, str]], key: str, target: float
) -> dict[str, str]:
    return min(rows, key=lambda row: abs(as_float(row, key) - target))


def select_dna(
    rows: list[dict[str, str]], count: int
) -> list[tuple[dict[str, str], str]]:
    choices: dict[str, tuple[dict[str, str], str]] = {}
    key = "run_id"
    add_choice(
        choices,
        min(rows, key=lambda row: as_float(row, "primary_mapping_percent")),
        key,
        "lowest_primary_mapping_percent",
    )
    gm = [row for row in rows if row.get("group") == "GM12878"]
    mix = [row for row in rows if row.get("group") == "cellMix"]
    if not gm or not mix:
        raise ValueError("DNA table must contain GM12878 and cellMix rows")
    gm_target = median(as_float(row, "expected_read_pairs") for row in gm)
    mix_target = median(
        as_float(row, "pig_percent_of_primary_mapped") for row in mix
    )
    add_choice(
        choices,
        closest(gm, "expected_read_pairs", gm_target),
        key,
        "GM12878_median_read_depth",
    )
    add_choice(
        choices,
        closest(mix, "pig_percent_of_primary_mapped", mix_target),
        key,
        "cellMix_median_pig_fraction",
    )
    add_choice(
        choices,
        max(rows, key=lambda row: as_float(row, "bam_bytes")),
        key,
        "largest_BAM",
    )
    ranked = sorted(
        rows, key=lambda row: as_float(row, "primary_mapping_percent")
    )
    quantile_rows = [
        (
            ranked[round(q * (len(ranked) - 1))],
            f"mapping_quantile_{q:.2f}",
        )
        for q in (index / max(1, count - 1) for index in range(count))
    ]
    for row, reason in quantile_rows:
        add_choice(choices, row, key, reason)
    for rank, row in enumerate(ranked):
        add_choice(choices, row, key, f"mapping_rank_{rank + 1}")
        if len(choices) >= count:
            break
    return list(choices.values())[:count]


def select_rna(
    rows: list[dict[str, str]], count: int
) -> list[tuple[dict[str, str], str]]:
    choices: dict[str, tuple[dict[str, str], str]] = {}
    key = "pair_id"
    criteria = (
        ("lowest_primary_mapping_percent", min, "primary_mapping_percent"),
        ("highest_unmapped_too_short", max, "unmapped_too_short_percent"),
        ("largest_input", max, "star_input_reads"),
        ("highest_multimapping", max, "multi_mapped_percent"),
        ("lowest_pig_fraction", min, "pig_percent_of_primary_mapped"),
        ("highest_pig_fraction", max, "pig_percent_of_primary_mapped"),
    )
    for reason, fn, column in criteria:
        add_choice(
            choices,
            fn(rows, key=lambda row, c=column: as_float(row, c)),
            key,
            reason,
        )
    too_short_ranked = sorted(
        rows,
        key=lambda row: as_float(row, "unmapped_too_short_percent"),
        reverse=True,
    )
    if len(too_short_ranked) >= 2:
        add_choice(
            choices,
            too_short_ranked[1],
            key,
            "second_highest_unmapped_too_short",
        )
    ranked = sorted(
        rows, key=lambda row: as_float(row, "primary_mapping_percent")
    )
    quantile_rows = [
        (
            ranked[round(q * (len(ranked) - 1))],
            f"mapping_quantile_{q:.2f}",
        )
        for q in (index / max(1, count - 1) for index in range(count))
    ]
    for row, reason in quantile_rows:
        add_choice(choices, row, key, reason)
    for rank, row in enumerate(ranked):
        add_choice(choices, row, key, f"mapping_rank_{rank + 1}")
        if len(choices) >= count:
            break
    return list(choices.values())[:count]


def atomic_tsv(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
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
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, default=PROJECT_DEFAULT)
    parser.add_argument("--dna-per-sample", type=Path)
    parser.add_argument("--rna-per-file", type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--dna-probe-samples", type=int, default=4)
    parser.add_argument("--rna-probe-files", type=int, default=8)
    args = parser.parse_args()

    project = args.project
    dna_path = args.dna_per_sample or (
        project
        / "08_qc_summary/12_dna_bwa_full/dna_bwa_full_per_sample.tsv"
    )
    rna_path = args.rna_per_file or (
        project
        / "08_qc_summary/09_rna_star_full/rna_star_full_per_file.tsv"
    )
    out_dir = args.out_dir or (
        project / "08_qc_summary/13_postmapping_probe"
    )
    if args.dna_probe_samples < 2 or args.rna_probe_files < 3:
        parser.error("Probe sizes are too small")

    dna_rows = read_tsv(dna_path)
    rna_rows = read_tsv(rna_path)
    if len(dna_rows) != 96:
        raise ValueError(f"expected 96 DNA rows, found {len(dna_rows)}")
    if len(rna_rows) != 96:
        raise ValueError(f"expected 96 RNA rows, found {len(rna_rows)}")
    if any(row.get("status") != "validated" for row in dna_rows):
        raise ValueError("all DNA rows must have status=validated")
    if any(as_float(row, "primary_mapping_percent") <= 0 for row in rna_rows):
        raise ValueError("RNA primary mapping percentages must be positive")
    if len({row["run_id"] for row in dna_rows}) != 96:
        raise ValueError("DNA run_id values are not unique")
    if len({row["pair_id"] for row in rna_rows}) != 96:
        raise ValueError("RNA pair_id values are not unique")

    dna_selected = select_dna(dna_rows, args.dna_probe_samples)
    rna_selected = select_rna(rna_rows, args.rna_probe_files)
    if len(dna_selected) != args.dna_probe_samples:
        raise ValueError("could not select requested number of DNA samples")
    if len(rna_selected) != args.rna_probe_files:
        raise ValueError("could not select requested number of RNA files")

    dna_manifest: list[dict[str, object]] = []
    for index, (row, reason) in enumerate(dna_selected, 1):
        dna_manifest.append(
            {
                "probe_index": index,
                "run_id": row["run_id"],
                "sample": row["sample"],
                "group": row["group"],
                "selection_reason": reason,
                "input_primary_records": row["primary_records"],
                "input_primary_mapped": row["primary_mapped"],
                "input_mapping_percent": row["primary_mapping_percent"],
                "input_bam_bytes": row["bam_bytes"],
                "bam": row["bam"],
            }
        )
    rna_manifest: list[dict[str, object]] = []
    for index, (row, reason) in enumerate(rna_selected, 1):
        rna_manifest.append(
            {
                "probe_index": index,
                "pair_id": row["pair_id"],
                "sample": row["sample"],
                "selection_reason": reason,
                "star_input_reads": row["star_input_reads"],
                "primary_mapped_reads": row["primary_mapped_reads"],
                "primary_mapping_percent": row["primary_mapping_percent"],
                "unmapped_too_short_percent": row[
                    "unmapped_too_short_percent"
                ],
                "bam": row["bam_path"],
            }
        )

    atomic_tsv(
        out_dir / "dna_cb_markdup_probe_manifest.tsv",
        list(dna_manifest[0]),
        dna_manifest,
    )
    atomic_tsv(
        out_dir / "rna_genefull_probe_manifest.tsv",
        list(rna_manifest[0]),
        rna_manifest,
    )

    review = [
        {"modality": "DNA", "metric": "validated_samples", "value": 96},
        {
            "modality": "DNA",
            "metric": "primary_mapping_percent_weighted",
            "value": (
                100
                * sum(as_float(row, "primary_mapped") for row in dna_rows)
                / sum(as_float(row, "primary_records") for row in dna_rows)
            ),
        },
        {
            "modality": "DNA",
            "metric": "lowest_sample_primary_mapping_percent",
            "value": min(
                as_float(row, "primary_mapping_percent") for row in dna_rows
            ),
        },
        {"modality": "RNA", "metric": "completed_files", "value": 96},
        {
            "modality": "RNA",
            "metric": "primary_mapping_percent_weighted",
            "value": (
                100
                * sum(
                    as_float(row, "primary_mapped_reads") for row in rna_rows
                )
                / sum(as_float(row, "star_input_reads") for row in rna_rows)
            ),
        },
        {
            "modality": "RNA",
            "metric": "lowest_file_primary_mapping_percent",
            "value": min(
                as_float(row, "primary_mapping_percent") for row in rna_rows
            ),
        },
    ]
    atomic_tsv(
        out_dir / "mapping_full_review.tsv",
        ["modality", "metric", "value"],
        review,
    )
    print(
        f"[OK] prepared {len(dna_manifest)} DNA and "
        f"{len(rna_manifest)} RNA Probe inputs in {out_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
