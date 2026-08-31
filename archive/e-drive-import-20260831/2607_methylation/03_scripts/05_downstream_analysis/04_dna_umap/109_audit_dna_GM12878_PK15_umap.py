#!/usr/bin/env python3
"""Step 109: audit D50 DNA UMAP against confirmed GM12878/PK15 identities.

The cell-line labels generated here are derived from human-versus-pig DNA
reference coverage.  They are useful for internal consistency checks but are
not independent validation of methylation-based cell-line separation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


PIPELINE_SIGNATURE = "snmar_dna_GM12878_PK15_umap_audit_20260806"
RUN_SIGNATURE = "snmar_dna_100kb_umap_20260806"
CELL_LINES = ("GM12878", "PK15", "ambiguous")


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
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    suffix = ".tmp.gz" if path.suffix == ".gz" else ".tmp"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=suffix, dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with open_text(temporary, "w") as handle:
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
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def atomic_metrics(path: Path, rows: list[tuple[str, object]]) -> None:
    atomic_tsv(
        path,
        ["metric", "value"],
        ({"metric": metric, "value": value} for metric, value in rows),
    )


def classify_cell_line(
    human_calls: int,
    pig_calls: int,
    minimum_major_fraction: float,
) -> tuple[str, float, float, float]:
    total = human_calls + pig_calls
    if total <= 0:
        return "ambiguous", 0.0, 0.0, 0.0
    human_fraction = human_calls / total
    pig_fraction = pig_calls / total
    major_fraction = max(human_fraction, pig_fraction)
    if major_fraction < minimum_major_fraction or human_calls == pig_calls:
        label = "ambiguous"
    elif human_calls > pig_calls:
        label = "GM12878"
    else:
        label = "PK15"
    return label, human_fraction, pig_fraction, major_fraction


def well_design(well_set: str) -> str:
    if well_set == "GM12878":
        return "GM12878_only_wells"
    if well_set == "cellMix":
        return "GM12878_PK15_mix_wells"
    raise ValueError(f"unexpected well_set: {well_set}")


def same_label_baseline(labels: list[str]) -> float:
    counts = Counter(labels)
    total = len(labels)
    if total < 2:
        return 0.0
    return sum(value * (value - 1) for value in counts.values()) / (
        total * (total - 1)
    )


def embedding_metrics(
    coordinates,
    labels: list[str],
    neighbors: int,
    seed: int,
    np,
    NearestNeighbors,
    silhouette_score,
) -> tuple[float, float]:
    if len(coordinates) < 3 or len(set(labels)) < 2:
        return 0.0, 0.0
    k = min(neighbors, len(coordinates) - 1)
    label_array = np.asarray(labels)
    indices = NearestNeighbors(n_neighbors=k + 1).fit(
        coordinates
    ).kneighbors(return_distance=False)[:, 1:]
    agreement = float(
        np.mean(label_array[indices] == label_array[:, None])
    )
    sample_size = min(5_000, len(coordinates))
    silhouette = float(
        silhouette_score(
            coordinates,
            label_array,
            sample_size=(
                sample_size if sample_size < len(coordinates) else None
            ),
            random_state=seed,
        )
    )
    return agreement, silhouette


def scatter_categories(ax, coordinates, labels, palette, order, size=5.0):
    for label in order:
        mask = [value == label for value in labels]
        if not any(mask):
            continue
        selected = coordinates[mask]
        ax.scatter(
            selected[:, 0],
            selected[:, 1],
            s=size,
            c=palette[label],
            label=label,
            alpha=0.72,
            linewidths=0,
            rasterized=True,
        )
    ax.legend(frameon=False, markerscale=2.2, fontsize=8)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("/share/org/YZWL/yzwl_wujx/Projects/2607_methylation"),
    )
    parser.add_argument("--threshold", default="D_total_50k")
    parser.add_argument("--minimum-major-fraction", type=float, default=0.90)
    parser.add_argument("--neighbors", type=int, default=30)
    parser.add_argument("--seed", type=int, default=2607)
    parser.add_argument("--input-summary", type=Path)
    parser.add_argument("--coordinates", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--qc-root", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 0.5 < args.minimum_major_fraction <= 1.0:
        raise ValueError("minimum major fraction must be in (0.5, 1]")
    if args.neighbors < 2:
        raise ValueError("neighbors must be at least 2")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from sklearn.metrics import silhouette_score
        from sklearn.neighbors import NearestNeighbors
    except ImportError as error:
        raise RuntimeError(
            "Step 109 requires numpy, scikit-learn and matplotlib"
        ) from error

    source_summary = args.input_summary or (
        args.project
        / "08_qc_summary/32_dna_100kb_umap"
        / args.threshold
        / "dna_100kb_umap_summary.tsv"
    )
    source = read_metrics(source_summary)
    if (
        source.get("status") != "validated"
        or source.get("pipeline_signature") != RUN_SIGNATURE
        or source.get("threshold_id") != args.threshold
    ):
        raise ValueError(f"invalid DNA UMAP summary: {source_summary}")
    coordinates_path = args.coordinates or Path(source["coordinates"])
    rows = read_tsv(coordinates_path)

    output_root = args.output_root or (
        args.project
        / "10_downstream_analysis/03_dna_cell_line_audit"
        / args.threshold
    )
    qc_root = args.qc_root or (
        args.project
        / "08_qc_summary/33_dna_cell_line_audit"
        / args.threshold
    )
    summary_path = qc_root / "dna_GM12878_PK15_umap_audit_summary.tsv"
    if summary_path.is_file() and not args.force:
        existing = read_metrics(summary_path)
        if (
            existing.get("status") == "validated"
            and existing.get("pipeline_signature") == PIPELINE_SIGNATURE
            and existing.get("threshold_id") == args.threshold
            and float(existing.get("minimum_major_fraction", "0"))
            == args.minimum_major_fraction
        ):
            print(f"[SKIP] validated audit exists: {summary_path}")
            return 0
        raise FileExistsError(
            f"non-current audit exists; inspect or use --force: {summary_path}"
        )

    required = {
        "cell_index",
        "well_index",
        "physical_well",
        "well_set",
        "cell_id",
        "CB",
        "human_signal_covered_calls",
        "pig_signal_covered_calls",
        "WCG_UMAP1",
        "WCG_UMAP2",
        "GCH_UMAP1",
        "GCH_UMAP2",
        "joint_UMAP1",
        "joint_UMAP2",
    }
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"coordinate columns missing: {sorted(missing)}")

    annotated: list[dict[str, object]] = []
    for row in rows:
        human_calls = int(row["human_signal_covered_calls"])
        pig_calls = int(row["pig_signal_covered_calls"])
        label, human_fraction, pig_fraction, major_fraction = (
            classify_cell_line(
                human_calls,
                pig_calls,
                args.minimum_major_fraction,
            )
        )
        annotated.append(
            {
                **row,
                "well_design": well_design(row["well_set"]),
                "DNA_cell_line_mapping_derived": label,
                "human_reference_fraction": human_fraction,
                "pig_reference_fraction": pig_fraction,
                "major_reference_fraction_recomputed": major_fraction,
                "cell_line_assignment_status": (
                    "high_confidence" if label != "ambiguous" else "ambiguous"
                ),
            }
        )

    if len({row["cell_id"] for row in annotated}) != len(annotated):
        raise ValueError("duplicate cell_id in UMAP coordinates")
    counts = Counter(
        str(row["DNA_cell_line_mapping_derived"]) for row in annotated
    )
    if sum(counts.values()) != len(annotated):
        raise AssertionError("cell-line counts do not reconcile")

    label_path = (
        output_root
        / f"dna_100kb_{args.threshold}_GM12878_PK15_labels.tsv.gz"
    )
    atomic_tsv(label_path, list(annotated[0]), annotated)

    count_rows: list[dict[str, object]] = []
    for scope, scope_rows in (
        ("all_cells", annotated),
        (
            "GM12878_only_wells",
            [row for row in annotated if row["well_set"] == "GM12878"],
        ),
        (
            "GM12878_PK15_mix_wells",
            [row for row in annotated if row["well_set"] == "cellMix"],
        ),
    ):
        scoped = Counter(
            str(row["DNA_cell_line_mapping_derived"]) for row in scope_rows
        )
        count_rows.append(
            {
                "scope": scope,
                "total_cells": len(scope_rows),
                "GM12878_cells": scoped["GM12878"],
                "PK15_cells": scoped["PK15"],
                "ambiguous_cells": scoped["ambiguous"],
                "ambiguous_fraction": (
                    scoped["ambiguous"] / len(scope_rows)
                    if scope_rows
                    else 0.0
                ),
            }
        )
    count_path = qc_root / "dna_GM12878_PK15_cell_line_counts.tsv"
    atomic_tsv(count_path, list(count_rows[0]), count_rows)

    per_well_rows: list[dict[str, object]] = []
    grouped: dict[tuple[int, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in annotated:
        key = (
            int(row["well_index"]),
            str(row["physical_well"]),
            str(row["well_design"]),
        )
        grouped[key].append(row)
    for (well_index, physical_well, design), values in sorted(grouped.items()):
        scoped = Counter(
            str(row["DNA_cell_line_mapping_derived"]) for row in values
        )
        per_well_rows.append(
            {
                "well_index": well_index,
                "physical_well": physical_well,
                "well_design": design,
                "total_cells": len(values),
                "GM12878_cells": scoped["GM12878"],
                "PK15_cells": scoped["PK15"],
                "ambiguous_cells": scoped["ambiguous"],
                "ambiguous_fraction": scoped["ambiguous"] / len(values),
            }
        )
    per_well_path = qc_root / "dna_GM12878_PK15_cell_line_per_well.tsv"
    atomic_tsv(per_well_path, list(per_well_rows[0]), per_well_rows)

    embeddings = {}
    for name in ("WCG", "GCH", "joint"):
        embeddings[name] = np.asarray(
            [
                [float(row[f"{name}_UMAP1"]), float(row[f"{name}_UMAP2"])]
                for row in annotated
            ],
            dtype=np.float32,
        )
    cell_labels = [
        str(row["DNA_cell_line_mapping_derived"]) for row in annotated
    ]
    confident_indices = np.asarray(
        [label != "ambiguous" for label in cell_labels]
    )
    confident_labels = [
        label for label in cell_labels if label != "ambiguous"
    ]
    gm_indices = np.asarray([label == "GM12878" for label in cell_labels])
    gm_sources = [
        str(row["well_design"])
        for row, keep in zip(annotated, gm_indices)
        if keep
    ]
    metric_rows: list[dict[str, object]] = []
    for name in ("WCG", "GCH", "joint"):
        agreement, silhouette = embedding_metrics(
            embeddings[name][confident_indices],
            confident_labels,
            args.neighbors,
            args.seed,
            np,
            NearestNeighbors,
            silhouette_score,
        )
        gm_agreement, gm_silhouette = embedding_metrics(
            embeddings[name][gm_indices],
            gm_sources,
            args.neighbors,
            args.seed,
            np,
            NearestNeighbors,
            silhouette_score,
        )
        baseline = same_label_baseline(gm_sources)
        metric_rows.append(
            {
                "embedding": name,
                "confident_cells": len(confident_labels),
                "GM12878_cells": confident_labels.count("GM12878"),
                "PK15_cells": confident_labels.count("PK15"),
                "cell_line_knn30_agreement": agreement,
                "cell_line_silhouette": silhouette,
                "GM12878_only_vs_mix_cells": len(gm_sources),
                "GM12878_well_set_knn30_agreement": gm_agreement,
                "GM12878_well_set_random_baseline": baseline,
                "GM12878_well_set_knn_excess": gm_agreement - baseline,
                "GM12878_well_set_silhouette": gm_silhouette,
            }
        )
    metric_path = qc_root / "dna_GM12878_PK15_embedding_metrics.tsv"
    atomic_tsv(metric_path, list(metric_rows[0]), metric_rows)

    label_palette = {
        "GM12878": "#377EB8",
        "PK15": "#E41A1C",
        "ambiguous": "#9E9E9E",
    }
    design_palette = {
        "GM12878_only_wells": "#1B9E77",
        "GM12878_PK15_mix_wells": "#BDBDBD",
    }
    designs = [str(row["well_design"]) for row in annotated]
    confidence = np.asarray(
        [float(row["major_reference_fraction_recomputed"]) for row in annotated]
    )
    figure, axes = plt.subplots(3, 3, figsize=(18, 16), constrained_layout=True)
    confidence_plot = None
    for row_index, name in enumerate(("WCG", "GCH", "joint")):
        scatter_categories(
            axes[row_index, 0],
            embeddings[name],
            cell_labels,
            label_palette,
            CELL_LINES,
        )
        axes[row_index, 0].set_title(
            f"{name}: DNA mapping-derived cell line"
        )
        scatter_categories(
            axes[row_index, 1],
            embeddings[name],
            designs,
            design_palette,
            ("GM12878_only_wells", "GM12878_PK15_mix_wells"),
        )
        axes[row_index, 1].set_title(f"{name}: physical-well design")
        confidence_plot = axes[row_index, 2].scatter(
            embeddings[name][:, 0],
            embeddings[name][:, 1],
            c=confidence,
            cmap="viridis",
            vmin=0.5,
            vmax=1.0,
            s=5,
            alpha=0.75,
            linewidths=0,
            rasterized=True,
        )
        axes[row_index, 2].set_title(f"{name}: mapping confidence")
        axes[row_index, 2].set_xlabel("UMAP1")
        axes[row_index, 2].set_ylabel("UMAP2")
    if confidence_plot is not None:
        figure.colorbar(
            confidence_plot,
            ax=axes[:, 2].tolist(),
            label="major human/pig reference fraction",
            shrink=0.78,
        )
    figure.suptitle(
        "snMAR DNA 100 kb GM12878/PK15 audit "
        f"({args.threshold}; mapping-derived labels, not independent validation)",
        fontsize=15,
    )
    png_path = output_root / "dna_GM12878_PK15_umap_audit.png"
    pdf_path = output_root / "dna_GM12878_PK15_umap_audit.pdf"
    output_root.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_path, dpi=220)
    figure.savefig(pdf_path)
    plt.close(figure)

    summary_rows: list[tuple[str, object]] = [
        ("status", "validated"),
        ("pipeline_signature", PIPELINE_SIGNATURE),
        ("threshold_id", args.threshold),
        ("confirmed_cell_lines", "GM12878,PK15"),
        ("minimum_major_fraction", args.minimum_major_fraction),
        ("total_cells", len(annotated)),
        ("GM12878_cells_mapping_derived", counts["GM12878"]),
        ("PK15_cells_mapping_derived", counts["PK15"]),
        ("ambiguous_cells", counts["ambiguous"]),
        ("ambiguous_fraction", counts["ambiguous"] / len(annotated)),
        ("interpretation_scope", "DNA_internal_consistency_only"),
        ("independent_validation", "false"),
        (
            "cell_line_assignment",
            "human_pig_reference_coverage_mapping_derived",
        ),
        (
            "claim_boundary",
            "does_not_prove_independent_methylation_cell_line_separation",
        ),
        ("source_umap_summary", source_summary),
        ("source_coordinates", coordinates_path),
        ("annotated_coordinates", label_path),
        ("cell_line_counts", count_path),
        ("per_well_counts", per_well_path),
        ("embedding_metrics", metric_path),
        ("overview_png", png_path),
        ("overview_pdf", pdf_path),
    ]
    atomic_metrics(summary_path, summary_rows)
    print(
        f"[OK] GM12878/PK15 DNA UMAP audit: cells={len(annotated)} "
        f"ambiguous={counts['ambiguous']} output={output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
