#!/usr/bin/env python3
"""Local tests for direct-full downstream code; no real BAM is required."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POST = ROOT / "03_scripts/04_postmapping"


def script(filename: str) -> Path:
    matches = list(POST.rglob(filename))
    if len(matches) != 1:
        raise RuntimeError(f"expected one script named {filename}: {matches}")
    return matches[0]


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, script(filename))
    if spec is None or spec.loader is None:
        raise RuntimeError(filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


prepare = load("downstream_full_prepare", "53_prepare_downstream_full.py")
taps = load("taps_full_shard", "54_call_dna_taps_cpg_full_shard.py")
merge = load("taps_full_merge", "55_merge_dna_taps_cpg_full_sample.py")
rna_filter = load(
    "rna_filtered_full", "60_filter_rna_reverse_matrix_full.py"
)


class FullPreparationTests(unittest.TestCase):
    def test_balances_all_contigs_once(self) -> None:
        contigs = [
            {
                "reference_order": index,
                "contig": f"human_{index}",
                "length": length,
                "species": "human",
            }
            for index, length in enumerate((100, 90, 80, 70), 1)
        ]
        shards = prepare.balance_contigs(contigs, 2)
        names = [row["contig"] for shard in shards for row in shard]
        self.assertEqual(sorted(names), sorted(row["contig"] for row in contigs))
        self.assertEqual(
            [sum(int(row["length"]) for row in shard) for shard in shards],
            [170, 170],
        )


class FullTapsTests(unittest.TestCase):
    def test_cigar_to_aligned_bases(self) -> None:
        observed = list(
            taps.aligned_bases(
                101,
                "2S3M1I2M1D1M",
                "AACGTTACG",
                "IIIIIIIII",
            )
        )
        self.assertEqual(
            [(position, base) for position, base, _quality in observed],
            [
                (100, "C"),
                (101, "G"),
                (102, "T"),
                (103, "A"),
                (104, "C"),
                (106, "G"),
            ],
        )

    def test_directional_calling(self) -> None:
        reference = "ACG"
        self.assertEqual(taps.taps_call(False, reference, 1, "T"), (1, 1))
        self.assertEqual(taps.taps_call(False, reference, 1, "C"), (1, 0))
        self.assertIsNone(taps.taps_call(False, reference, 2, "G"))
        self.assertEqual(taps.taps_call(True, reference, 2, "A"), (1, 1))
        self.assertEqual(taps.taps_call(True, reference, 2, "G"), (1, 0))
        self.assertIsNone(taps.taps_call(True, reference, 1, "C"))

    def test_coordinate_flush_keeps_reverse_boundary(self) -> None:
        key = ("sample|CB", "qname", 100)
        pending = {key: (30, 1)}
        heap = [(100, 1, key)]
        consumed = []
        conflicts = taps.flush_pending_before(
            100,
            heap,
            pending,
            lambda cell, site, modified: consumed.append(
                (cell, site, modified)
            ),
        )
        self.assertEqual(conflicts, 0)
        self.assertEqual(consumed, [])
        taps.flush_pending_before(
            101,
            heap,
            pending,
            lambda cell, site, modified: consumed.append(
                (cell, site, modified)
            ),
        )
        self.assertEqual(consumed, [("sample|CB", 100, 1)])

    def test_equal_quality_overlap_conflict_is_discarded(self) -> None:
        key = ("sample|CB", "qname", 100)
        pending = {}
        self.assertEqual(
            taps.update_template_call(pending, key, 30, 0), "new"
        )
        self.assertEqual(
            taps.update_template_call(pending, key, 30, 1), "conflict"
        )
        heap = [(100, 1, key)]
        consumed = []
        conflicts = taps.flush_pending_before(
            None,
            heap,
            pending,
            lambda cell, site, modified: consumed.append(
                (cell, site, modified)
            ),
        )
        self.assertEqual(conflicts, 1)
        self.assertEqual(consumed, [])


class FullDnaMergeTests(unittest.TestCase):
    def test_two_shard_merge_reconciles(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sample = "Sample-1"
            bam = root / "input.bam"
            bam.write_bytes(b"BAM")
            shard_root = root / sample
            output_root = root / "output"
            bin_fields = [
                "sample",
                "cell_id",
                "CB",
                "species",
                "chrom",
                "start",
                "end",
                "modified_template_cpg_calls",
                "covered_template_cpg_calls",
                "modified_rate",
            ]
            for shard, chrom, covered, modified in (
                (1, "human_1", 2, 1),
                (2, "pig_1", 3, 2),
            ):
                shard_dir = shard_root / f"shard_{shard:02d}"
                shard_dir.mkdir(parents=True)
                prefix = f"{sample}.shard_{shard:02d}"
                (shard_dir / f"{prefix}.summary.tsv").write_text(
                    "metric\tvalue\n"
                    "status\tvalidated\n"
                    f"sample\t{sample}\n"
                    f"shard\t{shard}\n"
                    f"input_bam\t{bam}\n"
                    "expected_contigs\t1\n"
                    "completed_contigs\t1\n"
                    "reference_bases\t100\n"
                    "bin_size\t100000\n"
                    "call_model\tforward_ref_C_CT_reverse_ref_G_GA\n"
                    "overlap_policy\tone_call_per_sample_CB_QNAME_CpG_highest_baseQ\n"
                    "missing_required_CB_or_XD\t0\n"
                    "nonzero_XD\t0\n"
                    f"covered_template_cpg_calls\t{covered}\n"
                    f"modified_template_cpg_calls\t{modified}\n"
                    f"human_covered_calls\t{covered if shard == 1 else 0}\n"
                    f"human_modified_calls\t{modified if shard == 1 else 0}\n"
                    f"pig_covered_calls\t{covered if shard == 2 else 0}\n"
                    f"pig_modified_calls\t{modified if shard == 2 else 0}\n",
                    encoding="utf-8",
                )
                (shard_dir / f"{prefix}.per_cell.tsv").write_text(
                    "sample\tcell_id\tCB\teligible_alignment_records\t"
                    "covered_template_cpg_calls\t"
                    "modified_template_cpg_calls\tmodified_rate\t"
                    "bins_with_calls\n"
                    f"{sample}\t{sample}|CB\tCB\t4\t{covered}\t"
                    f"{modified}\t{modified / covered}\t1\n",
                    encoding="utf-8",
                )
                with gzip.open(
                    shard_dir / f"{prefix}.bin_counts.tsv.gz",
                    "wt",
                    encoding="utf-8",
                    newline="",
                ) as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=bin_fields,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerow(
                        {
                            "sample": sample,
                            "cell_id": f"{sample}|CB",
                            "CB": "CB",
                            "species": "human" if shard == 1 else "pig",
                            "chrom": chrom,
                            "start": 0,
                            "end": 100000,
                            "modified_template_cpg_calls": modified,
                            "covered_template_cpg_calls": covered,
                            "modified_rate": modified / covered,
                        }
                    )
            old_argv = sys.argv
            try:
                sys.argv = [
                    "merge",
                    "--sample",
                    sample,
                    "--bam",
                    str(bam),
                    "--shard-root",
                    str(shard_root),
                    "--output-root",
                    str(output_root),
                    "--expected-shards",
                    "2",
                ]
                self.assertEqual(merge.main(), 0)
            finally:
                sys.argv = old_argv
            summary_path = (
                output_root / sample / f"{sample}.taps_cpg_full.summary.tsv"
            )
            with summary_path.open(
                "rt", encoding="utf-8", newline=""
            ) as handle:
                metrics = {
                    row["metric"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(metrics["covered_template_cpg_calls"], "5")
            self.assertEqual(metrics["modified_template_cpg_calls"], "3")


class FullRnaFilterTests(unittest.TestCase):
    def test_filters_complete_raw_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw"
            raw.mkdir()
            summary = root / "summary.tsv"
            cell_qc = root / "cell_qc.tsv.gz"
            output = root / "filtered"
            qc = root / "qc"
            summary.write_text(
                "metric\tvalue\n"
                "status\tvalidated\n"
                "cell_filtering\tnot_applied_raw_matrix\n"
                "raw_matrix_cells\t3\n"
                "raw_matrix_features\t2\n"
                "raw_matrix_nonzero_entries\t4\n"
                "raw_matrix_umi_sum\t9\n",
                encoding="utf-8",
            )
            with gzip.open(
                raw / "barcodes.tsv.gz", "wt", encoding="utf-8"
            ) as out:
                out.write("S1|A\nS1|B\nS2|C\n")
            with gzip.open(
                raw / "features.tsv.gz", "wt", encoding="utf-8"
            ) as out:
                out.write("G1\tG1\tGene Expression\nG2\tG2\tGene Expression\n")
            with gzip.open(
                raw / "feature_metadata.tsv.gz", "wt", encoding="utf-8"
            ) as out:
                out.write("gene_id\tgene_name\tspecies\nG1\tG1\thuman\nG2\tG2\tpig\n")
            with gzip.open(
                raw / "matrix.mtx.gz", "wt", encoding="utf-8"
            ) as out:
                out.write(
                    "%%MatrixMarket matrix coordinate integer general\n"
                    "% raw\n"
                    "2 3 4\n"
                    "1 1 1\n"
                    "2 1 2\n"
                    "1 2 1\n"
                    "2 3 5\n"
                )
            fields = [
                "pair_id",
                "sample",
                "cell_id",
                "CB",
                "CN",
                "tag_eligible_unique_primary",
                "gene_assigned_reads",
                "gene_assignment_fraction",
                "exact_umis",
                "genes_detected",
                "human_assigned_reads",
                "pig_assigned_reads",
                "human_exact_umis",
                "pig_exact_umis",
                "in_raw_matrix",
            ]
            with gzip.open(
                cell_qc, "wt", encoding="utf-8", newline=""
            ) as out:
                writer = csv.DictWriter(
                    out,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for cell_id, sample, umis, genes in (
                    ("S1|A", "S1", 3, 2),
                    ("S1|B", "S1", 1, 1),
                    ("S2|C", "S2", 5, 1),
                ):
                    writer.writerow(
                        {
                            "pair_id": sample,
                            "sample": sample,
                            "cell_id": cell_id,
                            "CB": cell_id.split("|")[1],
                            "CN": "CN",
                            "tag_eligible_unique_primary": umis,
                            "gene_assigned_reads": umis,
                            "gene_assignment_fraction": 1,
                            "exact_umis": umis,
                            "genes_detected": genes,
                            "human_assigned_reads": umis,
                            "pig_assigned_reads": 0,
                            "human_exact_umis": umis,
                            "pig_exact_umis": 0,
                            "in_raw_matrix": 1,
                        }
                    )
            old_argv = sys.argv
            try:
                sys.argv = [
                    "rna_filter",
                    "--raw-matrix-root",
                    str(raw),
                    "--cell-qc",
                    str(cell_qc),
                    "--summary",
                    str(summary),
                    "--output-root",
                    str(output),
                    "--qc-root",
                    str(qc),
                    "--min-exact-umis",
                    "2",
                    "--min-genes",
                    "2",
                ]
                self.assertEqual(rna_filter.main(), 0)
            finally:
                sys.argv = old_argv
            with gzip.open(
                output / "matrix.mtx.gz", "rt", encoding="utf-8"
            ) as handle:
                text = handle.read()
            self.assertIn("2 1 2", text)
            self.assertIn("1 1 1", text)
            self.assertIn("2 1 2", text)
            with (
                qc / "rna_filtered_matrix_full_summary.tsv"
            ).open("rt", encoding="utf-8", newline="") as handle:
                metrics = {
                    row["metric"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(metrics["filtered_cells"], "1")
            self.assertEqual(metrics["filtered_exact_umis"], "3")
            self.assertEqual(metrics["raw_exact_umis"], "9")


if __name__ == "__main__":
    unittest.main()
