#!/usr/bin/env python3
"""Local unit tests for v13 post-mapping helpers (no BAM files required)."""

from __future__ import annotations

import importlib.util
import csv
import gzip
import json
import sys
import tempfile
import types
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
    sys.modules.setdefault("pysam", types.SimpleNamespace())
    spec = importlib.util.spec_from_file_location(name, script(filename))
    if spec is None or spec.loader is None:
        raise RuntimeError(filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


prepare = load("prepare_postmapping", "33_prepare_postmapping_probe.py")
rna = load("rna_genefull", "39_count_rna_genefull_umi_probe.py")
dna_summary = load(
    "dna_cb_summary", "38_summarize_dna_cb_markdup_probe.py"
)
rna_full = load(
    "rna_reverse_full", "46_count_rna_reverse_genefull_one.py"
)
rna_merge = load(
    "rna_reverse_merge", "49_merge_rna_reverse_genefull_matrix.py"
)


class SelectionTests(unittest.TestCase):
    def test_add_choice_uses_row_identifier(self) -> None:
        choices = {}
        prepare.add_choice(choices, {"run_id": "A"}, "run_id", "first")
        prepare.add_choice(choices, {"run_id": "B"}, "run_id", "second")
        prepare.add_choice(choices, {"run_id": "A"}, "run_id", "duplicate")
        self.assertEqual(set(choices), {"A", "B"})
        self.assertEqual(choices["A"][1], "first")


class GeneAssignmentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.index = rna.GeneBins(bin_size=100)
        self.plus = rna.Gene(
            "HUMAN_G1", "G1", "human_1", 100, 300, "+", "human"
        )
        self.minus = rna.Gene(
            "HUMAN_G2", "G2", "human_1", 400, 600, "-", "human"
        )
        self.index.add(self.plus)
        self.index.add(self.minus)

    def test_half_open_interval_query(self) -> None:
        self.assertEqual(
            set(self.index.query("human_1", [(150, 180)])), {"HUMAN_G1"}
        )
        self.assertEqual(self.index.query("human_1", [(300, 350)]), {})

    def test_forward_reverse_assignment(self) -> None:
        candidates = {"HUMAN_G1": self.plus}
        status, gene = rna.assign_mode(candidates, "+", "forward")
        self.assertEqual(status, "assigned")
        self.assertEqual(gene, self.plus)
        status, gene = rna.assign_mode(candidates, "+", "reverse")
        self.assertEqual((status, gene), ("no_feature", None))

    def test_unstranded_overlap_is_ambiguous(self) -> None:
        candidates = {"HUMAN_G1": self.plus, "HUMAN_G2": self.minus}
        status, gene = rna.assign_mode(candidates, "+", "unstranded")
        self.assertEqual((status, gene), ("ambiguous", None))

    def test_spliced_cigar_blocks(self) -> None:
        self.assertEqual(
            rna.cigar_blocks(99, "10M90N5M2I3M"),
            [(99, 109), (199, 204), (204, 207)],
        )


class MarkdupStatsTests(unittest.TestCase):
    def test_markdup_stats_parser(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "stats.txt"
            path.write_text(
                "EXAMINED: 100\nDUPLICATE TOTAL: 25\nCOMMAND: samtools\n",
                encoding="utf-8",
            )
            parsed = dna_summary.read_markdup_stats(path)
        self.assertEqual(parsed["examined"], 100)
        self.assertEqual(parsed["duplicate_total"], 25)
        self.assertEqual(parsed["command"], "samtools")

    def test_dna_worker_uses_corrected_cb_template_markdup(self) -> None:
        text = script("35_run_dna_cb_markdup_one.sh").read_text(
            encoding="utf-8"
        )
        self.assertIn("--barcode-tag CB", text)
        self.assertIn("-c -s -m t", text)
        self.assertIn("-F 1024", text)
        self.assertIn(
            'PIPELINE_SIGNATURE="snmar_dna_cb_markdup_corrected_20260731"',
            text,
        )
        for block in text.split('"$SAMTOOLS_BIN" markdup \\')[1:]:
            self.assertNotIn("-l 1", "\n".join(block.splitlines()[:8]))

    def test_full_downstream_force_is_baked_into_rendered_jobs(self) -> None:
        template = (
            script("56_submit_dna_taps_cpg_full_template_32core.csub")
        ).read_text(encoding="utf-8")
        launcher = (
            script("57_submit_dna_taps_cpg_full_96jobs.sh")
        ).read_text(encoding="utf-8")
        self.assertIn('FORCE="${FORCE:-__FORCE__}"', template)
        self.assertIn('FORCE="${FORCE:-0}"', launcher)
        self.assertIn('s/__FORCE__/$FORCE/g', launcher)
        self.assertIn("if (( FORCE == 0 ))", launcher)


class FullRnaTests(unittest.TestCase):
    def test_reverse_assignment(self) -> None:
        plus = rna_full.Gene(
            "HUMAN_G1", "G1", "human_1", 0, 100, "+", "human"
        )
        status, gene = rna_full.reverse_gene({"HUMAN_G1": plus}, "-")
        self.assertEqual((status, gene), ("assigned", plus))
        status, gene = rna_full.reverse_gene({"HUMAN_G1": plus}, "+")
        self.assertEqual((status, gene), ("no_feature", None))

    def test_two_file_matrix_merge(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = root / "manifest.tsv"
            by_file = root / "by_file"
            output = root / "matrix"
            qc_root = root / "qc"
            manifest.write_text(
                "manifest_index\tshard\tshard_order\tpair_id\tsample\t"
                "star_input_reads\tprimary_mapped_reads\t"
                "primary_mapping_percent\tinput_bam\n"
                "1\t1\t1\tP1\tSample-1\t10\t8\t80\t/a.bam\n"
                "2\t1\t2\tP2\tSample-2\t20\t18\t90\t/b.bam\n",
                encoding="utf-8",
            )
            for pair_id, sample, input_reads, mapped, cell, gene, umis in (
                ("P1", "Sample-1", 10, 8, "AAA", "HUMAN_G1", 2),
                ("P2", "Sample-2", 20, 18, "BBB", "PIG_G2", 3),
            ):
                pair_root = by_file / pair_id
                pair_root.mkdir(parents=True)
                counts = pair_root / f"{pair_id}.cell_gene_counts.tsv.gz"
                cell_qc = pair_root / f"{pair_id}.cell_qc.tsv.gz"
                with gzip.open(counts, "wt", encoding="utf-8") as out:
                    out.write(
                        "pair_id\tsample\tcell_id\tCB\tCN\tgene_id\t"
                        "gene_name\tspecies\tassigned_reads\texact_umis\n"
                    )
                    species = "human" if gene.startswith("HUMAN_") else "pig"
                    out.write(
                        f"{pair_id}\t{sample}\t{sample}|{cell}\t{cell}\t"
                        f"CN\t{gene}\tGENE\t{species}\t5\t{umis}\n"
                    )
                with gzip.open(cell_qc, "wt", encoding="utf-8") as out:
                    out.write(
                        "pair_id\tsample\tcell_id\tCB\tCN\t"
                        "tag_eligible_unique_primary\tgene_assigned_reads\t"
                        "gene_assignment_fraction\texact_umis\t"
                        "genes_detected\thuman_assigned_reads\t"
                        "pig_assigned_reads\thuman_exact_umis\t"
                        "pig_exact_umis\n"
                    )
                    human_umi = umis if species == "human" else 0
                    pig_umi = umis if species == "pig" else 0
                    out.write(
                        f"{pair_id}\t{sample}\t{sample}|{cell}\t{cell}\t"
                        f"CN\t6\t5\t0.83333333\t{umis}\t1\t"
                        f"{5 if species == 'human' else 0}\t"
                        f"{5 if species == 'pig' else 0}\t"
                        f"{human_umi}\t{pig_umi}\n"
                    )
                counters = {
                    "primary_records": input_reads,
                    "primary_mapped": mapped,
                    "unique_primary": mapped,
                    "tag_eligible_unique_primary": 6,
                    "gene_assigned": 5,
                    "gene_ambiguous": 0,
                    "gene_no_feature": 1,
                    "missing_required_tags": 0,
                    "nonzero_XD": 0,
                }
                payload = {
                    "status": "validated",
                    "strand": "reverse",
                    "counting_feature": "GeneFull_gene_body",
                    "umi_deduplication": "exact_within_sample_CB_gene",
                    "counters": counters,
                    "exact_umis": umis,
                    "cells_tag_eligible": 1,
                    "cells_with_exact_umi": 1,
                    "detected_cell_gene_pairs": 1,
                    "cell_gene_counts": str(counts),
                    "cell_qc": str(cell_qc),
                }
                (pair_root / f"{pair_id}.summary.json").write_text(
                    json.dumps(payload), encoding="utf-8"
                )

            old_argv = sys.argv
            try:
                sys.argv = [
                    "merge",
                    "--manifest",
                    str(manifest),
                    "--by-file-root",
                    str(by_file),
                    "--output-root",
                    str(output),
                    "--qc-root",
                    str(qc_root),
                    "--expected-files",
                    "2",
                ]
                self.assertEqual(rna_merge.main(), 0)
            finally:
                sys.argv = old_argv
            with gzip.open(
                output / "matrix.mtx.gz", "rt", encoding="utf-8"
            ) as handle:
                matrix_text = handle.read()
            self.assertIn("2 2 2", matrix_text)
            self.assertIn("1 1 2", matrix_text)
            self.assertIn("2 2 3", matrix_text)
            with (
                qc_root / "rna_reverse_genefull_full_summary.tsv"
            ).open("rt", encoding="utf-8", newline="") as handle:
                summary = {
                    row["metric"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(summary["completed_files"], "2")
            self.assertEqual(summary["raw_matrix_umi_sum"], "5")


if __name__ == "__main__":
    unittest.main()
