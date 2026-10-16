#!/usr/bin/env python3
"""Tests for strict RNA exon/junction evidence recounting."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "03_scripts/05_downstream_analysis/05_rna_umap"
RECOUNT = SCRIPT_DIR / "113_recount_rna_expression_evidence_one.py"
TEMPLATE = SCRIPT_DIR / "114_submit_rna_expression_evidence_template_32core.csub"
LAUNCHER = SCRIPT_DIR / "115_submit_rna_expression_evidence_4jobs.sh"
SUMMARY = SCRIPT_DIR / "116_summarize_rna_expression_evidence.py"
SUMMARY_JOB = SCRIPT_DIR / "117_submit_rna_expression_evidence_summary_8core.csub"
MATRIX = SCRIPT_DIR / "118_build_rna_exon_salvage_matrix.py"
MATRIX_JOB = SCRIPT_DIR / "119_submit_rna_exon_salvage_matrix_16core.csub"
UMAP = SCRIPT_DIR / "120_build_rna_exon_salvage_umap.py"
UMAP_JOB = SCRIPT_DIR / "121_submit_rna_exon_salvage_umap_8core.csub"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class EvidenceDefinitionTests(unittest.TestCase):
    def test_cigar_reports_junction_and_aligned_blocks(self) -> None:
        module = load(RECOUNT, "rna_evidence_cigar_v26")
        blocks, spliced = module.cigar_blocks(100, "25M100N25M")
        self.assertEqual(blocks, [(100, 125), (225, 250)])
        self.assertTrue(spliced)
        blocks, spliced = module.cigar_blocks(100, "50M")
        self.assertEqual(blocks, [(100, 150)])
        self.assertFalse(spliced)

    def test_interval_is_checked_against_each_spliced_block(self) -> None:
        module = load(RECOUNT, "rna_evidence_intervals_v26")
        bins = module.IntervalBins(bin_size=1000)
        bins.add("human_1", 225, 250, "HUMAN_G1")
        self.assertEqual(
            bins.query("human_1", [(100, 125), (225, 250)]),
            {"HUMAN_G1"},
        )

    def test_reverse_strand_assignment_remains_explicit(self) -> None:
        module = load(RECOUNT, "rna_evidence_reverse_v26")
        gene = module.Gene(
            "HUMAN_G1", "G1", "human_1", 100, 300, "+", "human"
        )
        genes = {gene.gene_id: gene}
        self.assertEqual(
            module.reverse_unique_gene({gene.gene_id}, "-", genes)[0],
            "assigned",
        )
        self.assertEqual(
            module.reverse_unique_gene({gene.gene_id}, "+", genes)[0],
            "no_feature",
        )

    def test_joint_gtf_loads_gene_and_deduplicated_exons(self) -> None:
        module = load(RECOUNT, "rna_evidence_gtf_v26")
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "joint.gtf"
            path.write_text(
                "human_1\tx\tgene\t101\t300\t.\t+\t.\tgene_id \"HUMAN_G1\"; gene_name \"G1\";\n"
                "human_1\tx\texon\t101\t150\t.\t+\t.\tgene_id \"HUMAN_G1\"; transcript_id \"HUMAN_T1\";\n"
                "human_1\tx\texon\t101\t150\t.\t+\t.\tgene_id \"HUMAN_G1\"; transcript_id \"HUMAN_T2\";\n",
                encoding="utf-8",
            )
            models = module.load_gene_models(path)
            self.assertEqual(len(models.genes), 1)
            self.assertEqual(models.exon_intervals, 1)
            self.assertEqual(models.exons.query("human_1", [(100, 150)]), {"HUMAN_G1"})


class SubmissionTests(unittest.TestCase):
    def test_four_full_jobs_use_all_32_requested_cores(self) -> None:
        template = TEMPLATE.read_text(encoding="utf-8")
        launcher = LAUNCHER.read_text(encoding="utf-8")
        self.assertIn("#CSUB -n 32", template)
        self.assertIn("PARALLEL_FILES:-8", template)
        self.assertIn("SAMTOOLS_THREADS:-4", template)
        self.assertIn("PARALLEL_FILES * SAMTOOLS_THREADS <= 32", template)
        self.assertIn("for shard in 1 2 3 4", launcher)
        self.assertIn("CONFIRM_SUBMIT_RNA_EVIDENCE", launcher)

    def test_old_matrix_is_not_overwritten_and_summary_is_guarded(self) -> None:
        recount = RECOUNT.read_text(encoding="utf-8")
        summary = SUMMARY.read_text(encoding="utf-8")
        job = SUMMARY_JOB.read_text(encoding="utf-8")
        self.assertIn("salvage audit", recount.lower())
        self.assertIn("old_genefull_matrix_replaced", summary)
        self.assertIn("completed == 96", job)
        self.assertIn("RNA_UMAP_ENV:-snTAPS_analysis", job)

    def test_salvage_matrix_is_separate_and_umap_is_marked_exploratory(self) -> None:
        matrix = MATRIX.read_text(encoding="utf-8")
        matrix_job = MATRIX_JOB.read_text(encoding="utf-8")
        umap = UMAP.read_text(encoding="utf-8")
        umap_job = UMAP_JOB.read_text(encoding="utf-8")
        self.assertIn("06_rna_exon_salvage_matrix", matrix)
        self.assertIn("formal_expression_clustering_ready", matrix)
        self.assertIn("--min-umis 5", matrix_job)
        self.assertIn("--min-genes 3", matrix_job)
        self.assertIn("exploratory_low_depth", umap)
        self.assertIn("species-specific genes make this partly circular", umap)
        self.assertIn("RNA_UMAP_ENV:-snTAPS_analysis", umap_job)


class SalvageMatrixTests(unittest.TestCase):
    def test_species_label_keeps_ties_explicit(self) -> None:
        module = load(MATRIX, "rna_exon_matrix_v26")
        self.assertEqual(module.dominant_species(10, 2)[0], "human")
        self.assertEqual(module.dominant_species(2, 10)[0], "pig")
        self.assertEqual(module.dominant_species(5, 5)[0], "tie")
        self.assertEqual(module.dominant_species(0, 0)[0], "no_signal")

    def test_variable_features_are_balanced_by_reference_species(self) -> None:
        try:
            import numpy as np
        except ImportError:
            self.skipTest("numpy unavailable")
        module = load(UMAP, "rna_exon_umap_features_v26")
        variances = np.asarray([10.0, 9.0, 8.0, 7.0, 6.0])
        species = np.asarray(["human", "human", "human", "pig", "pig"])
        detected = np.asarray([4, 4, 4, 4, 4])
        selected = module.balanced_variable_features(
            variances, species, detected, maximum_per_species=2
        )
        self.assertEqual(set(selected.tolist()), {0, 1, 3, 4})


if __name__ == "__main__":
    unittest.main()
