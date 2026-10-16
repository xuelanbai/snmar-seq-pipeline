#!/usr/bin/env python3
"""Local tests for the DNA 100 kb UMAP stage."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "03_scripts/05_downstream_analysis/04_dna_umap"
BUILD = SCRIPT_DIR / "103_build_dna_100kb_umap.py"
D20 = SCRIPT_DIR / "104_submit_dna_100kb_umap_D20_32core.csub"
D50 = SCRIPT_DIR / "105_submit_dna_100kb_umap_D50_32core.csub"
LAUNCHER = SCRIPT_DIR / "106_submit_dna_100kb_umap_2jobs.sh"
SUMMARY = SCRIPT_DIR / "107_summarize_dna_100kb_umap.py"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FeatureSelectionTests(unittest.TestCase):
    def test_reference_feature_budget_is_balanced(self) -> None:
        module = load(BUILD, "dna_umap_v23")
        species = ["human"] * 6 + ["pig"] * 6
        selected = module.balanced_top_features(
            range(12), range(12), species, maximum=6
        )
        self.assertEqual(len(selected), 6)
        self.assertEqual(
            sum(species[index] == "human" for index in selected), 3
        )
        self.assertEqual(
            sum(species[index] == "pig" for index in selected), 3
        )

    def test_threshold_labels_are_stable(self) -> None:
        module = load(BUILD, "dna_umap_threshold_v23")
        self.assertEqual(module.threshold_id(20_000), "D_total_20k")
        self.assertEqual(module.threshold_id(50_000), "D_total_50k")
        with self.assertRaises(ValueError):
            module.threshold_id(0)

    def test_species_is_posthoc_not_a_filter_argument(self) -> None:
        text = BUILD.read_text(encoding="utf-8")
        self.assertIn(
            "unsupervised_embedding_species_posthoc_only", text
        )
        self.assertNotIn("--species", text)
        self.assertIn(
            "absent_cell_bin_is_NA_encoded_as_zero_residual", text
        )


class ClusterSubmissionTests(unittest.TestCase):
    def test_two_full_jobs_and_summary_are_wired(self) -> None:
        d20 = D20.read_text(encoding="utf-8")
        d50 = D50.read_text(encoding="utf-8")
        launcher = LAUNCHER.read_text(encoding="utf-8")
        summary = SUMMARY.read_text(encoding="utf-8")
        self.assertIn("#CSUB -n 32", d20)
        self.assertIn('DNA_UMAP_ENV:-snTAPS_analysis', d20)
        self.assertIn("--min-total-calls 20000", d20)
        self.assertIn("--min-total-calls 50000", d50)
        self.assertIn("CONFIRM_SUBMIT_DNA_UMAP", launcher)
        self.assertIn("[DRY RUN]", launcher)
        self.assertIn(
            'THRESHOLDS = ("D_total_20k", "D_total_50k")', summary
        )


if __name__ == "__main__":
    unittest.main()
