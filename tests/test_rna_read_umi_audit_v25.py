#!/usr/bin/env python3
"""Tests for RNA read-layer versus exact-UMI species auditing."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "03_scripts/05_downstream_analysis/05_rna_umap"
AUDIT = SCRIPT_DIR / "111_audit_rna_read_umi_cellline_labels.py"
SUBMIT = SCRIPT_DIR / "112_submit_rna_read_umi_cellline_audit_8core.csub"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class SpeciesLayerTests(unittest.TestCase):
    def test_high_confidence_species_has_ambiguous_and_no_signal(self) -> None:
        module = load(AUDIT, "rna_read_umi_species_v25")
        self.assertEqual(module.classify_species(90, 10, 0.90)[0], "GM12878")
        self.assertEqual(module.classify_species(10, 90, 0.90)[0], "PK15")
        self.assertEqual(module.classify_species(60, 40, 0.90)[0], "ambiguous")
        self.assertEqual(module.classify_species(0, 0, 0.90)[0], "no_signal")

    def test_read_and_umi_labels_are_not_collapsed(self) -> None:
        module = load(AUDIT, "rna_read_umi_status_v25")
        self.assertEqual(module.comparison_status("GM12878", "PK15"), "discordant")
        self.assertEqual(module.comparison_status("GM12878", "GM12878"), "concordant")
        self.assertEqual(module.comparison_status("GM12878", "ambiguous"), "read_only_confident")

    def test_recommendation_blocks_suspect_umi_labels(self) -> None:
        module = load(AUDIT, "rna_read_umi_recommendation_v25")
        metrics = {
            "read_GM12878_fraction": 0.90,
            "umi_GM12878_fraction": 0.10,
            "confident_discordant_fraction": 0.80,
        }
        self.assertEqual(
            module.recommendation_from_gm_group(metrics),
            "read_layer_provisional_only_do_not_use_UMI_species_labels",
        )

    def test_species_specific_collapsing_is_explicit(self) -> None:
        module = load(AUDIT, "rna_read_umi_collapsing_v25")
        self.assertEqual(module.reads_per_umi(100, 10), 10.0)
        self.assertIsNone(module.reads_per_umi(100, 0))


class SubmissionTests(unittest.TestCase):
    def test_submit_uses_analysis_environment_and_stated_filter(self) -> None:
        submit = SUBMIT.read_text(encoding="utf-8")
        audit = AUDIT.read_text(encoding="utf-8")
        self.assertIn("#CSUB -n 8", submit)
        self.assertIn("RNA_UMAP_ENV:-snTAPS_analysis", submit)
        self.assertIn("--min-umis 25", submit)
        self.assertIn("--min-genes 10", submit)
        self.assertIn("does not silently choose one label layer", audit)


if __name__ == "__main__":
    unittest.main()
