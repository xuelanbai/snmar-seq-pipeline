#!/usr/bin/env python3
"""Tests for the GM12878/PK15 DNA UMAP consistency audit."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "03_scripts/05_downstream_analysis/04_dna_umap"
AUDIT = SCRIPT_DIR / "109_audit_dna_GM12878_PK15_umap.py"
SUBMIT = SCRIPT_DIR / "110_submit_dna_GM12878_PK15_umap_audit_8core.csub"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class CellLineAssignmentTests(unittest.TestCase):
    def test_mapping_fraction_assignment_has_ambiguous_class(self) -> None:
        module = load(AUDIT, "dna_cellline_audit_v24")
        self.assertEqual(
            module.classify_cell_line(900, 100, 0.90)[0], "GM12878"
        )
        self.assertEqual(
            module.classify_cell_line(100, 900, 0.90)[0], "PK15"
        )
        self.assertEqual(
            module.classify_cell_line(600, 400, 0.90)[0], "ambiguous"
        )
        self.assertEqual(
            module.classify_cell_line(0, 0, 0.90)[0], "ambiguous"
        )

    def test_well_labels_are_not_cell_line_labels(self) -> None:
        module = load(AUDIT, "dna_cellline_wells_v24")
        self.assertEqual(
            module.well_design("GM12878"), "GM12878_only_wells"
        )
        self.assertEqual(
            module.well_design("cellMix"), "GM12878_PK15_mix_wells"
        )

    def test_random_neighbor_baseline_is_auditable(self) -> None:
        module = load(AUDIT, "dna_cellline_baseline_v24")
        self.assertAlmostEqual(
            module.same_label_baseline(["a", "a", "b"]), 1 / 3
        )


class SubmissionTests(unittest.TestCase):
    def test_full_D50_audit_job_is_wired(self) -> None:
        submit = SUBMIT.read_text(encoding="utf-8")
        audit = AUDIT.read_text(encoding="utf-8")
        self.assertIn("#CSUB -n 8", submit)
        self.assertIn("DNA_UMAP_ENV:-snTAPS_analysis", submit)
        self.assertIn("--threshold D_total_50k", submit)
        self.assertIn("--minimum-major-fraction 0.90", submit)
        self.assertIn("DNA_internal_consistency_only", audit)
        self.assertIn("does_not_prove_independent", audit)


if __name__ == "__main__":
    unittest.main()
