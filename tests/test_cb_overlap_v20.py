#!/usr/bin/env python3
"""Local tests for the exact-CB overlap audit across physical wells."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "03_scripts/05_downstream_analysis/02_cross_well_summary/"
    "94_audit_cb_overlap_96wells.py"
)


def load():
    spec = importlib.util.spec_from_file_location("cb_overlap_v20", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


audit = load()


class BarcodeValidationTests(unittest.TestCase):
    def test_exact_36bp_CB_is_required(self) -> None:
        audit.validate_cb("ACGT" * 9)
        with self.assertRaises(ValueError):
            audit.validate_cb("ACGT")
        with self.assertRaises(ValueError):
            audit.validate_cb("N" * 36)


class ScopeAuditTests(unittest.TestCase):
    def test_cross_well_reuse_does_not_collapse_cell_keys(self) -> None:
        well_sets = audit.make_well_sets()
        shared = "A" * 36
        well_sets[1] = {shared, "C" * 36}
        well_sets[2] = {shared, "G" * 36}
        well_sets[3] = {shared, "T" * 36}
        metadata = {
            index: (
                f"well-{index}",
                "GM12878" if index <= 48 else "cellMix",
            )
            for index in range(1, 97)
        }
        summary, per_well, pairwise, reused, distribution = audit.audit_scope(
            "test", well_sets, metadata
        )
        self.assertEqual(summary["well_CB_cell_keys"], 6)
        self.assertEqual(summary["unique_CB_sequences"], 4)
        self.assertEqual(summary["CB_sequences_in_multiple_wells"], 1)
        self.assertEqual(summary["well_CB_cell_keys_with_reused_CB"], 3)
        self.assertEqual(summary["well_pairs_with_nonzero_overlap"], 3)
        self.assertEqual(len(per_well), 96)
        self.assertEqual(len(pairwise), 4560)
        self.assertEqual(reused[0]["well_indices"], "1,2,3")
        distribution_by_wells = {
            int(row["number_of_wells_per_CB"]): int(row["CB_sequences"])
            for row in distribution
        }
        self.assertEqual(distribution_by_wells, {1: 3, 3: 1})


if __name__ == "__main__":
    unittest.main()
