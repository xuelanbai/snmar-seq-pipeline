#!/usr/bin/env python3
"""Local tests for the one-sample, 96-well species-aware summary."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "03_scripts/05_downstream_analysis/02_cross_well_summary/"
    "92_summarize_96well_dna_rna.py"
)


def load():
    spec = importlib.util.spec_from_file_location(
        "cross_well_species_v20", SCRIPT
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


summary = load()


class SpeciesClassificationTests(unittest.TestCase):
    def test_dominant_species(self) -> None:
        self.assertEqual(summary.dominant_species(10, 1, "no_signal"), "human")
        self.assertEqual(summary.dominant_species(1, 10, "no_signal"), "pig")
        self.assertEqual(summary.dominant_species(4, 4, "no_signal"), "tie")
        self.assertEqual(
            summary.dominant_species(0, 0, "no_signal"), "no_signal"
        )

    def test_negative_species_count_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            summary.dominant_species(-1, 0, "no_signal")

    def test_weighted_ratio_uses_pooled_counts(self) -> None:
        # Per-well arithmetic mean would be 50%; pooled count ratio is 10/11.
        numerator = 0 + 10
        denominator = 1 + 10
        self.assertAlmostEqual(summary.ratio(numerator, denominator), 10 / 11)


class OutputTests(unittest.TestCase):
    def test_atomic_gzip_tsv_roundtrip(self) -> None:
        rows = [
            {
                "cell_id": "well-1|AAAA",
                "species_label": "human",
            },
            {
                "cell_id": "well-2|AAAA",
                "species_label": "pig",
            },
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "labels.tsv.gz"
            summary.atomic_tsv(path, list(rows[0]), rows)
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                observed = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(observed, rows)
        self.assertNotEqual(observed[0]["cell_id"], observed[1]["cell_id"])

    def test_empty_cross_modal_table_keeps_header(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty.tsv.gz"
            summary.atomic_tsv(
                path,
                ["well_index", "CB"],
                [],
                allow_empty=True,
            )
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "well_index\tCB\n")


if __name__ == "__main__":
    unittest.main()
