#!/usr/bin/env python3
"""Local tests for cell-by-100 kb DNA matrix construction."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "03_scripts/05_downstream_analysis/01_dna_100kb_bins/"
    "87_build_dna_100kb_sample_matrix.py"
)


def load():
    spec = importlib.util.spec_from_file_location("dna_100kb_v20", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


dna_100kb = load()


class BinDefinitionTests(unittest.TestCase):
    def test_terminal_bin_is_clipped_and_retained(self) -> None:
        rows = [
            {
                "contig_index": "0",
                "species": "human",
                "contig": "human_1",
                "length": "250001",
            },
            {
                "contig_index": "1",
                "species": "pig",
                "contig": "pig_1",
                "length": "100000",
            },
        ]
        bins, offsets, lengths = dna_100kb.build_bins(rows, 100_000)
        self.assertEqual(offsets, [0, 3])
        self.assertEqual(lengths, [250001, 100000])
        self.assertEqual(len(bins), 4)
        self.assertEqual(
            (bins[2]["start0"], bins[2]["end0"], bins[2]["width"]),
            (200000, 250001, 50001),
        )
        self.assertEqual(bins[2]["terminal_partial_bin"], 1)
        self.assertEqual(bins[3]["terminal_partial_bin"], 0)
        self.assertTrue(
            all(
                int(row["end0"])
                <= lengths[int(row["contig_index"])]
                for row in bins
            )
        )


class TemplateRenderingTests(unittest.TestCase):
    def test_rendered_well_job_is_not_misclassified_as_template(self) -> None:
        template = (
            ROOT
            / "03_scripts/05_downstream_analysis/01_dna_100kb_bins/"
            "88_submit_dna_100kb_sample_template_8core.csub"
        ).read_text(encoding="utf-8")
        rendered = template.replace("__SAMPLE_INDEX__", "1").replace(
            "__FORCE__", "0"
        )
        self.assertNotIn("__SAMPLE_INDEX__", rendered)
        self.assertNotIn("__FORCE__", rendered)
        self.assertIn('SAMPLE_INDEX="1"', rendered)
        self.assertIn('FORCE="${FORCE:-0}"', rendered)
        self.assertIn('== *"__"*', rendered)
        self.assertNotIn('== "1" || "$FORCE" == "0"', rendered)


class BinAggregationTests(unittest.TestCase):
    def test_WCG_and_GCH_are_aggregated_separately(self) -> None:
        observed = dna_100kb.aggregate_arrays(
            cell_index=np.asarray([2, 2, 2, 3, 3], dtype=np.uint32),
            contig_index=np.asarray([0, 0, 0, 0, 1], dtype=np.uint32),
            position0=np.asarray(
                [0, 99999, 100000, 249999, 99999], dtype=np.uint32
            ),
            context_code=np.asarray([1, 1, 2, 1, 2], dtype=np.uint8),
            covered_calls=np.asarray([1, 2, 1, 3, 2], dtype=np.uint32),
            converted_calls=np.asarray([1, 0, 1, 1, 0], dtype=np.uint32),
            bin_offsets=[0, 3],
            contig_lengths=[250001, 100000],
            bin_size=100_000,
            total_bins=4,
        )
        rows = {
            (int(cell), int(bin_index), int(context)): (
                int(sites),
                int(covered),
                int(converted),
                float(fraction),
            )
            for (
                cell,
                bin_index,
                context,
                sites,
                covered,
                converted,
                fraction,
            ) in zip(
                observed["cell_index"],
                observed["bin_index"],
                observed["context_code"],
                observed["observed_sites"],
                observed["covered_calls"],
                observed["converted_calls"],
                observed["signal_fraction"],
            )
        }
        self.assertEqual(
            set(rows),
            {(2, 0, 1), (2, 1, 2), (3, 2, 1), (3, 3, 2)},
        )
        self.assertEqual(rows[(2, 0, 1)][:3], (2, 3, 1))
        self.assertAlmostEqual(rows[(2, 0, 1)][3], 1 / 3, places=6)
        self.assertEqual(rows[(2, 1, 2)], (1, 1, 1, 1.0))
        self.assertEqual(rows[(3, 2, 1)][:3], (1, 3, 1))
        self.assertEqual(rows[(3, 3, 2)], (1, 2, 0, 0.0))

    def test_position_outside_contig_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside its contig"):
            dna_100kb.aggregate_arrays(
                cell_index=np.asarray([0], dtype=np.uint32),
                contig_index=np.asarray([0], dtype=np.uint32),
                position0=np.asarray([100], dtype=np.uint32),
                context_code=np.asarray([1], dtype=np.uint8),
                covered_calls=np.asarray([1], dtype=np.uint32),
                converted_calls=np.asarray([0], dtype=np.uint32),
                bin_offsets=[0],
                contig_lengths=[100],
                bin_size=100,
                total_bins=1,
            )


if __name__ == "__main__":
    unittest.main()
