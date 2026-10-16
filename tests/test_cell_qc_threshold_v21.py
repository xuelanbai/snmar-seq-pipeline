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
    / "03_scripts/05_downstream_analysis/03_cell_qc_filtering/"
    "96_audit_cell_qc_thresholds.py"
)
CSUB = (
    ROOT
    / "03_scripts/05_downstream_analysis/03_cell_qc_filtering/"
    "97_submit_cell_qc_threshold_audit_4core.csub"
)


def load_module():
    spec = importlib.util.spec_from_file_location("cell_qc_threshold_audit", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CellQcThresholdAuditTests(unittest.TestCase):
    def test_threshold_grid_and_cross_modal_reconciliation(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dna_path = root / "dna.tsv.gz"
            rna_path = root / "rna.tsv.gz"
            dna_summary = root / "dna_summary.tsv"
            rna_summary = root / "rna_summary.tsv"
            output = root / "out"
            dna_rows = [
                {
                    "sample_index": 1,
                    "sample": "DNA-1",
                    "group": "GM12878",
                    "CB": "A" * 36,
                    "WCG_observed_sites": 2000,
                    "GCH_observed_sites": 10000,
                    "WCG_covered_calls": 4000,
                    "GCH_covered_calls": 8000,
                    "human_signal_covered_calls": 11000,
                    "pig_signal_covered_calls": 1000,
                },
                {
                    "sample_index": 1,
                    "sample": "DNA-1",
                    "group": "GM12878",
                    "CB": "C" * 36,
                    "WCG_observed_sites": 500,
                    "GCH_observed_sites": 2500,
                    "WCG_covered_calls": 1000,
                    "GCH_covered_calls": 2000,
                    "human_signal_covered_calls": 1500,
                    "pig_signal_covered_calls": 1500,
                },
                {
                    "sample_index": 2,
                    "sample": "DNA-2",
                    "group": "cellMix",
                    "CB": "G" * 36,
                    "WCG_observed_sites": 10000,
                    "GCH_observed_sites": 40000,
                    "WCG_covered_calls": 10000,
                    "GCH_covered_calls": 40000,
                    "human_signal_covered_calls": 1000,
                    "pig_signal_covered_calls": 49000,
                },
                {
                    "sample_index": 2,
                    "sample": "DNA-2",
                    "group": "cellMix",
                    "CB": "T" * 36,
                    "WCG_observed_sites": 0,
                    "GCH_observed_sites": 0,
                    "WCG_covered_calls": 0,
                    "GCH_covered_calls": 0,
                    "human_signal_covered_calls": 0,
                    "pig_signal_covered_calls": 0,
                },
            ]
            write_tsv(dna_path, list(dna_rows[0]), dna_rows)
            write_tsv(
                dna_summary,
                ["metric", "value"],
                [
                    {"metric": "status", "value": "validated"},
                    {"metric": "cells", "value": 4},
                ],
            )
            rna_rows = [
                {
                    "sample": "RNA-1",
                    "CB": "A" * 36,
                    "exact_umis": 600,
                    "genes_detected": 250,
                    "human_exact_umis": 580,
                    "pig_exact_umis": 20,
                    "in_raw_matrix": 1,
                },
                {
                    "sample": "RNA-1",
                    "CB": "C" * 36,
                    "exact_umis": 25,
                    "genes_detected": 10,
                    "human_exact_umis": 13,
                    "pig_exact_umis": 12,
                    "in_raw_matrix": 1,
                },
                {
                    "sample": "RNA-2",
                    "CB": "G" * 36,
                    "exact_umis": 50,
                    "genes_detected": 20,
                    "human_exact_umis": 0,
                    "pig_exact_umis": 50,
                    "in_raw_matrix": 1,
                },
                {
                    "sample": "RNA-2",
                    "CB": "T" * 36,
                    "exact_umis": 9,
                    "genes_detected": 5,
                    "human_exact_umis": 9,
                    "pig_exact_umis": 0,
                    "in_raw_matrix": 1,
                },
            ]
            write_tsv(rna_path, list(rna_rows[0]), rna_rows)
            write_tsv(
                rna_summary,
                ["metric", "value"],
                [
                    {"metric": "status", "value": "validated"},
                    {"metric": "raw_matrix_cells", "value": 4},
                    {"metric": "raw_matrix_umi_sum", "value": 684},
                ],
            )

            argv = sys.argv
            try:
                sys.argv = [
                    str(SCRIPT),
                    "--dna-cell-qc",
                    str(dna_path),
                    "--dna-summary",
                    str(dna_summary),
                    "--rna-cell-qc",
                    str(rna_path),
                    "--rna-summary",
                    str(rna_summary),
                    "--output-root",
                    str(output),
                    "--expected-wells",
                    "2",
                    "--approx-nuclei-per-well",
                    "2",
                ]
                self.assertEqual(module.main(), 0)
            finally:
                sys.argv = argv

            dna_grid = {
                row["threshold_id"]: row
                for row in read_tsv(output / "dna_threshold_grid.tsv")
            }
            rna_grid = {
                row["threshold_id"]: row
                for row in read_tsv(output / "rna_threshold_grid.tsv")
            }
            self.assertEqual(dna_grid["D_signal_any"]["retained_cells"], "3")
            self.assertEqual(
                dna_grid["D_total_10k_purity90"]["retained_cells"], "2"
            )
            self.assertEqual(rna_grid["R_UMI10_gene5"]["retained_cells"], "3")
            self.assertEqual(rna_grid["R_UMI500_gene200"]["retained_cells"], "1")

            intersections = read_tsv(output / "dna_rna_threshold_intersection.tsv")
            target = next(
                row
                for row in intersections
                if row["DNA_threshold_id"] == "D_total_10k_purity90"
                and row["RNA_threshold_id"] == "R_UMI25_gene10"
            )
            self.assertEqual(target["shared_well_CB"], "2")
            self.assertEqual(target["species_comparable_shared"], "2")
            self.assertEqual(target["species_concordant_shared"], "2")
            self.assertEqual(target["species_concordance_fraction"], "1.0")

            summary = {
                row["metric"]: row["value"]
                for row in read_tsv(output / "cell_qc_threshold_audit_summary.tsv")
            }
            self.assertEqual(summary["status"], "validated")
            self.assertEqual(summary["automatic_cell_filtering"], "none")

    def test_submitter_points_to_audit_only_script(self) -> None:
        text = CSUB.read_text(encoding="utf-8")
        self.assertIn("96_audit_cell_qc_thresholds.py", text)
        self.assertIn("#CSUB -n 4", text)
        self.assertNotIn("csub <", text)


if __name__ == "__main__":
    unittest.main()
