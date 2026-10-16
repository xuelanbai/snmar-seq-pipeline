import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = (
    ROOT
    / "03_scripts"
    / "05_downstream_analysis"
    / "06_dna_tss_umap"
)


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


umap_builder = load_module("dna_tss_umap_v29", "129_build_dna_tss_umap.py")
hvf_summary = load_module("dna_tss_hvf_summary_v29", "149_summarize_dna_tss_hvf.py")


class TestDnaTssHighVariablePromoters(unittest.TestCase):
    def test_run_tag_preserves_default_and_isolates_sensitivity_runs(self):
        defaults = SimpleNamespace(
            feature_space="ortholog_single_gene",
            min_total_calls=20_000,
            min_promoter_calls_per_context=25,
            run_tag="",
        )
        self.assertEqual(
            umap_builder.run_id(defaults),
            "ortholog_single_gene_D20k_TSSctx25",
        )
        defaults.min_promoter_calls_per_context = 100
        defaults.run_tag = "HVF2000_F20_C100"
        self.assertEqual(
            umap_builder.run_id(defaults),
            "ortholog_single_gene_D20k_TSSctx100_HVF2000_F20_C100",
        )

    def test_fixed_grid_has_four_unique_run_ids(self):
        names = [hvf_summary.run_id(depth, tag) for depth, _, tag in hvf_summary.RUNS]
        self.assertEqual(len(names), 4)
        self.assertEqual(len(set(names)), 4)
        self.assertIn(
            "ortholog_single_gene_D50k_TSSctx100_HVF1000_F20_C100",
            names,
        )

    def test_template_uses_predeclared_feature_thresholds(self):
        template = (
            SCRIPT_ROOT / "147_submit_dna_tss_hvf_template_32core.csub"
        ).read_text(encoding="utf-8")
        self.assertIn("--min-promoter-calls-per-context 100", template)
        self.assertIn("--min-feature-cells 20", template)
        self.assertIn("--min-feature-covered-calls 100", template)
        self.assertIn('--run-tag "$RUN_TAG"', template)
        self.assertIn("--drop-empty-selected-feature-cells", template)
        self.assertIn("--max-empty-selected-feature-cell-fraction 0.01", template)
        self.assertIn("submit through step 148", template)

    def test_empty_selected_feature_error_records_context_and_cells(self):
        error = umap_builder.EmptySelectedFeatureCells("WCG", [2, 7])
        self.assertEqual(error.context_name, "WCG")
        self.assertEqual(error.local_indices, (2, 7))
        self.assertIn("leave 2 cells empty", str(error))

    def test_summary_has_predeclared_wcg_gate(self):
        source = (
            SCRIPT_ROOT / "149_summarize_dna_tss_hvf.py"
        ).read_text(encoding="utf-8")
        self.assertIn('("predeclared_minimum_WCG_silhouette", 0.10)', source)
        self.assertIn(
            '("predeclared_minimum_WCG_knn_excess_over_random", 0.05)',
            source,
        )


if __name__ == "__main__":
    unittest.main()
