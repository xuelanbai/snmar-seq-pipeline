import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = ROOT / "03_scripts" / "05_downstream_analysis" / "06_dna_tss_umap"


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


builder = load_module("dna_tss_builder_v30", "129_build_dna_tss_umap.py")
audit = load_module(
    "dna_tss_reference_aware_audit_v30",
    "153_summarize_dna_tss_reference_aware.py",
)


class TestDnaReferenceAwarePromoter(unittest.TestCase):
    def test_run_ids_match_submitter(self):
        for depth, expected in (
            (20_000, audit.RUNS[0]),
            (50_000, audit.RUNS[1]),
        ):
            args = SimpleNamespace(
                feature_space="promoter_cluster",
                min_total_calls=depth,
                min_promoter_calls_per_context=25,
                run_tag="reference_aware",
            )
            self.assertEqual(builder.run_id(args), expected)

    def test_template_keeps_species_specific_promoter_axes(self):
        text = (
            SCRIPT_ROOT / "151_submit_dna_tss_reference_aware_template_32core.csub"
        ).read_text(encoding="utf-8")
        self.assertIn("--feature-space promoter_cluster", text)
        self.assertIn("--max-features-per-context 20000", text)
        self.assertIn("--run-tag reference_aware", text)
        self.assertIn("submit through step 152", text)

    def test_audit_requires_mixed_well_scope(self):
        source = (
            SCRIPT_ROOT / "153_summarize_dna_tss_reference_aware.py"
        ).read_text(encoding="utf-8")
        self.assertIn('row["well_set"] == "cellMix"', source)
        self.assertIn('("primary_scope", "mixed_wells_only")', source)
        self.assertIn("mapping_derived_posthoc_not_independent_validation", source)


if __name__ == "__main__":
    unittest.main()
