import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np


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


indexer = load_module("dna_tss_index_v27", "122_prepare_dna_tss_promoters.py")
builder = load_module("dna_tss_builder_v27", "123_build_dna_tss_well_matrix.py")
merger = load_module("dna_tss_merger_v27", "127_merge_dna_tss_complete_matrix.py")
umap_builder = load_module("dna_tss_umap_v27", "129_build_dna_tss_umap.py")
summarizer = load_module("dna_tss_summary_v27", "133_summarize_dna_tss_umap.py")


class TestDnaTssIndex(unittest.TestCase):
    def test_gtf_windows_are_strand_aware_and_clipped(self):
        contigs = [
            {
                "contig": "human_chr1",
                "species": "human",
                "length": 10000,
                "contig_index": 0,
            }
        ]
        text = (
            'human_chr1\ttest\tgene\t1\t100\t.\t+\t.\t'
            'gene_id "H1"; gene_name "EDGEPLUS";\n'
            'human_chr1\ttest\tgene\t9000\t10000\t.\t-\t.\t'
            'gene_id "H2"; gene_name "EDGEMINUS";\n'
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test.gtf"
            path.write_text(text, encoding="utf-8")
            rows = indexer.load_gene_windows(path, contigs, 1500)
        self.assertEqual(rows[0]["tss0"], 0)
        self.assertEqual((rows[0]["start0"], rows[0]["end0"]), (0, 1501))
        self.assertEqual(rows[0]["clipped_left"], 1)
        self.assertEqual(rows[1]["tss0"], 9999)
        self.assertEqual((rows[1]["start0"], rows[1]["end0"]), (8499, 10000))
        self.assertEqual(rows[1]["clipped_right"], 1)

    def test_overlaps_merge_and_ortholog_axes_are_conservative(self):
        def window(species, contig, start, end, gene_id, gene_name):
            return {
                "species": species,
                "contig_index": contig,
                "chrom": f"{species}_{contig}",
                "start0": start,
                "end0": end,
                "tss0": (start + end) // 2,
                "strand": "+",
                "gene_id": gene_id,
                "gene_name": gene_name,
                "gene_biotype": "protein_coding",
                "clipped_left": 0,
                "clipped_right": 0,
            }

        windows = [
            window("human", 0, 100, 200, "H_A", "ABC"),
            window("human", 0, 150, 250, "H_B", "MULTI"),
            window("human", 0, 300, 400, "H_C", "PAIR"),
            window("pig", 1, 500, 600, "P_C", "pair"),
            window("pig", 1, 700, 800, "P_D", "ONLYPIG"),
        ]
        clusters = indexer.build_promoter_clusters(windows)
        self.assertEqual(len(clusters), 4)
        self.assertEqual(clusters[0]["gene_count"], 2)
        pairs = indexer.add_conservative_ortholog_map(clusters)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["ortholog_key"], "PAIR")
        self.assertGreaterEqual(clusters[2]["ortholog_index"], 0)


class TestDnaTssAggregation(unittest.TestCase):
    def test_interval_assignment_and_count_reconciliation(self):
        feature_rows = [
            {"feature_index": "0", "contig_index": "0", "start0": "0", "end0": "10"},
            {"feature_index": "1", "contig_index": "0", "start0": "20", "end0": "30"},
        ]
        lookup = builder.build_interval_lookup(feature_rows, np)
        result, metrics = builder.aggregate_source_arrays(
            cell_index=np.asarray([0, 0, 0, 1]),
            contig_index=np.asarray([0, 0, 0, 0]),
            position0=np.asarray([5, 7, 15, 25]),
            context_code=np.asarray([1, 1, 1, 2]),
            covered_calls=np.asarray([1, 2, 1, 3]),
            converted_calls=np.asarray([1, 0, 0, 2]),
            lookup=lookup,
            total_features=2,
            np=np,
        )
        self.assertEqual(metrics["source_entries"], 4)
        self.assertEqual(metrics["selected_entries"], 3)
        self.assertEqual(result["observed_sites"].tolist(), [2, 1])
        self.assertEqual(result["covered_calls"].tolist(), [3, 3])
        self.assertEqual(result["converted_calls"].tolist(), [1, 2])
        combined = builder.combine_parts([result, result], 2, np)
        self.assertEqual(combined["observed_sites"].tolist(), [4, 2])
        self.assertEqual(combined["covered_calls"].tolist(), [6, 6])

    def test_overlapping_feature_index_is_rejected(self):
        rows = [
            {"feature_index": "0", "contig_index": "0", "start0": "0", "end0": "10"},
            {"feature_index": "1", "contig_index": "0", "start0": "9", "end0": "20"},
        ]
        with self.assertRaises(ValueError):
            builder.build_interval_lookup(rows, np)


class TestDnaTssFeatureSpace(unittest.TestCase):
    def test_human_and_pig_promoters_share_one_axis(self):
        promoters = [
            {"feature_index": "0", "feature_id": "H", "gene_names": "ABC", "species": "human"},
            {"feature_index": "1", "feature_id": "P", "gene_names": "ABC", "species": "pig"},
        ]
        orthologs = [
            {
                "ortholog_index": "0",
                "ortholog_key": "ABC",
                "human_feature_index": "0",
                "human_gene_id": "H1",
                "pig_feature_index": "1",
                "pig_gene_id": "P1",
                "pairing_rule": "test",
            }
        ]
        mapping, axes, species = umap_builder.load_feature_space(
            promoters, orthologs, "ortholog_single_gene", np
        )
        self.assertEqual(mapping.tolist(), [0, 0])
        self.assertEqual(len(axes), 1)
        self.assertEqual(species, ["shared_human_pig"])

    def test_mapping_label_requires_confidence(self):
        self.assertEqual(
            umap_builder.classify_cell_line(90, 10, 0.90),
            ("GM12878", 0.9),
        )
        self.assertEqual(
            umap_builder.classify_cell_line(89, 11, 0.90)[0],
            "ambiguous",
        )

    def test_feature_support_excludes_singleton_only_cells(self):
        arrays = {
            "cell": np.asarray([0, 1, 2], dtype=np.int32),
            "feature": np.asarray([0, 0, 1], dtype=np.int32),
            "covered": np.asarray([2.0, 3.0, 10.0]),
            "converted": np.asarray([1.0, 1.0, 5.0]),
        }
        mask = umap_builder.supported_cell_mask(
            arrays,
            n_cells=3,
            n_features=2,
            min_feature_cells=2,
            min_feature_covered_calls=5,
            np=np,
        )
        self.assertEqual(mask.tolist(), [True, True, False])

    def test_complete_stage_signatures_are_wired(self):
        self.assertEqual(
            merger.PIPELINE_SIGNATURE,
            umap_builder.SOURCE_SIGNATURE,
        )
        self.assertEqual(
            umap_builder.PIPELINE_SIGNATURE,
            summarizer.RUN_SIGNATURE,
        )


if __name__ == "__main__":
    unittest.main()
