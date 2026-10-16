import importlib.util
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import numpy as np

try:
    import scipy.sparse as sparse
except ImportError:  # Small semantic fallback for the bundled desktop runtime.
    class _DenseMatrix:
        def __init__(self, values):
            self.values = np.asarray(values)

        def tocsr(self):
            return self

        def astype(self, dtype):
            return _DenseMatrix(self.values.astype(dtype))

        def toarray(self):
            return self.values

    class _SparseFallback:
        @staticmethod
        def csr_matrix(values):
            return _DenseMatrix(values)

        @staticmethod
        def coo_matrix(specification, shape):
            data, (rows, columns) = specification
            values = np.zeros(shape, dtype=np.asarray(data).dtype)
            np.add.at(values, (rows, columns), data)
            return _DenseMatrix(values)

    sparse = _SparseFallback()


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = (
    ROOT
    / "03_scripts"
    / "05_downstream_analysis"
    / "07_dna_ortholog_100kb"
)


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


indexer = load_module(
    "dna_ortholog_100kb_index_v28",
    "135_prepare_dna_ortholog_100kb_index.py",
)
exporter = load_module(
    "dna_ortholog_100kb_export_v28",
    "140_export_dna_ortholog_100kb_seurat_matrices.py",
)
well_builder = load_module(
    "dna_ortholog_100kb_well_v28",
    "136_build_dna_ortholog_100kb_well_matrix.py",
)


class TestOrthologIndex(unittest.TestCase):
    def test_reference_wcg_count_matches_caller_contexts(self):
        self.assertEqual(indexer.count_wcg_in_sequence("ACG", 1, 1), 1)
        self.assertEqual(indexer.count_wcg_in_sequence("TCG", 1, 1), 1)
        self.assertEqual(indexer.count_wcg_in_sequence("GCG", 1, 1), 0)
        self.assertEqual(indexer.count_wcg_in_sequence("CGA", 1, 1), 1)
        self.assertEqual(indexer.count_wcg_in_sequence("CGT", 1, 1), 1)
        self.assertEqual(indexer.count_wcg_in_sequence("CGC", 1, 1), 0)
        # A CpG flanked by W on both sides can be called on both strands.
        self.assertEqual(indexer.count_wcg_in_sequence("ACGT", 1, 2), 2)

    def test_positive_chain_builds_two_nonoverlapping_orthobins(self):
        chain_text = "chain 100 chr1 1000 + 0 200 chr1 1000 + 100 300 1\n200\n"
        contigs = [
            {"contig": "human_1", "species": "human", "length": 1000, "contig_index": 0},
            {"contig": "pig_1", "species": "pig", "length": 1000, "contig_index": 1},
        ]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "test.chain"
            path.write_text(chain_text, encoding="utf-8")
            candidates, _metrics = indexer.collect_candidates(
                path,
                indexer.build_alias_lookup(contigs, "human"),
                indexer.build_alias_lookup(contigs, "pig"),
                100,
            )
        pairs, metrics = indexer.choose_pairs(candidates, 0.8, 0.7, 1.3)
        self.assertEqual(len(pairs), 2)
        self.assertEqual((pairs[0]["pig_start0"], pairs[0]["pig_end0"]), (100, 200))
        self.assertEqual((pairs[1]["pig_start0"], pairs[1]["pig_end0"]), (200, 300))
        self.assertEqual(metrics["accepted_before_WCG_filter"], 2)

    def test_negative_query_chain_uses_forward_pig_coordinates(self):
        header = indexer.parse_chain_header(
            "chain 100 chr1 1000 + 0 200 chr1 1000 - 100 300 7"
        )
        block = {
            **header,
            "t_block_start0": 0,
            "t_block_end0": 200,
            "q_block_start0": 700,
            "q_block_end0": 900,
        }
        self.assertEqual(indexer.project_segment(block, 0, 100), (800, 900))
        self.assertEqual(indexer.project_segment(block, 100, 200), (700, 800))


class TestOrthologMatrices(unittest.TestCase):
    def test_shared_interval_helper_is_loadable(self):
        helper = well_builder.load_interval_helpers()
        self.assertTrue(callable(helper.build_interval_lookup))
        self.assertTrue(callable(helper.aggregate_source_arrays))

    def test_two_species_physical_features_collapse_to_one_axis(self):
        features = [
            {"feature_index": "0", "orthobin_index": "0"},
            {"feature_index": "1", "orthobin_index": "0"},
        ]
        axes = [
            {
                "orthobin_index": "0",
                "human_feature_index": "0",
                "pig_feature_index": "1",
            }
        ]
        mapping = exporter.load_feature_axis(features, axes, np)
        self.assertEqual(mapping.tolist(), [0, 0])
        collapsed = exporter.collapse_arrays(
            cell=np.asarray([0, 0, 1]),
            feature=mapping[np.asarray([0, 1, 1])],
            covered=np.asarray([2, 3, 4]),
            converted=np.asarray([1, 1, 3]),
            n_features=1,
            np=np,
        )
        self.assertEqual(collapsed["covered"].tolist(), [5.0, 4.0])
        self.assertEqual(collapsed["converted"].tolist(), [2.0, 3.0])

    def test_missing_values_have_method_specific_neutral_meaning(self):
        arrays = {
            "cell": np.asarray([0, 1], dtype=np.int32),
            "feature": np.asarray([0, 0], dtype=np.int32),
            "covered": np.asarray([5.0, 5.0]),
            "converted": np.asarray([1.0, 4.0]),
        }
        args = Namespace(
            min_cell_bin_calls=5,
            prior_strength=20.0,
        )
        residual = np.asarray([-1.0, 1.0])
        matrices = exporter.build_matrices(
            arrays=arrays,
            selected=[0],
            pooled=np.asarray([0.5]),
            residual=residual,
            global_means=np.asarray([0.2, 0.8, 0.5]),
            args=args,
            np=np,
            sparse=sparse,
        )
        self.assertAlmostEqual(matrices["ratio_paper"].toarray()[2, 0], 0.5)
        self.assertAlmostEqual(matrices["beta_posterior"].toarray()[2, 0], 1.0)
        self.assertAlmostEqual(matrices["pearson_residual"].toarray()[2, 0], 0.0)
        self.assertAlmostEqual(matrices["coverage_control"].toarray()[2, 0], 0.0)

    def test_run_id_is_stable(self):
        args = Namespace(
            min_total_calls=20_000,
            min_orthobin_wcg_calls=100,
            min_cell_bin_calls=5,
            max_features=2_000,
        )
        self.assertEqual(
            exporter.run_id(args),
            "Dtotal20k_orthoWCG100_call5_top2000",
        )

    def test_feature_filter_requires_support_from_both_references(self):
        arrays = {
            "cell": np.asarray([0, 1, 0, 1, 0, 1, 0, 1], dtype=np.int32),
            "feature": np.asarray([0, 0, 1, 1, 2, 2, 3, 3], dtype=np.int32),
            "covered": np.asarray([10] * 8, dtype=float),
            "converted": np.asarray([2, 8, 2, 8, 3, 7, 4, 6], dtype=float),
        }
        side_support = {
            "cells": {
                "human": np.asarray([2, 2, 2, 2]),
                "pig": np.asarray([0, 2, 2, 2]),
            },
            "covered": {
                "human": np.asarray([20, 20, 20, 20], dtype=float),
                "pig": np.asarray([0, 20, 20, 20], dtype=float),
            },
        }
        args = Namespace(
            min_feature_cells=2,
            min_feature_covered_calls=20,
            min_feature_side_cells=1,
            min_feature_side_covered_calls=10,
            residual_clip=10.0,
            max_features=10,
        )
        result = exporter.feature_statistics(
            arrays, side_support, 4, args, np
        )
        self.assertEqual(result["eligible"].tolist(), [1, 2, 3])


if __name__ == "__main__":
    unittest.main()
