#!/usr/bin/env python3
"""Local tests for v19 context-aware DNA and complete RNA matrix code."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POST = ROOT / "03_scripts/04_postmapping"


def script(filename: str) -> Path:
    matches = list(POST.rglob(filename))
    if len(matches) != 1:
        raise RuntimeError(f"expected one script named {filename}: {matches}")
    return matches[0]


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, script(filename))
    if spec is None or spec.loader is None:
        raise RuntimeError(filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


caller = load(
    "single_cytosine_caller_v19",
    "65_call_dna_single_cytosine_task.py",
)
track = load(
    "single_cytosine_track_v19",
    "71_aggregate_single_cytosine_track_shard.py",
)
validator = load(
    "single_cytosine_validator_v19",
    "68_validate_dna_single_cytosine_sample.py",
)
rna = load("rna_complete_v19", "76_validate_rna_complete_matrix.py")
bigwig = load(
    "single_base_bigwig_v19", "81_build_dna_single_base_bigwigs.py"
)
cell_qc = load(
    "single_cell_coverage_qc_v19",
    "83_summarize_dna_single_cell_qc.py",
)


class SingleCytosineCallerTests(unittest.TestCase):
    def test_aligned_bases_retain_query_positions_for_end_exclusion(self) -> None:
        observed = list(
            caller.aligned_bases(11, "2S4M1I2M", "AACCGTACG", "I" * 9)
        )
        self.assertEqual(
            observed,
            [
                (10, 2, "C", 40),
                (11, 3, "C", 40),
                (12, 4, "G", 40),
                (13, 5, "T", 40),
                (14, 7, "C", 40),
                (15, 8, "G", 40),
            ],
        )

    def test_forward_contexts_distinguish_WCG_GCH_and_GCG(self) -> None:
        self.assertEqual(
            caller.context_call(False, "AACGTT", 2, "T"),
            (2, 1, caller.CONTEXT_WCG, 1),
        )
        self.assertEqual(
            caller.context_call(False, "AGCAT", 2, "C"),
            (2, 1, caller.CONTEXT_GCH, 0),
        )
        self.assertEqual(
            caller.context_call(False, "AGCGT", 2, "T"),
            (2, 1, caller.CONTEXT_GCG_AMBIGUOUS, 1),
        )

    def test_reverse_contexts_distinguish_WCG_GCH_and_GCG(self) -> None:
        self.assertEqual(
            caller.context_call(True, "CGA", 1, "A"),
            (1, -1, caller.CONTEXT_WCG, 1),
        )
        self.assertEqual(
            caller.context_call(True, "TGC", 1, "G"),
            (1, -1, caller.CONTEXT_GCH, 0),
        )
        self.assertEqual(
            caller.context_call(True, "CGC", 1, "A"),
            (1, -1, caller.CONTEXT_GCG_AMBIGUOUS, 1),
        )

    def test_non_context_and_wrong_read_base_are_rejected(self) -> None:
        self.assertIsNone(caller.context_call(False, "AACATT", 2, "T"))
        self.assertIsNone(caller.context_call(False, "AACGTT", 2, "A"))
        self.assertIsNone(caller.context_call(True, "CGA", 1, "T"))

    def test_overlap_policy_keeps_best_and_discards_equal_conflict(self) -> None:
        pending = {}
        key = (7, "read1", 100, 1, caller.CONTEXT_WCG)
        self.assertEqual(caller.update_template_call(pending, key, 20, 0), "new")
        self.assertEqual(
            caller.update_template_call(pending, key, 30, 1), "replaced"
        )
        self.assertEqual(
            caller.update_template_call(pending, key, 10, 0),
            "ignored_lower_quality",
        )
        self.assertEqual(
            caller.update_template_call(pending, key, 30, 0), "conflict"
        )
        self.assertEqual(pending[key], (30, None))


class TrackMergeTests(unittest.TestCase):
    def test_96_sorted_sample_tracks_merge_WCG_and_GCH_separately(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            input_root = (
                project / "08_qc_summary/20_single_cytosine_inputs"
            )
            input_root.mkdir(parents=True)
            tasks_path = input_root / "dna_single_cytosine_tasks.tsv"
            task_fields = [
                "task_index",
                "sample_index",
                "sample",
                "group",
                "shard",
                "input_bam",
                "reference_fasta",
                "contig_file",
                "cell_index_path",
                "contig_index_path",
                "output_dir",
            ]
            site_fields = [
                "contig_index",
                "species",
                "chrom",
                "start",
                "end",
                "strand",
                "context",
                "signal_type",
                "covered_cells",
                "covered_calls",
                "converted_calls",
                "signal_fraction",
            ]
            task_rows = []
            for index in range(1, 97):
                sample = f"S{index:03d}"
                output_dir = project / "shards" / sample
                output_dir.mkdir(parents=True)
                task_rows.append(
                    {
                        "task_index": index,
                        "sample_index": index,
                        "sample": sample,
                        "group": "test",
                        "shard": 1,
                        "input_bam": "x",
                        "reference_fasta": "x",
                        "contig_file": "x",
                        "cell_index_path": "x",
                        "contig_index_path": "x",
                        "output_dir": output_dir,
                    }
                )
                rows = [
                    {
                        "contig_index": 0,
                        "species": "human",
                        "chrom": "human_1",
                        "start": 10,
                        "end": 11,
                        "strand": "+",
                        "context": "WCG",
                        "signal_type": "endogenous_methylation",
                        "covered_cells": 1,
                        "covered_calls": 2,
                        "converted_calls": index % 2,
                        "signal_fraction": (index % 2) / 2,
                    },
                    {
                        "contig_index": 0,
                        "species": "human",
                        "chrom": "human_1",
                        "start": 20,
                        "end": 21,
                        "strand": "-",
                        "context": "GCH",
                        "signal_type": "GpC_accessibility",
                        "covered_cells": 1,
                        "covered_calls": 3,
                        "converted_calls": index % 3,
                        "signal_fraction": (index % 3) / 3,
                    },
                ]
                with gzip.open(
                    output_dir / f"{sample}.shard_01.site_summary.tsv.gz",
                    "wt",
                    encoding="utf-8",
                    newline="",
                ) as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=site_fields,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerows(rows)
                validation = (
                    project
                    / "08_qc_summary/21_dna_single_cytosine_validation"
                    / sample
                )
                validation.mkdir(parents=True)
                (
                    validation
                    / f"{sample}.single_cytosine.validation_summary.tsv"
                ).write_text(
                    "metric\tvalue\nstatus\tvalidated\n"
                    "pipeline_signature\t"
                    "snmar_single_cytosine_validation_20260801\n"
                    f"sample\t{sample}\n",
                    encoding="utf-8",
                )
            with tasks_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=task_fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(task_rows)
            (
                input_root / "dna_single_cytosine_contig_index.tsv"
            ).write_text(
                "contig_index\treference_order\tshard\tspecies\tcontig\tlength\n"
                "0\t1\t1\thuman\thuman_1\t100\n",
                encoding="utf-8",
            )
            old_argv = sys.argv
            try:
                sys.argv = [
                    "71_track",
                    "--project",
                    str(project),
                    "--shard",
                    "1",
                ]
                self.assertEqual(track.main(), 0)
            finally:
                sys.argv = old_argv
            result = (
                project
                / "07_postmapping/complete_matrices/dna_single_cytosine/tracks"
                / "dna_single_cytosine.pseudobulk.shard_01.tsv.gz"
            )
            with gzip.open(
                result, "rt", encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["context"], "WCG")
            self.assertEqual(int(rows[0]["covered_calls"]), 192)
            self.assertEqual(int(rows[0]["converted_calls"]), 48)
            self.assertEqual(rows[1]["context"], "GCH")
            self.assertEqual(int(rows[1]["covered_calls"]), 288)
            self.assertEqual(int(rows[1]["converted_calls"]), 96)
            for filename in (
                "dna_WCG_methylation.shard_01.bedGraph.gz",
                "dna_GCH_accessibility.shard_01.bedGraph.gz",
            ):
                self.assertTrue((result.parent / filename).is_file())


class SampleValidationTests(unittest.TestCase):
    def test_eight_shards_reconcile_context_and_GCG_audit_counts(self) -> None:
        try:
            import h5py
        except ImportError:
            self.skipTest("h5py is validated in the server snTAPS_align env")
        import numpy as np

        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            input_root = (
                project / "08_qc_summary/20_single_cytosine_inputs"
            )
            input_root.mkdir(parents=True)
            (
                input_root / "single_cytosine_prepare_summary.tsv"
            ).write_text(
                "metric\tvalue\nstatus\tvalidated\n",
                encoding="utf-8",
            )
            cell_fields = [
                "cell_index",
                "sample_index",
                "sample",
                "group",
                "cell_id",
                "CB",
                "legacy_covered_template_cpg_calls",
            ]
            with gzip.open(
                input_root / "dna_single_cytosine_cell_index.tsv.gz",
                "wt",
                encoding="utf-8",
                newline="",
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=cell_fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for cell_index in range(2):
                    writer.writerow(
                        {
                            "cell_index": cell_index,
                            "sample_index": 1,
                            "sample": "S001",
                            "group": "test",
                            "cell_id": f"S001|CB{cell_index}",
                            "CB": f"CB{cell_index}",
                            "legacy_covered_template_cpg_calls": 0,
                        }
                    )
            contig_fields = [
                "contig_index",
                "reference_order",
                "shard",
                "species",
                "contig",
                "length",
            ]
            with (
                input_root / "dna_single_cytosine_contig_index.tsv"
            ).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=contig_fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for shard in range(1, 9):
                    writer.writerow(
                        {
                            "contig_index": shard - 1,
                            "reference_order": shard,
                            "shard": shard,
                            "species": "human",
                            "contig": f"human_{shard}",
                            "length": 100,
                        }
                    )
            task_fields = [
                "task_index",
                "sample_index",
                "sample",
                "group",
                "shard",
                "input_bam",
                "reference_fasta",
                "contig_file",
                "cell_index_path",
                "contig_index_path",
                "output_dir",
                "minimum_mapq",
                "minimum_baseq",
                "excluded_bases_each_read_end",
            ]
            task_rows = []
            matrix_fields = caller.MATRIX_FIELDS
            for shard in range(1, 9):
                output_dir = project / "shards" / f"shard_{shard:02d}"
                output_dir.mkdir(parents=True)
                prefix = f"S001.shard_{shard:02d}"
                context = (
                    caller.CONTEXT_WCG
                    if shard % 2
                    else caller.CONTEXT_GCH
                )
                converted = int(shard % 2)
                with h5py.File(
                    output_dir / f"{prefix}.single_cytosine.h5", "w"
                ) as target:
                    group = target.create_group("matrix")
                    arrays = {
                        "cell_index": np.asarray([shard % 2], dtype=np.uint32),
                        "contig_index": np.asarray(
                            [shard - 1], dtype=np.uint32
                        ),
                        "cytosine_position0": np.asarray(
                            [10], dtype=np.uint32
                        ),
                        "strand": np.asarray(
                            [1 if shard % 2 else -1], dtype=np.int8
                        ),
                        "context_code": np.asarray(
                            [context], dtype=np.uint8
                        ),
                        "covered": np.asarray([1], dtype=np.uint32),
                        "converted": np.asarray(
                            [converted], dtype=np.uint32
                        ),
                        "signal_fraction": np.asarray(
                            [float(converted)], dtype=np.float32
                        ),
                    }
                    for field in matrix_fields:
                        group.create_dataset(field, data=arrays[field])
                    target.attrs["status"] = "validated"
                    target.attrs["sample"] = "S001"
                    target.attrs["shard"] = shard
                    target.attrs["minimum_mapq"] = 30
                    target.attrs["minimum_baseq"] = 20
                    target.attrs["excluded_bases_each_read_end"] = 10
                    target.attrs["ambiguous_context_policy"] = (
                        "GCG_excluded_from_both_signals"
                    )
                (
                    output_dir / f"{prefix}.site_summary.tsv.gz"
                ).write_bytes(b"x")
                wcg = int(context == caller.CONTEXT_WCG)
                gch = int(context == caller.CONTEXT_GCH)
                summary = [
                    ("status", "validated"),
                    ("sample", "S001"),
                    ("shard", shard),
                    ("minimum_mapq", 30),
                    ("minimum_baseq", 20),
                    ("excluded_bases_each_read_end", 10),
                    ("observed_cell_cytosine_entries", 1),
                    ("WCG_covered_calls", wcg),
                    ("WCG_converted_calls", wcg * converted),
                    ("GCH_covered_calls", gch),
                    ("GCH_converted_calls", gch * converted),
                    ("ambiguous_GCG_template_calls_excluded", 1),
                    (
                        "ambiguous_GCG_converted_template_calls_excluded",
                        converted,
                    ),
                ]
                (
                    output_dir / f"{prefix}.summary.tsv"
                ).write_text(
                    "metric\tvalue\n"
                    + "".join(f"{key}\t{value}\n" for key, value in summary),
                    encoding="utf-8",
                )
                task_rows.append(
                    {
                        "task_index": shard,
                        "sample_index": 1,
                        "sample": "S001",
                        "group": "test",
                        "shard": shard,
                        "input_bam": "x",
                        "reference_fasta": "x",
                        "contig_file": "x",
                        "cell_index_path": "x",
                        "contig_index_path": "x",
                        "output_dir": output_dir,
                        "minimum_mapq": 30,
                        "minimum_baseq": 20,
                        "excluded_bases_each_read_end": 10,
                    }
                )
            with (
                input_root / "dna_single_cytosine_tasks.tsv"
            ).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=task_fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(task_rows)
            old_argv = sys.argv
            try:
                sys.argv = [
                    "68_validate",
                    "--project",
                    str(project),
                    "--sample-index",
                    "1",
                ]
                self.assertEqual(validator.main(), 0)
            finally:
                sys.argv = old_argv
            summary_path = (
                project
                / "08_qc_summary/21_dna_single_cytosine_validation/S001"
                / "S001.single_cytosine.validation_summary.tsv"
            )
            with summary_path.open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                metrics = {
                    row["metric"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(metrics["WCG_covered_calls"], "4")
            self.assertEqual(metrics["GCH_covered_calls"], "4")
            self.assertEqual(
                metrics["pipeline_signature"],
                "snmar_single_cytosine_validation_20260801",
            )
            self.assertEqual(
                metrics["GCG_ambiguous_template_calls_excluded"], "8"
            )
            self.assertEqual(
                metrics[
                    "GCG_ambiguous_converted_template_calls_excluded"
                ],
                "4",
            )


class RnaCompleteMatrixTests(unittest.TestCase):
    def test_matrixmarket_conversion_preserves_gene_by_cell_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "matrix.mtx.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(
                    "%%MatrixMarket matrix coordinate integer general\n"
                    "% test\n"
                    "3 2 2\n"
                    "1 1 4\n"
                    "3 2 5\n"
                )
            self.assertEqual(
                rna.inspect_matrix(path),
                {
                    "features": 3,
                    "cells": 2,
                    "nonzero_entries": 2,
                    "exact_umis": 9,
                },
            )
            coo = rna.read_matrix_coo(path)
            self.assertEqual(coo["feature_index"].tolist(), [0, 2])
            self.assertEqual(coo["cell_index"].tolist(), [0, 1])
            self.assertEqual(coo["exact_umis"].tolist(), [4, 5])


class SingleBaseBigWigTests(unittest.TestCase):
    def test_reference_context_counter_matches_caller_definition(self) -> None:
        sequence = "AACGTTGCAGCGCGATGCTCGAGCCGC"
        expected = {"WCG": 0, "GCH": 0, "GCG": 0}
        names = {
            caller.CONTEXT_WCG: "WCG",
            caller.CONTEXT_GCH: "GCH",
            caller.CONTEXT_GCG_AMBIGUOUS: "GCG",
        }
        for position in range(1, len(sequence) - 1):
            for reverse, read_base in ((False, "C"), (True, "G")):
                result = caller.context_call(
                    reverse, sequence, position, read_base
                )
                if result is not None:
                    expected[names[result[2]]] += 1
        self.assertEqual(
            dict(bigwig.count_contexts_in_block(sequence.encode("ascii"))),
            expected,
        )

    def test_species_context_tracks_are_separated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            dedup_qc = project / "08_qc_summary/15_dna_cb_markdup_full"
            dedup_qc.mkdir(parents=True)
            dedup_summary = dedup_qc / "dna_cb_markdup_full_summary.tsv"
            dedup_summary.write_text(
                "metric\tvalue\n"
                "status\tvalidated\n"
                "pipeline_signature\t"
                "snmar_dna_cb_markdup_corrected_20260731\n",
                encoding="utf-8",
            )
            # Track outputs must be newer than the corrected dedup result.
            import os

            os.utime(dedup_summary, (1, 1))
            input_root = (
                project / "08_qc_summary/20_single_cytosine_inputs"
            )
            input_root.mkdir(parents=True)
            (
                input_root / "dna_single_cytosine_contig_index.tsv"
            ).write_text(
                "contig_index\treference_order\tshard\tspecies\tcontig\t"
                "length\n"
                "0\t1\t1\thuman\thuman_1\t100\n"
                "1\t2\t2\tpig\tpig_1\t100\n",
                encoding="utf-8",
            )
            reference = (
                project
                / "02_reference/rna_human_pig_joint/"
                "human_pig.prefixed.fa"
            )
            reference.parent.mkdir(parents=True)
            reference.write_text(
                ">human_1\n"
                + ("ACGGCA" * 16 + "ACGG")[:100]
                + "\n>pig_1\n"
                + ("TCGAGC" * 16 + "TCGA")[:100]
                + "\n",
                encoding="ascii",
            )
            tracks = (
                project
                / "07_postmapping/complete_matrices/"
                "dna_single_cytosine/tracks"
            )
            track_qc = (
                project / "08_qc_summary/22_dna_single_cytosine_tracks"
            )
            tracks.mkdir(parents=True)
            track_qc.mkdir(parents=True)
            rows_by_shard = {
                1: [
                    (
                        0,
                        "human",
                        "human_1",
                        10,
                        "+",
                        "WCG",
                        "endogenous_methylation",
                        2,
                        3,
                        2,
                    ),
                    (
                        0,
                        "human",
                        "human_1",
                        20,
                        "-",
                        "GCH",
                        "GpC_accessibility",
                        1,
                        2,
                        1,
                    ),
                ],
                2: [
                    (
                        1,
                        "pig",
                        "pig_1",
                        5,
                        "+",
                        "WCG",
                        "endogenous_methylation",
                        2,
                        4,
                        1,
                    ),
                    (
                        1,
                        "pig",
                        "pig_1",
                        15,
                        "-",
                        "GCH",
                        "GpC_accessibility",
                        1,
                        1,
                        1,
                    ),
                ],
            }
            for shard in range(1, 9):
                detailed = (
                    tracks
                    / "dna_single_cytosine.pseudobulk."
                    f"shard_{shard:02d}.tsv.gz"
                )
                with gzip.open(
                    detailed, "wt", encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=bigwig.FIELDS,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    for (
                        contig_index,
                        species,
                        chrom,
                        start,
                        strand,
                        context,
                        signal_type,
                        cells,
                        covered,
                        converted,
                    ) in rows_by_shard.get(shard, []):
                        writer.writerow(
                            {
                                "contig_index": contig_index,
                                "species": species,
                                "chrom": chrom,
                                "start": start,
                                "end": start + 1,
                                "strand": strand,
                                "context": context,
                                "signal_type": signal_type,
                                "covered_cells": cells,
                                "covered_calls": covered,
                                "converted_calls": converted,
                                "signal_fraction": converted / covered,
                            }
                        )
                (
                    track_qc
                    / "dna_single_cytosine."
                    f"track_shard_{shard:02d}.summary.tsv"
                ).write_text(
                    "metric\tvalue\n"
                    "status\tvalidated\n"
                    f"shard\t{shard}\n"
                    f"detailed_track\t{detailed}\n",
                    encoding="utf-8",
                )

            def fake_converter(command, check):
                self.assertTrue(check)
                Path(command[3]).write_bytes(b"fake-bigwig")
                return mock.Mock(returncode=0)

            old_argv = sys.argv
            try:
                sys.argv = ["81_bigwig", "--project", str(project)]
                with (
                    mock.patch.object(
                        bigwig.shutil, "which", return_value="fake"
                    ),
                    mock.patch.object(
                        bigwig.subprocess,
                        "run",
                        side_effect=fake_converter,
                    ),
                ):
                    self.assertEqual(bigwig.main(), 0)
            finally:
                sys.argv = old_argv
            output = (
                project
                / "07_postmapping/complete_matrices/"
                "dna_single_cytosine/bigwig"
            )
            for species in ("human", "pig"):
                for context, label in (
                    ("WCG", "methylation"),
                    ("GCH", "accessibility"),
                ):
                    self.assertTrue(
                        (
                            output
                            / f"{species}.{context}_{label}_fraction.bw"
                        ).is_file()
                    )
                    self.assertTrue(
                        (
                            output
                            / f"{species}.{context}_{label}_coverage.bw"
                        ).is_file()
                    )
                    self.assertTrue(
                        (
                            output
                            / f"{species}.{context}_{label}_covered_cells.bw"
                        ).is_file()
                    )
                    self.assertTrue(
                        (
                            output
                            / f"{species}.{context}_{label}."
                            "per_base.tsv.gz"
                        ).is_file()
                    )
            summary = bigwig.read_metrics(
                project
                / "08_qc_summary/24_dna_single_cytosine_bigwig/"
                "dna_single_base_bigwig_summary.tsv"
            )
            for species in ("human", "pig"):
                for context in ("WCG", "GCH"):
                    reference_sites = int(
                        summary[
                            f"{species}_{context}_reference_assayable_sites"
                        ]
                    )
                    self.assertGreater(reference_sites, 0)
                    self.assertAlmostEqual(
                        float(
                            summary[
                                f"{species}_{context}_"
                                "observed_site_coverage_fraction"
                            ]
                        ),
                        1 / reference_sites,
                    )


class SingleCellCoverageQCTests(unittest.TestCase):
    def test_species_reference_coverage_is_added_without_filtering(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            input_path = project / "cells.tsv.gz"
            fields = [
                "cell_id",
                "WCG_observed_sites",
                "WCG_covered_calls",
                "WCG_methylation_fraction",
                "GCH_observed_sites",
                "GCH_covered_calls",
                "GCH_accessibility_fraction",
                "human_signal_covered_calls",
                "pig_signal_covered_calls",
                "dominant_reference_signal",
                "human_WCG_observed_sites",
                "human_WCG_covered_calls",
                "human_GCH_observed_sites",
                "human_GCH_covered_calls",
                "pig_WCG_observed_sites",
                "pig_WCG_covered_calls",
                "pig_GCH_observed_sites",
                "pig_GCH_covered_calls",
            ]
            with gzip.open(
                input_path, "wt", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fields,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "cell_id": "S001|CB1",
                        "WCG_observed_sites": 2,
                        "WCG_covered_calls": 3,
                        "WCG_methylation_fraction": 2 / 3,
                        "GCH_observed_sites": 1,
                        "GCH_covered_calls": 1,
                        "GCH_accessibility_fraction": 1,
                        "human_signal_covered_calls": 4,
                        "pig_signal_covered_calls": 0,
                        "dominant_reference_signal": "human",
                        "human_WCG_observed_sites": 2,
                        "human_WCG_covered_calls": 3,
                        "human_GCH_observed_sites": 1,
                        "human_GCH_covered_calls": 1,
                        "pig_WCG_observed_sites": 0,
                        "pig_WCG_covered_calls": 0,
                        "pig_GCH_observed_sites": 0,
                        "pig_GCH_covered_calls": 0,
                    }
                )
            bigwig_summary = project / "bigwig.tsv"
            bigwig_summary.write_text(
                "metric\tvalue\n"
                "status\tvalidated\n"
                "human_WCG_reference_assayable_sites\t100\n"
                "human_GCH_reference_assayable_sites\t200\n"
                "pig_WCG_reference_assayable_sites\t50\n"
                "pig_GCH_reference_assayable_sites\t80\n",
                encoding="utf-8",
            )
            old_argv = sys.argv
            try:
                sys.argv = [
                    "83_cell_qc",
                    "--project",
                    str(project),
                    "--input-per-cell",
                    str(input_path),
                    "--bigwig-summary",
                    str(bigwig_summary),
                ]
                self.assertEqual(cell_qc.main(), 0)
            finally:
                sys.argv = old_argv
            output_root = (
                project / "08_qc_summary/25_dna_single_cell_qc"
            )
            with gzip.open(
                output_root / "dna_single_cell_coverage_qc.tsv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertAlmostEqual(
                float(row["human_WCG_reference_coverage_fraction"]),
                0.02,
            )
            self.assertAlmostEqual(
                float(row["human_GCH_reference_coverage_fraction"]),
                0.005,
            )
            self.assertEqual(row["major_reference_signal_fraction"], "1.0")
            metrics = cell_qc.read_metrics(
                output_root / "dna_single_cell_coverage_qc_summary.tsv"
            )
            self.assertEqual(metrics["cells"], "1")
            self.assertEqual(metrics["cell_filtering"], "none")


class ScriptInventoryTests(unittest.TestCase):
    def test_force_is_propagated_to_all_full_rerun_jobs(self) -> None:
        pairs = (
            (
                "66_submit_dna_single_cytosine_task_template_8core.csub",
                "67_submit_dna_single_cytosine_768jobs.sh",
            ),
            (
                "69_submit_dna_single_cytosine_validation_template_8core.csub",
                "70_submit_dna_single_cytosine_validation_96jobs.sh",
            ),
            (
                "72_submit_single_cytosine_track_template_8core.csub",
                "73_submit_single_cytosine_track_8jobs.sh",
            ),
        )
        for template_name, launcher_name in pairs:
            template = script(template_name).read_text(encoding="utf-8")
            launcher = script(launcher_name).read_text(encoding="utf-8")
            self.assertIn('FORCE="${FORCE:-__FORCE__}"', template)
            self.assertIn('FORCE="${FORCE:-0}"', launcher)
            self.assertIn('s/__FORCE__/$FORCE/g', launcher)
            self.assertIn("if (( FORCE == 0 ))", launcher)

    def test_bigwig_submission_defaults_to_foodie_environment(self) -> None:
        submission = (
            script("82_submit_dna_single_base_bigwigs_16core.csub")
        ).read_text(encoding="utf-8")
        self.assertIn('BIGWIG_ENV="${BIGWIG_ENV:-foodie_env}"', submission)
        self.assertIn("bedGraphToBigWig", submission)

    def test_v19_numbered_steps_are_complete(self) -> None:
        expected = {
            63: "63_prepare_single_cytosine_full.py",
            64: "64_submit_single_cytosine_prepare_8core.csub",
            65: "65_call_dna_single_cytosine_task.py",
            66: "66_submit_dna_single_cytosine_task_template_8core.csub",
            67: "67_submit_dna_single_cytosine_768jobs.sh",
            68: "68_validate_dna_single_cytosine_sample.py",
            69: "69_submit_dna_single_cytosine_validation_template_8core.csub",
            70: "70_submit_dna_single_cytosine_validation_96jobs.sh",
            71: "71_aggregate_single_cytosine_track_shard.py",
            72: "72_submit_single_cytosine_track_template_8core.csub",
            73: "73_submit_single_cytosine_track_8jobs.sh",
            74: "74_summarize_dna_single_cytosine_full.py",
            75: "75_submit_dna_single_cytosine_summary_8core.csub",
            76: "76_validate_rna_complete_matrix.py",
            77: "77_submit_rna_complete_matrix_8core.csub",
            78: "78_summarize_complete_modalities.py",
            79: "79_submit_complete_modalities_summary_4core.csub",
            80: "80_check_complete_matrix_progress.py",
            81: "81_build_dna_single_base_bigwigs.py",
            82: "82_submit_dna_single_base_bigwigs_16core.csub",
            83: "83_summarize_dna_single_cell_qc.py",
            84: "84_submit_dna_single_cell_qc_8core.csub",
        }
        for number, filename in expected.items():
            self.assertTrue(script(filename).is_file(), number)


if __name__ == "__main__":
    unittest.main()
