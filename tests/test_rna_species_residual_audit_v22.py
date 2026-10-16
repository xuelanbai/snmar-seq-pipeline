from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = (
    ROOT
    / "03_scripts/05_downstream_analysis/03_cell_qc_filtering"
)
ALIGN_SCRIPT = SCRIPT_DIR / "98_audit_rna_species_alignment_one.py"
SUMMARIZE_SCRIPT = SCRIPT_DIR / "101_summarize_rna_species_residual_audit.py"
TEMPLATE = SCRIPT_DIR / "99_submit_rna_species_alignment_audit_template_32core.csub"
LAUNCHER = SCRIPT_DIR / "100_submit_rna_species_alignment_audit_4jobs.sh"
SUMMARY_CSUB = SCRIPT_DIR / "102_submit_rna_species_residual_summary_4core.csub"


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
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


def sam_record(
    qname: str,
    flag: int,
    contig: str,
    cigar: str,
    nh: int = 1,
    include_tags: bool = True,
) -> str:
    fields = [
        qname,
        str(flag),
        contig,
        "1",
        "60",
        cigar,
        "*",
        "0",
        "0",
        "A" * 50,
        "I" * 50,
        f"NH:i:{nh}",
    ]
    if include_tags:
        fields.extend(
            [
                "CB:Z:" + "A" * 36,
                "CN:Z:cell",
                "UB:Z:ACGTACGT",
                "XD:i:0",
            ]
        )
    return "\t".join(fields) + "\n"


class AlignmentAuditTests(unittest.TestCase):
    def test_species_uniqueness_and_splice_counters(self) -> None:
        module = load(ALIGN_SCRIPT, "rna_species_alignment_audit")
        lines = [
            sam_record("h1", 0, "human_1", "25M100N25M"),
            sam_record("p1", 0, "pig_1", "50M"),
            sam_record("hmulti", 0, "human_2", "50M", nh=2),
            sam_record("unmapped", 4, "*", "*"),
            sam_record("secondary", 0x100, "pig_2", "50M"),
        ]
        counts = module.audit_sam_lines(lines)
        self.assertEqual(counts["primary_records"], 4)
        self.assertEqual(counts["primary_mapped"], 3)
        self.assertEqual(counts["unique_primary"], 2)
        self.assertEqual(counts["nonunique_primary"], 1)
        self.assertEqual(counts["eligible_spliced_human"], 1)
        self.assertEqual(counts["eligible_unspliced_pig"], 1)


class SummaryAuditTests(unittest.TestCase):
    def test_two_well_layer_reconciliation(self) -> None:
        module = load(SUMMARIZE_SCRIPT, "rna_species_residual_summary")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            per_file = root / "per_file.tsv"
            cell_qc = root / "cell_qc.tsv.gz"
            demux = root / "demux.tsv"
            alignment_root = root / "alignment"
            genefull_root = root / "genefull"
            output = root / "out"
            manifest_rows = [
                {
                    "manifest_index": 1,
                    "shard": 1,
                    "shard_order": 1,
                    "pair_id": "pair1",
                    "sample": "RNA-1",
                    "star_input_reads": 110,
                    "primary_mapped_reads": 100,
                    "primary_mapping_percent": 90.9,
                    "input_bam": "/tmp/one.bam",
                },
                {
                    "manifest_index": 2,
                    "shard": 1,
                    "shard_order": 2,
                    "pair_id": "pair2",
                    "sample": "RNA-2",
                    "star_input_reads": 105,
                    "primary_mapped_reads": 100,
                    "primary_mapping_percent": 95.2,
                    "input_bam": "/tmp/two.bam",
                },
            ]
            write_tsv(manifest, list(manifest_rows[0]), manifest_rows)
            per_file_rows = [
                {
                    "pair_id": "pair1",
                    "sample": "RNA-1",
                    "star_input_reads": 110,
                    "primary_mapped_reads": 100,
                    "unique_primary_reads": 80,
                    "tag_eligible_unique_primary": 80,
                    "exact_umis": 55,
                    "human_exact_umis": 30,
                    "pig_exact_umis": 25,
                },
                {
                    "pair_id": "pair2",
                    "sample": "RNA-2",
                    "star_input_reads": 105,
                    "primary_mapped_reads": 100,
                    "unique_primary_reads": 80,
                    "tag_eligible_unique_primary": 80,
                    "exact_umis": 55,
                    "human_exact_umis": 5,
                    "pig_exact_umis": 50,
                },
            ]
            write_tsv(per_file, list(per_file_rows[0]), per_file_rows)

            audit_rows = [
                {
                    "status": "validated",
                    "pipeline_signature": module.ALIGNMENT_SIGNATURE,
                    "manifest_index": 1,
                    "shard": 1,
                    "pair_id": "pair1",
                    "sample": "RNA-1",
                    "well_index": 1,
                    "expected_well_set": "GM12878",
                    "primary_mapped": 100,
                    "primary_human": 60,
                    "primary_pig": 40,
                    "primary_human_fraction": 0.6,
                    "unique_primary": 80,
                    "unique_human": 50,
                    "unique_pig": 30,
                    "unique_human_fraction": 0.625,
                    "tag_eligible_unique_primary": 80,
                    "eligible_human": 50,
                    "eligible_pig": 30,
                    "eligible_spliced_human": 30,
                    "eligible_spliced_pig": 3,
                    "eligible_human_spliced_fraction": 0.6,
                    "eligible_pig_spliced_fraction": 0.1,
                },
                {
                    "status": "validated",
                    "pipeline_signature": module.ALIGNMENT_SIGNATURE,
                    "manifest_index": 2,
                    "shard": 1,
                    "pair_id": "pair2",
                    "sample": "RNA-2",
                    "well_index": 2,
                    "expected_well_set": "cellMix",
                    "primary_mapped": 100,
                    "primary_human": 20,
                    "primary_pig": 80,
                    "primary_human_fraction": 0.2,
                    "unique_primary": 80,
                    "unique_human": 10,
                    "unique_pig": 70,
                    "unique_human_fraction": 0.125,
                    "tag_eligible_unique_primary": 80,
                    "eligible_human": 10,
                    "eligible_pig": 70,
                    "eligible_spliced_human": 5,
                    "eligible_spliced_pig": 35,
                    "eligible_human_spliced_fraction": 0.5,
                    "eligible_pig_spliced_fraction": 0.5,
                },
            ]
            for row in audit_rows:
                write_tsv(
                    alignment_root
                    / f"{row['pair_id']}.rna_species_alignment_audit.tsv",
                    list(row),
                    [row],
                )

            for pair_id, human_reads, pig_reads in (
                ("pair1", 20, 40),
                ("pair2", 5, 50),
            ):
                path = genefull_root / pair_id / f"{pair_id}.summary.json"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(
                    json.dumps(
                        {
                            "status": "validated",
                            "pair_id": pair_id,
                            "counters": {
                                "gene_assigned": human_reads + pig_reads,
                                "human_assigned_reads": human_reads,
                                "pig_assigned_reads": pig_reads,
                            },
                        }
                    ),
                    encoding="utf-8",
                )

            demux_rows = []
            for pair_id, sample in (("pair1", "RNA-1"), ("pair2", "RNA-2")):
                for metric, count in (
                    ("total_pairs", 1000),
                    ("strict_assigned", 400),
                    ("DNA_residual_GTCA_linker_given_T7", 300),
                ):
                    demux_rows.append(
                        {
                            "pair_id": pair_id,
                            "sample": sample,
                            "metric": metric,
                            "count": count,
                            "percent": count / 10,
                        }
                    )
            write_tsv(demux, list(demux_rows[0]), demux_rows)

            cell_rows = [
                {
                    "sample": "RNA-1",
                    "CB": "A" * 36,
                    "exact_umis": 30,
                    "genes_detected": 12,
                    "human_exact_umis": 5,
                    "pig_exact_umis": 25,
                    "in_raw_matrix": 1,
                },
                {
                    "sample": "RNA-1",
                    "CB": "C" * 36,
                    "exact_umis": 25,
                    "genes_detected": 10,
                    "human_exact_umis": 25,
                    "pig_exact_umis": 0,
                    "in_raw_matrix": 1,
                },
                {
                    "sample": "RNA-2",
                    "CB": "G" * 36,
                    "exact_umis": 55,
                    "genes_detected": 20,
                    "human_exact_umis": 5,
                    "pig_exact_umis": 50,
                    "in_raw_matrix": 1,
                },
            ]
            write_tsv(cell_qc, list(cell_rows[0]), cell_rows)

            old_argv = sys.argv
            try:
                sys.argv = [
                    str(SUMMARIZE_SCRIPT),
                    "--manifest",
                    str(manifest),
                    "--alignment-audit-root",
                    str(alignment_root),
                    "--rna-per-file",
                    str(per_file),
                    "--rna-cell-qc",
                    str(cell_qc),
                    "--demux-per-pair",
                    str(demux),
                    "--genefull-by-file-root",
                    str(genefull_root),
                    "--output-root",
                    str(output),
                    "--expected-wells",
                    "2",
                    "--human-only-wells",
                    "1",
                ]
                self.assertEqual(module.main(), 0)
            finally:
                sys.argv = old_argv

            per_set = {
                row["expected_well_set"]: row
                for row in read_tsv(output / "rna_species_residual_audit_per_set.tsv")
            }
            self.assertEqual(per_set["GM12878"]["unique_human_fraction"], "0.625")
            self.assertEqual(
                per_set["GM12878"]["UMI20_gene10_human_cells"], "1"
            )
            self.assertEqual(
                per_set["GM12878"]["UMI20_gene10_pig_cells"], "1"
            )
            summary = {
                row["metric"]: row["value"]
                for row in read_tsv(output / "rna_species_residual_audit_summary.tsv")
            }
            self.assertEqual(summary["status"], "validated")
            self.assertEqual(summary["biological_species_validation"], "requires_review")
            self.assertEqual(
                summary["known_DNA_residual_GTCA_linker_fraction_of_raw_pairs"],
                "0.3",
            )

    def test_cluster_scripts_are_guarded_and_parallel(self) -> None:
        template = TEMPLATE.read_text(encoding="utf-8")
        launcher = LAUNCHER.read_text(encoding="utf-8")
        summary = SUMMARY_CSUB.read_text(encoding="utf-8")
        self.assertIn("#CSUB -n 32", template)
        self.assertIn("__SHARD_INDEX__", template)
        self.assertIn("xargs -r -P", template)
        self.assertIn("CONFIRM_SUBMIT_RNA_SPECIES_AUDIT", launcher)
        self.assertIn("[DRY RUN]", launcher)
        self.assertIn("expected 96 alignment audit files", summary)


if __name__ == "__main__":
    unittest.main()
