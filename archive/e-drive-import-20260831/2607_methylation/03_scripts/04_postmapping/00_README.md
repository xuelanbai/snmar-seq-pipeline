# 04 post-mapping and direct full processing

This directory starts from the validated full BAMs. It does not remap reads.

## Why a Probe is required

- DNA BAMs contain many cells. Ordinary whole-BAM duplicate marking would
  compare molecules from different cells. The DNA Probe therefore runs
  `samtools markdup --barcode-tag CB`.
- RNA already has `CB` and `UB` tags in the BAM. The RNA Probe counts exact
  `(sample, CB, gene, UB)` molecules from unique primary alignments.
- snMAR RNA is single-nucleus RNA. Gene assignment therefore uses complete
  gene bodies (exons plus introns; STARsolo `GeneFull` semantics).
- RNA strandedness is not assumed. The Probe reports forward, reverse, and
  unstranded assignments separately.

## Steps

```text
00_cleanup_obsolete_matrix_code.sh
33_prepare_postmapping_probe.py
34_check_postmapping_requirements.sh
35_run_dna_cb_markdup_one.sh
36_submit_dna_cb_markdup_probe_template_64core.csub
37_submit_dna_cb_markdup_probe_4jobs.sh
38_summarize_dna_cb_markdup_probe.py
39_count_rna_genefull_umi_probe.py
40_summarize_rna_genefull_umi_probe.py
41_submit_rna_genefull_umi_probe_32core.csub
42_prepare_postmapping_full_manifests.py
43_submit_dna_cb_markdup_full_template_64core.csub
44_submit_dna_cb_markdup_full_96jobs.sh
45_summarize_dna_cb_markdup_full.py
46_count_rna_reverse_genefull_one.py
47_submit_rna_reverse_genefull_full_template_32core.csub
48_submit_rna_reverse_genefull_full_4jobs.sh
49_merge_rna_reverse_genefull_matrix.py
50_submit_rna_reverse_matrix_merge_8core.csub
51_submit_dna_cb_markdup_full_summary_8core.csub
52_prepare_and_submit_full_dedup.sh
53_prepare_downstream_full.py
54_call_dna_taps_cpg_full_shard.py
55_merge_dna_taps_cpg_full_sample.py
56_submit_dna_taps_cpg_full_template_32core.csub
57_submit_dna_taps_cpg_full_96jobs.sh
58_summarize_dna_taps_cpg_full.py
59_submit_dna_taps_cpg_full_summary_8core.csub
60_filter_rna_reverse_matrix_full.py
61_submit_rna_filtered_matrix_full_8core.csub
62_prepare_and_submit_downstream_full.sh
63_prepare_single_cytosine_full.py
64_submit_single_cytosine_prepare_8core.csub
65_call_dna_single_cytosine_task.py
66_submit_dna_single_cytosine_task_template_8core.csub
67_submit_dna_single_cytosine_768jobs.sh
68_validate_dna_single_cytosine_sample.py
69_submit_dna_single_cytosine_validation_template_8core.csub
70_submit_dna_single_cytosine_validation_96jobs.sh
71_aggregate_single_cytosine_track_shard.py
72_submit_single_cytosine_track_template_8core.csub
73_submit_single_cytosine_track_8jobs.sh
74_summarize_dna_single_cytosine_full.py
75_submit_dna_single_cytosine_summary_8core.csub
76_validate_rna_complete_matrix.py
77_submit_rna_complete_matrix_8core.csub
78_summarize_complete_modalities.py
79_submit_complete_modalities_summary_4core.csub
80_check_complete_matrix_progress.py
81_build_dna_single_base_bigwigs.py
82_submit_dna_single_base_bigwigs_16core.csub
83_summarize_dna_single_cell_qc.py
84_submit_dna_single_cell_qc_8core.csub
```

## Prepare

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate snTAPS_align

python 03_scripts/04_postmapping/33_prepare_postmapping_probe.py
bash 03_scripts/04_postmapping/34_check_postmapping_requirements.sh
```

Step 33 validates the 96 DNA and 96 RNA full mapping rows and writes
deterministic Probe manifests under:

```text
08_qc_summary/13_postmapping_probe/
```

The Python workers read BAM records through `samtools view`; no `pysam`
installation and no new conda environment are required.

## Submit DNA Probe

Preview the four independent 64-CPU jobs:

```bash
bash 03_scripts/04_postmapping/37_submit_dna_cb_markdup_probe_4jobs.sh
```

Submit:

```bash
CONFIRM_SUBMIT_DNA_CB_MARKDUP=1 \
  bash 03_scripts/04_postmapping/37_submit_dna_cb_markdup_probe_4jobs.sh
```

After all four jobs finish:

```bash
python 03_scripts/04_postmapping/38_summarize_dna_cb_markdup_probe.py
```

The worker performs:

```text
query-name sort
  -> samtools fixmate -m
  -> coordinate sort
  -> samtools markdup -c -m t --barcode-tag CB
  -> filter records carrying duplicate flag 0x400
  -> BAM index and validation
```

DNA has no molecular UMI. Duplicate identity is
`CB + reference + unclipped paired-template endpoints + orientation`; read
length is not used as an independent duplicate key. The corrected command
does not pass `markdup -l 1`: for `markdup`, `-l` is expected maximum read
length, not BAM compression level. The Probe sets `KEEP_MARKED_BAM=1` so its
pre/post-dedup BAMs can be compared. Full jobs retain only the deduplicated
BAM to avoid duplicating terabytes of storage.

## Submit RNA Probe

```bash
csub < 03_scripts/04_postmapping/41_submit_rna_genefull_umi_probe_32core.csub
```

The Probe processes eight representative BAMs and reports three modes:

```text
forward
reverse
unstranded
```

Only unique primary alignments (`NH:i:1`) with strict `XD:i:0`, `CB`, and
`UB` are eligible. Exact UMIs are deduplicated within sample, cell, and gene.
No barcode correction or one-base UMI rescue is performed in this Probe.

## Decision gate before full processing

Do not write or submit full DNA deduplication or full RNA count-matrix jobs
until both Probe summaries have been reviewed. In particular:

- inspect the DNA CB-aware duplicate fraction and tag retention;
- choose RNA forward or reverse only if it clearly outperforms the other;
- retain unstranded only if the protocol evidence and Probe results justify it.

## Corrected DNA run mode

The earlier 97.36% DNA result was produced with the invalid
`samtools markdup -l 1` option and is not an accepted decision. The current
decision is to skip a new Probe and run the corrected worker directly on all
96 mapped samples. Step 44 defaults to `DIRECT_DNA_FULL=1`. The worker
requires pipeline signature `snmar_dna_cb_markdup_corrected_20260731`; old
outputs do not satisfy this signature and are replaced only after a new BAM
has passed `samtools quickcheck`.

RNA:

```text
strand                             reverse
feature                            GeneFull gene body
UMI deduplication                  exact within sample+CB+gene
reverse assigned fraction         24.651184%
forward assigned fraction         18.729669%
reverse exact UMIs                48,971
forward exact UMIs                15,424
```

## Steps 42-52: full DNA/RNA deduplication

Step 42 validates both Probe summaries and prepares:

```text
08_qc_summary/14_postmapping_full_inputs/
  dna_cb_markdup_full_manifest.tsv
  rna_reverse_genefull_full_manifest.tsv
  postmapping_full_decisions.tsv
```

The four RNA shards each contain 24 files and are greedily balanced by STAR
input-read count.

### Combined dry run and first submission

The combined launcher defaults to DNA samples 1-24 plus all four RNA shards:

```bash
bash 03_scripts/04_postmapping/52_prepare_and_submit_full_dedup.sh

CONFIRM_SUBMIT_FULL_DEDUP=1 \
  bash 03_scripts/04_postmapping/52_prepare_and_submit_full_dedup.sh
```

DNA uses one independent 64-CPU job per sample. RNA uses four independent
32-CPU jobs, each processing 24 files with eight files in parallel.

### Remaining DNA ranges

Do not rerun step 52 merely to submit later DNA ranges, because step 52 also
submits the RNA shards. Use step 44 directly:

Submit all 96 corrected DNA jobs directly, without a Probe:

```bash
START_INDEX=1 END_INDEX=96 DIRECT_DNA_FULL=1 \
CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=1 \
  bash 03_scripts/04_postmapping/44_submit_dna_cb_markdup_full_96jobs.sh
```

This submits 96 independent 64-CPU jobs. After step 51 validates the
corrected signature, do not call step 62 because it also submits RNA.
Regenerate the DNA manifest with step 53 and force one DNA-only downstream
refresh:

```bash
python 03_scripts/04_postmapping/53_prepare_downstream_full.py

START_INDEX=1 END_INDEX=96 FORCE=1 \
CONFIRM_SUBMIT_DNA_TAPS_FULL=1 \
  bash 03_scripts/04_postmapping/57_submit_dna_taps_cpg_full_96jobs.sh
```

`FORCE=1` is rendered into each CSUB file, so it remains effective even if
the cluster submission command does not export the caller environment. This
one forced refresh is required because the corrected BAM deliberately
replaces the old BAM at the same pathname.

```bash
START_INDEX=25 END_INDEX=48 \
CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=1 \
  bash 03_scripts/04_postmapping/44_submit_dna_cb_markdup_full_96jobs.sh

START_INDEX=49 END_INDEX=72 \
CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=1 \
  bash 03_scripts/04_postmapping/44_submit_dna_cb_markdup_full_96jobs.sh

START_INDEX=73 END_INDEX=96 \
CONFIRM_SUBMIT_DNA_CB_MARKDUP_FULL=1 \
  bash 03_scripts/04_postmapping/44_submit_dna_cb_markdup_full_96jobs.sh
```

The DNA worker skips an existing completed BAM only when its BAM, BAI, stats,
and done manifest all exist and `samtools quickcheck` passes.

### Final summaries

After all four RNA shards finish:

```bash
csub < 03_scripts/04_postmapping/50_submit_rna_reverse_matrix_merge_8core.csub
```

This creates an unfiltered raw matrix:

```text
07_postmapping/rna_reverse_genefull_full/raw_feature_bc_matrix/
  matrix.mtx.gz
  barcodes.tsv.gz
  features.tsv.gz
  feature_metadata.tsv.gz
08_qc_summary/16_rna_reverse_genefull_full/
  rna_reverse_cell_qc.tsv.gz
  rna_reverse_genefull_full_per_file.tsv
  rna_reverse_genefull_full_summary.tsv
```

No cell filtering is applied at this stage.

After all 96 DNA jobs finish:

```bash
csub < 03_scripts/04_postmapping/51_submit_dna_cb_markdup_full_summary_8core.csub
```

DNA full outputs:

```text
07_postmapping/dna_cb_dedup_full/<sample>/<sample>.cb_dedup.bam
08_qc_summary/15_dna_cb_markdup_full/
  dna_cb_markdup_full_per_sample.tsv
  dna_cb_markdup_full_summary.tsv
```

## Steps 53-62: direct full DNA/RNA processing

Step 53 revalidates all 96 full DNA dedup outputs and all 96 full RNA count
outputs. It requires exact total reconciliation before writing:

```text
08_qc_summary/17_downstream_full_inputs/
  dna_taps_cpg_full_manifest.tsv
  dna_contig_shards.tsv
  dna_contig_shards/shard_01..08.contigs.txt
  downstream_full_decisions.tsv
```

### DNA full TAPS-CpG calling

Step 54 implements strand-directional TAPS calls:

```text
forward alignment + reference C in CpG: C/T
reverse alignment + reference G in CpG: G/A
```

The reference G observation is mapped back to the canonical reference-C
coordinate. Paired reads overlapping the same CpG are collapsed within
`sample + CB + QNAME + CpG`; the higher-base-quality observation wins and an
equal-quality conflict is discarded. This is deliberately different from
scanning both C and G positions on every alignment, which would double count
one CpG and bias the rate downward.

Every prefixed human and pig reference contig is assigned exactly once.
Step 57 submits 96 sample jobs. Each job requests 32 CPUs and runs eight
contig shards in parallel with four samtools threads per shard. Intermediate
shards are deleted only after their counts reconcile and the final sample
files have been written.

### RNA full filtered matrix

Step 60 reads:

```text
08_qc_summary/16_rna_reverse_genefull_full/rna_reverse_cell_qc.tsv.gz
```

It filters every raw-matrix column using the default reversible rule:

```text
exact UMIs >= 10 AND genes detected >= 5
```

The raw matrix is retained. `MIN_EXACT_UMIS` and `MIN_GENES` can be changed
when step 61 is resubmitted.

### Preview and submit

```bash
bash 03_scripts/04_postmapping/62_prepare_and_submit_downstream_full.sh

CONFIRM_SUBMIT_DOWNSTREAM_FULL=1 \
  bash 03_scripts/04_postmapping/62_prepare_and_submit_downstream_full.sh
```

Review:

```text
08_qc_summary/18_dna_taps_cpg_full/dna_taps_cpg_full_summary.tsv
08_qc_summary/18_dna_taps_cpg_full/dna_taps_cpg_full_per_sample.tsv
08_qc_summary/19_rna_filtered_matrix_full/minUMI10_minGene5/
  rna_filtered_matrix_full_summary.tsv
  rna_filtered_matrix_full_per_sample.tsv
```

After all 96 DNA sample jobs finish, submit step 59 to validate and summarize
the full DNA outputs.

## Steps 63-84: complete matrices, single-base tracks and per-cell QC

The upstream DNA product is an observed single-cytosine sparse matrix with
explicit strand and sequence context. A 100-kb bin is not an upstream axis
and no genomic bases are discarded to force fixed-width windows. Downstream
100-kb, TSS, gene-body and other matrices are derived separately from the
WCG or GCH signal, never by mixing both.

The context model is:

```text
WCG (W=A/T): endogenous methylation
GCH (H=A/C/T): GpC-enzyme accessibility
GCG: ambiguous CpG/GpC overlap, counted for audit and excluded
```

Forward-strand calls store the reference C position and `C/T`; reverse-strand
calls store the reference G position, `strand=-1`, and `G/A`. Converted
`T/A` fractions are labelled WCG methylation or GCH accessibility according
to context. Calls require MAPQ >= 30 and baseQ >= 20. Ten query bases are
excluded at each read end. Overlapping mate observations are collapsed
within `cell + QNAME + coordinate + strand + context`; highest baseQ wins
and an equal-quality conflict is discarded. Missing cell-site coordinates
are `NA/no coverage`, never 0%.

Verify the matrix dependencies already installed in the existing environment:

```bash
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda activate snTAPS_align
python -c 'import numpy, h5py; print(numpy.__version__, h5py.__version__)'
```

Prepare the cell axis from all 768 formal demultiplex barcode-count files,
build the reference-contig index, and create 768 task rows:

```bash
csub < 03_scripts/04_postmapping/64_submit_single_cytosine_prepare_8core.csub
```

After step 64 finishes, preview and submit the 96 samples x 8 contig shards.
Submission ranges are resumable; validated tasks are skipped:

```bash
bash 03_scripts/04_postmapping/67_submit_dna_single_cytosine_768jobs.sh

START_TASK=1 END_TASK=192 CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE=1 \
  bash 03_scripts/04_postmapping/67_submit_dna_single_cytosine_768jobs.sh
```

When rebuilding from the corrected deduplication outputs, add `FORCE=1` to
all four task ranges, for example:

```bash
START_TASK=1 END_TASK=192 FORCE=1 \
CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE=1 \
  bash 03_scripts/04_postmapping/67_submit_dna_single_cytosine_768jobs.sh
```

Repeat ranges 193-384, 385-576, and 577-768. Each job requests 8 CPUs and
uses four samtools threads. After all 768 jobs finish, submit 96 sample
validators:

```bash
CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_VALIDATION=1 \
  bash 03_scripts/04_postmapping/70_submit_dna_single_cytosine_validation_96jobs.sh
```

Use `FORCE=1` here as well for the first corrected rebuild.

After all 96 validations finish, submit eight pseudo-bulk track mergers:

```bash
CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_TRACKS=1 \
  bash 03_scripts/04_postmapping/73_submit_single_cytosine_track_8jobs.sh
```

Use `FORCE=1` here as well for the first corrected rebuild. Later resumptions
should return to the default `FORCE=0`, which skips validated outputs.

Step 75 is not a replacement for steps 70 and 73. If it reports a missing
`*.single_cytosine.validation_summary.tsv`, inspect progress first:

```bash
python 03_scripts/04_postmapping/80_check_complete_matrix_progress.py
```

If validation is below `96/96`, regenerate all current-version validation
tables from the completed caller shards:

```bash
FORCE=1 CONFIRM_SUBMIT_DNA_SINGLE_CYTOSINE_VALIDATION=1 \
  bash 03_scripts/04_postmapping/70_submit_dna_single_cytosine_validation_96jobs.sh
```

Wait for `96/96`, then run step 73. After the track stage reaches `8/8`,
steps 75 and 82 may be submitted in parallel. Step 75 writes
`dna_single_cytosine_validation_preflight.tsv` and reports all missing or
stale samples together instead of stopping at the first path.

Then build the small HDF5 virtual dataset spanning all 768 physical shards:

```bash
csub < 03_scripts/04_postmapping/75_submit_dna_single_cytosine_summary_8core.csub
```

Primary DNA products:

```text
07_postmapping/complete_matrices/dna_single_cytosine/
  shards/<sample>/shard_01..08/*.single_cytosine.h5
  tracks/dna_single_cytosine.pseudobulk.shard_*.tsv.gz
  tracks/dna_WCG_methylation.shard_*.bedGraph.gz
  tracks/dna_GCH_accessibility.shard_*.bedGraph.gz
  dna_single_cytosine_complete_matrix.vds.h5
  dna_single_cytosine_shard_manifest.tsv
  bigwig/human.WCG_methylation.per_base.tsv.gz
  bigwig/human.WCG_methylation_fraction.bw
  bigwig/human.WCG_methylation_coverage.bw
  bigwig/human.WCG_methylation_covered_cells.bw
  bigwig/human.GCH_accessibility.per_base.tsv.gz
  bigwig/human.GCH_accessibility_fraction.bw
  bigwig/human.GCH_accessibility_coverage.bw
  bigwig/human.GCH_accessibility_covered_cells.bw
  bigwig/pig.WCG_methylation.per_base.tsv.gz
  bigwig/pig.WCG_methylation_fraction.bw
  bigwig/pig.WCG_methylation_coverage.bw
  bigwig/pig.WCG_methylation_covered_cells.bw
  bigwig/pig.GCH_accessibility.per_base.tsv.gz
  bigwig/pig.GCH_accessibility_fraction.bw
  bigwig/pig.GCH_accessibility_coverage.bw
  bigwig/pig.GCH_accessibility_covered_cells.bw
```

The complete cell-level QC table is:

```text
08_qc_summary/23_dna_single_cytosine_complete/
  dna_single_cytosine_complete_per_cell.tsv.gz
```

The HDF5 COO datasets are:

```text
cell_index
contig_index
cytosine_position0
strand
context_code
covered
converted
signal_fraction
```

`context_code=1` is WCG and `context_code=2` is GCH. GCG is intentionally
absent from the primary matrix and signal tracks, but its excluded call count
is retained in task, sample and complete summaries.

RNA is independent and retains both complete raw and filtered GeneFull
exact-UMI gene-by-cell matrices. MatrixMarket remains authoritative; step 76
also writes a convenient zero-based sparse HDF5 and records cell metadata:

```bash
csub < 03_scripts/04_postmapping/77_submit_rna_complete_matrix_8core.csub
```

RNA primary products:

```text
07_postmapping/complete_matrices/rna_expression/
  rna_expression_complete_sparse.h5
  rna_expression_matrix_manifest.tsv
```

After DNA step 75 and RNA step 77 both finish:

```bash
csub < 03_scripts/04_postmapping/79_submit_complete_modalities_summary_4core.csub
```

The final summary deliberately does not assume that DNA and RNA `cell_id`
values are one-to-one paired. Biological interpretation is not part of
steps 63-84.

Progress can be checked at any time without reading BAM/HDF5 payloads:

```bash
python 03_scripts/04_postmapping/80_check_complete_matrix_progress.py
```

After the corrected DNA dedup result and rerun steps 53-75 are complete,
build species-separated per-base tables and pseudo-bulk BigWigs:

```bash
csub < 03_scripts/04_postmapping/82_submit_dna_single_base_bigwigs_16core.csub
```

Step 82 activates `foodie_env` by default because it contains
`bedGraphToBigWig`. Override with `BIGWIG_ENV=<environment>` only when another
environment provides both `numpy` and `bedGraphToBigWig`.

`*_fraction.bw` stores a value from 0 to 1 at each observed cytosine.
`*_coverage.bw` stores the number of PCR-deduplicated template calls. WCG is
endogenous methylation, GCH is GpC accessibility, and ambiguous GCG is
excluded. BigWig cannot store a cell dimension, so these are pseudo-bulk
tracks; the HDF5 matrix remains the complete single-cell result. Human and
pig are emitted separately using the prefixed joint-reference contigs.
`*_covered_cells.bw` stores the number of distinct cells covering each site.
The QC summary also reports reference-assayable WCG/GCH sites and
`observed_sites / reference_assayable_sites` for each species and context.

After steps 74-75 and 81-82 are both validated, create the complete per-cell
coverage QC table:

```bash
csub < 03_scripts/04_postmapping/84_submit_dna_single_cell_qc_8core.csub
```

This report keeps human and pig WCG/GCH observations separate. It reports
per-cell observed sites, PCR-deduplicated template calls, signal fractions,
species-specific reference coverage fractions, distribution quantiles, and
counts above 1k/10k/100k observed-site thresholds. It does not impose a cell
filter or a species-doublet cutoff.
