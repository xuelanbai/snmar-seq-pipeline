# 03 mapping and mapping preparation

Project root:

```text
/share/org/YZWL/yzwl_wujx/Projects/2607_methylation
```

## DNA steps 01-07

DNA is the TAPS branch of snMAR. The current strict analysis uses only
`XD:i:0` reads.

```text
01_stream_strict_fastq.py
02_discover_dna_assigned_pairs.py
03_prepare_dna_xd0_front_trim.py
04_summarize_dna_fastp_probe.py
05_submit_dna_fastp_probe_64core.csub
06_summarize_dna_fastp_full.py
07_submit_dna_fastp_full_64core.csub
```

Current DNA structure trimming:

```text
R1: 30 bp
R2: 43 bp
```

DNA fastp outputs:

```text
05_fastp/dna_strict/
08_qc_summary/04_dna_fastp_full/
```

## RNA steps 08-11: R1-only workflow

RNA library structure:

```text
R1: transcript read
R2: T7(8) + linker(22-25) + T5(8) + UMI(8) + dT(31) + read
```

Barcode and UMI are extracted during strict10 demultiplex and copied to the
R1 FASTQ comment. RNA fastp and mapping therefore use only R1. R2 is not read
or written by the following steps.

```text
08_discover_rna_r1_assigned.py
09_summarize_rna_r1_fastp_full.py
10_submit_rna_r1_fastp_full_64core.csub
11_finalize_remove_old_rna_pe_fastp.sh
```

Run:

```bash
cd /share/org/YZWL/yzwl_wujx/Projects/2607_methylation
unset FORCE
csub < 03_scripts/03_mapping/10_submit_rna_r1_fastp_full_64core.csub
```

Input:

```text
04_fastq_split/rna_strict10/<sample>/<pair_id>.assigned_R1.fq.gz
```

Output:

```text
05_fastp/rna_r1_strict/<sample>/<pair_id>.clean_R1.fq.gz
08_qc_summary/07_rna_r1_fastp_full/rna_r1_fastp_full_summary.tsv
08_qc_summary/07_rna_r1_fastp_full/rna_r1_fastp_full_per_file.tsv
08_qc_summary/07_rna_r1_fastp_full/fastp_json/
09_logs/03_mapping/rna_r1_fastp_full_by_file/
```

The expected total fastp input is exactly the total `assigned_XD0` count in:

```text
08_qc_summary/02_rna_demultiplex/write_strict10_per_pair.tsv
```

For the current dataset this is:

```text
58,287,579 R1 reads
```

R1 header tags required after fastp:

```text
CR CB CN UR UB XD XE XL XI
```

## Safe cleanup

The old paired-end RNA fastp code/results are removed only after the new
R1-only summary, 96 clean files, demultiplex input total, and zero failed-file
checks all pass.

Dry run:

```bash
bash 03_scripts/03_mapping/11_finalize_remove_old_rna_pe_fastp.sh
```

Confirmed cleanup:

```bash
CONFIRM_REMOVE_OLD_RNA_PE=1 \
  bash 03_scripts/03_mapping/11_finalize_remove_old_rna_pe_fastp.sh
```

Strict split R2 source FASTQ files are preserved.

## RNA steps 12-20: joint human+pig STAR mapping

The current RNA fastp result contains:

```text
96 clean R1 files
24,737,553 clean R1 reads
```

Scripts:

```text
12_prepare_joint_rna_reference.py
13_submit_joint_rna_star_index_64core.csub
14_discover_rna_star_inputs.py
15_rna_fastq_sam_tags.py
16_run_rna_star_one.sh
17_summarize_rna_star.py
18_submit_rna_star_probe_64core.csub
19_submit_rna_star_full_64core.csub
20_finalize_remove_rna_star_probe.sh
21_submit_joint_dna_bwa_index_8core.csub
22_discover_dna_bwa_inputs.py
23_prepare_dna_bwa_probe_fastq.py
24_validate_dna_bwa_bam.py
25_run_dna_bwa_one.sh
26_summarize_dna_bwa.py
27_submit_dna_bwa_prepare_8core.csub
28_submit_dna_bwa_probe_64core.csub
29_submit_dna_bwa_sample_template_64core.csub
30_submit_dna_bwa_full_96jobs.sh
31_submit_dna_bwa_full_summary_8core.csub
32_finalize_remove_dna_bwa_probe.sh
```

### One-time STAR environment setup

The current `snTAPS_align` environment contains samtools but not STAR. Install
STAR directly into that existing alignment environment:

```bash
source /share/apps/anaconda3/etc/profile.d/conda.sh
conda install -n snTAPS_align \
  -c conda-forge -c bioconda \
  star=2.7.11b
```

Steps 13, 18, and 19 activate `snTAPS_align` automatically. A different
environment can be selected with `RNA_ALIGN_ENV=<name>`.

### Reference and index

Step 12 creates:

```text
02_reference/rna_human_pig_joint/human_pig.prefixed.fa
02_reference/rna_human_pig_joint/human_pig.prefixed.gtf
02_reference/rna_human_pig_joint/reference_manifest.tsv
```

Human contigs use `human_`; pig contigs use `pig_`. Human stable feature IDs
use `HUMAN_`; pig stable feature IDs use `PIG_`. This prevents chromosome and
feature-ID collisions in the union reference.

Build the STAR index:

```bash
csub < 03_scripts/03_mapping/13_submit_joint_rna_star_index_64core.csub
```

Index:

```text
02_reference/rna_human_pig_joint/star_index_sjdb149/
```

### Why step 15 is required

The RNA `CR/CB/CN/UR/UB/XD/XE/XL/XI` tags currently reside after whitespace
in the FASTQ header. A standard FASTQ STAR run does not convert those comment
fields into BAM auxiliary tags. Step 15 temporarily encodes them into a
SAM-compatible QNAME, and step 16 restores the original QNAME and BAM tags
after mapping. The encoded FASTQ and encoded BAM are temporary and removed
after validation.

### Probe

Step 14 deterministically selects 8 files spanning low, median, and high
fastp retention plus small/large clean files. Run:

```bash
csub < 03_scripts/03_mapping/18_submit_rna_star_probe_64core.csub
```

Outputs:

```text
06_mapping/rna_star_joint_probe/
08_qc_summary/08_rna_star_probe/rna_star_probe_summary.tsv
08_qc_summary/08_rna_star_probe/rna_star_probe_per_file.tsv
```

The probe is diagnostic. Review unique, multi-mapped, total primary mapping,
STAR unmapped-too-short, and human/pig percentages before step 19.

### Full mapping

```bash
csub < 03_scripts/03_mapping/19_submit_rna_star_full_64core.csub
```

Outputs:

```text
06_mapping/rna_star_joint/
08_qc_summary/09_rna_star_full/rna_star_full_summary.tsv
08_qc_summary/09_rna_star_full/rna_star_full_per_file.tsv
09_logs/03_mapping/rna_star_full_by_file/
```

The full job runs 8 files in parallel with 8 STAR threads per file. It loads
one shared-memory index per allocated node. Probe and full jobs should not be
run simultaneously on the same node when automatic shared-index removal is
enabled.

### Probe cleanup after full validation

Dry run:

```bash
bash 03_scripts/03_mapping/20_finalize_remove_rna_star_probe.sh
```

Confirmed:

```bash
CONFIRM_REMOVE_RNA_STAR_PROBE=1 \
  bash 03_scripts/03_mapping/20_finalize_remove_rna_star_probe.sh
```

Only probe BAM/log outputs are removed. Probe QC tables and all full outputs
are retained.

## DNA paired-end BWA mapping

DNA mapping reuses the `human_`/`pig_` joint FASTA from step 12. Build a
separate BWA index and generate validated input manifests:

```bash
csub < 03_scripts/03_mapping/21_submit_joint_dna_bwa_index_8core.csub
csub < 03_scripts/03_mapping/27_submit_dna_bwa_prepare_8core.csub
```

The two preparation jobs are independent. After both succeed, run:

```bash
csub < 03_scripts/03_mapping/28_submit_dna_bwa_probe_64core.csub
```

The probe maps 8 representative lane pairs, 1 million read pairs per lane.
It uses paired-end BWA-MEM followed by coordinate sorting and indexing.
Global `markdup` is deliberately not used: one BAM contains many `CB`-tagged
cells, and coordinate-only duplicate marking across cells is invalid.

Probe QC:

```text
08_qc_summary/11_dna_bwa_probe/dna_bwa_probe_summary.tsv
08_qc_summary/11_dna_bwa_probe/dna_bwa_probe_per_run.tsv
```

Full mapping is split into 96 independent 64-CPU jobs, one sample and eight
ordered lane pairs per job. Preview any sample-index range:

```bash
START_INDEX=1 END_INDEX=24 \
  bash 03_scripts/03_mapping/30_submit_dna_bwa_full_96jobs.sh
```

Submit the reviewed range:

```bash
START_INDEX=1 END_INDEX=24 CONFIRM_SUBMIT_DNA_BWA=1 \
  bash 03_scripts/03_mapping/30_submit_dna_bwa_full_96jobs.sh
```

After all 96 jobs have finished:

```bash
csub < 03_scripts/03_mapping/31_submit_dna_bwa_full_summary_8core.csub
```

The source demultiplex tag `XE` is the four-segment mismatch vector. It is
renamed to `ZE` before `bwa mem -C` so it cannot collide with a BWA-owned
`XE` tag. BAM validation requires `CR/CB/CN/XD/ZE`, requires `XD:i:0`, checks
primary-record count closure, and reports human/pig primary mapping.

If duplicate removal is required later, it must be done within each `CB`, not
globally across the sample BAM.
