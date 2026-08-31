# snMAR 2607 data retention for the server move

Source project: `/share/org/YZWL/yzwl_wujx/Projects/2607_methylation`

The original FASTQ data already exists on external storage and is therefore
not a P0 copy to the new compute node. No old-server data should be deleted
until the destination audit passes.

## P0: time-expensive restart points and downstream state

Retain these project units:

- `README_snMAR_2607.md`
- `02_reference/`
- `03_scripts/`
- `06_mapping/`
- `07_postmapping/`
- `08_qc_summary/`
- `09_logs/`
- `10_downstream_analysis/`

Within `06_mapping/`, the DNA BWA and RNA STAR BAM/BAI outputs are restart
points for the unfinished post-mapping and downstream work. Within
`07_postmapping/`, keep CB-aware deduplicated BAMs, single-cytosine caller
outputs, matrices, BigWigs, manifests, summaries, and validation tables.

The mapping BAMs remain necessary even when a deduplicated BAM exists: the
current analysis was still being corrected and may need deduplication or
calling rerun without returning to FASTQ processing.

## P1: expensive but reproducible intermediates

If destination capacity permits, also retain:

- `05_fastp/`
- selected completed assigned FASTQs under `04_fastq_split/`

These are lower priority than BAMs and later products because the raw reads are
already backed up elsewhere.

## Not copied to the compute node in the first pass

- `00raw_data/`
- probe/test intermediates and temporary sort files
- failed or incomplete task scratch directories that have no unique output

The first-pass exclusion is a capacity decision, not permission to delete the
old-server copy.

## Verification

Record source path, byte size, mtime, and symlink target before transfer. Audit
file counts and bytes after each directory finishes. Validate all BAMs with
`samtools quickcheck` on the destination and verify every expected BAI. Hash
manifests, matrices, HDF5 files, final summaries, and figures with SHA-256.
