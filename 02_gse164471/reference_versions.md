# Reference resources used for GSE164471

This file records the reference resources used in the GSE164471 skeletal-muscle RNA-seq analysis.

## Reference genome

- **Genome build:** GRCh38
- **Assembly scope:** primary assembly
- **File:** `GRCh38.primary_assembly.genome.fa`
- **Use:** STAR genome index construction and RNA-seq alignment
- **Source:** GENCODE release 26 download location used in the project pipeline

## Gene annotation

- **Annotation:** GENCODE v26
- **Assembly:** GRCh38 primary assembly
- **File:** `gencode.v26.primary_assembly.annotation.gtf`
- **Use:** STAR genome index construction, gene-level annotation, TEcount gene annotation, and downstream genomic-context analysis with GenomicRanges

## TEcount transposable-element annotation

- **File:** `GRCh38_GENCODE_rmsk_TE.gtf`
- **Coordinate system:** GRCh38
- **Use:** TEcount transposable-element quantification together with the GENCODE v26 gene annotation
- **Primary ERV definition used downstream:** RepeatMasker families `ERV1`, `ERVK`, and `ERVL`
- **Excluded from the strict ERV set:** `ERVL-MaLR`

The repository records the exact TE GTF filename used by the TEcount pipeline. The full annotation file is not stored in GitHub because reference annotation files are treated as external inputs.

## Telescope HERV annotation

- **Annotation:** `HERV_rmsk.hg38.v2`
- **File:** `HERV_rmsk.hg38.v2/transcripts.gtf`
- **Assembly:** hg38 / GRCh38
- **Use:** locus-specific HERV quantification with Telescope and downstream locus-to-subfamily mapping
- **Annotated HERV loci:** 14,968

The Telescope annotation provides locus identifiers such as `MER41_5q14.3a` and an `intModel` field used in this project as the ERV subfamily/model annotation.

## Coordinate compatibility

All genomic analyses in the GSE164471 workflow were performed in the GRCh38/hg38 coordinate system. This allowed the STAR-aligned reads, GENCODE v26 gene annotation, TEcount annotation, Telescope HERV annotation, and downstream GenomicRanges overlap analysis to be compared without cross-build coordinate conversion.

## Files intentionally not stored in the repository

Large external reference resources are not version-controlled, including:

- genome FASTA files
- GENCODE GTF files
- TE annotation GTF files
- Telescope HERV GTF files
- STAR genome indices

Instead, this document and the pipeline scripts record the filenames, versions, and roles required to reconstruct the analysis environment.
