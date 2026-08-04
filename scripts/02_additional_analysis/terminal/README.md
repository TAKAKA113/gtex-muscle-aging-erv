# Terminal preprocessing scripts

This directory is reserved for the command-line and Slurm scripts used for the GSE164471 skeletal-muscle RNA-seq analysis.

## Planned structure

```text
terminal/
├── shell/
├── slurm/
├── R/
└── environment_specs/
```

The scripts should be copied from the BlueBEAR project directory without renumbering so that script names remain consistent with the original logs and output summaries.

Main steps represented by these scripts are:

1. FASTQ quality control.
2. GRCh38 STAR index construction.
3. Repeat-aware STAR alignment.
4. BAM quality control and indexing.
5. Empirical strandedness assessment.
6. TEcount quantification and matrix merging.
7. Gene and TE feature annotation.
8. ERV subfamily whitelist construction.
9. Telescope environment installation.
10. Telescope HERV locus quantification.

Large input data, BAM files, STAR indices, Conda environments and temporary files are excluded from version control.
