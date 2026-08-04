# Terminal preprocessing pipeline

This directory contains the command-line and Slurm scripts used for the GSE164471 skeletal-muscle RNA-seq analysis on BlueBEAR.

## Pipeline overview

```text
metadata取得
    ↓
FASTQ download
    ↓
FastQC / MultiQC
    ↓
STAR index
    ↓
STAR alignment
    ↓
BAM QC
    ↓
strandedness確認
    ↓
TEcount
    ↓
Telescope
```

TEcount and Telescope both use the STAR-aligned BAM files. TEcount provides Gene and TE-subfamily counts, whereas Telescope provides HERV locus-level counts.

## Slurm scripts

| Step | Script | Purpose |
|---|---|---|
| Pilot download | `slurm/01_download_pilot.sbatch` | Download and convert one SRA run to FASTQ for testing |
| Pilot FastQC | `slurm/02_fastqc_pilot.sbatch` | Run FastQC on the pilot FASTQ |
| FASTQ download | `slurm/03_download_all.sbatch` | Download and convert all 53 SRA runs |
| FastQC | `slurm/04_fastqc_all.sbatch` | Run FastQC on all FASTQ files |
| FASTQ MultiQC | `slurm/05_multiqc_fastqc.sbatch` | Combine FastQC results across samples |
| Reference download | `slurm/06_download_gencode_v26.sbatch` | Download GRCh38 and GENCODE v26 reference files |
| STAR index | `slurm/07_build_star_index.sbatch` | Build the STAR genome index |
| STAR alignment | `slurm/08_star_align_all.sbatch` | Align all samples while retaining multimapping reads |
| STAR MultiQC | `slurm/09_multiqc_star.sbatch` | Combine STAR alignment summaries |
| BAM QC | `slurm/10_samtools_bam_qc.sbatch` | Check, index and summarise aligned BAM files |
| Strandedness | `slurm/11_featurecounts_all_strand_test.sbatch` | Compare strandedness settings using featureCounts |
| TEcount | `slurm/12_tecount_all_samples.sbatch` | Quantify Gene and TE-subfamily expression |
| Feature annotation | `slurm/14_build_feature_annotation.sbatch` | Build Gene and TE feature annotations |
| Telescope installation | `slurm/16_install_telescope.sbatch` | Create the Telescope software environment |
| Telescope | `slurm/17_telescope_all_samples.sbatch` | Quantify HERV expression at individual loci |

## Main outputs

```text
53 FASTQ files
    ↓
53 coordinate-sorted STAR BAM files
    ├── TEcount: 58,278 Gene features + 1,330 TE subfamilies
    └── Telescope: 14,968 annotated HERV loci per sample
```

The TEcount branch produced a raw Gene–TE matrix with 59,608 features across 53 samples. Telescope completed successfully for 53 of 53 samples and produced one locus-level report per sample.

## Repository policy

FASTQ files, BAM files, STAR indices, Conda environments and temporary files are not included in GitHub. The scripts are retained with their original names so that they remain consistent with the BlueBEAR logs and output files.
