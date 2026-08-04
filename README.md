# GTEx Skeletal Muscle Ageing and ERV Analysis

This repository contains an MSc Bioinformatics dissertation project investigating ageing-related gene co-expression networks and endogenous retroviruses (ERVs) in human skeletal muscle.

## Research aims

1. Identify gene co-expression modules associated with age in GTEx skeletal muscle using WGCNA.
2. Characterise ageing-related modules using hub genes and ageing-related gene sets.
3. Quantify ERV expression at subfamily and locus levels in an additional skeletal-muscle RNA-seq dataset.
4. Compare ERV-associated gene programmes between the additional dataset and GTEx.

## Datasets

### GTEx skeletal muscle

- **Project:** GTEx v10
- **Tissue:** Muscle - Skeletal
- **Samples:** 818
- **Input:** Gene-level RNA-seq read counts

### Additional skeletal-muscle RNA-seq dataset

- **Study:** GSE164471 / SRP300916
- **Samples:** 53
- **Input:** Single-end RNA-seq reads
- **Reference:** GRCh38 with GENCODE v26

The original sequencing data and large intermediate files are not included in this repository.

## Analysis workflow

### GTEx branch

1. Select GTEx skeletal-muscle samples.
2. Remove low-count genes.
3. Apply variance-stabilising transformation with DESeq2.
4. Select the 20,000 most variable genes by median absolute deviation.
5. Construct a signed WGCNA network.
6. Test module associations with age and other donor traits.
7. Examine hub genes and ageing-related gene enrichment.

### Additional RNA-seq branch

1. Perform FASTQ quality control.
2. Align reads to GRCh38 with STAR while retaining multimapping reads.
3. Quantify genes and TE subfamilies with TEcount.
4. Define strict and broad ERV subfamily sets.
5. Quantify individual HERV loci with Telescope.
6. Continue with age regression, WGCNA and pathway analysis in R.

## Current status

- GTEx WGCNA, module-trait analysis and initial hub-gene analysis are complete.
- Terminal preprocessing and ERV quantification for all 53 additional samples are complete.
- The additional dataset contains 58,278 gene features and 1,330 TE subfamilies from TEcount.
- The strict ERV set contains 478 subfamilies.
- Telescope completed for 53 of 53 samples using an annotation of 14,968 HERV loci.
- Downstream statistical analysis in R is the next stage.

## Repository structure

- `scripts/01_gtex/` — GTEx analysis scripts
- `scripts/02_additional_analysis/` — additional RNA-seq analysis scripts
- `scripts/03_integration/` — cross-dataset integration scripts
- `results/01_gtex/` — GTEx figures and tables
- `results/02_additional_analysis/` — additional dataset figures and tables
- `results/03_integration/` — integrated outputs
- `reports/` — concise workflow and progress reports
- `data_description.md` — dataset and preprocessing information

Detailed command-line workflow:

- `reports/02_additional_analysis_terminal_pipeline.md`

## Note

This is an ongoing MSc dissertation project. Scripts, documentation and selected results will be added progressively.
