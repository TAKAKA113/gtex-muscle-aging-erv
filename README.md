# GTEx Skeletal Muscle Ageing and ERV Analysis

This repository contains an MSc Bioinformatics dissertation project investigating ageing-related gene co-expression networks and endogenous retroviruses (ERVs) in human skeletal muscle.

## Research aims

1. Identify gene co-expression modules associated with age in GTEx skeletal muscle using WGCNA.
2. Characterise ageing-related modules using hub genes and ageing-related gene sets.
3. Quantify ERV expression at subfamily and locus levels in GSE164471 skeletal-muscle RNA-seq data.
4. Integrate ERV-associated genes, modules and pathways across GSE164471 and GTEx.

## Datasets

### GTEx skeletal muscle

- **Project:** GTEx v10
- **Tissue:** Muscle - Skeletal
- **Samples:** 818
- **Input:** Gene-level RNA-seq read counts

### GSE164471 skeletal-muscle RNA-seq

- **Study:** GSE164471 / SRP300916
- **Samples:** 53
- **Input:** Single-end RNA-seq reads
- **Reference:** GRCh38 with GENCODE v26

The original sequencing data and large intermediate files are not included in this repository.

## Analysis components

### 1. GTEx analysis

1. Select GTEx skeletal-muscle samples.
2. Remove low-count genes.
3. Apply variance-stabilising transformation with DESeq2.
4. Select the 20,000 most variable genes by median absolute deviation.
5. Construct a signed WGCNA network.
6. Test module associations with age and donor traits.
7. Examine hub genes, GenAge enrichment and ERV genomic context.

### 2. GSE164471 analysis

1. Retrieve metadata and download FASTQ files.
2. Perform FastQC and MultiQC.
3. Build the STAR index and align reads while retaining multimapping reads.
4. Perform BAM QC and confirm library strandedness.
5. Quantify genes and TE subfamilies with TEcount.
6. Quantify individual HERV loci with Telescope.
7. Continue with filtering, normalisation, age association, WGCNA, enrichment and genomic annotation in R.

### 3. Cross-dataset integration

1. Compare age-associated genes and modules.
2. Compare hub genes and ageing-related gene sets.
3. Compare enriched biological pathways.
4. Integrate GSE164471 ERV candidates with GTEx gene-network evidence.

## Current status

- GTEx WGCNA, module-trait analysis, hub-gene analysis and initial ERV genomic-context analysis are complete.
- Terminal preprocessing and ERV quantification are complete for all 53 GSE164471 samples.
- TEcount produced 58,278 gene features and 1,330 TE subfamilies.
- The strict ERV input set contains 478 subfamilies.
- Telescope completed for 53 of 53 samples using an annotation of 14,968 HERV loci.
- Downstream statistical analysis of GSE164471 in R is the next stage.

## Repository structure

```text
gtex-muscle-aging-erv/
├── README.md
├── 01_gtex/
│   ├── data_description.md
│   ├── scripts/
│   ├── results/
│   └── reports/
├── 02_gse164471/
│   ├── scripts/
│   ├── results/
│   └── reports/
└── 03_integration/
    ├── scripts/
    ├── results/
    └── reports/
```

- `01_gtex/` contains the GTEx analysis, selected outputs and progress report.
- `02_gse164471/` contains the BlueBEAR preprocessing pipeline and downstream R analysis.
- `03_integration/` is reserved for analyses that directly compare or combine the two datasets.

Detailed GSE164471 command-line workflow:

- `02_gse164471/scripts/terminal/README.md`
- `02_gse164471/reports/terminal_pipeline.md`

## Repository policy

FASTQ files, BAM files, genome indices, Conda environments, large R objects and temporary files are excluded from version control. Scripts, compact tables, selected figures and reports are retained for reproducibility.

## Note

This is an ongoing MSc dissertation project. Scripts, documentation and selected results will be added progressively.
