# GTEx Skeletal Muscle Ageing and ERV Analysis

This repository contains an MSc Bioinformatics dissertation project investigating ageing-related gene co-expression networks and endogenous retroviruses (ERVs) in human skeletal muscle.

## Research aims

1. Identify gene co-expression modules associated with age in GTEx skeletal muscle using WGCNA.
2. Characterise ageing-related modules using hub genes and the GenAge database.
3. Investigate associations between ERV expression and ageing-related gene modules.

## Dataset

- **Project:** GTEx v8
- **Tissue:** Muscle - Skeletal
- **Samples:** 818
- **Input:** Gene-level RNA-seq read counts

The original GTEx data and large intermediate files are not included in this repository.

## Analysis workflow

1. Select GTEx skeletal muscle samples.
2. Remove low-count genes.
3. Apply variance-stabilising transformation with DESeq2.
4. Select the 20,000 most variable genes by median absolute deviation (MAD).
5. Construct a signed weighted gene co-expression network with WGCNA.
6. Test associations between module eigengenes and donor traits.
7. Examine hub genes and GenAge enrichment.
8. Assess postmortem covariates and integrate ERV expression.

## Current status

The WGCNA module detection, module-trait analysis and initial GenAge analysis have been completed. Postmortem covariate assessment and ERV integration are in progress.

## Repository structure

- `scripts/` — R scripts used in the analysis
- `results/` — selected figures and summary tables
- `data_description.md` — dataset and preprocessing information

## Note

This is an ongoing MSc dissertation project. Scripts, documentation and selected results will be added progressively.
