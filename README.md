# Endogenous Retroviral Expression in Human Skeletal Muscle Ageing

This repository contains the reproducible analysis for an MSc Bioinformatics dissertation examining endogenous retroviruses (ERVs) in human skeletal muscle at three levels: genomic context, ERV subfamily expression and individual ERV-locus expression.

The project uses two complementary datasets. They were analysed independently because direct repeat-aware ERV quantification was available only for GSE164471.

## Study design

| Study | Data | Samples | Primary analysis |
|---|---|---:|---|
| GTEx skeletal muscle (v10) | Gene-level counts | 818 | Gene WGCNA followed by genomic ERV-context analysis |
| GSE164471 / SRP300916 | Single-end total RNA-seq | 53 | Repeat-aware STAR alignment, TEcount, gene–ERV WGCNA and Telescope |

```text
GTEx gene counts                         GSE164471 FASTQ
       |                                       |
       v                                       v
gene-only WGCNA                     repeat-aware STAR alignment
       |                                  /            \
       v                                 v              v
age-associated candidate modules      TEcount       Telescope
       |                                 |              |
       v                                 v              v
promoter overlap, nearest ERV       subfamily       locus-level age and
and ERV-family composition          WGCNA/DESeq2    ERV–gene analyses
       \                                  /
        \________________________________/
                 interpretation across studies
```

## Main findings

- GTEx identified seven candidate age-associated gene modules. ERV fragments were common near skeletal-muscle genes, but promoter overlap, nearest-ERV distance and family composition did not show a consistent broad enrichment across these modules.
- In GSE164471, 461 of 465 filtered ERV subfamilies (99.1%) were assigned to the turquoise or yellow WGCNA modules. This was a descriptive network concentration, not a formal enrichment result.
- No WGCNA module or TEcount ERV subfamily was significantly associated with age after Benjamini–Hochberg correction.
- Telescope retained 11,503 expressed ERV loci. **MER41_5q14.3a** was the only locus significantly associated with age (FDR = 0.0212; estimated +25.3% expression per decade).
- An independent exploratory analysis identified nine strongly co-expressed ERV–gene pairs with direct gene-body overlap: four exonic and five intronic.

These results do not support uniform ERV activation with chronological age in skeletal muscle. They instead nominate locus-specific candidates. Co-expression, proximity and overlap do not demonstrate causal regulation.

## Repository structure

```text
gtex-muscle-aging-erv/
├── README.md
├── REPRODUCIBILITY.md
├── 01_gtex/
│   ├── data_description.md
│   ├── scripts/
│   ├── results/
│   └── reports/
└── 02_gse164471/
    ├── reference_versions.md
    ├── scripts/
    │   ├── terminal/
    │   └── R/
    ├── results/
    └── reports/
```

Start here:

- [GTEx analysis summary](01_gtex/reports/GTEx_summary.md)
- [GSE164471 terminal pipeline](02_gse164471/reports/terminal_pipeline.md)
- [GSE164471 final analysis summary](02_gse164471/reports/final_analysis_summary.md)
- [Reference genome and annotation versions](02_gse164471/reference_versions.md)
- [Reproducibility guide](REPRODUCIBILITY.md)

## Analysis outline

### 1. GTEx genomic-context study

1. Filter GTEx skeletal-muscle gene counts.
2. Apply DESeq2 variance-stabilising transformation.
3. Select the 20,000 genes with highest median absolute deviation.
4. Construct a signed WGCNA network and identify candidate age-associated modules.
5. Rank hub-gene candidates and compare modules with GenAge.
6. Use GRCh38 RepeatMasker annotations to test promoter ERV overlap, nearest-ERV distance and ERV-family composition.
7. Add gene-body and Roadmap Epigenomics context for prioritised ERV fragments.

ERV expression was not quantified in GTEx; its WGCNA matrix was gene-only.

### 2. GSE164471 repeat-aware expression study

1. Run FastQC/MultiQC and align reads to GRCh38 with STAR while retaining multimapping reads.
2. Quantify genes and TE subfamilies with TEcount.
3. Retain 465 ERV1, ERVK and ERVL subfamilies and construct a combined gene–ERV WGCNA network.
4. Test module and feature associations with age while adjusting for sex.
5. Perform functional enrichment for the ERV-rich turquoise and yellow modules.
6. Quantify individual ERV loci with Telescope.
7. Test locus-level age association with DESeq2.
8. Prioritise module-constrained ERV-locus–gene correlations and annotate their genomic relationships.

## Data and repository policy

Raw FASTQ files, BAM files, STAR indices, reference annotations, Conda environments, large R objects and temporary files are not version-controlled. The repository retains scripts, compact result tables, selected figures and reports needed to understand and reconstruct the workflow.

The GSE164471 raw sequencing data are available through GEO/SRA under **GSE164471 / SRP300916**. GTEx data access and use remain subject to GTEx policies.

## Project status

The dissertation analysis is complete. The repository represents the final analytical workflow and distinguishes confirmatory statistical results from exploratory candidate prioritisation.

