# Skeletal Muscle Ageing and ERV Analysis

This repository contains an MSc Bioinformatics dissertation project investigating endogenous retroviruses (ERVs) in human skeletal muscle using two complementary datasets.

The project is organised as two related but analytically distinct studies rather than a formal cross-dataset integration. GTEx provides a large-scale skeletal-muscle gene co-expression and genomic ERV-context analysis, whereas GSE164471 enables direct ERV expression quantification at subfamily and individual-locus resolution.

## Research aims

1. Characterise age-associated gene co-expression structure in GTEx skeletal muscle and examine the genomic distribution of annotated ERV fragments around candidate genes and modules.
2. Quantify ERV expression directly from GSE164471 skeletal-muscle RNA-seq data at subfamily and locus levels.
3. Test whether ERV expression shows broad age association or preferential co-expression-network structure.
4. Resolve individual ERV loci and identify locus-specific age associations and candidate ERV–gene relationships using co-expression and genomic overlap.

## Datasets

### Study 1: GTEx skeletal muscle

- **Project:** GTEx v10
- **Tissue:** Muscle - Skeletal
- **Samples:** 818
- **Input:** Gene-level RNA-seq read counts
- **Main role:** Large-scale gene co-expression and genomic ERV-context analysis

### Study 2: GSE164471 skeletal-muscle RNA-seq

- **Study:** GSE164471 / SRP300916
- **Samples:** 53 healthy skeletal-muscle samples
- **Age range:** 22–83 years
- **Input:** Single-end total RNA-seq reads
- **Reference:** GRCh38 with GENCODE v26
- **Main role:** Direct ERV expression analysis at subfamily and individual-locus resolution

The original sequencing data and large intermediate files are not included in this repository.

## Analysis components

### 1. GTEx analysis

1. Select GTEx skeletal-muscle samples.
2. Remove low-count genes.
3. Apply variance-stabilising transformation with DESeq2.
4. Select the 20,000 most variable genes by median absolute deviation.
5. Construct a signed WGCNA network.
6. Identify age-associated candidate modules and candidate hub genes.
7. Compare candidate genes and modules with GenAge annotations.
8. Map annotated ERV fragments to candidate-gene genomic context.
9. Evaluate promoter ERV overlap, nearest-ERV distance and ERV-family-specific promoter overlap.

The GTEx analysis showed that ERV fragments are common around skeletal-muscle genes, including genes in age-associated candidate modules, but did not provide strong evidence for broad ERV enrichment specifically within ageing-related modules.

### 2. GSE164471 analysis

1. Retrieve metadata and download FASTQ files.
2. Perform FastQC and MultiQC.
3. Build the STAR index and align reads while retaining multimapping reads.
4. Perform alignment and BAM quality control.
5. Quantify genes and TE subfamilies with TEcount.
6. Restrict the TE analysis to ERV1, ERVK and ERVL subfamilies.
7. Construct a combined gene + ERV WGCNA network.
8. Test module and feature associations with age.
9. Perform functional enrichment for ERV-rich co-expression modules.
10. Quantify individual ERV loci with Telescope.
11. Test locus-level age association.
12. Evaluate module-constrained ERV-locus–gene co-expression using biweight midcorrelation.
13. Annotate genomic relationships between prioritised ERV loci and co-expressed genes.

Key GSE164471 findings include:

- After filtering, 465 ERV subfamilies entered the WGCNA analysis.
- 461 of 465 ERV subfamilies were concentrated in the turquoise and yellow modules.
- No WGCNA module or TEcount ERV subfamily showed significant age association after FDR correction.
- Telescope retained 11,503 expressed ERV loci across 60 ERV subfamilies.
- **MER41_5q14.3a** was the only Telescope locus significantly associated with age after FDR correction (FDR = 0.0212; estimated +25.3% expression per decade).
- MER41_5q14.3a did not overlap an annotated gene body; the nearest annotated TSS was 17,013 bp away.
- An independent co-expression and genomic-overlap analysis identified nine ERV–gene candidate associations, comprising four exonic and five intronic ERV loci.

## Study structure

The final dissertation treats GTEx and GSE164471 as two complementary studies rather than forcing direct statistical integration:

```text
Study 1: GTEx
large-scale skeletal-muscle gene network
+ genomic ERV annotation

Study 2: GSE164471
repeat-aware RNA-seq
→ ERV subfamily expression
→ WGCNA
→ Telescope locus resolution
→ locus-specific age and ERV–gene analyses

Final interpretation
→ compare conclusions and limitations across the two studies in the Discussion
```

This distinction is important because ERV expression was directly quantified in GSE164471 but not in GTEx. Therefore, GTEx is not treated as a direct ERV-expression replication cohort.

## Repository structure

```text
gtex-muscle-aging-erv/
├── README.md
├── 01_gtex/
│   ├── data_description.md
│   ├── scripts/
│   ├── results/
│   └── reports/
└── 02_gse164471/
    ├── scripts/
    ├── results/
    └── reports/
```

- `01_gtex/` contains the GTEx WGCNA, GenAge and genomic ERV-context analyses.
- `02_gse164471/` contains the BlueBEAR preprocessing workflow, TEcount/Telescope quantification and downstream R analyses.
- The previously planned `03_integration/` component was removed because no formal cross-dataset integration was performed.

Detailed GSE164471 command-line workflow:

- `02_gse164471/scripts/terminal/README.md`
- `02_gse164471/reports/terminal_pipeline.md`

## Repository policy

FASTQ files, BAM files, genome indices, Conda environments, large R objects and temporary files are excluded from version control. Scripts, compact tables, selected figures and reports are retained for reproducibility.

## Project status

The core analytical workflow is complete. Remaining work is focused on final figure/table curation and dissertation writing rather than additional cross-dataset analysis.
