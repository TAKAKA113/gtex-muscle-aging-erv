# Reproducibility guide

This guide records the intended execution order and the boundaries of what is stored in the repository.

## 1. External inputs

The analysis requires external data and reference files that are not committed to GitHub:

- GTEx v10 skeletal-muscle gene counts and phenotype metadata
- GSE164471 / SRP300916 FASTQ files for 53 samples
- GRCh38 primary-assembly FASTA
- GENCODE v39 annotation for the GTEx genomic-context analysis
- GENCODE v26 primary-assembly annotation for GSE164471
- TEtranscripts-compatible GRCh38 RepeatMasker GTF
- Telescope `HERV_rmsk.hg38.v2` annotation
- UCSC GRCh38 RepeatMasker annotation
- Roadmap Epigenomics skeletal-muscle ChromHMM annotation
- GenAge human gene annotation

See `01_gtex/data_description.md` and `02_gse164471/reference_versions.md` for study-specific details.

## 2. GTEx execution order

Run the numbered scripts in `01_gtex/scripts/` in ascending order. The main stages are:

```text
sample and phenotype preparation
-> expression-matrix preparation
-> low-expression filtering and DESeq2 VST
-> WGCNA
-> module-trait and GenAge analyses
-> gene/ERV coordinate preparation
-> promoter overlap and nearest-distance tests
-> hub-gene and epigenomic-context analyses
```

The GTEx workflow does not quantify ERV expression. ERV annotation is added after construction of the gene-only network.

## 3. GSE164471 execution order

First run the numbered shell/SLURM workflow documented in:

- `02_gse164471/scripts/terminal/README.md`
- `02_gse164471/reports/terminal_pipeline.md`

This produces gene/TE-subfamily counts and sample-level Telescope reports. Then run the numbered scripts in `02_gse164471/scripts/R/` in ascending order.

The main downstream stages are:

```text
TEcount matrix preparation
-> ERV filtering and gene-biotype selection
-> combined gene–ERV WGCNA
-> DESeq2 age-association tests
-> functional enrichment
-> Telescope locus-matrix preparation
-> locus-level DESeq2 analysis
-> module-constrained ERV-locus–gene correlation
-> genomic annotation and final figures
```

## 4. Statistical interpretation

- Age models include sex as a covariate where specified in the analysis reports.
- Multiple-testing correction uses the Benjamini–Hochberg method.
- FDR < 0.05 is the significance threshold for the primary age-association analyses.
- The top 0.1% absolute bicor threshold is an exploratory prioritisation rule, not a statistical-significance threshold.
- ERV subfamily concentration in the turquoise and yellow WGCNA modules is descriptive and should not be reported as formal enrichment.
- Co-expression and genomic overlap do not establish causal regulation.

## 5. Expected key checkpoints

| Checkpoint | Expected value |
|---|---:|
| GTEx samples | 818 |
| GTEx WGCNA input genes | 20,000 |
| GSE164471 samples | 53 |
| GSE164471 WGCNA input features | 18,230 |
| ERV subfamilies in WGCNA | 465 |
| Expressed Telescope loci | 11,503 |
| FDR-significant ERV subfamilies for age | 0 |
| FDR-significant ERV loci for age | 1 |
| Primary overlapping ERV–gene candidates | 9 |

If these checkpoints differ, first verify input versions, sample matching, filtering rules and annotation compatibility.

