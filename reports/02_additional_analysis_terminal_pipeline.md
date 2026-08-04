# Additional RNA-seq terminal pipeline

## Dataset

- **Study:** GSE164471 / SRP300916
- **Tissue:** Human skeletal muscle
- **Samples:** 53
- **Library type:** Single-end, unstranded RNA-seq
- **Reference genome:** GRCh38
- **Gene annotation:** GENCODE v26

## Workflow

```text
FASTQ
  ↓
FastQC / MultiQC
  ↓
STAR alignment retaining multimapping reads
  ↓
Coordinate-sorted BAM
  ├── TEcount → Gene + TE-subfamily counts
  └── samtools name sorting → Telescope → HERV-locus counts
```

## Processing summary

| Stage | Result |
|---|---:|
| RNA-seq samples | 53 |
| Samples retained after QC | 53 |
| STAR overall mapping rate | Median 93.42% |
| STAR unique mapping rate | Median 89.71% |
| Gene + TE features from TEcount | 59,608 |
| Gene features | 58,278 |
| TE subfamilies | 1,330 |
| Protein-coding genes | 19,846 |
| Strict ERV subfamilies | 478 |
| Broad ERV/LTR set | 564 |
| Telescope HERV loci | 14,968 |
| Telescope samples completed | 53 / 53 |

## TEcount branch

TEcount was run from the STAR coordinate-sorted BAM files using an unstranded library setting. The merged raw count matrix contains 58,278 gene features and 1,330 TE subfamilies across 53 samples.

Gene biotype was retained to define the gene population used in the primary downstream analysis. The planned primary feature set consists of protein-coding genes and strictly defined ERV subfamilies. Other gene biotypes remain in the original count matrix.

The strict ERV set contains RepeatMasker families **ERV1, ERVK and ERVL**. The broader sensitivity set additionally contains **ERVL-MaLR**.

## Telescope branch

Telescope 1.0.4.1 was run separately for each sample using the official `HERV_rmsk.hg38.v2/transcripts.gtf` annotation. This annotation contains 14,968 HERV loci.

The coordinate-sorted STAR BAM files were temporarily sorted by read name with SAMtools. The temporary BAM files were created in compute-node scratch storage and removed after Telescope completed. The final sample-level outputs are the `*-telescope_report.tsv` files.

## Main outputs

```text
05_counts/tecount/gene_te_counts_raw.tsv
05_counts/tecount/gene_te_feature_annotation.tsv
03_reference/te_annotation/erv_subfamily_whitelist.tsv
05_counts/telescope/<sample>/<sample>-telescope_report.tsv
```

Large sequencing files, BAM files, software environments and intermediate files are not stored in this repository.

## Downstream analysis

The following steps will be performed in R:

1. Build the Telescope locus-by-sample matrix.
2. Select protein-coding genes and ERV features.
3. Apply low-expression filtering and transformation.
4. Perform feature-wise age regression.
5. Construct the primary age-blind WGCNA network.
6. Identify age-associated modules, hub genes and ERVs.
7. Perform pathway enrichment and locus annotation.
8. Compare the additional dataset with the GTEx analysis.
