# Data description

## Source data

This project uses gene-level RNA-seq read counts from GTEx v8.

- **Tissue:** Muscle - Skeletal
- **Number of samples:** 818
- **Initial number of genes:** approximately 59,000

GTEx source data are not redistributed in this repository. Large input and intermediate files remain in the original research computing environment.

## Preprocessing summary

1. Skeletal muscle samples were selected from the GTEx dataset.
2. Genes with fewer than 50 total reads across all samples were removed.
3. Variance-stabilising transformation was performed with DESeq2 using a design of `~ 1`.
4. The expression matrix was transposed into samples by genes for WGCNA.
5. The 20,000 genes with the highest median absolute deviation (MAD) were retained.

## WGCNA input

- **Samples:** 818
- **Genes:** 20,000
- **Network type:** Signed
- **Soft-thresholding power:** 6
- **Minimum module size:** 30
- **Module merge cut height:** 0.25

## Files not included

The following large files are intentionally excluded:

- Raw GTEx expression matrices
- FASTQ and BAM files
- Transformed full expression matrices
- Large R objects such as `.rds` files
- Temporary and intermediate analysis files
