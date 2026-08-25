# GSE164471 repeat-aware RNA-seq terminal pipeline

## Dataset

- **Study:** GSE164471 / SRP300916
- **Tissue:** Human left vastus lateralis skeletal muscle
- **Samples:** 53 healthy individuals, ages 22–83 years
- **Library:** Single-end total RNA-seq; treated as unstranded after strandedness assessment
- **Reference genome:** GRCh38 primary assembly
- **Gene annotation:** GENCODE v26

## Workflow

```text
SRA accessions
  -> prefetch / fasterq-dump
  -> FASTQ
  -> FastQC / MultiQC
  -> STAR alignment retaining multimapping reads
  -> coordinate-sorted BAM
       |-> TEcount -> gene + TE-subfamily counts
       `-> samtools name sort -> Telescope -> HERV-locus counts
```

No adapter trimming, quality trimming or deduplication was applied after QC.

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
| Protein-coding genes before expression filtering | 19,846 |
| Strict ERV subfamilies before expression filtering | 478 |
| Strict ERV subfamilies after expression filtering | 465 |
| Telescope annotation loci | 14,968 |
| Telescope loci reported before expression filtering | 14,177 |
| Telescope loci after expression filtering | 11,503 |
| Telescope samples completed | 53 / 53 |

## STAR alignment

STAR v2.7.11b was used with GENCODE v26. Multimapping reads were retained for repeat-aware quantification using:

```text
--outFilterMultimapNmax 100
--winAnchorMultimapNmax 200
```

Alignments were stored as coordinate-sorted BAM files.

## TEcount branch

TEcount was run on the coordinate-sorted STAR BAM files in unstranded, multimapping mode. The merged raw matrix contained 58,278 gene features and 1,330 TE subfamilies across 53 samples.

The primary downstream ERV definition retained RepeatMasker families **ERV1, ERVK and ERVL**. **ERVL-MaLR** was excluded from this strict set. Protein-coding genes and ERVs were retained when they had at least 10 raw counts in at least 6 samples.

After filtering, the matrix contained 17,826 protein-coding genes and 465 ERV subfamilies. DESeq2 variance-stabilising transformation and removal of zero-MAD features produced the final WGCNA input of 17,765 genes plus 465 ERV subfamilies (18,230 features).

## Telescope branch

Telescope v1.0.4.1 was run separately for each sample using `HERV_rmsk.hg38.v2/transcripts.gtf`, which contains 14,968 annotated HERV loci. Telescope was run in unstranded mode with a minimum overlap threshold of 0.2.

The coordinate-sorted STAR BAM files were temporarily sorted by read name with SAMtools. Temporary BAM files were created in compute-node scratch storage and removed after Telescope completed. Telescope reported 14,177 loci; applying the same count threshold retained 11,503 expressed loci for downstream testing.

## Main intermediate outputs

```text
05_counts/tecount/gene_te_counts_raw.tsv
05_counts/tecount/gene_te_feature_annotation.tsv
03_reference/te_annotation/erv_subfamily_whitelist.tsv
05_counts/telescope/<sample>/<sample>-telescope_report.tsv
```

Large sequencing files, BAM files, reference resources, software environments and temporary files are not stored in the repository.

## Completed downstream analysis

The terminal outputs were analysed with the numbered scripts in `02_gse164471/scripts/R/`. Completed stages were:

1. TEcount matrix construction and strict ERV filtering.
2. Combined gene–ERV-subfamily WGCNA.
3. Module- and feature-level age association adjusted for sex.
4. Functional enrichment of the ERV-rich turquoise and yellow modules.
5. Telescope locus-matrix construction and locus-level age association.
6. Module-constrained ERV-locus–gene biweight midcorrelation.
7. Genomic annotation and prioritisation of nine gene-overlapping candidates.

The complete numerical results are recorded in `02_gse164471/reports/final_analysis_summary.md`.

