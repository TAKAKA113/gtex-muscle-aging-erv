# GSE164471 ERV analysis — final summary

## Study overview

This analysis used 53 healthy human skeletal-muscle RNA-seq samples (ages 22–83 years) to test whether endogenous retrovirus (ERV) expression is associated with ageing and to characterise ERV–gene relationships at subfamily and individual-locus resolution.

The analysis proceeded from repeat-aware read alignment to TE subfamily quantification, gene–ERV co-expression network analysis, locus-specific Telescope quantification, age association, functional enrichment, ERV–gene co-expression, and genomic annotation.

---

## 1. Repeat-aware RNA-seq processing

The 53 samples were processed from raw sequencing reads using a repeat-aware STAR alignment strategy that retained multimapping reads for downstream TE analysis.

Main downstream quantification layers:

- Gene-level counts
- TE/ERV subfamily-level counts with TEcount
- Individual HERV/ERV loci with Telescope

A strict ERV definition was used for the TEcount analysis: LTR elements belonging to ERV1, ERVK, or ERVL. ERVL-MaLR was excluded from this strict ERV set.

---

## 2. TEcount preprocessing and WGCNA1

TEcount initially produced 1,330 TE subfamilies, of which 478 met the strict ERV definition.

```text
Initial combined matrix
20,324 features
= 19,846 genes + 478 ERV subfamilies

        ↓ low-expression filter
          count ≥ 10 in ≥ 6 samples

18,291 features
= 17,826 genes + 465 ERV subfamilies

        ↓ DESeq2 VST
        ↓ remove MAD = 0 features

WGCNA1 input
18,230 features
= 17,765 genes + 465 ERV subfamilies
```

WGCNA settings:

| Parameter | Value |
|---|---:|
| Samples | 53 |
| Correlation | bicor |
| Network | signed |
| TOM | signed |
| Soft-threshold power | 22 |
| Minimum module size | 30 |
| Merge cut height | 0.25 |

### Module composition

| Module | Genes | ERVs | Total |
|---|---:|---:|---:|
| turquoise | 5,599 | 397 | 5,996 |
| blue | 3,994 | 0 | 3,994 |
| grey | 2,740 | 4 | 2,744 |
| brown | 1,667 | 0 | 1,667 |
| yellow | 1,437 | 64 | 1,501 |
| green | 783 | 0 | 783 |
| red | 520 | 0 | 520 |
| black | 502 | 0 | 502 |
| pink | 410 | 0 | 410 |
| magenta | 113 | 0 | 113 |
| **Total** | **17,765** | **465** | **18,230** |

Of the 465 ERV subfamilies, 461 (99.1%) were assigned to the turquoise or yellow module. This is described as a strong concentration of ERV subfamilies in these two modules, not as formal statistical enrichment.

---

## 3. Association with age

### Module-level association

Module eigengenes were tested using:

```text
module eigengene ~ age + sex
```

No module was significantly associated with age after multiple-testing correction (FDR < 0.05).

The lowest raw age p-value was observed for the red module:

- age beta = −0.00209
- p = 0.0349
- FDR = 0.1819

The ERV-rich turquoise and yellow modules were not significantly associated with age.

### Feature-level age association

All low-expression-filtered features were tested with DESeq2 using:

```text
~ sex + age
```

Results:

| Result | Genes | ERV subfamilies |
|---|---:|---:|
| FDR < 0.05 | 67 | 0 |
| nominal p < 0.05 | 1,114 | 1 |

The only ERV subfamily with nominal p < 0.05 was LTR12B:ERV1:LTR (p = 0.0425, FDR = 0.543), therefore no ERV subfamily was considered significantly age-associated.

Together, the module- and subfamily-level analyses did not support a broad ageing-associated ERV-expression programme in this dataset.

---

## 4. Functional context of ERV-rich WGCNA modules

Functional enrichment was performed for the gene components of the turquoise and yellow ERV-rich modules using clusterProfiler with the 17,765 protein-coding WGCNA genes as the background universe.

### Turquoise module

The turquoise module showed a heterogeneous set of enriched biological processes, including:

- immune and defence-related processes
- keratin/intermediate-filament biology
- neuronal/synaptic categories
- cilium-related processes
- cell adhesion
- sensory/olfactory-receptor-related categories

The olfactory/sensory terms are interpreted as annotation/pathway categories driven by receptor-gene families rather than evidence of olfactory function in skeletal muscle.

### Yellow module

The yellow module showed a more coherent cilium/centrosome-related signal, including:

- cilium organisation
- cilium assembly
- ciliary localisation
- intraflagellar transport
- centriole replication

No KEGG pathway passed the enrichment criteria for the yellow module.

---

## 5. Telescope locus-level ERV quantification

Telescope produced locus-specific ERV counts using the HERV reference annotation.

- Raw Telescope loci: 14,177
- After the same low-expression rule used above (count ≥ 10 in ≥ 6 samples): 11,503 loci
- All 11,503 filtered loci were found in the Telescope HERV GTF
- Telescope loci represented 60 ERV subfamilies (`intModel` annotations)

Mapping the 60 ERV subfamilies back to WGCNA1 showed:

- 8,975 expressed loci belonging to turquoise-assigned ERV subfamilies
- 2,528 expressed loci belonging to yellow-assigned ERV subfamilies

The Telescope loci themselves were not used to construct WGCNA1; module labels refer to the WGCNA assignment of the corresponding ERV subfamily.

---

## 6. Telescope locus-level age association

A joint DESeq2 model was fitted to host genes and Telescope loci using:

```text
~ sex + age
```

Multiple-testing correction for the locus analysis was then applied across the 11,503 Telescope loci.

### Main result

Only one ERV locus was significantly associated with age after FDR correction:

**MER41_5q14.3a**

- ERV subfamily: MER41-int
- chromosome: chr5
- coordinates: 89,917,678–89,919,657 (GRCh38)
- strand: −
- WGCNA subfamily context: yellow
- baseMean: 158.06
- log2 fold change per year: 0.03258
- p = 1.85 × 10⁻⁶
- FDR = 0.0212
- estimated change per decade: +25.3%

### Genomic context of MER41_5q14.3a

GENCODE v26 annotation showed that MER41_5q14.3a does not directly overlap an annotated gene body.

Its nearest annotated gene was:

- **AC113167.1** (ENSG00000214942.4)
- nearest TSS distance: 17,013 bp

Therefore, MER41_5q14.3a is reported as an age-associated locus but not as a direct gene-overlapping or promoter-proximal ERV under the proximity definitions used in this project.

---

## 7. Module-constrained ERV locus–gene co-expression

Locus–gene co-expression was assessed using VST-transformed expression and bicor within the WGCNA-defined module context.

This was not a second WGCNA analysis. Instead, Telescope loci were correlated with genes belonging to the module associated with the corresponding ERV subfamily.

### Yellow module

- 2,528 ERV loci × 1,437 genes
- 3,632,736 locus–gene pairs
- top 0.1% by absolute bicor: 3,633 pairs
- top pair: MER4_4q31.3b – ASPM, bicor = 0.9366

### Turquoise module

- 8,975 ERV loci × 5,599 genes
- 50,251,025 locus–gene pairs
- top 0.1% by absolute bicor: 50,252 pairs
- top pair: HERVK11D_5p13.2 – UNC80, bicor = 0.9748

The top 0.1% threshold was used for exploratory candidate prioritisation and is not interpreted as a statistical-significance threshold.

---

## 8. Genomic annotation of strong ERV–gene associations

The 53,885 strong co-expression pairs were annotated with GenomicRanges using the same GENCODE v26 gene annotation used in the RNA-seq pipeline.

For each pair, the analysis assessed:

- chromosome concordance
- direct gene-body overlap
- exon versus intron overlap
- distance to the paired gene TSS when non-overlapping
- whether the paired gene was also the nearest annotated gene

Summary:

| Metric | Number |
|---|---:|
| Strong co-expression pairs | 53,885 |
| Same chromosome | 3,307 |
| Direct gene-body overlap | 9 |
| Exonic overlap | 4 |
| Intronic overlap | 5 |
| Non-overlapping pairs within 10 kb of paired-gene TSS | 0 |
| Paired gene was nearest annotated gene | 31 |
| Nearest but non-overlapping pairs | 22 |

The 22 nearest but non-overlapping pairs were all more than 10 kb from the paired-gene TSS (minimum 24,233 bp), so they were not treated as promoter-proximal candidates.

---

## 9. Final ERV–gene candidate associations

Nine ERV–gene candidate associations combined strong expression co-variation with direct gene-body overlap.

| ERV locus | ERV subfamily | Overlapping gene | bicor | Genomic context |
|---|---|---|---:|---|
| HERVH_5q32b | HERVH-int | SCGB3A2 | 0.944 | Exonic |
| HERVL_18p11.21a | HERVL-int | PIEZO2 | 0.922 | Intronic |
| HERVH_5q35.1a | HERVH-int | GABRP | 0.914 | Exonic |
| HERVH_19p13.3b | HERVH-int | FUT3 | 0.912 | Exonic |
| ERV316A3_8p11.22b | ERV3-16A3_I-int | ADAM18 | 0.904 | Intronic |
| HERVL_3q13.31 | HERVL-int | LSAMP | 0.904 | Intronic |
| HERVH_15q26.3b | HERVH-int | LRRK1 | 0.903 | Exonic |
| HERVL40_12q21.31 | HERVL40-int | LRRIQ1 | 0.896 | Intronic |
| PABLA_7q21.13 | PABL_A-int | CFAP69 | 0.873 | Intronic |

These are candidate associations, not evidence of causal regulation. Direct genomic overlap plus co-expression increases biological plausibility but does not establish that the ERV controls the overlapping gene.

---

## 10. Overall interpretation

The analysis did not identify a broad age-associated ERV programme at the WGCNA-module or TEcount-subfamily level. However, ERV subfamilies showed a highly structured network distribution, with 99.1% concentrated in two WGCNA modules with distinct biological annotations.

Locus-resolved Telescope analysis added two complementary findings:

1. **MER41_5q14.3a** was the only individual ERV locus significantly associated with age after FDR correction.
2. An independent co-expression-plus-genomic-overlap strategy identified **nine ERV–gene candidate associations** with strong bicor values and direct exonic or intronic overlap.

The age-associated MER41 locus was not one of the nine gene-overlapping candidates. These analyses therefore capture different properties of ERV loci: direct association with age versus strong local ERV–gene coupling.

---

## 11. Relationship to the GTEx analysis

The GTEx and GSE164471 analyses are retained as two complementary studies rather than formally integrated datasets.

- GTEx provides a large-sample gene-expression network and genomic ERV-context analysis.
- GSE164471 provides direct ERV expression quantification at subfamily and individual-locus resolution.

GTEx was not used as an ERV-expression replication cohort because repeat-aware ERV expression could not be quantified from the available GTEx input used in this project. Cross-study interpretation is therefore performed at the dissertation Discussion level rather than through a formal statistical integration pipeline.

---

## 12. Main limitations

- GSE164471 contains 53 samples, limiting power for age-association testing.
- TEcount subfamily-level quantification may obscure locus-specific behaviour.
- Telescope locus assignment depends on the selected HERV reference and multimapping-read model.
- The top 0.1% bicor threshold is an exploratory prioritisation criterion, not an FDR-controlled significance threshold.
- Co-expression and genomic proximity/overlap do not demonstrate causal ERV regulation of host genes.
- GENCODE v26 was retained for annotation consistency with the original RNA-seq pipeline, although it is not the most recent annotation release.
