# GSE164471 WGCNA1 summary

TEcount出力
1,330 TE subfamilies
        ↓ strict ERV filter
478 strict ERV subfamilies
        ↓ low-expression filter
465 ERV subfamilies
        ↓ VST + MAD > 0
465 ERV subfamilies
        ↓
WGCNA1






## Dataset and input

- Healthy human skeletal muscle RNA-seq: **53 samples**
- Age range: **22–83 years**
- Quantification: protein-coding genes + strict ERV **subfamilies**

```text
Initial TEcount matrix
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

## WGCNA settings

| Parameter | Value |
|---|---:|
| Samples | 53 |
| Correlation | bicor |
| Network | signed |
| TOM | signed |
| Soft-threshold power | 22 |
| Minimum module size | 30 |
| Merge cut height | 0.25 |

## Module composition

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

## ERV distribution

| Module | ERV subfamilies | Percentage of all ERVs |
|---|---:|---:|
| turquoise | 397 | 85.4% |
| yellow | 64 | 13.8% |
| grey | 4 | 0.9% |

**461 of 465 ERV subfamilies (99.1%) were assigned to the turquoise or yellow module.**

## Association with age

Module eigengenes were tested using:

```text
module eigengene ~ age + sex
```

- Age-associated modules at **FDR < 0.05: 0**
- Lowest raw age p-value: **red module**
  - age beta = **−0.00209**
  - p = **0.0349**
  - FDR = **0.1819**
- The ERV-rich turquoise and yellow modules were **not significantly associated with age**.

## Feature-level age regression

All low-expression-filtered features were tested using:

```text
raw count ~ sex + age
```

| Result | Genes | ERV subfamilies |
|---|---:|---:|
| FDR < 0.05 | 67 | 0 |
| Non-grey WGCNA features with FDR < 0.05 | 27 | 0 |
| FDR < 0.05 and \|MM\| ≥ 0.8 | 4 | 0 |

## Main result

ERV subfamily expression formed a highly structured co-expression network in skeletal muscle, with **99.1% of ERV subfamilies concentrated in two modules**. However, neither ERV-rich module showed a significant association with age, and no individual ERV subfamily was significant in feature-level age regression after multiple-testing correction.

## Limitation

TEcount aggregates reads at the **ERV subfamily level**. Opposing age-related changes among individual loci within the same subfamily may therefore be masked. Telescope locus-level analysis is used as a follow-up.