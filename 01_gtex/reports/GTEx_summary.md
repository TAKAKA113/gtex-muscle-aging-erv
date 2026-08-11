# GTEx skeletal muscle summary

## Dataset and input

- GTEx v10 skeletal muscle RNA-seq
- **818 samples**
- Gene-level RNA-seq counts

```text
Raw gene counts
        ↓ low-expression filtering
        ↓ DESeq2 VST
        ↓ select top 20,000 genes by MAD

WGCNA input
818 samples × 20,000 genes
```

## WGCNA settings

| Parameter | Value |
|---|---:|
| Samples | 818 |
| Genes | 20,000 |
| Network | signed |
| Soft-threshold power | 6 |
| Minimum module size | 30 |
| Merge cut height | 0.25 |

## Module detection

- **33 module eigengenes** were obtained.
- Grey genes: **4,658**

Major modules:

| Module | Genes |
|---|---:|
| turquoise | 2,608 |
| blue | 2,499 |
| brown | 1,506 |
| yellow | 1,205 |
| green | 740 |
| red | 630 |
| black | 629 |
| pink | 606 |
| greenyellow | 445 |
| midnightblue | 328 |
| darkgreen | 119 |

## Age-associated candidate modules

Module eigengenes were compared with AGE and SEX.

Seven modules were selected as AGE-associated candidates:

```text
greenyellow
red
darkgreen
turquoise
midnightblue
blue
pink
```

The observed AGE correlations were modest, with the largest absolute correlations around **0.35–0.39**.

## Hub-gene candidates

For each gene:

```text
MM
= correlation with its module eigengene

GS_AGE
= correlation between gene expression and AGE
```

Genes with high module membership and stronger age association were prioritised.

```text
7 target modules × top 30 genes
= 210 hub-gene candidates
```

## GenAge annotation

- WGCNA genes: **20,000**
- GenAge genes present in WGCNA: **191**
- Hub candidates: **210**
- GenAge genes among hub candidates: **3**

The three GenAge genes found among the hub candidates were:

| Module | Gene | MM | GS_AGE | Rank in module |
|---|---|---:|---:|---:|
| red | BSCL2 | 0.842 | -0.323 | 3 |
| darkgreen | S100B | 0.738 | 0.089 | 24 |
| pink | PML | 0.831 | 0.233 | 12 |

## ERV genomic-context analysis

Important: **ERV expression was not included directly in the GTEx WGCNA matrix.**
The WGCNA network was gene-only, and ERV annotation was added afterwards to examine the genomic ERV context of genes in the seven AGE-associated candidate modules.

The following analyses were performed:

```text
AGE-associated candidate module genes
        ↓
promoter ERV overlap
        ↓
ERV-family-specific promoter overlap
        ↓
distance from gene TSS to nearest ERV fragment
        ↓
module-level ERV enrichment / depletion
```

### Promoter ERV overlap

Among all eligible WGCNA genes:

```text
promoter ERV overlap = 23.86%
```

Approximate module-level promoter ERV overlap:

| Module | Genes with promoter ERV overlap |
|---|---:|
| greenyellow | ~24–25% |
| blue | ~24% |
| turquoise | ~23% |
| red | ~20–21% |
| midnightblue | ~19–20% |
| darkgreen | ~19% |
| pink | ~17% |

These values were also evaluated as enrichment/depletion using odds ratios and 95% confidence intervals.

### ERV family-specific overlap

Promoter overlap was further separated into:

```text
ERV1
ERVK
ERVL
ERVL-MaLR
```

For each module and ERV family, enrichment/depletion was summarised using log2 odds ratios, with FDR correction.

### Distance to nearest ERV

For each module gene, the distance from the annotated transcription start site (TSS) to the nearest ERV fragment was calculated.

```text
Median distance across all eligible genes = 3,803 bp
```

This was used to compare whether genes in specific AGE-associated candidate modules tended to lie closer to ERV fragments than the overall WGCNA background.

### Why this ERV information matters for integration

GTEx provides **gene-network + genomic ERV-context evidence**, whereas GSE164471 provides **ERV expression evidence**.

```text
GTEx
AGE-associated gene modules
+ promoter / nearby ERV annotation

GSE164471
TEcount ERV subfamily expression
+ Telescope locus expression
+ ERV–gene co-expression

        ↓
Cross-dataset integration
```

The most useful integration targets will be genes that are supported by both datasets, for example:

```text
GTEx: AGE-associated module gene with nearby/promoter ERV
+
GSE164471: gene strongly co-expressed with the same ERV family/subfamily
+
Telescope: corresponding ERV locus near that gene
```

A Japanese explanation of the Step10–12 ERV figures is available here:

```text
01_gtex/reports/ERV_Figure_guide_JP.md
```

## Current GTEx workflow

```text
GTEx skeletal muscle gene counts
        ↓
Filtering + VST + MAD selection
        ↓
WGCNA
        ↓
AGE / SEX module association
        ↓
7 AGE-associated candidate modules
        ↓
MM + GS_AGE hub ranking
        ↓
210 hub candidates
        ↓
GenAge annotation
        ↓
ERV genomic-context analysis
        ↓
Integration with GSE164471 ERV expression / Telescope loci
```

## Main numerical summary

| Result | Number |
|---|---:|
| Samples | 818 |
| WGCNA genes | 20,000 |
| Module eigengenes | 33 |
| AGE-associated candidate modules | 7 |
| Hub-gene candidates | 210 |
| GenAge genes in WGCNA | 191 |
| GenAge genes among hub candidates | 3 |
| Promoter ERV overlap in eligible WGCNA genes | 23.86% |
| Median nearest ERV distance | 3,803 bp |
