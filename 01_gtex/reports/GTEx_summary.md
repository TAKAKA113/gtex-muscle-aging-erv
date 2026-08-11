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

ERV annotations were compared with genes in the seven AGE-associated candidate modules.

Analyses included:

```text
Promoter ERV overlap
ERV-family-specific promoter overlap
Distance from gene TSS to nearest ERV fragment
Promoter ERV enrichment/depletion by module
```

Reference values from the genomic-context analysis:

- Promoter ERV overlap among all eligible WGCNA genes: **23.86%**
- Median distance from annotated TSS to nearest ERV fragment: **3,803 bp**

The ERV families examined were:

```text
ERV1
ERVK
ERVL
ERVL-MaLR
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
