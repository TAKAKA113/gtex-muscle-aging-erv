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
The WGCNA network was gene-only, and ERV annotation was added afterwards to examine the genomic ERV context of genes in the seven AGE-associated candidate modules and their hub genes.

The downstream ERV analysis was performed in stages:

```text
Step10
AGE-associated candidate module genes
→ promoter ERV overlap
→ ERV-family-specific enrichment/depletion
→ distance from TSS to nearest ERV

Step11
210 hub-gene candidates
→ MM / GS_AGE together with ERV context
→ promoter ERV overlap among hub genes

Step12A
hub-associated ERVs
→ genic context
→ comparison across WGCNA modules

Step12B
hub-associated ERVs
→ Roadmap Epigenomics annotation
→ enhancer / regulatory chromatin context
```

### Step10 reference values

Among all eligible WGCNA genes:

```text
promoter ERV overlap = 23.86%
median distance from TSS to nearest ERV = 3,803 bp
```

ERV families examined:

```text
ERV1
ERVK
ERVL
ERVL-MaLR
```

### Why Step11–12 matter

Step10 asks whether ERV annotation is associated with an entire AGE-associated candidate module.
Step11 then asks whether the same ERV context is also present among the **network hub genes**.
Step12 further divides those hub-associated ERVs by genomic context and adds Roadmap epigenomic information.

Thus the GTEx ERV evidence progresses from:

```text
AGE-associated module
        ↓
hub gene
        ↓
nearby/promoter ERV
        ↓
ERV genomic context
        ↓
Roadmap regulatory annotation
```

This is genomic/annotation evidence rather than direct ERV-expression evidence.

### Relevance for GSE164471 integration

GTEx and GSE164471 provide complementary evidence:

```text
GTEx
Gene expression WGCNA
+ AGE-associated modules/hubs
+ ERV genomic position
+ Roadmap epigenomic annotation

GSE164471
Gene + ERV subfamily expression
+ ERV–gene co-expression
+ Telescope locus-level expression
```

The strongest cross-dataset candidates will be genes supported by both sides, for example:

```text
GTEx:
AGE-associated hub gene
+ promoter/nearby ERV
+ regulatory epigenomic annotation

AND

GSE164471:
strong co-expression with an ERV subfamily
+ Telescope locus from that subfamily near the same gene
```

Detailed Japanese explanations of each GTEx figure are available directly beside the figures:

```text
01_gtex/results/figures/FIGURE_GUIDE_JP.md
```

The guide covers:

```text
06_GenAge_enrichment_light.pdf
Step10_ERV_family_heatmap.pdf
Step10_nearest_ERV_distance_boxplot.pdf
Step10_promoter_ERV_forest.pdf
Step10_promoter_ERV_percent.pdf
Step11_hub_MM_GS_ERV_scatter.pdf
Step11_hub_promoter_ERV_percent.pdf
Step12A_hub_ERV_context_by_module.pdf
Step12A_hub_ERV_genic_context.pdf
Step12B_Roadmap_ERV_family_by_epigenome.pdf
Step12B_Roadmap_hub_enhanc...pdf
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
Step10 module-level ERV genomic context
        ↓
Step11 hub-level ERV genomic context
        ↓
Step12A genic-context annotation
        ↓
Step12B Roadmap epigenomic annotation
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
