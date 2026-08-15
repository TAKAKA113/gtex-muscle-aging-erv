# 06_module_enrichment.R
#
# 目的：
#   WGCNA1でstrict ERV subfamilyが集中したmoduleについて、
#   module内のprotein-coding geneがどのような
#   biological functionに関与しているかを調べる。

# turquoise
# 5,599 genes + 397 ERV subfamilies
# 
# yellow
# 1,437 genes + 64 ERV subfamilies
# にERV-rich moduleが得られ


#しかし「ERVが集まっている」と分かっただけでは、
#そのERV-rich networkが筋収縮なのか、ミトコンドリアなのか、免疫なのか、RNA processingなのかわからない
# 
# そのためmoduleに所属するprotein-coding genesを使ってGO enrichmentを行う
# 
# GO enrichmentとは
# Gene Ontology（GO）はGeneを生物学的機能で分類したannotation
# 
# 今回はまず Biological Process (BP) を中心に見る

# 1. Set project and load packages ----

project_root <- "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
setwd(project_root)
getwd()

# GO enrichmentを行うpackage。
library(clusterProfiler)

# Human gene annotation database。
# Ensembl ID / Entrez ID / gene symbolの変換にも使用する。
library(org.Hs.eg.db)

library(enrichplot)
library(ggplot2)
options(stringsAsFactors = FALSE)

dir.create(
  file.path(project_root, "results", "objects"),
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Read input objects ----

wgcna1_membership <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_all_feature_module_membership.rds"
  )
)

feature_annotation <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_DESeq2_feature_annotation.rds"
  )
)


# 3. Extract genes from ERV-rich modules ----

# WGCNA1でERV subfamilyが集中していた
# turquoise / yellow moduleからprotein-coding geneのみを抽出する。

turquoise_genes <- subset(
  wgcna1_membership,
  module == "turquoise" &
    feature_type == "GENE" &
    gene_biotype == "protein_coding"
)

yellow_genes <- subset(
  wgcna1_membership,
  module == "yellow" &
    feature_type == "GENE" &
    gene_biotype == "protein_coding"
)


# Gene数を確認する。
cat("\nTurquoise protein-coding genes:",
    nrow(turquoise_genes), "\n")

cat("Yellow protein-coding genes:",
    nrow(yellow_genes), "\n")

# Turquoise protein-coding genes: 5599
# Yellow protein-coding genes: 1437


cat("\nTurquoise gene IDs:\n")
print(head(
  turquoise_genes[
    ,
    c(
      "gene_id_unversioned",
      "gene_symbol",
      "module",
      "MM"
    )
  ]
))


cat("\nYellow gene IDs:\n")
print(head(
  yellow_genes[
    ,
    c(
      "gene_id_unversioned",
      "gene_symbol",
      "module",
      "MM"
    )
  ]
))


# 4. Prepare WGCNA1 background genes ----

# WGCNA1に実際に投入されたfeatureのannotation。
#
# 低発現filter後の18,291 featuresから、
# MAD = 0の61 genesを除外した後の
# WGCNA1 input featureを含むことを確認する。
wgcna1_feature_annotation <- readRDS(
  "results/objects/TEcount_WGCNA_feature_annotation.rds"
)

# objectの構造を確認する。
cat("\n===== WGCNA1 feature annotation =====\n")

cat("Dimensions:\n")
print(dim(wgcna1_feature_annotation))

cat("\nColumn names:\n")
print(colnames(wgcna1_feature_annotation))

cat("\nFeature type counts:\n")
print(table(wgcna1_feature_annotation$feature_type))

cat("\nGene biotype counts:\n")
print(table(wgcna1_feature_annotation$gene_biotype))

dim(wgcna1_feature_annotation)
# [1] 18230    10



# WGCNA1に投入されたprotein-coding geneのみを
# enrichment解析のbackgroundとして使用する。
wgcna1_background_genes <- subset(
  wgcna1_feature_annotation,
  feature_type == "GENE" &
    gene_biotype == "protein_coding"
)


# Background gene数を確認する。
cat(
  "\nWGCNA1 background protein-coding genes:",
  nrow(wgcna1_background_genes),
  "\n"
)


dim(wgcna1_feature_annotation)
# [1] 18230    10
table(wgcna1_feature_annotation$feature_type)
#GENE    TE 
#17765   465 
table(wgcna1_feature_annotation$gene_biotype)
#protein_coding 
#465          17765 



# 5. Convert Ensembl IDs to Entrez IDs ----

# WGCNA1 background genesについて、
# Ensembl gene IDからEntrez gene IDへ変換する。
#
# GO enrichmentではEnsembl IDを直接使用可能だが、
# KEGG enrichmentではEntrez IDを使用するため、
# ここで対応表を作成する。

background_id_map <- bitr(
  unique(wgcna1_background_genes$gene_id_unversioned),
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)


turquoise_id_map <- bitr(
  unique(turquoise_genes$gene_id_unversioned),
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)


yellow_id_map <- bitr(
  unique(yellow_genes$gene_id_unversioned),
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)


#Check ID conversion ----

cat("\n===== ID conversion =====\n")

cat(
  "Background:",
  length(unique(background_id_map$ENSEMBL)),
  "/",
  length(unique(wgcna1_background_genes$gene_id_unversioned)),
  "\n"
)

cat(
  "Turquoise:",
  length(unique(turquoise_id_map$ENSEMBL)),
  "/",
  length(unique(turquoise_genes$gene_id_unversioned)),
  "\n"
)

cat(
  "Yellow:",
  length(unique(yellow_id_map$ENSEMBL)),
  "/",
  length(unique(yellow_genes$gene_id_unversioned)),
  "\n"
)


#重複mappingを確認
yellow_mapping_count <- table(yellow_id_map$ENSEMBL)

yellow_mapping_count[
  yellow_mapping_count > 1
]


turquoise_mapping_count <- table(turquoise_id_map$ENSEMBL)

turquoise_mapping_count[
  turquoise_mapping_count > 1
]



length(unique(background_id_map$ENSEMBL))
length(unique(turquoise_id_map$ENSEMBL))
length(unique(yellow_id_map$ENSEMBL))


# 6. GO Biological Process enrichment ----

# Turquoise module
go_turquoise <- enrichGO(
  gene = unique(turquoise_genes$gene_id_unversioned),
  universe = unique(wgcna1_background_genes$gene_id_unversioned),
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# Yellow module
go_yellow <- enrichGO(
  gene = unique(yellow_genes$gene_id_unversioned),
  universe = unique(wgcna1_background_genes$gene_id_unversioned),
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


#check the result
cat("\n===== Turquoise GO BP =====\n")
print(head(as.data.frame(go_turquoise), 10))

cat("\nNumber of significant GO BP terms:",
    nrow(as.data.frame(go_turquoise)), "\n")


cat("\n===== Yellow GO BP =====\n")
print(head(as.data.frame(go_yellow), 10))

cat("\nNumber of significant GO BP terms:",
    nrow(as.data.frame(go_yellow)), "\n")




# 7. Reduce redundant GO terms ----

# GO ontologyには類似した親子termが多数存在するため、
# semantic similarityを用いて冗長なtermを整理する。
#
# cutoff = 0.7:
#   similarity > 0.7のGO term群から、
#   p.adjustが最も小さい代表termを残す。
#cutoff = 0.7 はclusterProfilerで示されている標準的な例
go_turquoise_simple <- simplify(
  go_turquoise,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

go_yellow_simple <- simplify(
  go_yellow,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)



cat("\n===== Simplified Turquoise GO BP =====\n")
print(
  head(
    as.data.frame(go_turquoise_simple)[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "FoldEnrichment",
        "p.adjust",
        "Count"
      )
    ],
    20
  )
)

cat(
  "\nNumber of simplified turquoise GO terms:",
  nrow(as.data.frame(go_turquoise_simple)),
  "\n"
)


cat("\n===== Simplified Yellow GO BP =====\n")
print(
  head(
    as.data.frame(go_yellow_simple)[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "FoldEnrichment",
        "p.adjust",
        "Count"
      )
    ],
    20
  )
)

cat(
  "\nNumber of simplified yellow GO terms:",
  nrow(as.data.frame(go_yellow_simple)),
  "\n"
)


# 8. Save simplified GO results ----

turquoise_go_df <- as.data.frame(go_turquoise_simple)
yellow_go_df    <- as.data.frame(go_yellow_simple)

write.table(
  turquoise_go_df,
  "results/tables/06_turquoise_GO_BP.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  yellow_go_df,
  "results/tables/06_yellow_GO_BP.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 9. Visualize GO BP enrichment ----

p_go_turquoise <- dotplot(
  go_turquoise_simple,
  showCategory = 15,
  x = "GeneRatio",
  color = "p.adjust"
) +
  ggtitle("Turquoise module: GO Biological Process")


p_go_yellow <- dotplot(
  go_yellow_simple,
  showCategory = 9,
  x = "GeneRatio",
  color = "p.adjust"
) +
  ggtitle("Yellow module: GO Biological Process")


print(p_go_turquoise)
print(p_go_yellow)


ggsave(
  "results/figures/06_turquoise_GO_BP.pdf",
  p_go_turquoise,
  width = 9,
  height = 7
)

ggsave(
  "results/figures/06_yellow_GO_BP.pdf",
  p_go_yellow,
  width = 9,
  height = 6
)


# 10. Prepare Entrez IDs for KEGG ----

background_entrez <- unique(background_id_map$ENTREZID)
turquoise_entrez  <- unique(turquoise_id_map$ENTREZID)
yellow_entrez     <- unique(yellow_id_map$ENTREZID)


cat("\n===== KEGG input =====\n")

cat(
  "Background Entrez genes:",
  length(background_entrez),
  "\n"
)

cat(
  "Turquoise Entrez genes:",
  length(turquoise_entrez),
  "\n"
)

cat(
  "Yellow Entrez genes:",
  length(yellow_entrez),
  "\n"
)


# 11. KEGG pathway enrichment ----

kegg_turquoise <- enrichKEGG(
  gene = turquoise_entrez,
  organism = "hsa",
  universe = background_entrez,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)


kegg_yellow <- enrichKEGG(
  gene = yellow_entrez,
  organism = "hsa",
  universe = background_entrez,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)


cat("\n===== Turquoise KEGG =====\n")

print(
  head(
    as.data.frame(kegg_turquoise)[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ],
    15
  )
)

cat(
  "\nNumber of significant KEGG pathways:",
  nrow(as.data.frame(kegg_turquoise)),
  "\n"
)


cat("\n===== Yellow KEGG =====\n")

print(
  head(
    as.data.frame(kegg_yellow)[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ],
    15
  )
)

cat(
  "\nNumber of significant KEGG pathways:",
  nrow(as.data.frame(kegg_yellow)),
  "\n"
)


# 12. Save KEGG results ----

turquoise_kegg_df <- as.data.frame(kegg_turquoise)
yellow_kegg_df    <- as.data.frame(kegg_yellow)

write.table(
  turquoise_kegg_df,
  "results/tables/06_turquoise_KEGG.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  yellow_kegg_df,
  "results/tables/06_yellow_KEGG.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 13. Visualize KEGG enrichment ----

p_kegg_turquoise <- dotplot(
  kegg_turquoise,
  showCategory = 15,
  x = "GeneRatio",
  color = "p.adjust"
) +
  ggtitle("Turquoise module: KEGG pathway enrichment")

print(p_kegg_turquoise)

ggsave(
  "results/figures/06_turquoise_KEGG.pdf",
  p_kegg_turquoise,
  width = 9,
  height = 7
)



