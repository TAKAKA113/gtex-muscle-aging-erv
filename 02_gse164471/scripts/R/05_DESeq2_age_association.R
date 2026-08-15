# 05_DESeq2_age_association.R
#
# 目的：
#   低発現filter後のgene・strict ERV subfamilyについて、
#   性別を調整した年齢関連解析を行う。
#
# model：
#   raw count ~ sex + age
#
# さらに、
#   DESeq2 age association
#   ＋ WGCNA1 module membership
# を統合し、年齢関連network hub候補を確認する。
#
# 最後に、WGCNA2で使用可能な複数の
# age-associated feature setを保存する。


# 1. Set project and load package ----

project_root <- paste0(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/",
  "gse164471_erv_aging"
)

setwd(project_root)

library(DESeq2)

options(stringsAsFactors = FALSE)


# 出力directoryが存在しない場合は作成する。
dir.create(
  file.path(project_root, "results", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(project_root, "results", "objects"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(project_root, "results", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Read input objects ----

# 低発現filter後のraw countを含むDESeq2 object。
dds_age <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_dds_low_expression_filtered.rds"
  )
)


# DESeq2対象featureのannotation。
feature_annotation <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_DESeq2_feature_annotation.rds"
  )
)


# WGCNA1で計算した全featureのmodule membership。
wgcna1_membership <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_all_feature_module_membership.rds"
  )
)


# objectの大きさを確認する。
cat(
  "DDS dimensions:",
  dim(dds_age),
  "\n"
)

cat(
  "Annotation dimensions:",
  dim(feature_annotation),
  "\n"
)

cat(
  "WGCNA1 membership dimensions:",
  dim(wgcna1_membership),
  "\n"
)


# 3. Prepare DESeq2 design ----

# ageを数値として明示する。
dds_age$age <- as.numeric(
  as.character(dds_age$age)
)


# Femaleをsexの基準にする。
dds_age$sex <- relevel(
  factor(dds_age$sex),
  ref = "F"
)


# 性別を調整しながら年齢効果を解析する。
design(dds_age) <- ~ sex + age


cat(
  "DESeq2 design:",
  deparse(design(dds_age)),
  "\n"
)


cat(
  "Sex levels:",
  levels(dds_age$sex),
  "\n"
)


cat(
  "Age range:",
  range(dds_age$age),
  "\n"
)


# 4. Run DESeq2 age association ----

# 各gene・ERVにnegative binomial modelを適合する。
dds_age <- DESeq(
  dds_age
)


# 利用可能な係数名を確認する。
cat(
  "\nDESeq2 coefficient names:\n"
)

print(
  resultsNames(dds_age)
)


# 性別調整後のage coefficientを取り出す。
age_results <- results(
  dds_age,
  name = "age",
  alpha = 0.05
)


# data frameへ変換する。
age_results_df <- as.data.frame(
  age_results
)


# feature IDを通常の列として追加する。
age_results_df$feature_id <- row.names(
  age_results_df
)


# 列順を整理する。
age_results_df <- age_results_df[
  ,
  c(
    "feature_id",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )
]


# 5. Add interpretable age-effect columns ----

# DESeq2のage coefficientは、
# 年齢が1歳増加したときのlog2 fold changeである。
age_results_df$log2FC_per_10_years <- (
  age_results_df$log2FoldChange * 10
)


# 10年間のfold changeへ変換する。
age_results_df$fold_change_per_10_years <- (
  2 ^ age_results_df$log2FC_per_10_years
)


# 10年間の発現変化率へ変換する。
#
# 例：
#   20  → 10年間で20%上昇
#  -20  → 10年間で20%低下
age_results_df$percent_change_per_10_years <- (
  (
    age_results_df$fold_change_per_10_years - 1
  ) * 100
)


# 6. Add feature annotation ----

annotation_index <- match(
  age_results_df$feature_id,
  feature_annotation$feature_id
)


# annotationに対応しないfeatureがないか確認する。
cat(
  "\nUnmatched annotations:",
  sum(is.na(annotation_index)),
  "\n"
)


age_results_annotated <- cbind(
  age_results_df,
  feature_annotation[
    annotation_index,
    setdiff(
      colnames(feature_annotation),
      "feature_id"
    ),
    drop = FALSE
  ]
)


# FDRとp値が小さい順に並べる。
age_results_annotated <- age_results_annotated[
  order(
    age_results_annotated$padj,
    age_results_annotated$pvalue,
    na.last = TRUE
  ),
]


# 7. Classify age-association results ----

age_results_annotated$age_direction <- "not_significant"


age_results_annotated$age_direction[
  !is.na(age_results_annotated$padj) &
    age_results_annotated$padj < 0.05 &
    age_results_annotated$log2FoldChange > 0
] <- "age_up"


age_results_annotated$age_direction[
  !is.na(age_results_annotated$padj) &
    age_results_annotated$padj < 0.05 &
    age_results_annotated$log2FoldChange < 0
] <- "age_down"


age_results_annotated$age_direction[
  is.na(age_results_annotated$pvalue)
] <- "not_testable"


# FDR < 0.05。
age_significant <- age_results_annotated[
  !is.na(age_results_annotated$padj) &
    age_results_annotated$padj < 0.05,
]


# 年齢とともに上昇するfeature。
age_up <- age_significant[
  age_significant$log2FoldChange > 0,
]


# 年齢とともに低下するfeature。
age_down <- age_significant[
  age_significant$log2FoldChange < 0,
]


# nominal p < 0.05。
#
# WGCNA2候補としては使えるが、
# 正式な有意featureとは扱わない。
age_nominal <- age_results_annotated[
  !is.na(age_results_annotated$pvalue) &
    age_results_annotated$pvalue < 0.05,
]


# 8. Summarise DESeq2 results ----

result_summary <- data.frame(
  metric = c(
    "Total tested features",
    "Features with non-NA p-value",
    "FDR < 0.05 features",
    "FDR < 0.05 age-up features",
    "FDR < 0.05 age-down features",
    "Nominal p < 0.05 features",
    "FDR < 0.05 genes",
    "FDR < 0.05 ERVs",
    "Nominal p < 0.05 genes",
    "Nominal p < 0.05 ERVs"
  ),
  value = c(
    nrow(age_results_annotated),
    
    sum(
      !is.na(age_results_annotated$pvalue)
    ),
    
    nrow(age_significant),
    
    nrow(age_up),
    
    nrow(age_down),
    
    nrow(age_nominal),
    
    sum(
      age_significant$feature_type == "GENE",
      na.rm = TRUE
    ),
    
    sum(
      age_significant$feature_type == "TE",
      na.rm = TRUE
    ),
    
    sum(
      age_nominal$feature_type == "GENE",
      na.rm = TRUE
    ),
    
    sum(
      age_nominal$feature_type == "TE",
      na.rm = TRUE
    )
  )
)


cat(
  "\nDESeq2 age-association summary:\n"
)

print(
  result_summary,
  row.names = FALSE
)


# 9. Extract ERV-specific results ----

# 全ERVのage association結果。
erv_age_results <- age_results_annotated[
  age_results_annotated$feature_type == "TE",
]


# p値が小さい順に並べる。
erv_age_results <- erv_age_results[
  order(
    erv_age_results$pvalue,
    na.last = TRUE
  ),
]


# FDR < 0.05のERV。
erv_age_significant <- erv_age_results[
  !is.na(erv_age_results$padj) &
    erv_age_results$padj < 0.05,
]


# nominal p < 0.05のERV。
erv_age_nominal <- erv_age_results[
  !is.na(erv_age_results$pvalue) &
    erv_age_results$pvalue < 0.05,
]


cat(
  "\nStrict ERV age-association summary:\n"
)

cat(
  "Total tested ERVs:",
  nrow(erv_age_results),
  "\n"
)

cat(
  "ERVs with FDR < 0.05:",
  nrow(erv_age_significant),
  "\n"
)

cat(
  "ERVs with nominal p < 0.05:",
  nrow(erv_age_nominal),
  "\n"
)


# 上位20 ERVsをconsoleに表示する。
cat(
  "\nTop 20 ERVs ordered by age p-value:\n"
)

print(
  head(
    erv_age_results[
      ,
      c(
        "feature_id",
        "baseMean",
        "log2FC_per_10_years",
        "percent_change_per_10_years",
        "pvalue",
        "padj"
      )
    ],
    20
  ),
  row.names = FALSE
)


# 10. Save DESeq2 results ----

write.table(
  age_results_annotated,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_all_features.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  age_significant,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_FDR05_features.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  age_nominal,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_nominal_p05_features.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  erv_age_results,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_all_strict_ERVs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  erv_age_significant,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_FDR05_strict_ERVs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  erv_age_nominal,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_nominal_p05_strict_ERVs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  result_summary,
  file.path(
    project_root,
    "results",
    "tables",
    "DESeq2_age_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


saveRDS(
  dds_age,
  file.path(
    project_root,
    "results",
    "objects",
    "DESeq2_age_dds.rds"
  )
)


saveRDS(
  age_results_annotated,
  file.path(
    project_root,
    "results",
    "objects",
    "DESeq2_age_all_features.rds"
  )
)


# 11. Create DESeq2 MA plot ----

pdf(
  file.path(
    project_root,
    "results",
    "figures",
    "DESeq2_age_MA_plot.pdf"
  ),
  width = 7,
  height = 6
)


plotMA(
  age_results,
  alpha = 0.05,
  main = "Age association: gene and strict ERV features"
)


dev.off()


# 12. Create volcano plot ----

volcano_data <- age_results_annotated[
  !is.na(age_results_annotated$pvalue) &
    !is.na(age_results_annotated$log2FC_per_10_years),
]


volcano_colour <- ifelse(
  !is.na(volcano_data$padj) &
    volcano_data$padj < 0.05,
  "red",
  "grey70"
)


pdf(
  file.path(
    project_root,
    "results",
    "figures",
    "DESeq2_age_volcano_plot.pdf"
  ),
  width = 7,
  height = 6
)


plot(
  volcano_data$log2FC_per_10_years,
  -log10(volcano_data$pvalue),
  pch = 16,
  cex = 0.6,
  col = volcano_colour,
  xlab = "log2 fold change per 10 years",
  ylab = "-log10(p-value)",
  main = "Feature-wise age association"
)


abline(
  h = -log10(0.05),
  lty = 2
)


abline(
  v = 0,
  lty = 2
)


legend(
  "topright",
  legend = c(
    "FDR < 0.05",
    "Not significant"
  ),
  col = c(
    "red",
    "grey70"
  ),
  pch = 16,
  bty = "n"
)


dev.off()


# 13. Integrate DESeq2 results with WGCNA1 membership ----

age_index <- match(
  wgcna1_membership$feature_id,
  age_results_annotated$feature_id
)


cat(
  "\nWGCNA1 features without DESeq2 result:",
  sum(is.na(age_index)),
  "\n"
)


# WGCNA1 tableにはannotationがすでに含まれているため、
# DESeq2結果の主要列だけ追加する。
wgcna1_age_integrated <- cbind(
  wgcna1_membership,
  
  age_results_annotated[
    age_index,
    c(
      "baseMean",
      "log2FoldChange",
      "log2FC_per_10_years",
      "fold_change_per_10_years",
      "percent_change_per_10_years",
      "stat",
      "pvalue",
      "padj",
      "age_direction"
    ),
    drop = FALSE
  ]
)


# FDR < 0.05かつgrey以外。
wgcna1_age_significant <- wgcna1_age_integrated[
  wgcna1_age_integrated$module != "grey" &
    !is.na(wgcna1_age_integrated$padj) &
    wgcna1_age_integrated$padj < 0.05,
]


# FDR < 0.05かつ|MM| >= 0.8。
#
# 現段階では最終hubではなく、
# age-associated network hub candidateとする。
wgcna1_age_high_mm <- wgcna1_age_integrated[
  wgcna1_age_integrated$module != "grey" &
    !is.na(wgcna1_age_integrated$padj) &
    wgcna1_age_integrated$padj < 0.05 &
    abs(wgcna1_age_integrated$MM) >= 0.8,
]


# FDR < 0.05かつ|MM| >= 0.8のERV。
wgcna1_age_high_mm_erv <- wgcna1_age_high_mm[
  wgcna1_age_high_mm$feature_type == "TE",
]


# 候補をFDR、次にMMの順で並べる。
wgcna1_age_high_mm <- wgcna1_age_high_mm[
  order(
    wgcna1_age_high_mm$padj,
    -abs(wgcna1_age_high_mm$MM),
    na.last = TRUE
  ),
]


wgcna1_age_high_mm_erv <- wgcna1_age_high_mm_erv[
  order(
    wgcna1_age_high_mm_erv$padj,
    -abs(wgcna1_age_high_mm_erv$MM),
    na.last = TRUE
  ),
]


cat(
  "\nWGCNA1 and DESeq2 integration:\n"
)

cat(
  "FDR < 0.05 features in non-grey modules:",
  nrow(wgcna1_age_significant),
  "\n"
)

cat(
  "FDR < 0.05 and |MM| >= 0.8:",
  nrow(wgcna1_age_high_mm),
  "\n"
)

cat(
  "FDR < 0.05 and |MM| >= 0.8 ERVs:",
  nrow(wgcna1_age_high_mm_erv),
  "\n"
)


# 14. Save WGCNA1–DESeq2 integrated results ----

write.table(
  wgcna1_age_integrated,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_DESeq2_age_integrated.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  wgcna1_age_high_mm,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_age_associated_high_MM_candidates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  wgcna1_age_high_mm_erv,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_age_associated_high_MM_ERV_candidates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


saveRDS(
  wgcna1_age_integrated,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_DESeq2_age_integrated.rds"
  )
)


# 15. Prepare candidate feature sets for WGCNA2 ----

# FDR < 0.05。
feature_set_fdr05 <- age_results_annotated$feature_id[
  !is.na(age_results_annotated$padj) &
    age_results_annotated$padj < 0.05
]


# FDR < 0.10。
feature_set_fdr10 <- age_results_annotated$feature_id[
  !is.na(age_results_annotated$padj) &
    age_results_annotated$padj < 0.10
]


# nominal p < 0.05。
feature_set_p05 <- age_results_annotated$feature_id[
  !is.na(age_results_annotated$pvalue) &
    age_results_annotated$pvalue < 0.05
]


# p値上位1,000 features。
valid_pvalue_results <- age_results_annotated[
  !is.na(age_results_annotated$pvalue),
]


valid_pvalue_results <- valid_pvalue_results[
  order(valid_pvalue_results$pvalue),
]


feature_set_top1000 <- head(
  valid_pvalue_results$feature_id,
  1000
)


# p値上位2,000 features。
feature_set_top2000 <- head(
  valid_pvalue_results$feature_id,
  2000
)


wgcna2_feature_sets <- list(
  FDR_0.05 = feature_set_fdr05,
  FDR_0.10 = feature_set_fdr10,
  nominal_p_0.05 = feature_set_p05,
  top_1000_by_pvalue = feature_set_top1000,
  top_2000_by_pvalue = feature_set_top2000
)


# 各feature setに含まれるGene・ERV数を集計する。
feature_set_summary <- do.call(
  rbind,
  lapply(
    names(wgcna2_feature_sets),
    function(set_name) {
      
      feature_ids <- wgcna2_feature_sets[[set_name]]
      
      feature_index <- match(
        feature_ids,
        feature_annotation$feature_id
      )
      
      set_annotation <- feature_annotation[
        feature_index,
      ]
      
      data.frame(
        feature_set = set_name,
        total_features = length(feature_ids),
        
        genes = sum(
          set_annotation$feature_type == "GENE",
          na.rm = TRUE
        ),
        
        ERVs = sum(
          set_annotation$feature_type == "TE",
          na.rm = TRUE
        )
      )
    }
  )
)


cat(
  "\nCandidate feature sets for WGCNA2:\n"
)

print(
  feature_set_summary,
  row.names = FALSE
)


saveRDS(
  wgcna2_feature_sets,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA2_candidate_feature_sets.rds"
  )
)


write.table(
  feature_set_summary,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA2_candidate_feature_set_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 16. Final message ----

cat(
  "\n05_DESeq2_age_association.R completed successfully.\n"
)

cat(
  "Full DESeq2 results:\n",
  "results/tables/DESeq2_age_all_features.tsv\n"
)

cat(
  "ERV results:\n",
  "results/tables/DESeq2_age_all_strict_ERVs.tsv\n"
)

cat(
  "WGCNA1 integration:\n",
  "results/tables/WGCNA1_DESeq2_age_integrated.tsv\n"
)

cat(
  "WGCNA2 candidate sets:\n",
  "results/tables/WGCNA2_candidate_feature_set_summary.tsv\n"
)