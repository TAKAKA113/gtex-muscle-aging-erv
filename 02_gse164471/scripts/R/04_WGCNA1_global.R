# 04_WGCNA1_global.R
#
# Purpose:
#   Protein-coding genesとstrict ERV subfamiliesを統合した
#   global co-expression networkをWGCNAで構築する。
#
#   Network構築にはageを使用せず、
#   module作成後にmodule eigengeneとage・sexの関連を評価する。


# 1. Project directory and packages ----

project_root <- "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
setwd(project_root)

library(WGCNA)

# 文字列を自動的にfactorへ変換しない。
options(stringsAsFactors = FALSE)


# 2. Read WGCNA input data ----

datExpr <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_WGCNA_input.rds"
  )
)

sample_metadata <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_sample_metadata.rds"
  )
)

feature_annotation <- readRDS(
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_WGCNA_feature_annotation.rds"
  )
)


# Data dimensions
dim(datExpr)
dim(sample_metadata)
dim(feature_annotation)


# 3. Check WGCNA input ----

# goodSamplesGenes()は、
# missing valueが多いsampleやfeature、
# 変動を計算できないfeatureがないか確認する。
wgcna_check <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

wgcna_check$allOK




# 4. Select soft-threshold power ----

# R²がpower 20でも上昇中だったため、
# power 30まで候補を広げて確認する。
powers <- c(
  1:10,
  seq(12, 30, by = 2)
)

soft_power_result <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.1),
  verbose = 5
)

soft_power_result$fitIndices[, c(
  "Power",
  "SFT.R.sq",
  "slope",
  "mean.k."
)]



# 5. Plot soft-threshold results ----

soft_power <- 22

soft_power_file <- file.path(
  project_root,
  "results",
  "figures",
  "TEcount_WGCNA_soft_threshold.pdf"
)

pdf(
  soft_power_file,
  width = 10,
  height = 5
)

par(mfrow = c(1, 2))

# Scale-free topology fit
plot(
  soft_power_result$fitIndices$Power,
  -sign(soft_power_result$fitIndices$slope) *
    soft_power_result$fitIndices$SFT.R.sq,
  xlab = "Soft-threshold power",
  ylab = "Signed scale-free topology fit (R²)",
  type = "n",
  main = "Scale independence"
)

text(
  soft_power_result$fitIndices$Power,
  -sign(soft_power_result$fitIndices$slope) *
    soft_power_result$fitIndices$SFT.R.sq,
  labels = soft_power_result$fitIndices$Power
)

abline(h = 0.8, lty = 2)
abline(v = soft_power, lty = 2)


# Mean connectivity
plot(
  soft_power_result$fitIndices$Power,
  soft_power_result$fitIndices$mean.k.,
  xlab = "Soft-threshold power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  soft_power_result$fitIndices$Power,
  soft_power_result$fitIndices$mean.k.,
  labels = soft_power_result$fitIndices$Power
)

abline(v = soft_power, lty = 2)

dev.off()


#soft-threshold power = 22

# 6. Detect co-expression modules ----

# blockwiseModules()は、
# feature間の相関計算からmodule検出までを一括して行う。
#
# maxBlockSize = 5000：
#   18,230 featuresを複数blockに分けて計算し、
#   memory使用量を抑える。
#   featureを除外する設定ではない。
#
# minModuleSize = 30：
#   30 features未満の小さなgroupはmoduleとして扱わない。
#
# mergeCutHeight = 0.25：
#   発現パターンが非常に似たmodule同士を統合する。
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  corType = "bicor",
  maxPOutliers = 0.1,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  maxBlockSize = 5000,
  numericLabels = TRUE,
  randomSeed = 12345,
  verbose = 3
)


# module番号をWGCNAの色名へ変換する。
module_colors <- labels2colors(net$colors)

# 各moduleに含まれるfeature数を確認する。
sort(
  table(module_colors),
  decreasing = TRUE
)


# module_colors
# turquoise      blue      grey     brown    yellow     green       red     black      pink 
# 5996      3994      2744      1667      1501       783       520       502       410 
# magenta 
# 113 

# Moduleに割り当て済み：15,486 features（約84.9%）
# Grey                  ： 2,744 features（約15.1%）




# 7. Summarise genes and ERVs in each module ----

# 各featureにmodule情報を付ける。
module_table <- data.frame(
  feature_id = colnames(datExpr),
  module = module_colors
)

# feature annotationからGENE / TE分類を対応させる。
annotation_order <- match(
  module_table$feature_id,
  feature_annotation$feature_id
)

module_table$feature_type <- feature_annotation$feature_type[
  annotation_order
]

# 各moduleのGene数とERV数を確認する。
module_composition <- table(
  module_table$module,
  module_table$feature_type
)

module_composition

# GENE   TE
# black      502    0
# blue      3994    0
# brown     1667    0
# green      783    0
# grey      2740    4
# magenta    113    0
# pink       410    0
# red        520    0
# turquoise 5599  397
# yellow    1437   64

#465個のstrict ERVは、ほぼ2つのmoduleに集中
#greyはmoduleに配置されなかったやつ


# 8. Calculate module eigengenes ----

# 各moduleの代表的な発現パターンを計算する。
#
# excludeGrey = TRUE：
#   greyは共発現moduleではないため除外する。
module_eigengenes <- moduleEigengenes(
  datExpr,
  colors = module_colors,
  excludeGrey = TRUE
)$eigengenes


# module eigengeneの列順を整理する。
module_eigengenes <- orderMEs(module_eigengenes)


# 大きさとmodule名を確認する。
dim(module_eigengenes)

colnames(module_eigengenes)


# sample metadataとの行順が一致しているか確認する。
identical(
  row.names(module_eigengenes),
  row.names(sample_metadata)
)


# 9. Test module associations with age and sex ----
#ME(odule eigengene)とはかくmoduleが2つの変化する形質と関係しているかを示す


# 線形モデルで使用するsample情報を準備する。
trait_data <- data.frame(
  age = as.numeric(sample_metadata$age),
  sex = factor(
    sample_metadata$sex,
    levels = c("F", "M")
  ),
  row.names = row.names(sample_metadata)
)


# 各moduleの解析結果を入れる空のdata frameを作る。
module_trait_results <- data.frame()


# 9つのmoduleを1つずつ解析する。
for (me_name in colnames(module_eigengenes)) {
  
  # 1つのmodule eigengeneとmetadataをまとめる。
  model_data <- data.frame(
    eigengene = module_eigengenes[, me_name],
    age = trait_data$age,
    sex = trait_data$sex
  )
  
  # 性別を考慮しながら、年齢との関連を調べる。
  model <- lm(
    eigengene ~ age + sex,
    data = model_data
  )
  
  # 係数とp値を取り出す。
  model_summary <- summary(model)$coefficients
  
  # moduleごとの結果を追加する。
  module_trait_results <- rbind(
    module_trait_results,
    data.frame(
      module = sub("^ME", "", me_name),
      age_beta = model_summary["age", "Estimate"],
      age_p = model_summary["age", "Pr(>|t|)"],
      sex_M_beta = model_summary["sexM", "Estimate"],
      sex_p = model_summary["sexM", "Pr(>|t|)"]
    )
  )
}


# 9 modulesのage p値を多重検定補正する。
module_trait_results$age_FDR <- p.adjust(
  module_trait_results$age_p,
  method = "BH"
)


# 年齢との関連が強い順に並べる。
module_trait_results <- module_trait_results[
  order(module_trait_results$age_p),
]


module_trait_results


# 10. Save global WGCNA results ----

write.table(
  module_trait_results,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_module_age_sex_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  net,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_network.rds"
  )
)

saveRDS(
  module_table,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_module_assignment.rds"
  )
)

saveRDS(
  module_eigengenes,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_module_eigengenes.rds"
  )
)


# 11. Create module–trait correlation heatmap ----

# Heatmap表示用にtraitを数値化する。
#
# Age：
#   実年齢
#
# Sex_M：
#   Female = 0
#   Male   = 1
trait_numeric <- data.frame(
  Age = as.numeric(sample_metadata$age),
  Sex_M = ifelse(sample_metadata$sex == "M", 1, 0),
  row.names = row.names(sample_metadata)
)


# module eigengeneとage・sexのbicorおよびp値を計算する。
module_trait_test <- bicorAndPvalue(
  module_eigengenes,
  trait_numeric,
  maxPOutliers = 0.1
)

module_trait_cor <- module_trait_test$bicor
module_trait_p <- module_trait_test$p


# 各cellに、
# 相関係数
# (p値)
# を表示する。
heatmap_text <- matrix(
  paste0(
    sprintf("%.2f", module_trait_cor),
    "\n(p=",
    formatC(
      module_trait_p,
      format = "e",
      digits = 1
    ),
    ")"
  ),
  nrow = nrow(module_trait_cor),
  ncol = ncol(module_trait_cor)
)


heatmap_file <- file.path(
  project_root,
  "results",
  "figures",
  "WGCNA1_module_trait_heatmap.pdf"
)

pdf(
  heatmap_file,
  width = 7,
  height = 7
)

par(
  mar = c(6, 8, 3, 3)
)

labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = colnames(module_trait_cor),
  yLabels = sub(
    "^ME",
    "",
    row.names(module_trait_cor)
  ),
  ySymbols = sub(
    "^ME",
    "",
    row.names(module_trait_cor)
  ),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = heatmap_text,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module–trait correlations"
)

dev.off()


saveRDS(
  module_trait_cor,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_module_trait_correlation.rds"
  )
)

saveRDS(
  module_trait_p,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_module_trait_pvalue.rds"
  )
)



# 12. Create FDR-adjusted module–trait heatmap ----

# AgeとSexについて、それぞれ9 modules間で
# BH法による多重検定補正を行う。
module_trait_fdr <- module_trait_p

module_trait_fdr[, "Age"] <- p.adjust(
  module_trait_p[, "Age"],
  method = "BH"
)

module_trait_fdr[, "Sex_M"] <- p.adjust(
  module_trait_p[, "Sex_M"],
  method = "BH"
)


# FDR < 0.05のcellに*を付ける。
significance_mark <- ifelse(
  module_trait_fdr < 0.05,
  "*",
  ""
)


# 相関係数とFDRを表示する文字列を作る。
heatmap_text_fdr <- matrix(
  paste0(
    sprintf("%.2f", module_trait_cor),
    significance_mark,
    "\n(FDR=",
    formatC(
      module_trait_fdr,
      format = "e",
      digits = 1
    ),
    ")"
  ),
  nrow = nrow(module_trait_cor),
  ncol = ncol(module_trait_cor)
)


fdr_heatmap_file <- file.path(
  project_root,
  "results",
  "figures",
  "WGCNA1_module_trait_heatmap_FDR.pdf"
)

pdf(
  fdr_heatmap_file,
  width = 7,
  height = 7
)

par(
  mar = c(6, 8, 3, 3)
)

labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = colnames(module_trait_cor),
  yLabels = sub("^ME", "", row.names(module_trait_cor)),
  ySymbols = sub("^ME", "", row.names(module_trait_cor)),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = heatmap_text_fdr,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module–trait correlations"
)

dev.off()


# 13. Calculate module membership and extract network hubs ----

# 各featureと9つのmodule eigengeneのbicorを計算する。
#
# この相関がModule Membership（MM）である。
# MMが高いfeatureほど、そのmoduleの代表的な発現パターンを示す。
mm_test <- bicorAndPvalue(
  datExpr,
  module_eigengenes,
  maxPOutliers = 0.1
)

mm_matrix <- mm_test$bicor
mm_p_matrix <- mm_test$p


# module情報をhub解析用tableへコピーする。
module_hub_table <- module_table


# 各module名に対応するeigengene列名を作る。
#
# red       → MEred
# turquoise → MEturquoise
own_me_name <- paste0(
  "ME",
  module_hub_table$module
)


# 対応するmodule eigengeneの列番号を取得する。
#
# greyにはeigengeneを作っていないためNAになる。
own_me_index <- match(
  own_me_name,
  colnames(mm_matrix)
)


# 結果を入れる列を作る。
module_hub_table$MM <- NA_real_
module_hub_table$MM_p <- NA_real_


# grey以外のfeatureだけMMを取り出す。
valid_feature <- !is.na(own_me_index)

module_hub_table$MM[valid_feature] <- mm_matrix[
  cbind(
    which(valid_feature),
    own_me_index[valid_feature]
  )
]

module_hub_table$MM_p[valid_feature] <- mm_p_matrix[
  cbind(
    which(valid_feature),
    own_me_index[valid_feature]
  )
]


# feature annotationの行をmodule tableに合わせる。
annotation_index <- match(
  module_hub_table$feature_id,
  feature_annotation$feature_id
)


# すでに存在するfeature_idとfeature_type以外のannotationを追加する。
annotation_columns <- setdiff(
  colnames(feature_annotation),
  c("feature_id", "feature_type")
)

module_hub_table <- cbind(
  module_hub_table,
  feature_annotation[
    annotation_index,
    annotation_columns,
    drop = FALSE
  ]
)


# greyを除外する。
module_hub_ranked <- module_hub_table[
  module_hub_table$module != "grey" &
    !is.na(module_hub_table$MM),
]


# moduleごとにMMが高い順に並べる。
module_hub_ranked <- module_hub_ranked[
  order(
    module_hub_ranked$module,
    -module_hub_ranked$MM
  ),
]


# 各module内でMM順位を付ける。
module_hub_ranked$rank_in_module <- ave(
  module_hub_ranked$MM,
  module_hub_ranked$module,
  FUN = function(x) {
    rank(
      -x,
      ties.method = "first"
    )
  }
)


# 各moduleのtop 30 network hub featuresを抽出する。
module_hub_top30 <- module_hub_ranked[
  module_hub_ranked$rank_in_module <= 30,
]


# moduleごとの件数を確認する。
table(module_hub_top30$module)


# ERVだけを抽出する。
erv_hub_ranked <- module_hub_ranked[
  module_hub_ranked$feature_type == "TE",
]


# 各module内でERVだけのMM順位を付ける。
erv_hub_ranked$rank_ERV_in_module <- ave(
  erv_hub_ranked$MM,
  erv_hub_ranked$module,
  FUN = function(x) {
    rank(
      -x,
      ties.method = "first"
    )
  }
)


# 各ERV-containing moduleからtop 30 ERVsを抽出する。
erv_hub_top30 <- erv_hub_ranked[
  erv_hub_ranked$rank_ERV_in_module <= 30,
]


table(erv_hub_top30$module)


# 14. Save WGCNA1 network hub results ----

write.table(
  module_hub_ranked,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_all_feature_module_membership.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  module_hub_top30,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_top30_network_hub_features.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  erv_hub_top30,
  file.path(
    project_root,
    "results",
    "tables",
    "WGCNA1_top30_ERV_network_hubs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


saveRDS(
  module_hub_ranked,
  file.path(
    project_root,
    "results",
    "objects",
    "WGCNA1_all_feature_module_membership.rds"
  )
)



