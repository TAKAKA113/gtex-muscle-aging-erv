# ============================================================
# 10_module_ERV_enrichment.R
#
# 目的：
# GTEx V10骨格筋WGCNAで選んだ7つのAGE関連候補moduleについて、
# promoter近傍のERV fragmentが多いか、少ないかを検定し、
# module間の違いを可視化する。
#
# 【TSSとは】
# TSS = Transcription Start Site（転写開始点）。
# RNAへの転写が始まるゲノム上の位置を指す。
#
# + strand：TSS = gene_start
# - strand：TSS = gene_end
#
# 【今回のpromoter定義】
# TSS上流2,000 bpから下流500 bpまでを
# promoter-proximal regionとして操作的に定義した。
# これは絶対的なpromoter境界ではない。
# 下流500 bpはfirst exon、5' UTR、first intronなどと
# 重なる可能性があるため、promoter overlapとexon overlapは同義ではない。
#
# 【解釈上の注意】
# promoter近傍にERV fragmentが存在することは、
# そのERVが実際に遺伝子発現を制御していることを証明しない。
# 今回はreference genome上の静的な位置関係を解析する。
#
# 【主解析】
# 1. 各moduleのpromoter ERV overlap率
# 2. ERV全体のFisher's exact test
# 3. ERV family別のFisher's exact test
# 4. nearest ERV distanceのWilcoxon rank-sum test
# 5. 棒グラフ、forest plot、heatmap、boxplot
#
# 入力：
# clean/geneModuleTable_with_ERV_overlap_mad20k_power6.rds
# clean/promoter_ERV_pairs_mad20k_power6.rds
#
# 出力：
# results/target_module_promoter_ERV_fisher.csv
# results/target_module_promoter_ERV_family_fisher.csv
# results/target_module_nearest_ERV_distance_wilcoxon.csv
# figures/Step10_promoter_ERV_percent.pdf
# figures/Step10_promoter_ERV_forest.pdf
# figures/Step10_ERV_family_heatmap.pdf
# figures/Step10_nearest_ERV_distance_boxplot.pdf
# ============================================================


## Step 1
# 出力フォルダを作る

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


## Step 2
# packageを準備する

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)


## Step 3
# Step 09の結果を読む

geneModuleTable <- readRDS(
  "clean/geneModuleTable_with_ERV_overlap_mad20k_power6.rds"
)

promoter_ERV_pairs <- readRDS(
  "clean/promoter_ERV_pairs_mad20k_power6.rds"
)

cat("WGCNA gene数:", nrow(geneModuleTable), "\n")
cat("promoter-ERV pair数:", nrow(promoter_ERV_pairs), "\n")


## Step 4
# 必要な列を確認する

required_gene_columns <- c(
  "gene",
  "symbol",
  "module",
  "ERV_position_eligible",
  "promoter_ERV_any",
  "nearest_ERV_distance"
)

required_pair_columns <- c(
  "gene",
  "erv_fragment_id",
  "repFamily"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  colnames(geneModuleTable)
)

missing_pair_columns <- setdiff(
  required_pair_columns,
  colnames(promoter_ERV_pairs)
)

if (length(missing_gene_columns) > 0) {
  stop(
    paste(
      "geneModuleTableに不足している列:",
      paste(missing_gene_columns, collapse = ", ")
    )
  )
}

if (length(missing_pair_columns) > 0) {
  stop(
    paste(
      "promoter_ERV_pairsに不足している列:",
      paste(missing_pair_columns, collapse = ", ")
    )
  )
}


## Step 5
# AGE関連候補moduleを指定する

target_modules <- c(
  "greenyellow",
  "red",
  "darkgreen",
  "turquoise",
  "midnightblue",
  "blue",
  "pink"
)

missing_modules <- setdiff(
  target_modules,
  unique(geneModuleTable$module)
)

if (length(missing_modules) > 0) {
  stop(
    paste(
      "geneModuleTableに存在しないmodule:",
      paste(missing_modules, collapse = ", ")
    )
  )
}

print(
  table(
    geneModuleTable$module[
      geneModuleTable$module %in% target_modules
    ]
  )
)


## Step 6
# ERV位置解析が可能な遺伝子だけをbackgroundにする

eligible_df <- geneModuleTable[
  geneModuleTable$ERV_position_eligible == TRUE &
    !is.na(geneModuleTable$promoter_ERV_any),
]

# gene IDの重複を確認
if (sum(duplicated(eligible_df$gene)) > 0) {
  stop("eligible_df内でgene IDが重複しています。")
}

background_gene_n <- nrow(eligible_df)

background_promoter_ERV_n <- sum(
  eligible_df$promoter_ERV_any
)

background_promoter_ERV_percent <-
  background_promoter_ERV_n /
  background_gene_n * 100

background_median_distance <- median(
  eligible_df$nearest_ERV_distance,
  na.rm = TRUE
)

cat("解析対象gene数:", background_gene_n, "\n")
cat("promoter ERVありgene数:", background_promoter_ERV_n, "\n")
cat(
  "背景promoter ERV率:",
  round(background_promoter_ERV_percent, 2),
  "%\n"
)
cat(
  "背景nearest ERV distance中央値:",
  background_median_distance,
  "bp\n"
)


## Step 7
# Fisher's exact testを行う関数を定義する
#
# positive_vector：
# 各geneについて、対象ERVがpromoterにある場合TRUE

run_fisher_test <- function(module_name, positive_vector, data_table) {
  
  in_module <- data_table$module == module_name
  
  # 対象module内：ERVあり
  a <- sum(in_module & positive_vector)
  
  # 対象module内：ERVなし
  b <- sum(in_module & !positive_vector)
  
  # 対象module以外：ERVあり
  c <- sum(!in_module & positive_vector)
  
  # 対象module以外：ERVなし
  d <- sum(!in_module & !positive_vector)
  
  contingency_table <- matrix(
    c(a, b, c, d),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(contingency_table) <- c(
    module_name,
    "other_genes"
  )
  
  colnames(contingency_table) <- c(
    "ERV_yes",
    "ERV_no"
  )
  
  fisher_result <- fisher.test(
    contingency_table,
    alternative = "two.sided"
  )
  
  data.frame(
    module = module_name,
    module_gene_n = a + b,
    positive_gene_n = a,
    positive_percent = a / (a + b) * 100,
    other_gene_n = c + d,
    other_positive_gene_n = c,
    other_positive_percent = c / (c + d) * 100,
    odds_ratio = unname(fisher_result$estimate),
    confidence_interval_low = fisher_result$conf.int[1],
    confidence_interval_high = fisher_result$conf.int[2],
    p_value = fisher_result$p.value,
    stringsAsFactors = FALSE
  )
}


## Step 8
# ERV全体についてmodule別Fisher検定を行う

overall_results_list <- lapply(
  target_modules,
  function(mod) {
    run_fisher_test(
      module_name = mod,
      positive_vector = eligible_df$promoter_ERV_any,
      data_table = eligible_df
    )
  }
)

overall_ERV_results <- do.call(
  rbind,
  overall_results_list
)

# 列名を分かりやすく変更
colnames(overall_ERV_results)[
  colnames(overall_ERV_results) == "positive_gene_n"
] <- "promoter_ERV_gene_n"

colnames(overall_ERV_results)[
  colnames(overall_ERV_results) == "positive_percent"
] <- "promoter_ERV_percent"

colnames(overall_ERV_results)[
  colnames(overall_ERV_results) == "other_positive_gene_n"
] <- "other_promoter_ERV_gene_n"

colnames(overall_ERV_results)[
  colnames(overall_ERV_results) == "other_positive_percent"
] <- "other_promoter_ERV_percent"


## Step 9
# ERV全体のp値をFDR補正する

overall_ERV_results$FDR <- p.adjust(
  overall_ERV_results$p_value,
  method = "BH"
)

overall_ERV_results$direction <- ifelse(
  overall_ERV_results$odds_ratio > 1,
  "higher",
  ifelse(
    overall_ERV_results$odds_ratio < 1,
    "lower",
    "equal"
  )
)

overall_ERV_results$result <- ifelse(
  overall_ERV_results$FDR < 0.05 &
    overall_ERV_results$odds_ratio > 1,
  "enriched",
  ifelse(
    overall_ERV_results$FDR < 0.05 &
      overall_ERV_results$odds_ratio < 1,
    "depleted",
    "not_significant"
  )
)

overall_ERV_results <- overall_ERV_results[
  order(overall_ERV_results$FDR),
]

print(overall_ERV_results)

write.csv(
  overall_ERV_results,
  "results/target_module_promoter_ERV_fisher.csv",
  row.names = FALSE
)


## Step 10
# module別promoter ERV率の棒グラフを作る

module_colors <- c(
  "greenyellow" = "greenyellow",
  "red" = "red",
  "darkgreen" = "darkgreen",
  "turquoise" = "turquoise",
  "midnightblue" = "midnightblue",
  "blue" = "blue",
  "pink" = "pink"
)

percent_plot_df <- overall_ERV_results

percent_plot_df$module_plot <- factor(
  percent_plot_df$module,
  levels = percent_plot_df$module[
    order(percent_plot_df$promoter_ERV_percent)
  ]
)

promoter_percent_plot <- ggplot(
  percent_plot_df,
  aes(
    x = module_plot,
    y = promoter_ERV_percent,
    fill = module
  )
) +
  geom_col(width = 0.7) +
  geom_hline(
    yintercept = background_promoter_ERV_percent,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = module_colors,
    guide = "none"
  ) +
  coord_flip() +
  labs(
    x = "WGCNA module",
    y = "Genes with promoter ERV overlap (%)",
    title = "Promoter ERV overlap in age-associated candidate modules",
    subtitle = paste0(
      "Dashed line: all eligible WGCNA genes (",
      round(background_promoter_ERV_percent, 2),
      "%)"
    )
  ) +
  theme_classic(base_size = 12)

print(promoter_percent_plot)

ggsave(
  filename = "figures/Step10_promoter_ERV_percent.pdf",
  plot = promoter_percent_plot,
  width = 8,
  height = 5
)


## Step 11
# ERV全体のodds ratioをforest plotにする

forest_df <- overall_ERV_results[
  is.finite(overall_ERV_results$odds_ratio) &
    overall_ERV_results$odds_ratio > 0 &
    is.finite(overall_ERV_results$confidence_interval_low) &
    overall_ERV_results$confidence_interval_low > 0 &
    is.finite(overall_ERV_results$confidence_interval_high) &
    overall_ERV_results$confidence_interval_high > 0,
]

forest_df$module_forest <- factor(
  forest_df$module,
  levels = forest_df$module[
    order(forest_df$odds_ratio)
  ]
)

forest_plot <- ggplot(
  forest_df,
  aes(
    x = module_forest,
    y = odds_ratio,
    colour = module
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = confidence_interval_low,
      ymax = confidence_interval_high
    ),
    width = 0.2
  ) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = module_colors,
    guide = "none"
  ) +
  scale_y_log10() +
  coord_flip() +
  labs(
    x = "WGCNA module",
    y = "Odds ratio (log scale)",
    title = "Promoter ERV enrichment or depletion",
    subtitle = "Points: odds ratio; bars: 95% confidence interval"
  ) +
  theme_classic(base_size = 12)

print(forest_plot)

ggsave(
  filename = "figures/Step10_promoter_ERV_forest.pdf",
  plot = forest_plot,
  width = 8,
  height = 5
)


## Step 12
# ERV family別のpromoter-positive gene setを作る
#
# 同じgeneのpromoterに同じfamilyのERV fragmentが複数あっても、
# gene-levelでは1回だけ数える。

erv_families <- c(
  "ERV1",
  "ERVK",
  "ERVL",
  "ERVL-MaLR"
)

family_gene_sets <- setNames(
  lapply(
    erv_families,
    function(fam) {
      unique(
        promoter_ERV_pairs$gene[
          promoter_ERV_pairs$repFamily == fam
        ]
      )
    }
  ),
  erv_families
)

for (fam in erv_families) {
  cat(
    fam,
    ":",
    sum(eligible_df$gene %in% family_gene_sets[[fam]]),
    "genes\n"
  )
}


## Step 13
# module × ERV familyのFisher検定を行う

family_results_list <- list()
result_counter <- 1

for (mod in target_modules) {
  
  for (fam in erv_families) {
    
    family_positive <- eligible_df$gene %in%
      family_gene_sets[[fam]]
    
    family_result <- run_fisher_test(
      module_name = mod,
      positive_vector = family_positive,
      data_table = eligible_df
    )
    
    family_result$ERV_family <- fam
    
    family_results_list[[result_counter]] <- family_result
    
    result_counter <- result_counter + 1
  }
}

family_ERV_results <- do.call(
  rbind,
  family_results_list
)

# family解析用に列名を変更
colnames(family_ERV_results)[
  colnames(family_ERV_results) == "positive_gene_n"
] <- "family_positive_gene_n"

colnames(family_ERV_results)[
  colnames(family_ERV_results) == "positive_percent"
] <- "family_positive_percent"

colnames(family_ERV_results)[
  colnames(family_ERV_results) == "other_positive_gene_n"
] <- "other_family_positive_gene_n"

colnames(family_ERV_results)[
  colnames(family_ERV_results) == "other_positive_percent"
] <- "other_family_positive_percent"

# 列順を整理
family_ERV_results <- family_ERV_results[
  , c(
    "module",
    "ERV_family",
    "module_gene_n",
    "family_positive_gene_n",
    "family_positive_percent",
    "other_gene_n",
    "other_family_positive_gene_n",
    "other_family_positive_percent",
    "odds_ratio",
    "confidence_interval_low",
    "confidence_interval_high",
    "p_value"
  )
]


## Step 14
# family解析のp値をFDR補正する
# 7 modules × 4 families = 28検定

family_ERV_results$FDR <- p.adjust(
  family_ERV_results$p_value,
  method = "BH"
)

family_ERV_results$direction <- ifelse(
  family_ERV_results$odds_ratio > 1,
  "higher",
  ifelse(
    family_ERV_results$odds_ratio < 1,
    "lower",
    "equal"
  )
)

family_ERV_results$result <- ifelse(
  family_ERV_results$FDR < 0.05 &
    family_ERV_results$odds_ratio > 1,
  "enriched",
  ifelse(
    family_ERV_results$FDR < 0.05 &
      family_ERV_results$odds_ratio < 1,
    "depleted",
    "not_significant"
  )
)

family_ERV_results <- family_ERV_results[
  order(family_ERV_results$FDR),
]

print(family_ERV_results)

write.csv(
  family_ERV_results,
  "results/target_module_promoter_ERV_family_fisher.csv",
  row.names = FALSE
)


## Step 15
# module × ERV familyのheatmapを作る

family_heatmap_df <- family_ERV_results

family_heatmap_df$log2_odds_ratio <- log2(
  family_heatmap_df$odds_ratio
)

# 表示上のみlog2 ORを-3から+3に制限する
family_heatmap_df$log2_OR_plot <- pmax(
  pmin(family_heatmap_df$log2_odds_ratio, 3),
  -3
)

family_heatmap_df$significance_label <- ifelse(
  family_heatmap_df$FDR < 0.05,
  "*",
  ""
)

family_heatmap_df$module_heatmap <- factor(
  family_heatmap_df$module,
  levels = rev(target_modules)
)

family_heatmap_df$ERV_family <- factor(
  family_heatmap_df$ERV_family,
  levels = erv_families
)

family_heatmap <- ggplot(
  family_heatmap_df,
  aes(
    x = ERV_family,
    y = module_heatmap,
    fill = log2_OR_plot
  )
) +
  geom_tile(colour = "white") +
  geom_text(
    aes(label = significance_label),
    size = 5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    name = "log2 OR"
  ) +
  labs(
    x = "ERV family",
    y = "WGCNA module",
    title = "ERV-family-specific promoter overlap",
    subtitle = "* FDR < 0.05; colour scale limited to log2 OR ±3"
  ) +
  theme_classic(base_size = 12)

print(family_heatmap)

ggsave(
  filename = "figures/Step10_ERV_family_heatmap.pdf",
  plot = family_heatmap,
  width = 7,
  height = 5
)


## Step 16
# nearest ERV distanceをmoduleとその他のgeneで比較する

distance_results_list <- lapply(
  target_modules,
  function(mod) {
    
    module_distance <- eligible_df$nearest_ERV_distance[
      eligible_df$module == mod
    ]
    
    other_distance <- eligible_df$nearest_ERV_distance[
      eligible_df$module != mod
    ]
    
    module_distance <- module_distance[
      !is.na(module_distance)
    ]
    
    other_distance <- other_distance[
      !is.na(other_distance)
    ]
    
    wilcox_result <- wilcox.test(
      module_distance,
      other_distance,
      alternative = "two.sided",
      exact = FALSE
    )
    
    data.frame(
      module = mod,
      module_gene_n = length(module_distance),
      other_gene_n = length(other_distance),
      median_distance_module = median(module_distance),
      median_distance_other = median(other_distance),
      median_difference =
        median(module_distance) - median(other_distance),
      p_value = wilcox_result$p.value,
      stringsAsFactors = FALSE
    )
  }
)

distance_results <- do.call(
  rbind,
  distance_results_list
)


## Step 17
# distance解析のp値をFDR補正する

distance_results$FDR <- p.adjust(
  distance_results$p_value,
  method = "BH"
)

distance_results$direction <- ifelse(
  distance_results$median_difference < 0,
  "closer",
  ifelse(
    distance_results$median_difference > 0,
    "farther",
    "equal"
  )
)

distance_results$result <- ifelse(
  distance_results$FDR < 0.05,
  distance_results$direction,
  "not_significant"
)

distance_results <- distance_results[
  order(distance_results$FDR),
]

print(distance_results)

write.csv(
  distance_results,
  "results/target_module_nearest_ERV_distance_wilcoxon.csv",
  row.names = FALSE
)


## Step 18
# 7 moduleのnearest ERV distanceをboxplotにする

target_distance_df <- eligible_df[
  eligible_df$module %in% target_modules &
    !is.na(eligible_df$nearest_ERV_distance),
]

target_distance_df$module <- factor(
  target_distance_df$module,
  levels = target_modules
)

distance_boxplot <- ggplot(
  target_distance_df,
  aes(
    x = module,
    y = log10(nearest_ERV_distance + 1),
    fill = module
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(
    yintercept = log10(background_median_distance + 1),
    linetype = "dashed"
  ) +
  scale_fill_manual(
    values = module_colors,
    guide = "none"
  ) +
  coord_flip() +
  labs(
    x = "WGCNA module",
    y = "log10(nearest ERV distance + 1 bp)",
    title = "Distance from annotated TSS to nearest ERV fragment",
    subtitle = paste0(
      "Dashed line: median of all eligible genes (",
      background_median_distance,
      " bp)"
    )
  ) +
  theme_classic(base_size = 12)

print(distance_boxplot)

ggsave(
  filename = "figures/Step10_nearest_ERV_distance_boxplot.pdf",
  plot = distance_boxplot,
  width = 8,
  height = 5
)


## Step 19
# 最終確認

cat("\n==============================\n")
cat("背景gene数:", background_gene_n, "\n")
cat(
  "背景promoter ERV率:",
  round(background_promoter_ERV_percent, 2),
  "%\n"
)

cat("\nERV全体のmodule検定結果:\n")
print(
  overall_ERV_results[
    , c(
      "module",
      "module_gene_n",
      "promoter_ERV_percent",
      "odds_ratio",
      "p_value",
      "FDR",
      "result"
    )
  ]
)

cat("\nFamily別でFDR < 0.05の結果:\n")
print(
  family_ERV_results[
    family_ERV_results$FDR < 0.05,
    c(
      "module",
      "ERV_family",
      "family_positive_percent",
      "odds_ratio",
      "p_value",
      "FDR",
      "result"
    )
  ]
)

cat("\nNearest ERV distance検定結果:\n")
print(
  distance_results[
    , c(
      "module",
      "median_distance_module",
      "median_distance_other",
      "median_difference",
      "p_value",
      "FDR",
      "result"
    )
  ]
)

cat("==============================\n")


# Step 10の統計結果を読み込む

overall_ERV_results <- read.csv(
  "results/target_module_promoter_ERV_fisher.csv"
)

family_ERV_results <- read.csv(
  "results/target_module_promoter_ERV_family_fisher.csv"
)

distance_results <- read.csv(
  "results/target_module_nearest_ERV_distance_wilcoxon.csv"
)


#確認
# ERV全体のpromoter overlap
overall_ERV_results[
  order(overall_ERV_results$FDR),
  c(
    "module",
    "module_gene_n",
    "promoter_ERV_percent",
    "odds_ratio",
    "p_value",
    "FDR",
    "result"
  )
]

# ERV family別
family_ERV_results[
  order(family_ERV_results$FDR),
  c(
    "module",
    "ERV_family",
    "family_positive_gene_n",
    "family_positive_percent",
    "odds_ratio",
    "p_value",
    "FDR",
    "result"
  )
]

# nearest ERV distance
distance_results[
  order(distance_results$FDR),
  c(
    "module",
    "median_distance_module",
    "median_distance_other",
    "median_difference",
    "p_value",
    "FDR",
    "result"
  )
]
