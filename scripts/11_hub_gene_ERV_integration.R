# ============================================================
# 11_hub_gene_ERV_integration.R
#
# 目的：
# 7つのAGE関連候補moduleから抽出したtop 30 hub genes
# （7 modules × 30 genes = 210 genes）に、
# GenAge情報とERV位置情報を統合する。
#
# 【このStepの位置づけ】
#
# Step 10ではmodule全体を対象に統計検定を行った。
# その結果、pink moduleではpromoter近傍ERV overlapが
# 有意に少なく、複数moduleでnearest ERVまでの距離が
# 有意に長いことが分かった。
#
# Step 11では、module全体の統計検定を繰り返すのではなく、
# 210 hub genesを個別にannotateし、候補gene–ERV pairを整理する。
#
# 【MMとは】
#
# MM = Module Membership
# gene発現と所属module eigengeneの相関であり、
# module内での中心性を表す。
# |MM|が大きいほど、そのmoduleを代表するgeneと考えられる。
#
# 【GS_AGEとは】
#
# GS_AGE = Gene Significance for AGE
# gene発現と年齢カテゴリーの相関である。
#
# GS_AGE > 0：年齢とともに発現が高くなる傾向
# GS_AGE < 0：年齢とともに発現が低くなる傾向
#
# 【promoter ERVとは】
#
# annotated TSSの上流2,000 bpから下流500 bpまでに
# RepeatMasker ERV fragmentが重なる場合をpromoter ERVありとした。
# これは位置的な関連であり、ERVによる発現制御を証明しない。
#
# 【候補順位の考え方】
#
# promoter ERV overlap、GenAge登録、|GS_AGE|、|MM|を
# 個別の情報として保持し、透明な規則で並べる。
# 統合スコアを恣意的に作ることはしない。
#
# 入力：
# clean/top30_hub_candidates_summary_mad20k_power6.rds
# clean/geneModuleTable_with_ERV_overlap_mad20k_power6.rds
# clean/promoter_ERV_pairs_mad20k_power6.rds
#
# 出力：
# clean/hub_gene_ERV_annotation_all.rds
# results/hub_gene_ERV_annotation_all.csv
# results/hub_promoter_ERV_pairs.csv
# results/hub_ERV_summary_by_module.csv
# results/hub_gene_ERV_priority.csv
# results/hub_gene_ERV_priority_top20.csv
# results/GenAge_hub_ERV_annotation.csv
# results/pink_hub_ERV_annotation.csv
# figures/Step11_hub_promoter_ERV_percent.pdf
# figures/Step11_hub_MM_GS_ERV_scatter.pdf
# ============================================================


## Step 1
# 出力フォルダを作る

dir.create(
  "clean",
  showWarnings = FALSE
)

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)


## Step 2
# 必要なpackageを準備する

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)


## Step 3
# 入力データを読む

hub_summary_all <- readRDS(
  "clean/top30_hub_candidates_summary_mad20k_power6.rds"
)

geneModuleTable <- readRDS(
  "clean/geneModuleTable_with_ERV_overlap_mad20k_power6.rds"
)

promoter_ERV_pairs <- readRDS(
  "clean/promoter_ERV_pairs_mad20k_power6.rds"
)


## Step 4
# 入力データの基本確認

cat(
  "hub候補gene数:",
  nrow(hub_summary_all),
  "\n"
)

cat(
  "WGCNA全gene数:",
  nrow(geneModuleTable),
  "\n"
)

cat(
  "promoter-ERV pair数:",
  nrow(promoter_ERV_pairs),
  "\n"
)

# hub表に必須の列
required_hub_columns <- c(
  "gene",
  "module",
  "MM",
  "GS_AGE"
)

missing_hub_columns <- setdiff(
  required_hub_columns,
  colnames(hub_summary_all)
)

if (length(missing_hub_columns) > 0) {
  stop(
    paste0(
      "hub_summary_allに必要な列がありません: ",
      paste(missing_hub_columns, collapse = ", "),
      "\n実際の列: ",
      paste(colnames(hub_summary_all), collapse = ", ")
    )
  )
}

# gene–ERV統合表に必須の列
required_gene_columns <- c(
  "gene",
  "symbol",
  "module",
  "is_GenAge",
  "chromosome",
  "strand",
  "TSS",
  "ERV_position_eligible",
  "promoter_ERV_n",
  "promoter_ERV_any",
  "promoter_ERV_families",
  "nearest_ERV_distance",
  "nearest_ERV_fragment_id",
  "nearest_ERV_repName",
  "nearest_ERV_family",
  "nearest_ERV_part",
  "nearest_ERV_direction"
)

missing_gene_columns <- setdiff(
  required_gene_columns,
  colnames(geneModuleTable)
)

if (length(missing_gene_columns) > 0) {
  stop(
    paste0(
      "geneModuleTableに必要な列がありません: ",
      paste(missing_gene_columns, collapse = ", ")
    )
  )
}

# pair表に必須の列
required_pair_columns <- c(
  "gene",
  "erv_fragment_id",
  "repeat_start",
  "repeat_end",
  "erv_strand",
  "repName",
  "repFamily",
  "repeat_part",
  "percent_divergence"
)

missing_pair_columns <- setdiff(
  required_pair_columns,
  colnames(promoter_ERV_pairs)
)

if (length(missing_pair_columns) > 0) {
  stop(
    paste0(
      "promoter_ERV_pairsに必要な列がありません: ",
      paste(missing_pair_columns, collapse = ", ")
    )
  )
}


## Step 5
# target moduleを指定する

target_modules <- c(
  "greenyellow",
  "red",
  "darkgreen",
  "turquoise",
  "midnightblue",
  "blue",
  "pink"
)

# 7 modules × 30 genesになっているか確認
hub_module_counts <- table(
  hub_summary_all$module
)

hub_module_counts

unexpected_modules <- setdiff(
  unique(hub_summary_all$module),
  target_modules
)

if (length(unexpected_modules) > 0) {
  warning(
    paste0(
      "想定外のmoduleがあります: ",
      paste(unexpected_modules, collapse = ", ")
    )
  )
}


## Step 6
# rank_in_module列がなければ作る

# 元ファイルにrank_in_moduleがある場合はそのまま使用する。
# ない場合は、各module内で|MM|、|GS_AGE|の順に並べて順位を作る。

if (!"rank_in_module" %in% colnames(hub_summary_all)) {
  
  hub_summary_all$rank_in_module <- NA_integer_
  
  for (mod in unique(hub_summary_all$module)) {
    
    module_rows <- which(
      hub_summary_all$module == mod
    )
    
    module_order <- order(
      abs(hub_summary_all$MM[module_rows]),
      abs(hub_summary_all$GS_AGE[module_rows]),
      decreasing = TRUE
    )
    
    hub_summary_all$rank_in_module[
      module_rows[module_order]
    ] <- seq_along(module_rows)
  }
}


## Step 7
# hub geneをWGCNA・ERV統合表と照合する

hub_match_index <- match(
  hub_summary_all$gene,
  geneModuleTable$gene
)

cat(
  "ERV統合表に一致したhub gene数:",
  sum(!is.na(hub_match_index)),
  "\n"
)

cat(
  "一致しなかったhub gene数:",
  sum(is.na(hub_match_index)),
  "\n"
)

if (any(is.na(hub_match_index))) {
  
  unmatched_hubs <- hub_summary_all[
    is.na(hub_match_index),
    c(
      "gene",
      "module",
      "MM",
      "GS_AGE"
    )
  ]
  
  write.csv(
    unmatched_hubs,
    "results/unmatched_hub_genes_ERV_integration.csv",
    row.names = FALSE
  )
  
  stop(
    "一部hub geneがgeneModuleTableに一致しません。unmatched_hub_genes_ERV_integration.csvを確認してください。"
  )
}

# moduleが一致するか確認
module_match <- hub_summary_all$module ==
  geneModuleTable$module[hub_match_index]

if (!all(module_match)) {
  stop("hub表とgeneModuleTableでmodule assignmentが一致しません。")
}


## Step 8
# hub表にGenAge・座標・ERV情報を追加する

hub_annotated <- hub_summary_all

columns_to_add <- c(
  "symbol",
  "is_GenAge",
  "gene_type",
  "chromosome",
  "gene_start",
  "gene_end",
  "strand",
  "TSS",
  "gene_length",
  "ERV_position_eligible",
  "promoter_ERV_n",
  "promoter_ERV_any",
  "promoter_ERV_families",
  "nearest_ERV_distance",
  "nearest_ERV_signed_distance",
  "nearest_ERV_fragment_id",
  "nearest_ERV_repName",
  "nearest_ERV_family",
  "nearest_ERV_part",
  "nearest_ERV_start",
  "nearest_ERV_end",
  "nearest_ERV_strand",
  "nearest_ERV_direction"
)

# geneModuleTableに実際に存在する列だけを追加する
columns_to_add <- intersect(
  columns_to_add,
  colnames(geneModuleTable)
)

for (column_name in columns_to_add) {
  hub_annotated[[column_name]] <-
    geneModuleTable[[column_name]][hub_match_index]
}

# is_GenAgeのNAはFALSEとして扱う
hub_annotated$is_GenAge[
  is.na(hub_annotated$is_GenAge)
] <- FALSE

# 年齢との相関方向を明示する
hub_annotated$AGE_expression_direction <- ifelse(
  hub_annotated$GS_AGE > 0,
  "increases_with_age",
  ifelse(
    hub_annotated$GS_AGE < 0,
    "decreases_with_age",
    "no_direction"
  )
)

# 並べ替えに使う絶対値
hub_annotated$abs_MM <- abs(
  hub_annotated$MM
)

hub_annotated$abs_GS_AGE <- abs(
  hub_annotated$GS_AGE
)


## Step 9
# promoter ERVとGenAgeの組合せを分類する

hub_annotated$candidate_group <- ifelse(
  
  hub_annotated$is_GenAge &
    hub_annotated$promoter_ERV_any,
  
  "GenAge_and_promoter_ERV",
  
  ifelse(
    
    hub_annotated$promoter_ERV_any,
    
    "promoter_ERV",
    
    ifelse(
      hub_annotated$is_GenAge,
      "GenAge_without_promoter_ERV",
      "other_hub"
    )
  )
)

# candidate_groupの優先順
candidate_group_order <- c(
  "GenAge_and_promoter_ERV",
  "promoter_ERV",
  "GenAge_without_promoter_ERV",
  "other_hub"
)

hub_annotated$candidate_group <- factor(
  hub_annotated$candidate_group,
  levels = candidate_group_order
)


## Step 10
# hub annotation表を並べ替える

# 統合スコアは作らず、以下の順に透明に並べる。
# 1. candidate group
# 2. |GS_AGE|
# 3. |MM|
# 4. promoter ERV fragment数

hub_priority_order <- order(
  hub_annotated$candidate_group,
  -hub_annotated$abs_GS_AGE,
  -hub_annotated$abs_MM,
  -hub_annotated$promoter_ERV_n,
  hub_annotated$module,
  hub_annotated$rank_in_module
)

hub_priority <- hub_annotated[
  hub_priority_order,
]

# 優先順位番号を追加
hub_priority$priority_rank <- seq_len(
  nrow(hub_priority)
)

# 行名を通常の連番に戻す
rownames(hub_priority) <- NULL


## Step 11
# hub promoter–ERV pair表を作る

hub_gene_ids <- unique(
  hub_annotated$gene
)

hub_promoter_ERV_pairs <- promoter_ERV_pairs[
  promoter_ERV_pairs$gene %in% hub_gene_ids,
]

# 各pairへhub情報を追加する
pair_hub_index <- match(
  hub_promoter_ERV_pairs$gene,
  hub_annotated$gene
)

hub_promoter_ERV_pairs$hub_module <-
  hub_annotated$module[pair_hub_index]

hub_promoter_ERV_pairs$rank_in_module <-
  hub_annotated$rank_in_module[pair_hub_index]

hub_promoter_ERV_pairs$MM <-
  hub_annotated$MM[pair_hub_index]

hub_promoter_ERV_pairs$GS_AGE <-
  hub_annotated$GS_AGE[pair_hub_index]

hub_promoter_ERV_pairs$is_GenAge <-
  hub_annotated$is_GenAge[pair_hub_index]

# 見やすい順に並べる
hub_promoter_ERV_pairs <- hub_promoter_ERV_pairs[
  order(
    hub_promoter_ERV_pairs$hub_module,
    hub_promoter_ERV_pairs$rank_in_module,
    hub_promoter_ERV_pairs$gene,
    hub_promoter_ERV_pairs$repeat_start
  ),
]

rownames(hub_promoter_ERV_pairs) <- NULL

# gene–ERV fragment pairの重複確認
hub_pair_duplicates <- sum(
  duplicated(
    hub_promoter_ERV_pairs[
      , c(
        "gene",
        "erv_fragment_id"
      )
    ]
  )
)


## Step 12
# module別のhub ERV summaryを作る

background_eligible <- geneModuleTable[
  geneModuleTable$ERV_position_eligible == TRUE &
    !is.na(geneModuleTable$promoter_ERV_any),
]

background_promoter_ERV_percent <-
  mean(background_eligible$promoter_ERV_any) * 100

hub_module_summary_list <- vector(
  "list",
  length(target_modules)
)

for (i in seq_along(target_modules)) {
  
  mod <- target_modules[i]
  
  module_hubs <- hub_annotated[
    hub_annotated$module == mod,
  ]
  
  eligible_hubs <- module_hubs[
    module_hubs$ERV_position_eligible == TRUE &
      !is.na(module_hubs$promoter_ERV_any),
  ]
  
  promoter_positive_n <- sum(
    eligible_hubs$promoter_ERV_any
  )
  
  hub_module_summary_list[[i]] <- data.frame(
    
    module = mod,
    
    hub_gene_n = nrow(module_hubs),
    
    ERV_position_eligible_hub_n = nrow(eligible_hubs),
    
    promoter_ERV_hub_n = promoter_positive_n,
    
    promoter_ERV_hub_percent = ifelse(
      nrow(eligible_hubs) > 0,
      promoter_positive_n / nrow(eligible_hubs) * 100,
      NA_real_
    ),
    
    promoter_ERV_fragment_n = sum(
      eligible_hubs$promoter_ERV_n,
      na.rm = TRUE
    ),
    
    GenAge_hub_n = sum(
      module_hubs$is_GenAge,
      na.rm = TRUE
    ),
    
    GenAge_and_promoter_ERV_hub_n = sum(
      module_hubs$is_GenAge &
        module_hubs$promoter_ERV_any,
      na.rm = TRUE
    ),
    
    median_nearest_ERV_distance = median(
      eligible_hubs$nearest_ERV_distance,
      na.rm = TRUE
    ),
    
    background_promoter_ERV_percent =
      background_promoter_ERV_percent,
    
    stringsAsFactors = FALSE
  )
}

hub_module_summary <- do.call(
  rbind,
  hub_module_summary_list
)


## Step 13
# 特に確認するsubsetsを作る

# promoterにERVが重なるhub genes
promoter_positive_hubs <- hub_priority[
  hub_priority$promoter_ERV_any %in% TRUE,
]

# GenAgeに登録されたhub genes
GenAge_hubs <- hub_priority[
  hub_priority$is_GenAge %in% TRUE,
]

# Step 10で最も明確なdepletionを示したpink moduleのhub genes
pink_hubs <- hub_annotated[
  hub_annotated$module == "pink",
]

pink_hubs <- pink_hubs[
  order(
    pink_hubs$rank_in_module
  ),
]


## Step 14
# 結果を保存する

saveRDS(
  hub_annotated,
  "clean/hub_gene_ERV_annotation_all.rds"
)

write.csv(
  hub_annotated,
  "results/hub_gene_ERV_annotation_all.csv",
  row.names = FALSE
)

write.csv(
  hub_promoter_ERV_pairs,
  "results/hub_promoter_ERV_pairs.csv",
  row.names = FALSE
)

write.csv(
  hub_module_summary,
  "results/hub_ERV_summary_by_module.csv",
  row.names = FALSE
)

write.csv(
  hub_priority,
  "results/hub_gene_ERV_priority.csv",
  row.names = FALSE
)

write.csv(
  head(hub_priority, 20),
  "results/hub_gene_ERV_priority_top20.csv",
  row.names = FALSE
)

write.csv(
  promoter_positive_hubs,
  "results/promoter_positive_hub_genes.csv",
  row.names = FALSE
)

write.csv(
  GenAge_hubs,
  "results/GenAge_hub_ERV_annotation.csv",
  row.names = FALSE
)

write.csv(
  pink_hubs,
  "results/pink_hub_ERV_annotation.csv",
  row.names = FALSE
)


## Step 15
# Figure 1：module別hub promoter ERV率

module_colors <- c(
  "greenyellow" = "greenyellow",
  "red" = "red",
  "darkgreen" = "darkgreen",
  "turquoise" = "turquoise",
  "midnightblue" = "midnightblue",
  "blue" = "blue",
  "pink" = "pink"
)

hub_module_summary$module_plot <- factor(
  hub_module_summary$module,
  levels = rev(target_modules)
)

hub_percent_plot <- ggplot(
  hub_module_summary,
  aes(
    x = module_plot,
    y = promoter_ERV_hub_percent,
    fill = module
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
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
    y = "Hub genes with promoter ERV overlap (%)",
    title = "Promoter ERV overlap among top 30 hub genes",
    subtitle = paste0(
      "Descriptive analysis; dashed line: all eligible WGCNA genes (",
      round(background_promoter_ERV_percent, 2),
      "%)"
    )
  ) +
  
  theme_classic(
    base_size = 12
  )

# 画面表示は行わず、直接PDFへ保存する

ggsave(
  filename = "figures/Step11_hub_promoter_ERV_percent.pdf",
  plot = hub_percent_plot,
  device = "pdf",
  width = 8,
  height = 5
)


## Step 16
# Figure 2：MM・GS_AGEとpromoter ERVの関係

hub_annotated$promoter_status <- ifelse(
  hub_annotated$promoter_ERV_any %in% TRUE,
  "Promoter ERV",
  "No promoter ERV"
)

hub_annotated$module_plot <- factor(
  hub_annotated$module,
  levels = target_modules
)

hub_scatter_plot <- ggplot(
  hub_annotated,
  aes(
    x = abs_MM,
    y = abs_GS_AGE,
    colour = module,
    shape = promoter_status
  )
) +
  
  geom_point(
    size = 2.6,
    alpha = 0.85
  ) +
  
  geom_text(
    data = hub_annotated[
      hub_annotated$is_GenAge %in% TRUE,
    ],
    aes(
      label = symbol
    ),
    vjust = -0.7,
    size = 3,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  
  scale_colour_manual(
    values = module_colors,
    name = "WGCNA module"
  ) +
  
  scale_shape_manual(
    values = c(
      "No promoter ERV" = 1,
      "Promoter ERV" = 16
    ),
    name = "Promoter overlap"
  ) +
  
  labs(
    x = "Absolute module membership |MM|",
    y = "Absolute gene significance for age |GS_AGE|",
    title = "ERV annotation of age-associated hub gene candidates",
    subtitle = "Labels indicate GenAge hub genes"
  ) +
  
  theme_classic(
    base_size = 12
  )

# 画面表示は行わず、直接PDFへ保存する

ggsave(
  filename = "figures/Step11_hub_MM_GS_ERV_scatter.pdf",
  plot = hub_scatter_plot,
  device = "pdf",
  width = 9,
  height = 6
)


## Step 17
# 最終確認

cat(
  "\n==============================\n"
)

cat(
  "hub gene総数:",
  nrow(hub_annotated),
  "\n"
)

cat(
  "module数:",
  length(unique(hub_annotated$module)),
  "\n"
)

cat(
  "ERV位置解析対象hub数:",
  sum(hub_annotated$ERV_position_eligible),
  "\n"
)

cat(
  "promoter ERVありhub数:",
  sum(hub_annotated$promoter_ERV_any, na.rm = TRUE),
  "\n"
)

cat(
  "hub promoter-ERV pair数:",
  nrow(hub_promoter_ERV_pairs),
  "\n"
)

cat(
  "hub promoter-ERV pair重複数:",
  hub_pair_duplicates,
  "\n"
)

cat(
  "GenAge hub数:",
  sum(hub_annotated$is_GenAge),
  "\n"
)

cat(
  "GenAgeかつpromoter ERVありhub数:",
  sum(
    hub_annotated$is_GenAge &
      hub_annotated$promoter_ERV_any,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nmodule別summary:\n"
)

print(
  hub_module_summary
)

cat(
  "\nGenAge hub genes:\n"
)

print(
  GenAge_hubs[
    , c(
      "gene",
      "symbol",
      "module",
      "rank_in_module",
      "MM",
      "GS_AGE",
      "promoter_ERV_any",
      "promoter_ERV_n",
      "promoter_ERV_families",
      "nearest_ERV_repName",
      "nearest_ERV_family",
      "nearest_ERV_distance"
    )
  ]
)

cat(
  "\n優先順位上位20 hub genes:\n"
)

print(
  head(
    hub_priority[
      , c(
        "priority_rank",
        "gene",
        "symbol",
        "module",
        "rank_in_module",
        "MM",
        "GS_AGE",
        "is_GenAge",
        "promoter_ERV_any",
        "promoter_ERV_n",
        "promoter_ERV_families",
        "nearest_ERV_repName",
        "nearest_ERV_family",
        "nearest_ERV_distance",
        "candidate_group"
      )
    ],
    20
  )
)

cat(
  "==============================\n"
)