# 目的：
# Step 11で得られた33個のpromoter-ERV陽性hub genesと
# 63個のhub–ERV fragment pairsについて、
# ERV fragmentが遺伝子構造のどこに位置するかを分類する。
#
# 【今回分類する領域】
#
# ・TSS_overlap
#   ERV fragmentが転写開始点そのものを含む
#
# ・first_exon
#   ERV fragmentがTSS側の最初のexonに重なる
#
# ・other_exon
#   ERV fragmentがfirst exon以外のexonに重なる
#
# ・intron
#   ERV fragmentがexon間のintronに重なる
#
# ・upstream_promoter
#   ERV fragmentがTSS上流のpromoter領域にあり、
#   TSS・exon・intronには重ならない
#
# ・downstream_promoter
#   ERV fragmentがTSS下流側にあるが、
#   exon・intronには重ならない
#
# 【TSSとは】
#
# TSS = Transcription Start Site
# RNAへの転写が始まるゲノム上の位置である。
#
# + strand：TSS = gene_start
# - strand：TSS = gene_end
#
# 【first exonとは】
#
# + strandでは、TSS側にある最も座標の小さいexonをfirst exonとする。
# - strandでは、TSS側にある最も座標の大きいexonをfirst exonとする。
#
# 今回使用するGTEx V10用GENCODE v39 gene modelは、
# geneごとに統合されたtranscript modelである。
# したがって、ここでのfirst exonはこのannotation model上のfirst exonであり、
# 骨格筋で実際に使用されるisoformのfirst exonを直接証明するものではない。
#
# 【重要な注意】
#
# 1つのERV fragmentが境界をまたぐ場合、
# first exonとintronなど複数のbinary flagがTRUEになる可能性がある。
# main_contextでは、次の優先順で1つの代表分類を付ける。
#
# TSS_overlap > first_exon > other_exon > intron
# > upstream_promoter / downstream_promoter
#
# また、repeat_partのLTR_or_otherはsolo LTRを意味しない。
# 単にrepNameが「-int」で終わらないfragmentを表す。
#
# 【このStepで分からないこと】
#
# ・ERV fragmentが実際に転写されているか
# ・遺伝子発現を調節しているか
# ・enhancerとして機能するか
#
# enhancer overlapは次のStep 12Bで別annotationを用いて解析する。
#
# 入力：
# clean/hub_gene_ERV_annotation_all.rds
# clean/promoter_ERV_pairs_mad20k_power6.rds
# /rds/homes/t/txk567/gtex_wgcna/raw/gencode.v39.GRCh38.genes.gtf
#
# 出力：
# clean/hub_promoter_ERV_pairs_genic_context.rds
# results/hub_promoter_ERV_pairs_genic_context.csv
# results/hub_promoter_ERV_genic_context_summary.csv
# results/hub_promoter_ERV_context_by_module.csv
# results/hub_promoter_ERV_candidate_gene_summary.csv
# results/hub_promoter_ERV_context_priority.csv
# figures/Step12A_hub_ERV_genic_context.pdf
# figures/Step12A_hub_ERV_context_by_module.pdf
# ============================================================


## Step 1
# 出力フォルダを作る

dir.create("clean", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


## Step 2
# 必要なpackageを準備する

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(rtracklayer)
library(GenomicRanges)
library(ggplot2)


## Step 3
# 入力ファイルを読む

hub_annotated <- readRDS(
  "clean/hub_gene_ERV_annotation_all.rds"
)

promoter_ERV_pairs <- readRDS(
  "clean/promoter_ERV_pairs_mad20k_power6.rds"
)

gtf_file <- paste0(
  "/rds/homes/t/txk567/gtex_wgcna/raw/",
  "gencode.v39.GRCh38.genes.gtf"
)

if (!file.exists(gtf_file)) {
  stop(
    paste0(
      "GTFファイルが見つかりません：\n",
      gtf_file
    )
  )
}


## Step 4
# 入力データを確認する

required_hub_columns <- c(
  "gene",
  "symbol",
  "module",
  "rank_in_module",
  "MM",
  "GS_AGE",
  "is_GenAge",
  "gene_start",
  "gene_end",
  "strand",
  "TSS",
  "promoter_ERV_any"
)

required_pair_columns <- c(
  "gene",
  "symbol",
  "module",
  "chromosome",
  "gene_strand",
  "TSS",
  "promoter_start",
  "promoter_end",
  "erv_fragment_id",
  "repeat_start",
  "repeat_end",
  "erv_strand",
  "repName",
  "repFamily",
  "repeat_part",
  "percent_divergence"
)

missing_hub_columns <- setdiff(
  required_hub_columns,
  colnames(hub_annotated)
)

missing_pair_columns <- setdiff(
  required_pair_columns,
  colnames(promoter_ERV_pairs)
)

if (length(missing_hub_columns) > 0) {
  stop(
    paste0(
      "hub_annotatedに必要な列がありません：",
      paste(missing_hub_columns, collapse = ", ")
    )
  )
}

if (length(missing_pair_columns) > 0) {
  stop(
    paste0(
      "promoter_ERV_pairsに必要な列がありません：",
      paste(missing_pair_columns, collapse = ", ")
    )
  )
}

cat("hub gene総数:", nrow(hub_annotated), "\n")
cat(
  "promoter ERVありhub数:",
  sum(hub_annotated$promoter_ERV_any, na.rm = TRUE),
  "\n"
)


## Step 5
# hub geneに対応するpromoter–ERV pairsだけを取り出す

hub_gene_ids <- unique(hub_annotated$gene)

hub_pairs <- promoter_ERV_pairs[
  promoter_ERV_pairs$gene %in% hub_gene_ids,
]

rownames(hub_pairs) <- NULL

cat("hub promoter–ERV pair数:", nrow(hub_pairs), "\n")
cat(
  "pairに含まれるhub gene数:",
  length(unique(hub_pairs$gene)),
  "\n"
)

if (nrow(hub_pairs) == 0) {
  stop("hub promoter–ERV pairがありません。Step 11の入力を確認してください。")
}


## Step 6
# hub情報をpair表へ追加する

hub_index <- match(
  hub_pairs$gene,
  hub_annotated$gene
)

if (any(is.na(hub_index))) {
  stop("一部のhub pairがhub annotation表に一致しません。")
}

hub_pairs$hub_module <- hub_annotated$module[hub_index]
hub_pairs$rank_in_module <- hub_annotated$rank_in_module[hub_index]
hub_pairs$MM <- hub_annotated$MM[hub_index]
hub_pairs$GS_AGE <- hub_annotated$GS_AGE[hub_index]
hub_pairs$is_GenAge <- hub_annotated$is_GenAge[hub_index]
hub_pairs$gene_start <- hub_annotated$gene_start[hub_index]
hub_pairs$gene_end <- hub_annotated$gene_end[hub_index]


## Step 7
# 各ERV fragmentのTSSに対する方向と距離を計算する

hub_pairs$overlaps_TSS <-
  hub_pairs$repeat_start <= hub_pairs$TSS &
  hub_pairs$repeat_end >= hub_pairs$TSS

hub_pairs$relative_to_TSS <- ifelse(
  hub_pairs$overlaps_TSS,
  "overlaps_TSS",
  ifelse(
    hub_pairs$gene_strand == "+",
    ifelse(
      hub_pairs$repeat_end < hub_pairs$TSS,
      "upstream",
      "downstream"
    ),
    ifelse(
      hub_pairs$repeat_start > hub_pairs$TSS,
      "upstream",
      "downstream"
    )
  )
)

# TSSとERV fragmentの最も近い端との距離
hub_pairs$distance_to_TSS <- ifelse(
  hub_pairs$overlaps_TSS,
  0L,
  ifelse(
    hub_pairs$repeat_end < hub_pairs$TSS,
    hub_pairs$TSS - hub_pairs$repeat_end,
    hub_pairs$repeat_start - hub_pairs$TSS
  )
)

# upstreamを負、downstreamを正として表す
hub_pairs$signed_distance_to_TSS <- hub_pairs$distance_to_TSS

hub_pairs$signed_distance_to_TSS[
  hub_pairs$relative_to_TSS == "upstream"
] <- -hub_pairs$distance_to_TSS[
  hub_pairs$relative_to_TSS == "upstream"
]

hub_pairs$signed_distance_to_TSS[
  hub_pairs$relative_to_TSS == "overlaps_TSS"
] <- 0L

# promoterとERVが実際に重なる塩基数
hub_pairs$promoter_overlap_bp <- pmax(
  0L,
  pmin(hub_pairs$repeat_end, hub_pairs$promoter_end) -
    pmax(hub_pairs$repeat_start, hub_pairs$promoter_start) + 1L
)

# geneとERVの向きが同じか反対か
hub_pairs$strand_relation <- ifelse(
  hub_pairs$gene_strand == hub_pairs$erv_strand,
  "same_strand",
  "opposite_strand"
)


## Step 8
# GENCODE v39 GTFからexonを読み込む

gtf_v39 <- import(gtf_file)

exons_v39 <- gtf_v39[
  gtf_v39$type == "exon"
]

# hub genesに対応するexonだけを残す
hub_exons_v39 <- exons_v39[
  exons_v39$gene_id %in% hub_gene_ids
]

exon_table <- data.frame(
  gene = hub_exons_v39$gene_id,
  transcript_id = hub_exons_v39$transcript_id,
  chromosome = as.character(seqnames(hub_exons_v39)),
  exon_start = start(hub_exons_v39),
  exon_end = end(hub_exons_v39),
  gene_strand = as.character(strand(hub_exons_v39)),
  stringsAsFactors = FALSE
)

cat("hub genesに対応したexon数:", nrow(exon_table), "\n")
cat(
  "exonが見つかったhub gene数:",
  length(unique(exon_table$gene)),
  "\n"
)

missing_exon_genes <- setdiff(
  hub_gene_ids,
  unique(exon_table$gene)
)

if (length(missing_exon_genes) > 0) {
  write.csv(
    data.frame(gene = missing_exon_genes),
    "results/hub_genes_without_GENCODE_exons.csv",
    row.names = FALSE
  )
  
  warning(
    paste0(
      length(missing_exon_genes),
      " hub genesでexonが見つかりませんでした。"
    )
  )
}


## Step 9
# annotation model上のfirst exonを判定する

exon_table$exon_rank_from_TSS <- NA_integer_

for (gene_id in unique(exon_table$gene)) {
  
  exon_rows <- which(
    exon_table$gene == gene_id
  )
  
  current_strand <- exon_table$gene_strand[exon_rows[1]]
  
  if (current_strand == "+") {
    exon_order <- order(
      exon_table$exon_start[exon_rows],
      exon_table$exon_end[exon_rows]
    )
  } else {
    exon_order <- order(
      -exon_table$exon_end[exon_rows],
      -exon_table$exon_start[exon_rows]
    )
  }
  
  exon_table$exon_rank_from_TSS[
    exon_rows[exon_order]
  ] <- seq_along(exon_rows)
}

exon_table$is_first_exon <-
  exon_table$exon_rank_from_TSS == 1L

cat(
  "first exonとして判定されたexon数:",
  sum(exon_table$is_first_exon),
  "\n"
)


## Step 10
# hub–ERV pairsとexonの重なりを調べる

pair_gr <- GRanges(
  seqnames = hub_pairs$chromosome,
  ranges = IRanges(
    start = hub_pairs$repeat_start,
    end = hub_pairs$repeat_end
  ),
  strand = hub_pairs$erv_strand
)

exon_gr <- GRanges(
  seqnames = exon_table$chromosome,
  ranges = IRanges(
    start = exon_table$exon_start,
    end = exon_table$exon_end
  ),
  strand = exon_table$gene_strand
)

exon_hits <- findOverlaps(
  pair_gr,
  exon_gr,
  ignore.strand = TRUE
)

exon_query <- queryHits(exon_hits)
exon_subject <- subjectHits(exon_hits)

# 別遺伝子のexonとの偶然の重なりは除外する
same_gene_exon_hit <-
  hub_pairs$gene[exon_query] ==
  exon_table$gene[exon_subject]

exon_query <- exon_query[same_gene_exon_hit]
exon_subject <- exon_subject[same_gene_exon_hit]

hub_pairs$overlaps_first_exon <- FALSE
hub_pairs$overlaps_other_exon <- FALSE

first_exon_queries <- unique(
  exon_query[
    exon_table$is_first_exon[exon_subject]
  ]
)

other_exon_queries <- unique(
  exon_query[
    !exon_table$is_first_exon[exon_subject]
  ]
)

hub_pairs$overlaps_first_exon[first_exon_queries] <- TRUE
hub_pairs$overlaps_other_exon[other_exon_queries] <- TRUE
hub_pairs$overlaps_any_exon <-
  hub_pairs$overlaps_first_exon |
  hub_pairs$overlaps_other_exon


## Step 11
# exon間の領域からintron座標を作る

intron_list <- list()
intron_counter <- 1L

for (gene_id in unique(exon_table$gene)) {
  
  exon_rows <- which(
    exon_table$gene == gene_id
  )
  
  reduced_exons <- reduce(
    IRanges(
      start = exon_table$exon_start[exon_rows],
      end = exon_table$exon_end[exon_rows]
    )
  )
  
  if (length(reduced_exons) < 2L) {
    next
  }
  
  intron_start <- end(reduced_exons)[
    seq_len(length(reduced_exons) - 1L)
  ] + 1L
  
  intron_end <- start(reduced_exons)[
    seq.int(2L, length(reduced_exons))
  ] - 1L
  
  valid_introns <- intron_start <= intron_end
  
  if (!any(valid_introns)) {
    next
  }
  
  first_exon_row <- exon_rows[1]
  
  intron_list[[intron_counter]] <- data.frame(
    gene = gene_id,
    chromosome = exon_table$chromosome[first_exon_row],
    intron_start = intron_start[valid_introns],
    intron_end = intron_end[valid_introns],
    gene_strand = exon_table$gene_strand[first_exon_row],
    stringsAsFactors = FALSE
  )
  
  intron_counter <- intron_counter + 1L
}

if (length(intron_list) > 0) {
  intron_table <- do.call(rbind, intron_list)
} else {
  intron_table <- data.frame(
    gene = character(0),
    chromosome = character(0),
    intron_start = integer(0),
    intron_end = integer(0),
    gene_strand = character(0),
    stringsAsFactors = FALSE
  )
}

cat("作成したintron区間数:", nrow(intron_table), "\n")


## Step 12
# hub–ERV pairsとintronの重なりを調べる

hub_pairs$overlaps_intron <- FALSE

if (nrow(intron_table) > 0) {
  
  intron_gr <- GRanges(
    seqnames = intron_table$chromosome,
    ranges = IRanges(
      start = intron_table$intron_start,
      end = intron_table$intron_end
    ),
    strand = intron_table$gene_strand
  )
  
  intron_hits <- findOverlaps(
    pair_gr,
    intron_gr,
    ignore.strand = TRUE
  )
  
  intron_query <- queryHits(intron_hits)
  intron_subject <- subjectHits(intron_hits)
  
  same_gene_intron_hit <-
    hub_pairs$gene[intron_query] ==
    intron_table$gene[intron_subject]
  
  intron_query <- intron_query[same_gene_intron_hit]
  
  hub_pairs$overlaps_intron[
    unique(intron_query)
  ] <- TRUE
}


## Step 13
# pairごとに代表的なgenic contextを付ける

context_levels <- c(
  "TSS_overlap",
  "first_exon",
  "other_exon",
  "intron",
  "upstream_promoter",
  "downstream_promoter"
)

hub_pairs$main_context <- ifelse(
  hub_pairs$relative_to_TSS == "upstream",
  "upstream_promoter",
  "downstream_promoter"
)

hub_pairs$main_context[
  hub_pairs$overlaps_intron
] <- "intron"

hub_pairs$main_context[
  hub_pairs$overlaps_other_exon
] <- "other_exon"

hub_pairs$main_context[
  hub_pairs$overlaps_first_exon
] <- "first_exon"

hub_pairs$main_context[
  hub_pairs$overlaps_TSS
] <- "TSS_overlap"

hub_pairs$main_context <- factor(
  hub_pairs$main_context,
  levels = context_levels
)


## Step 14
# pair表を見やすい順に並べる

context_rank <- match(
  hub_pairs$main_context,
  context_levels
)

hub_pairs <- hub_pairs[
  order(
    context_rank,
    -abs(hub_pairs$GS_AGE),
    -abs(hub_pairs$MM),
    hub_pairs$hub_module,
    hub_pairs$gene,
    hub_pairs$repeat_start
  ),
]

rownames(hub_pairs) <- NULL

hub_pairs$context_priority_rank <- seq_len(
  nrow(hub_pairs)
)


## Step 15
# genic context全体のsummaryを作る

context_summary <- as.data.frame(
  table(hub_pairs$main_context),
  stringsAsFactors = FALSE
)

colnames(context_summary) <- c(
  "main_context",
  "pair_n"
)

context_summary$pair_percent <-
  context_summary$pair_n /
  sum(context_summary$pair_n) * 100

context_summary


## Step 16
# module × genic contextのsummaryを作る

module_context_summary <- as.data.frame(
  table(
    hub_pairs$hub_module,
    hub_pairs$main_context
  ),
  stringsAsFactors = FALSE
)

colnames(module_context_summary) <- c(
  "module",
  "main_context",
  "pair_n"
)

module_context_summary$module_pair_total <- ave(
  module_context_summary$pair_n,
  module_context_summary$module,
  FUN = sum
)

module_context_summary$within_module_percent <- ifelse(
  module_context_summary$module_pair_total > 0,
  module_context_summary$pair_n /
    module_context_summary$module_pair_total * 100,
  0
)


## Step 17
# promoter-positive hub geneごとのsummaryを作る

pair_split <- split(
  hub_pairs,
  hub_pairs$gene
)

gene_summary_list <- vector(
  "list",
  length(pair_split)
)

gene_names <- names(pair_split)

for (i in seq_along(gene_names)) {
  
  gene_id <- gene_names[i]
  pair_df <- pair_split[[gene_id]]
  hub_row <- hub_annotated[
    match(gene_id, hub_annotated$gene),
  ]
  
  gene_summary_list[[i]] <- data.frame(
    gene = gene_id,
    symbol = hub_row$symbol,
    module = hub_row$module,
    rank_in_module = hub_row$rank_in_module,
    MM = hub_row$MM,
    GS_AGE = hub_row$GS_AGE,
    is_GenAge = hub_row$is_GenAge,
    promoter_ERV_pair_n = nrow(pair_df),
    unique_ERV_fragment_n = length(unique(pair_df$erv_fragment_id)),
    TSS_overlap_n = sum(pair_df$main_context == "TSS_overlap"),
    first_exon_n = sum(pair_df$main_context == "first_exon"),
    other_exon_n = sum(pair_df$main_context == "other_exon"),
    intron_n = sum(pair_df$main_context == "intron"),
    upstream_promoter_n = sum(pair_df$main_context == "upstream_promoter"),
    downstream_promoter_n = sum(pair_df$main_context == "downstream_promoter"),
    same_strand_n = sum(pair_df$strand_relation == "same_strand"),
    opposite_strand_n = sum(pair_df$strand_relation == "opposite_strand"),
    minimum_TSS_distance = min(pair_df$distance_to_TSS),
    ERV_families = paste(
      sort(unique(pair_df$repFamily)),
      collapse = ";"
    ),
    ERV_names = paste(
      sort(unique(pair_df$repName)),
      collapse = ";"
    ),
    stringsAsFactors = FALSE
  )
}

gene_context_summary <- do.call(
  rbind,
  gene_summary_list
)

gene_context_summary <- gene_context_summary[
  order(
    -gene_context_summary$TSS_overlap_n,
    -gene_context_summary$first_exon_n,
    -gene_context_summary$upstream_promoter_n,
    -abs(gene_context_summary$GS_AGE),
    -abs(gene_context_summary$MM)
  ),
]

rownames(gene_context_summary) <- NULL


## Step 18
# contextを優先したpair表を保存用に整理する

priority_columns <- c(
  "context_priority_rank",
  "gene",
  "symbol",
  "hub_module",
  "rank_in_module",
  "MM",
  "GS_AGE",
  "is_GenAge",
  "chromosome",
  "gene_strand",
  "TSS",
  "promoter_start",
  "promoter_end",
  "erv_fragment_id",
  "repeat_start",
  "repeat_end",
  "erv_strand",
  "strand_relation",
  "repName",
  "repFamily",
  "repeat_part",
  "percent_divergence",
  "promoter_overlap_bp",
  "distance_to_TSS",
  "signed_distance_to_TSS",
  "relative_to_TSS",
  "overlaps_TSS",
  "overlaps_first_exon",
  "overlaps_other_exon",
  "overlaps_intron",
  "main_context"
)

hub_pairs_priority <- hub_pairs[
  , priority_columns
]


## Step 19
# 結果を保存する

saveRDS(
  hub_pairs,
  "clean/hub_promoter_ERV_pairs_genic_context.rds"
)

write.csv(
  hub_pairs,
  "results/hub_promoter_ERV_pairs_genic_context.csv",
  row.names = FALSE
)

write.csv(
  context_summary,
  "results/hub_promoter_ERV_genic_context_summary.csv",
  row.names = FALSE
)

write.csv(
  module_context_summary,
  "results/hub_promoter_ERV_context_by_module.csv",
  row.names = FALSE
)

write.csv(
  gene_context_summary,
  "results/hub_promoter_ERV_candidate_gene_summary.csv",
  row.names = FALSE
)

write.csv(
  hub_pairs_priority,
  "results/hub_promoter_ERV_context_priority.csv",
  row.names = FALSE
)


## Step 20
# Figure 1：genic context別pair数

context_plot_data <- context_summary[
  context_summary$pair_n > 0,
]

context_plot_data$main_context <- factor(
  context_plot_data$main_context,
  levels = rev(context_levels)
)

context_plot <- ggplot(
  context_plot_data,
  aes(
    x = main_context,
    y = pair_n
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = pair_n),
    hjust = -0.15,
    size = 3.5
  ) +
  coord_flip() +
  expand_limits(
    y = max(context_plot_data$pair_n) * 1.12
  ) +
  labs(
    x = "Genic context",
    y = "Hub–ERV fragment pairs",
    title = "Genic context of promoter-overlapping ERV fragments",
    subtitle = "Main context uses the priority: TSS > first exon > other exon > intron > promoter"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = "figures/Step12A_hub_ERV_genic_context.pdf",
  plot = context_plot,
  device = "pdf",
  width = 8,
  height = 5
)


## Step 21
# Figure 2：module別のgenic context

module_colors <- c(
  "greenyellow" = "greenyellow",
  "red" = "red",
  "darkgreen" = "darkgreen",
  "turquoise" = "turquoise",
  "midnightblue" = "midnightblue",
  "blue" = "blue",
  "pink" = "pink"
)

module_context_plot_data <- module_context_summary[
  module_context_summary$pair_n > 0,
]

module_context_plot_data$main_context <- factor(
  module_context_plot_data$main_context,
  levels = context_levels
)

module_context_plot <- ggplot(
  module_context_plot_data,
  aes(
    x = main_context,
    y = pair_n,
    fill = module
  )
) +
  geom_col(position = "stack") +
  scale_fill_manual(values = module_colors) +
  labs(
    x = "Genic context",
    y = "Hub–ERV fragment pairs",
    fill = "WGCNA module",
    title = "Genic context of hub-associated ERV fragments by module"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

ggsave(
  filename = "figures/Step12A_hub_ERV_context_by_module.pdf",
  plot = module_context_plot,
  device = "pdf",
  width = 9,
  height = 5.5
)


## Step 22
# 最終確認

cat("\n==============================\n")
cat("hub promoter–ERV pair数:", nrow(hub_pairs), "\n")
cat("promoter-positive hub gene数:", length(unique(hub_pairs$gene)), "\n")
cat(
  "gene–ERV pair重複数:",
  sum(
    duplicated(
      hub_pairs[, c("gene", "erv_fragment_id")]
    )
  ),
  "\n"
)

cat("\nMain context summary:\n")
print(context_summary)

cat("\nBinary overlap flags:\n")
cat("TSS overlap:", sum(hub_pairs$overlaps_TSS), "\n")
cat("first exon overlap:", sum(hub_pairs$overlaps_first_exon), "\n")
cat("other exon overlap:", sum(hub_pairs$overlaps_other_exon), "\n")
cat("intron overlap:", sum(hub_pairs$overlaps_intron), "\n")

cat("\nStrand relation:\n")
print(table(hub_pairs$strand_relation))

cat("\n候補gene summary上位20:\n")
print(head(gene_context_summary, 20))
cat("==============================\n")