# 目的：
# Roadmap Epigenomicsのskeletal muscle chromatin-state annotationを用いて、
# ERV fragmentが骨格筋でenhancer-likeなクロマチン状態と重なるか調べる。
#
# その後、enhancer-associated ERV fragmentが、210 hub genesの
# annotated TSSから±100 kb以内に存在するかを調べる。
#
#
# 【前に使用したENCODE ENCSR530NDCを使わない理由】
#
# ENCSR530NDCのBEDでは、skeletal muscleに対する分類が
#
#   High-H3K4me3
#   Unclassified
#
# のみで、全行が
#
#   Missing-data/Partial-classification
#
# と記録されていた。
#
# pELS・dELSの判定にはDNase accessibilityとH3K27acが必要であるため、
# このファイルからskeletal-muscle enhancerを抽出することはできない。
#
# そこで、このStepではH3K27acを含むRoadmap Epigenomicsの
# 18-state ChromHMM modelを使用する。
#
#
# 【使用する骨格筋reference epigenome】
#
# E100：Psoas muscle
# E108：Skeletal muscle, female
#
# 両方について、Roadmapが提供するGRCh38 lift-over済みBEDを使用する。
#
# 1つのreferenceだけに依存しないように、各ERV fragmentについて、
#
#   support_n = 1：E100またはE108の片方でenhancer state
#   support_n = 2：E100とE108の両方でenhancer state
#
# を記録する。
#
# support_n = 2を、より保守的なconsensus evidenceとして優先する。
#
#
# 【ChromHMM enhancer state】
#
# 主解析では次のactive enhancer statesを使用する。
#
#   EnhG1：Genic enhancer 1
#   EnhG2：Genic enhancer 2
#   EnhA1：Active enhancer 1
#   EnhA2：Active enhancer 2
#
# EnhWkはweak enhancerなので、件数確認用に読み込むが、
# 主解析のactive enhancer setには含めない。
#
#
# 【Enhancerとは】
#
# Enhancerは遺伝子転写を促進し得るcis-regulatory DNA elementである。
# Promoterとは異なり、TSSから離れた位置やintron内にも存在する。
# また、最も近い遺伝子を調節するとは限らない。
#
#
# 【±100 kb window】
#
# 各hub geneのannotated TSSを中心に±100 kbを探索する。
# これはcandidate prioritisation用の操作的な範囲であり、
# enhancer-gene regulationを証明する境界ではない。
#
#
# 【この解析が証明しないこと】
#
# ・ERVがGTEx sampleで発現している
# ・ERVが加齢により活性化している
# ・ERVが実際にenhancerとして機能する
# ・近傍hub geneがそのenhancerの標的である
#
# 今回示すのは、reference genomeとreference epigenomeを用いた
# expression-informed genomic candidate associationである。
#
#
# 入力：
# clean/RepeatMasker_hg38_ERV_primary_GRCh38.rds
# clean/hub_gene_ERV_annotation_all.rds
# clean/hub_promoter_ERV_pairs_genic_context.rds
#
# 外部データ：
# Roadmap Epigenomics 18-state ChromHMM, GRCh38 lift-over
# E100 / E108
#
# 主な出力：
# raw/E100_18_core_K27ac_hg38lift_mnemonics.bed.gz
# raw/E108_18_core_K27ac_hg38lift_mnemonics.bed.gz
# results/Roadmap_muscle_chromatin_state_counts.csv
# results/Roadmap_muscle_active_enhancer_ERV_overlap_pairs.csv
# clean/Roadmap_muscle_active_enhancer_ERV_summary.rds
# results/hub_Roadmap_muscle_enhancer_ERV_pairs_100kb.csv
# results/hub_Roadmap_muscle_enhancer_ERV_gene_summary_100kb.csv
# results/GenAge_hub_Roadmap_muscle_enhancer_ERV_summary_100kb.csv
# results/hub_Roadmap_muscle_enhancer_ERV_priority_100kb.csv
# figures/Step12B_Roadmap_ERV_family_by_epigenome.pdf
# figures/Step12B_Roadmap_hub_enhancer_ERV_percent_by_module.pdf
# ============================================================


## Step 1
# 出力フォルダを作る

dir.create("raw", showWarnings = FALSE)
dir.create("clean", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


## Step 2
# 必要なpackageを準備する

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(GenomicRanges)
library(ggplot2)


## Step 3
# 解析済みデータを読む

erv_annotation <- readRDS(
  "clean/RepeatMasker_hg38_ERV_primary_GRCh38.rds"
)

hub_annotated <- readRDS(
  "clean/hub_gene_ERV_annotation_all.rds"
)

hub_promoter_context <- readRDS(
  "clean/hub_promoter_ERV_pairs_genic_context.rds"
)

cat("主解析用ERV fragment数:", nrow(erv_annotation), "\n")
cat("hub gene数:", nrow(hub_annotated), "\n")
cat("hub promoter-ERV context pair数:", nrow(hub_promoter_context), "\n")


## Step 4
# 必要な列を確認する

required_erv_columns <- c(
  "erv_fragment_id",
  "chromosome",
  "repeat_start",
  "repeat_end",
  "repeat_length",
  "strand",
  "repName",
  "repFamily",
  "repeat_part",
  "percent_divergence"
)

required_hub_columns <- c(
  "gene",
  "symbol",
  "module",
  "rank_in_module",
  "MM",
  "GS_AGE",
  "is_GenAge",
  "chromosome",
  "strand",
  "TSS"
)

required_context_columns <- c(
  "gene",
  "erv_fragment_id",
  "main_context",
  "overlaps_TSS",
  "overlaps_first_exon",
  "overlaps_other_exon",
  "overlaps_intron"
)

missing_erv_columns <- setdiff(
  required_erv_columns,
  colnames(erv_annotation)
)

missing_hub_columns <- setdiff(
  required_hub_columns,
  colnames(hub_annotated)
)

missing_context_columns <- setdiff(
  required_context_columns,
  colnames(hub_promoter_context)
)

if (length(missing_erv_columns) > 0) {
  stop(
    paste0(
      "ERV annotationに必要な列がありません: ",
      paste(missing_erv_columns, collapse = ", ")
    )
  )
}

if (length(missing_hub_columns) > 0) {
  stop(
    paste0(
      "hub annotationに必要な列がありません: ",
      paste(missing_hub_columns, collapse = ", ")
    )
  )
}

if (length(missing_context_columns) > 0) {
  stop(
    paste0(
      "Step 12A context表に必要な列がありません: ",
      paste(missing_context_columns, collapse = ", ")
    )
  )
}

if (anyDuplicated(hub_annotated$gene) > 0) {
  stop("hub annotationのgene IDに重複があります。")
}


## Step 5
# Roadmap ChromHMMファイルの場所を指定する

roadmap_base_url <- paste0(
  "https://egg2.wustl.edu/roadmap/data/byFileType/",
  "chromhmmSegmentations/ChmmModels/core_K27ac/",
  "jointModel/final/"
)

epigenome_metadata <- data.frame(
  epigenome = c("E100", "E108"),
  tissue = c(
    "Psoas muscle",
    "Skeletal muscle female"
  ),
  filename = c(
    "E100_18_core_K27ac_hg38lift_mnemonics.bed.gz",
    "E108_18_core_K27ac_hg38lift_mnemonics.bed.gz"
  ),
  stringsAsFactors = FALSE
)

epigenome_metadata$download_url <- paste0(
  roadmap_base_url,
  epigenome_metadata$filename
)

epigenome_metadata$local_file <- file.path(
  "raw",
  epigenome_metadata$filename
)

write.csv(
  epigenome_metadata,
  "results/Roadmap_muscle_epigenome_download_metadata.csv",
  row.names = FALSE
)


## Step 6
# Roadmapファイルをダウンロードする

for (i in seq_len(nrow(epigenome_metadata))) {
  
  local_file <- epigenome_metadata$local_file[i]
  download_url <- epigenome_metadata$download_url[i]
  
  if (!file.exists(local_file)) {
    download.file(
      url = download_url,
      destfile = local_file,
      mode = "wb",
      method = "libcurl"
    )
  }
  
  if (!file.exists(local_file)) {
    stop(
      paste0(
        epigenome_metadata$epigenome[i],
        "のRoadmap BEDを取得できませんでした。"
      )
    )
  }
  
  if (file.info(local_file)$size == 0) {
    stop(
      paste0(
        epigenome_metadata$epigenome[i],
        "のRoadmap BEDのサイズが0です。"
      )
    )
  }
  
  cat(
    epigenome_metadata$epigenome[i],
    "file:",
    local_file,
    "size:",
    file.info(local_file)$size,
    "bytes\n"
  )
}


## Step 7
# Roadmap 4-column BEDを読む関数

read_roadmap_chromhmm <- function(
    bed_file,
    epigenome_id,
    tissue_name
) {
  
  chromhmm_raw <- read.delim(
    bed_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "#",
    quote = ""
  )
  
  if (ncol(chromhmm_raw) < 4) {
    stop(
      paste0(
        epigenome_id,
        "のChromHMM BEDが4列未満です。"
      )
    )
  }
  
  chromhmm_raw <- chromhmm_raw[, 1:4]
  
  colnames(chromhmm_raw) <- c(
    "chromosome",
    "start_0based",
    "end",
    "chromatin_state"
  )
  
  chromhmm_raw$chromosome <- as.character(
    chromhmm_raw$chromosome
  )
  
  chromhmm_raw$chromatin_state <- as.character(
    chromhmm_raw$chromatin_state
  )
  
  # Roadmapの番号付きstate名を通常のmnemonicへ変換
  # 例：7_EnhG1 → EnhG1
  chromhmm_raw$chromatin_state <- sub(
    "^[0-9]+_",
    "",
    chromhmm_raw$chromatin_state
  )
  # BED startは0-basedなので1を足す
  chromhmm_raw$state_start <- as.integer(
    chromhmm_raw$start_0based
  ) + 1L
  
  # BED endはhalf-openなので、その値をinclusive endとして使う
  chromhmm_raw$state_end <- as.integer(
    chromhmm_raw$end
  )
  
  chromhmm_raw$epigenome <- epigenome_id
  chromhmm_raw$tissue <- tissue_name
  
  standard_chromosomes <- c(
    paste0("chr", 1:22),
    "chrX",
    "chrY"
  )
  
  valid_rows <-
    chromhmm_raw$chromosome %in% standard_chromosomes &
    !is.na(chromhmm_raw$state_start) &
    !is.na(chromhmm_raw$state_end) &
    chromhmm_raw$state_start >= 1 &
    chromhmm_raw$state_end >= chromhmm_raw$state_start
  
  chromhmm_raw <- chromhmm_raw[
    valid_rows,
    c(
      "epigenome",
      "tissue",
      "chromosome",
      "state_start",
      "state_end",
      "chromatin_state"
    )
  ]
  
  rownames(chromhmm_raw) <- NULL
  
  chromhmm_raw
}


## Step 8
# E100とE108を読み込む

chromhmm_list <- vector(
  mode = "list",
  length = nrow(epigenome_metadata)
)

for (i in seq_len(nrow(epigenome_metadata))) {
  
  chromhmm_list[[i]] <- read_roadmap_chromhmm(
    bed_file = epigenome_metadata$local_file[i],
    epigenome_id = epigenome_metadata$epigenome[i],
    tissue_name = epigenome_metadata$tissue[i]
  )
  
  cat(
    epigenome_metadata$epigenome[i],
    "ChromHMM interval数:",
    nrow(chromhmm_list[[i]]),
    "\n"
  )
}

names(chromhmm_list) <- epigenome_metadata$epigenome


## Step 9
# ChromHMM stateの内訳を保存する

chromatin_state_summary <- do.call(
  rbind,
  lapply(
    names(chromhmm_list),
    function(epigenome_id) {
      
      state_table <- as.data.frame(
        table(
          chromhmm_list[[epigenome_id]]$chromatin_state
        ),
        stringsAsFactors = FALSE
      )
      
      colnames(state_table) <- c(
        "chromatin_state",
        "interval_n"
      )
      
      state_table$epigenome <- epigenome_id
      
      state_table[, c(
        "epigenome",
        "chromatin_state",
        "interval_n"
      )]
    }
  )
)

write.csv(
  chromatin_state_summary,
  "results/Roadmap_muscle_chromatin_state_counts.csv",
  row.names = FALSE
)

cat("\nChromHMM state summary:\n")
print(chromatin_state_summary)


## Step 10
# Active enhancer statesを抽出する

active_enhancer_states <- c(
  "EnhG1",
  "EnhG2",
  "EnhA1",
  "EnhA2"
)

weak_enhancer_state <- "EnhWk"

active_enhancer_list <- lapply(
  chromhmm_list,
  function(chromhmm_df) {
    chromhmm_df[
      chromhmm_df$chromatin_state %in%
        active_enhancer_states,
      ,
      drop = FALSE
    ]
  }
)

weak_enhancer_list <- lapply(
  chromhmm_list,
  function(chromhmm_df) {
    chromhmm_df[
      chromhmm_df$chromatin_state ==
        weak_enhancer_state,
      ,
      drop = FALSE
    ]
  }
)

for (epigenome_id in names(active_enhancer_list)) {
  cat(
    epigenome_id,
    "active enhancer interval数:",
    nrow(active_enhancer_list[[epigenome_id]]),
    "\n"
  )
  
  cat(
    epigenome_id,
    "weak enhancer interval数:",
    nrow(weak_enhancer_list[[epigenome_id]]),
    "\n"
  )
  
  if (nrow(active_enhancer_list[[epigenome_id]]) == 0) {
    stop(
      paste0(
        epigenome_id,
        "でactive enhancer stateが見つかりません。"
      )
    )
  }
}

saveRDS(
  active_enhancer_list,
  "clean/Roadmap_muscle_active_enhancer_intervals_GRCh38.rds"
)

saveRDS(
  weak_enhancer_list,
  "clean/Roadmap_muscle_weak_enhancer_intervals_GRCh38.rds"
)


## Step 11
# ERV annotationをGRangesへ変換する

erv_gr <- GRanges(
  seqnames = erv_annotation$chromosome,
  ranges = IRanges(
    start = erv_annotation$repeat_start,
    end = erv_annotation$repeat_end
  ),
  strand = erv_annotation$strand
)


## Step 12
# 各epigenomeでERVとactive enhancer stateを重ねる

enhancer_erv_pair_list <- list()

for (epigenome_id in names(active_enhancer_list)) {
  
  enhancer_df <- active_enhancer_list[[epigenome_id]]
  
  enhancer_gr <- GRanges(
    seqnames = enhancer_df$chromosome,
    ranges = IRanges(
      start = enhancer_df$state_start,
      end = enhancer_df$state_end
    ),
    strand = "*"
  )
  
  overlap_hits <- findOverlaps(
    erv_gr,
    enhancer_gr,
    ignore.strand = TRUE
  )
  
  erv_index <- queryHits(overlap_hits)
  enhancer_index <- subjectHits(overlap_hits)
  
  cat(
    epigenome_id,
    "ERV-enhancer overlap pair数:",
    length(overlap_hits),
    "\n"
  )
  
  cat(
    epigenome_id,
    "enhancerと重なるunique ERV fragment数:",
    length(unique(erv_index)),
    "\n"
  )
  
  enhancer_erv_pairs <- data.frame(
    epigenome = epigenome_id,
    tissue = enhancer_df$tissue[enhancer_index],
    erv_fragment_id = erv_annotation$erv_fragment_id[erv_index],
    chromosome = erv_annotation$chromosome[erv_index],
    repeat_start = erv_annotation$repeat_start[erv_index],
    repeat_end = erv_annotation$repeat_end[erv_index],
    repeat_length = erv_annotation$repeat_length[erv_index],
    erv_strand = erv_annotation$strand[erv_index],
    repName = erv_annotation$repName[erv_index],
    repFamily = erv_annotation$repFamily[erv_index],
    repeat_part = erv_annotation$repeat_part[erv_index],
    percent_divergence = erv_annotation$percent_divergence[erv_index],
    enhancer_start = enhancer_df$state_start[enhancer_index],
    enhancer_end = enhancer_df$state_end[enhancer_index],
    enhancer_state = enhancer_df$chromatin_state[enhancer_index],
    stringsAsFactors = FALSE
  )
  
  enhancer_erv_pairs <- enhancer_erv_pairs[
    !duplicated(
      enhancer_erv_pairs[, c(
        "epigenome",
        "erv_fragment_id",
        "enhancer_start",
        "enhancer_end",
        "enhancer_state"
      )]
    ),
    ,
    drop = FALSE
  ]
  
  rownames(enhancer_erv_pairs) <- NULL
  
  enhancer_erv_pair_list[[epigenome_id]] <- enhancer_erv_pairs
}

all_enhancer_ERV_pairs <- do.call(
  rbind,
  enhancer_erv_pair_list
)

rownames(all_enhancer_ERV_pairs) <- NULL

write.csv(
  all_enhancer_ERV_pairs,
  "results/Roadmap_muscle_active_enhancer_ERV_overlap_pairs.csv",
  row.names = FALSE
)

saveRDS(
  all_enhancer_ERV_pairs,
  "clean/Roadmap_muscle_active_enhancer_ERV_overlap_pairs.rds"
)


## Step 13
# ERV fragmentごとにE100・E108のsupportをまとめる

erv_pair_rows <- split(
  seq_len(nrow(all_enhancer_ERV_pairs)),
  all_enhancer_ERV_pairs$erv_fragment_id
)

active_enhancer_ERV_summary <- do.call(
  rbind,
  lapply(
    names(erv_pair_rows),
    function(fragment_id) {
      
      rows <- erv_pair_rows[[fragment_id]]
      first_row <- rows[1]
      
      supporting_epigenomes <- sort(
        unique(
          all_enhancer_ERV_pairs$epigenome[rows]
        )
      )
      
      data.frame(
        erv_fragment_id = fragment_id,
        chromosome = all_enhancer_ERV_pairs$chromosome[first_row],
        repeat_start = all_enhancer_ERV_pairs$repeat_start[first_row],
        repeat_end = all_enhancer_ERV_pairs$repeat_end[first_row],
        repeat_length = all_enhancer_ERV_pairs$repeat_length[first_row],
        erv_strand = all_enhancer_ERV_pairs$erv_strand[first_row],
        repName = all_enhancer_ERV_pairs$repName[first_row],
        repFamily = all_enhancer_ERV_pairs$repFamily[first_row],
        repeat_part = all_enhancer_ERV_pairs$repeat_part[first_row],
        percent_divergence = all_enhancer_ERV_pairs$percent_divergence[first_row],
        epigenome_support_n = length(supporting_epigenomes),
        supporting_epigenomes = paste(
          supporting_epigenomes,
          collapse = ";"
        ),
        supported_in_E100 = "E100" %in% supporting_epigenomes,
        supported_in_E108 = "E108" %in% supporting_epigenomes,
        enhancer_states = paste(
          sort(
            unique(
              all_enhancer_ERV_pairs$enhancer_state[rows]
            )
          ),
          collapse = ";"
        ),
        enhancer_overlap_pair_n = length(rows),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(active_enhancer_ERV_summary) <- NULL

active_enhancer_ERV_summary$consensus_support <-
  active_enhancer_ERV_summary$epigenome_support_n == 2

saveRDS(
  active_enhancer_ERV_summary,
  "clean/Roadmap_muscle_active_enhancer_ERV_summary.rds"
)

write.csv(
  active_enhancer_ERV_summary,
  "results/Roadmap_muscle_active_enhancer_ERV_summary.csv",
  row.names = FALSE
)

cat(
  "active enhancer-associated unique ERV fragment数:",
  nrow(active_enhancer_ERV_summary),
  "\n"
)

cat(
  "E100とE108の両方でsupportされたERV fragment数:",
  sum(active_enhancer_ERV_summary$consensus_support),
  "\n"
)


## Step 14
# 210 hub genesのTSS ±100 kb windowを作る

hub_window_size <- 100000L

hub_window_start <- pmax(
  as.integer(hub_annotated$TSS) - hub_window_size,
  1L
)

hub_window_end <-
  as.integer(hub_annotated$TSS) + hub_window_size

hub_window_gr <- GRanges(
  seqnames = hub_annotated$chromosome,
  ranges = IRanges(
    start = hub_window_start,
    end = hub_window_end
  ),
  strand = hub_annotated$strand
)

enhancer_erv_gr <- GRanges(
  seqnames = active_enhancer_ERV_summary$chromosome,
  ranges = IRanges(
    start = active_enhancer_ERV_summary$repeat_start,
    end = active_enhancer_ERV_summary$repeat_end
  ),
  strand = active_enhancer_ERV_summary$erv_strand
)


## Step 15
# Hub TSS ±100 kbにあるenhancer-associated ERVを探す

hub_enhancer_hits <- findOverlaps(
  hub_window_gr,
  enhancer_erv_gr,
  ignore.strand = TRUE
)

hub_index <- queryHits(hub_enhancer_hits)
erv_index <- subjectHits(hub_enhancer_hits)

cat(
  "hub-enhancer ERV pair数（±100 kb）:",
  length(hub_enhancer_hits),
  "\n"
)

cat(
  "enhancer-associated ERVを持つhub gene数:",
  length(unique(hub_index)),
  "\n"
)

if (length(hub_enhancer_hits) == 0) {
  stop("Hub TSS ±100 kbにenhancer-associated ERVがありません。")
}


## Step 16
# Hub-enhancer ERV pair表を作る

hub_enhancer_pairs <- data.frame(
  gene = hub_annotated$gene[hub_index],
  symbol = hub_annotated$symbol[hub_index],
  module = hub_annotated$module[hub_index],
  rank_in_module = hub_annotated$rank_in_module[hub_index],
  MM = hub_annotated$MM[hub_index],
  GS_AGE = hub_annotated$GS_AGE[hub_index],
  is_GenAge = hub_annotated$is_GenAge[hub_index],
  chromosome = hub_annotated$chromosome[hub_index],
  gene_strand = hub_annotated$strand[hub_index],
  TSS = hub_annotated$TSS[hub_index],
  window_start = hub_window_start[hub_index],
  window_end = hub_window_end[hub_index],
  erv_fragment_id = active_enhancer_ERV_summary$erv_fragment_id[erv_index],
  repeat_start = active_enhancer_ERV_summary$repeat_start[erv_index],
  repeat_end = active_enhancer_ERV_summary$repeat_end[erv_index],
  erv_strand = active_enhancer_ERV_summary$erv_strand[erv_index],
  repName = active_enhancer_ERV_summary$repName[erv_index],
  repFamily = active_enhancer_ERV_summary$repFamily[erv_index],
  repeat_part = active_enhancer_ERV_summary$repeat_part[erv_index],
  percent_divergence = active_enhancer_ERV_summary$percent_divergence[erv_index],
  epigenome_support_n = active_enhancer_ERV_summary$epigenome_support_n[erv_index],
  supporting_epigenomes = active_enhancer_ERV_summary$supporting_epigenomes[erv_index],
  supported_in_E100 = active_enhancer_ERV_summary$supported_in_E100[erv_index],
  supported_in_E108 = active_enhancer_ERV_summary$supported_in_E108[erv_index],
  consensus_support = active_enhancer_ERV_summary$consensus_support[erv_index],
  enhancer_states = active_enhancer_ERV_summary$enhancer_states[erv_index],
  stringsAsFactors = FALSE
)

hub_enhancer_pairs <- hub_enhancer_pairs[
  !duplicated(
    hub_enhancer_pairs[, c(
      "gene",
      "erv_fragment_id"
    )]
  ),
  ,
  drop = FALSE
]

rownames(hub_enhancer_pairs) <- NULL


## Step 17
# TSSからERVまでの距離と方向を計算する

hub_enhancer_pairs$TSS_distance <- ifelse(
  hub_enhancer_pairs$repeat_start <= hub_enhancer_pairs$TSS &
    hub_enhancer_pairs$repeat_end >= hub_enhancer_pairs$TSS,
  0L,
  pmin(
    abs(
      hub_enhancer_pairs$TSS -
        hub_enhancer_pairs$repeat_start
    ),
    abs(
      hub_enhancer_pairs$TSS -
        hub_enhancer_pairs$repeat_end
    )
  )
)

hub_enhancer_pairs$TSS_direction <- ifelse(
  hub_enhancer_pairs$TSS_distance == 0,
  "overlaps_TSS",
  ifelse(
    hub_enhancer_pairs$gene_strand == "+",
    ifelse(
      hub_enhancer_pairs$repeat_end <
        hub_enhancer_pairs$TSS,
      "upstream",
      "downstream"
    ),
    ifelse(
      hub_enhancer_pairs$repeat_start >
        hub_enhancer_pairs$TSS,
      "upstream",
      "downstream"
    )
  )
)

hub_enhancer_pairs$signed_TSS_distance <-
  hub_enhancer_pairs$TSS_distance

hub_enhancer_pairs$signed_TSS_distance[
  hub_enhancer_pairs$TSS_direction == "upstream"
] <- -hub_enhancer_pairs$TSS_distance[
  hub_enhancer_pairs$TSS_direction == "upstream"
]

hub_enhancer_pairs$same_strand <-
  hub_enhancer_pairs$gene_strand ==
  hub_enhancer_pairs$erv_strand

hub_enhancer_pairs$proximity_group <- cut(
  hub_enhancer_pairs$TSS_distance,
  breaks = c(-1, 0, 10000, 50000, 100000),
  labels = c(
    "TSS_overlap",
    "within_10kb",
    "10_to_50kb",
    "50_to_100kb"
  )
)


## Step 18
# Step 12Aのpromoter・exon・intron contextを追加する

context_key <- paste(
  hub_promoter_context$gene,
  hub_promoter_context$erv_fragment_id,
  sep = "||"
)

pair_key <- paste(
  hub_enhancer_pairs$gene,
  hub_enhancer_pairs$erv_fragment_id,
  sep = "||"
)

context_index <- match(
  pair_key,
  context_key
)

hub_enhancer_pairs$is_promoter_ERV_pair <-
  !is.na(context_index)

hub_enhancer_pairs$promoter_main_context <-
  hub_promoter_context$main_context[context_index]

hub_enhancer_pairs$overlaps_TSS <-
  hub_promoter_context$overlaps_TSS[context_index]

hub_enhancer_pairs$overlaps_first_exon <-
  hub_promoter_context$overlaps_first_exon[context_index]

hub_enhancer_pairs$overlaps_other_exon <-
  hub_promoter_context$overlaps_other_exon[context_index]

hub_enhancer_pairs$overlaps_intron <-
  hub_promoter_context$overlaps_intron[context_index]

binary_context_columns <- c(
  "overlaps_TSS",
  "overlaps_first_exon",
  "overlaps_other_exon",
  "overlaps_intron"
)

for (column_name in binary_context_columns) {
  hub_enhancer_pairs[[column_name]][
    is.na(hub_enhancer_pairs[[column_name]])
  ] <- FALSE
}


## Step 19
# Candidate groupと優先順位を付ける

hub_enhancer_pairs$candidate_group <- ifelse(
  hub_enhancer_pairs$is_GenAge &
    hub_enhancer_pairs$consensus_support,
  "GenAge_consensus_enhancer_ERV",
  ifelse(
    hub_enhancer_pairs$is_GenAge,
    "GenAge_enhancer_ERV",
    ifelse(
      hub_enhancer_pairs$is_promoter_ERV_pair &
        hub_enhancer_pairs$consensus_support,
      "promoter_consensus_enhancer_ERV",
      ifelse(
        hub_enhancer_pairs$consensus_support,
        "consensus_enhancer_ERV",
        "single_epigenome_enhancer_ERV"
      )
    )
  )
)

candidate_group_order <- c(
  "GenAge_consensus_enhancer_ERV",
  "GenAge_enhancer_ERV",
  "promoter_consensus_enhancer_ERV",
  "consensus_enhancer_ERV",
  "single_epigenome_enhancer_ERV"
)

hub_enhancer_pairs$candidate_group_rank <- match(
  hub_enhancer_pairs$candidate_group,
  candidate_group_order
)

priority_order <- order(
  hub_enhancer_pairs$candidate_group_rank,
  -abs(hub_enhancer_pairs$GS_AGE),
  -abs(hub_enhancer_pairs$MM),
  hub_enhancer_pairs$TSS_distance,
  hub_enhancer_pairs$gene,
  hub_enhancer_pairs$erv_fragment_id
)

hub_enhancer_priority <- hub_enhancer_pairs[
  priority_order,
  ,
  drop = FALSE
]

hub_enhancer_priority$priority_rank <- seq_len(
  nrow(hub_enhancer_priority)
)

hub_enhancer_priority <- hub_enhancer_priority[, c(
  "priority_rank",
  setdiff(
    colnames(hub_enhancer_priority),
    "priority_rank"
  )
)]


## Step 20
# Hub geneごとのsummaryを作る

pair_rows_by_gene <- split(
  seq_len(nrow(hub_enhancer_pairs)),
  hub_enhancer_pairs$gene
)

hub_enhancer_gene_summary <- hub_annotated

hub_enhancer_gene_summary$enhancer_ERV_any <- FALSE
hub_enhancer_gene_summary$enhancer_ERV_n <- 0L
hub_enhancer_gene_summary$consensus_enhancer_ERV_n <- 0L
hub_enhancer_gene_summary$promoter_and_enhancer_ERV_n <- 0L
hub_enhancer_gene_summary$minimum_enhancer_ERV_distance <- NA_integer_
hub_enhancer_gene_summary$enhancer_ERV_families <- NA_character_
hub_enhancer_gene_summary$enhancer_ERV_names <- NA_character_
hub_enhancer_gene_summary$enhancer_states <- NA_character_
hub_enhancer_gene_summary$supporting_epigenomes <- NA_character_

for (gene_id in names(pair_rows_by_gene)) {
  
  rows <- pair_rows_by_gene[[gene_id]]
  hub_row <- match(
    gene_id,
    hub_enhancer_gene_summary$gene
  )
  
  hub_enhancer_gene_summary$enhancer_ERV_any[hub_row] <- TRUE
  
  hub_enhancer_gene_summary$enhancer_ERV_n[hub_row] <-
    length(
      unique(
        hub_enhancer_pairs$erv_fragment_id[rows]
      )
    )
  
  hub_enhancer_gene_summary$consensus_enhancer_ERV_n[hub_row] <-
    length(
      unique(
        hub_enhancer_pairs$erv_fragment_id[
          rows[
            hub_enhancer_pairs$consensus_support[rows]
          ]
        ]
      )
    )
  
  hub_enhancer_gene_summary$promoter_and_enhancer_ERV_n[hub_row] <-
    length(
      unique(
        hub_enhancer_pairs$erv_fragment_id[
          rows[
            hub_enhancer_pairs$is_promoter_ERV_pair[rows]
          ]
        ]
      )
    )
  
  hub_enhancer_gene_summary$minimum_enhancer_ERV_distance[hub_row] <-
    min(
      hub_enhancer_pairs$TSS_distance[rows],
      na.rm = TRUE
    )
  
  hub_enhancer_gene_summary$enhancer_ERV_families[hub_row] <-
    paste(
      sort(
        unique(
          hub_enhancer_pairs$repFamily[rows]
        )
      ),
      collapse = ";"
    )
  
  hub_enhancer_gene_summary$enhancer_ERV_names[hub_row] <-
    paste(
      sort(
        unique(
          hub_enhancer_pairs$repName[rows]
        )
      ),
      collapse = ";"
    )
  
  hub_enhancer_gene_summary$enhancer_states[hub_row] <-
    paste(
      sort(
        unique(
          hub_enhancer_pairs$enhancer_states[rows]
        )
      ),
      collapse = ";"
    )
  
  hub_enhancer_gene_summary$supporting_epigenomes[hub_row] <-
    paste(
      sort(
        unique(
          unlist(
            strsplit(
              hub_enhancer_pairs$supporting_epigenomes[rows],
              ";",
              fixed = TRUE
            )
          )
        )
      ),
      collapse = ";"
    )
}


## Step 21
# Module別summaryを作る

module_list <- unique(hub_enhancer_gene_summary$module)

module_summary <- do.call(
  rbind,
  lapply(
    module_list,
    function(module_name) {
      
      module_rows <- hub_enhancer_gene_summary$module ==
        module_name
      
      data.frame(
        module = module_name,
        hub_gene_n = sum(module_rows),
        enhancer_ERV_hub_n = sum(
          hub_enhancer_gene_summary$enhancer_ERV_any[
            module_rows
          ]
        ),
        enhancer_ERV_hub_percent = mean(
          hub_enhancer_gene_summary$enhancer_ERV_any[
            module_rows
          ]
        ) * 100,
        consensus_enhancer_ERV_hub_n = sum(
          hub_enhancer_gene_summary$consensus_enhancer_ERV_n[
            module_rows
          ] > 0
        ),
        consensus_enhancer_ERV_hub_percent = mean(
          hub_enhancer_gene_summary$consensus_enhancer_ERV_n[
            module_rows
          ] > 0
        ) * 100,
        promoter_and_enhancer_hub_n = sum(
          hub_enhancer_gene_summary$promoter_and_enhancer_ERV_n[
            module_rows
          ] > 0
        ),
        median_minimum_distance = median(
          hub_enhancer_gene_summary$minimum_enhancer_ERV_distance[
            module_rows
          ],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(module_summary) <- NULL


## Step 22
# GenAge hubsだけを抽出する

GenAge_hub_enhancer_summary <- hub_enhancer_gene_summary[
  hub_enhancer_gene_summary$is_GenAge,
  ,
  drop = FALSE
]


## Step 23
# 結果を保存する

saveRDS(
  hub_enhancer_pairs,
  "clean/hub_Roadmap_muscle_enhancer_ERV_pairs_100kb.rds"
)

write.csv(
  hub_enhancer_pairs,
  "results/hub_Roadmap_muscle_enhancer_ERV_pairs_100kb.csv",
  row.names = FALSE
)

write.csv(
  hub_enhancer_gene_summary,
  "results/hub_Roadmap_muscle_enhancer_ERV_gene_summary_100kb.csv",
  row.names = FALSE
)

write.csv(
  GenAge_hub_enhancer_summary,
  "results/GenAge_hub_Roadmap_muscle_enhancer_ERV_summary_100kb.csv",
  row.names = FALSE
)

write.csv(
  hub_enhancer_priority,
  "results/hub_Roadmap_muscle_enhancer_ERV_priority_100kb.csv",
  row.names = FALSE
)

write.csv(
  head(hub_enhancer_priority, 50),
  "results/hub_Roadmap_muscle_enhancer_ERV_priority_top50_100kb.csv",
  row.names = FALSE
)

write.csv(
  module_summary,
  "results/hub_Roadmap_muscle_enhancer_ERV_summary_by_module.csv",
  row.names = FALSE
)


## Step 24
# Figure 1：epigenome別・ERV family別のoverlap数

family_epigenome_summary <- aggregate(
  erv_fragment_id ~ epigenome + repFamily,
  data = unique(
    all_enhancer_ERV_pairs[, c(
      "epigenome",
      "repFamily",
      "erv_fragment_id"
    )]
  ),
  FUN = length
)

colnames(family_epigenome_summary)[3] <-
  "unique_ERV_fragment_n"

family_plot <- ggplot(
  family_epigenome_summary,
  aes(
    x = repFamily,
    y = unique_ERV_fragment_n,
    fill = epigenome
  )
) +
  geom_col(
    position = "dodge"
  ) +
  labs(
    x = "ERV family",
    y = "Unique ERV fragments overlapping active enhancer states",
    title = "Roadmap skeletal-muscle enhancer-associated ERV fragments",
    subtitle = "E100: psoas muscle; E108: skeletal muscle female"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = "figures/Step12B_Roadmap_ERV_family_by_epigenome.pdf",
  plot = family_plot,
  device = grDevices::pdf,
  width = 8,
  height = 5
)


## Step 25
# Figure 2：module別hub gene割合

module_plot_df <- rbind(
  data.frame(
    module = module_summary$module,
    evidence = "Any active enhancer support",
    hub_percent = module_summary$enhancer_ERV_hub_percent,
    stringsAsFactors = FALSE
  ),
  data.frame(
    module = module_summary$module,
    evidence = "Consensus E100 and E108",
    hub_percent = module_summary$consensus_enhancer_ERV_hub_percent,
    stringsAsFactors = FALSE
  )
)

module_colors <- c(
  "greenyellow" = "greenyellow",
  "red" = "red",
  "darkgreen" = "darkgreen",
  "turquoise" = "turquoise",
  "midnightblue" = "midnightblue",
  "blue" = "blue",
  "pink" = "pink"
)

module_plot_df$module <- factor(
  module_plot_df$module,
  levels = c(
    "greenyellow",
    "red",
    "darkgreen",
    "turquoise",
    "midnightblue",
    "blue",
    "pink"
  )
)

module_plot <- ggplot(
  module_plot_df,
  aes(
    x = module,
    y = hub_percent,
    fill = module,
    alpha = evidence
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  scale_fill_manual(
    values = module_colors,
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(
      "Any active enhancer support" = 0.55,
      "Consensus E100 and E108" = 1
    ),
    name = "Evidence"
  ) +
  coord_flip() +
  labs(
    x = "WGCNA module",
    y = "Hub genes with enhancer-associated ERV within 100 kb (%)",
    title = "Skeletal-muscle enhancer-associated ERVs near hub genes"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  filename = "figures/Step12B_Roadmap_hub_enhancer_ERV_percent_by_module.pdf",
  plot = module_plot,
  device = grDevices::pdf,
  width = 9,
  height = 5.5
)


## Step 26
# 最終確認

cat("\n==============================\n")

cat(
  "active enhancer-associated unique ERV fragment数:",
  nrow(active_enhancer_ERV_summary),
  "\n"
)

cat(
  "consensus support ERV fragment数:",
  sum(active_enhancer_ERV_summary$consensus_support),
  "\n"
)

cat(
  "hub-gene enhancer-associated ERV pair数:",
  nrow(hub_enhancer_pairs),
  "\n"
)

cat(
  "enhancer-associated ERVを100 kb以内に持つhub gene数:",
  sum(hub_enhancer_gene_summary$enhancer_ERV_any),
  "\n"
)

cat(
  "consensus enhancer-associated ERVを持つhub gene数:",
  sum(
    hub_enhancer_gene_summary$consensus_enhancer_ERV_n > 0
  ),
  "\n"
)

cat(
  "promoterかつenhancer-associated ERVを持つhub gene数:",
  sum(
    hub_enhancer_gene_summary$promoter_and_enhancer_ERV_n > 0
  ),
  "\n"
)

cat(
  "GenAge hubのうちenhancer-associated ERV candidateあり:",
  sum(GenAge_hub_enhancer_summary$enhancer_ERV_any),
  "/",
  nrow(GenAge_hub_enhancer_summary),
  "\n"
)

cat("\nModule別summary:\n")
print(module_summary)

cat("\nGenAge hub summary:\n")
print(
  GenAge_hub_enhancer_summary[, c(
    "gene",
    "symbol",
    "module",
    "MM",
    "GS_AGE",
    "enhancer_ERV_any",
    "enhancer_ERV_n",
    "consensus_enhancer_ERV_n",
    "promoter_and_enhancer_ERV_n",
    "minimum_enhancer_ERV_distance",
    "enhancer_ERV_families",
    "enhancer_ERV_names",
    "supporting_epigenomes"
  )]
)

cat("\n優先順位上位20 pairs:\n")
print(
  head(
    hub_enhancer_priority[, c(
      "priority_rank",
      "gene",
      "symbol",
      "module",
      "MM",
      "GS_AGE",
      "is_GenAge",
      "erv_fragment_id",
      "repName",
      "repFamily",
      "epigenome_support_n",
      "supporting_epigenomes",
      "enhancer_states",
      "TSS_distance",
      "is_promoter_ERV_pair",
      "candidate_group"
    )],
    20
  )
)

cat("==============================\n")