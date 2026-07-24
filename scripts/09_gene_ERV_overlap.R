# ============================================================
# 09_gene_ERV_overlap.R
#
# 目的：
# GTEx V10骨格筋WGCNA遺伝子と、
# GRCh38上のRepeatMasker ERV fragmentの位置関係を調べる
#
#
# 【TSSとは】
#
# TSSはTranscription Start Siteの略で、
# 遺伝子の転写が始まるゲノム上の位置を指す。
#
# + strandの遺伝子：
#   TSS = gene_start
#
# - strandの遺伝子：
#   TSS = gene_end
#
# TSSはタンパク質翻訳の開始点ではない。
#
# TSS        = RNA転写の開始位置
# start codon = タンパク質翻訳の開始位置
#
#
# 【なぜTSS周辺のERVを調べるのか】
#
# ERV由来LTRは、promoterやenhancerなどの
# cis-regulatory elementとして働く可能性がある。
#
# そのため、AGE関連moduleの遺伝子について、
#
# ・promoter領域にERVが重なるか
# ・TSSから最も近いERVまでの距離
#
# を調べる。
#
#
# 【今回のpromoter定義】
#
# TSSの上流2,000 bpから下流500 bpまでを
# promoter領域として定義する。
#
# + strand：
#   TSSより座標が小さい方向が上流
#
# - strand：
#   TSSより座標が大きい方向が上流
#
# promoter座標は遺伝子のstrandを考慮して作成する。
#
#
# 【ERVのstrandについて】
#
# promoterとERVの重なりでは、
# 遺伝子とERVのstrandが同じかどうかは要求しない。
#
# ERV由来配列は、遺伝子と反対向きに存在していても
# regulatory elementとして影響する可能性があるため、
# findOverlapsではignore.strand = TRUEを使用する。
#
#
# 【重要な注意】
#
# 今回のERV annotationの1行は、
# 完全なprovirusではなくRepeatMasker annotation fragmentである。
#
# したがって、以下では
# 「ERV locus」ではなく「ERV fragment」と表現する。
#
#
# 入力：
# clean/geneModuleTable_with_GenAge_coordinates_mad20k_power6.rds
# clean/RepeatMasker_hg38_ERV_primary_GRCh38.rds
#
# 出力：
# clean/geneModuleTable_with_ERV_overlap_mad20k_power6.rds
# results/geneModuleTable_with_ERV_overlap_mad20k_power6.csv
# clean/promoter_ERV_pairs_mad20k_power6.rds
# results/promoter_ERV_pairs_mad20k_power6.csv
# results/module_ERV_overlap_descriptive_summary.csv
# ============================================================


## Step 1
# 必要なフォルダを確認する

dir.create(
  "clean",
  showWarnings = FALSE
)

dir.create(
  "results",
  showWarnings = FALSE
)


## Step 2
# 必要なpackageを準備する

# GenomicRangesは、
# ゲノム上の区間同士の重なりや距離を計算するpackage

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges")
}

library(GenomicRanges)


## Step 3
# 座標付きWGCNA遺伝子表を読む

geneModuleTable <- readRDS(
  "clean/geneModuleTable_with_GenAge_coordinates_mad20k_power6.rds"
)

# gene数
nrow(geneModuleTable)

# 必要な列を確認
colnames(geneModuleTable)

head(
  geneModuleTable[
    , c(
      "gene",
      "symbol",
      "module",
      "chromosome",
      "gene_start",
      "gene_end",
      "strand",
      "TSS"
    )
  ]
)


## Step 4
# 主解析用ERV annotationを読む

erv_annotation_primary <- readRDS(
  "clean/RepeatMasker_hg38_ERV_primary_GRCh38.rds"
)

# ERV fragment数
nrow(erv_annotation_primary)

# family別件数
table(
  erv_annotation_primary$repFamily
)

# 内容確認
head(erv_annotation_primary)


## Step 5
# 遺伝子とERVの共通染色体を確認する

common_chromosomes <- intersect(
  unique(geneModuleTable$chromosome),
  unique(erv_annotation_primary$chromosome)
)

sort(common_chromosomes)

# ERV位置解析に使用できるgeneを判定する
#
# chrMには今回のERV annotationがないため、
# chrM geneは位置解析の対象外となる

geneModuleTable$ERV_position_eligible <-
  geneModuleTable$chromosome %in% common_chromosomes &
  !is.na(geneModuleTable$TSS)


# 解析可能なgene数
table(
  geneModuleTable$ERV_position_eligible
)


## Step 6
# ERV位置解析に使用するgeneだけを取り出す

eligible_rows <- which(
  geneModuleTable$ERV_position_eligible
)

gene_position_table <- geneModuleTable[
  eligible_rows,
]

# 元のgeneModuleTableの何行目か記録する
gene_position_table$original_row <- eligible_rows

# gene数確認
nrow(gene_position_table)

# chrMなどの除外gene数
nrow(geneModuleTable) -
  nrow(gene_position_table)


## Step 7
# 遺伝子をGRangesへ変換する

# GRangesは、
# chromosome・start・end・strandをまとめて管理する形式

gene_gr <- GRanges(
  seqnames = gene_position_table$chromosome,
  
  ranges = IRanges(
    start = gene_position_table$gene_start,
    end = gene_position_table$gene_end
  ),
  
  strand = gene_position_table$strand
)

# gene情報をGRangesへ追加
mcols(gene_gr)$original_row <-
  gene_position_table$original_row

mcols(gene_gr)$gene <-
  gene_position_table$gene

mcols(gene_gr)$symbol <-
  gene_position_table$symbol

mcols(gene_gr)$module <-
  gene_position_table$module


## Step 8
# TSSを1 bpのGRangesとして作る

# TSSはゲノム上の1つの位置なので、
# startとendに同じTSS座標を指定する

tss_gr <- GRanges(
  seqnames = gene_position_table$chromosome,
  
  ranges = IRanges(
    start = gene_position_table$TSS,
    end = gene_position_table$TSS
  ),
  
  strand = gene_position_table$strand
)

mcols(tss_gr)$original_row <-
  gene_position_table$original_row

mcols(tss_gr)$gene <-
  gene_position_table$gene

mcols(tss_gr)$symbol <-
  gene_position_table$symbol

mcols(tss_gr)$module <-
  gene_position_table$module


## Step 9
# TSSからpromoter領域を作る

# promoter：
# TSS上流2,000 bp ～ TSS下流500 bp
#
# promoters()はstrandを考慮して、
# 上流・下流を自動的に判断する

promoter_gr <- promoters(
  tss_gr,
  upstream = 2000,
  downstream = 500
)

# 染色体の先頭付近でstartが0以下になる場合に備えて、
# promoter startの最小値を1にする

start(promoter_gr) <- pmax(
  start(promoter_gr),
  1
)

# promoter情報を確認
promoter_gr[1:6]

# promoter幅を確認
width(promoter_gr)[1:6]


## Step 10
# ERV annotationをGRangesへ変換する

erv_gr <- GRanges(
  seqnames = erv_annotation_primary$chromosome,
  
  ranges = IRanges(
    start = erv_annotation_primary$repeat_start,
    end = erv_annotation_primary$repeat_end
  ),
  
  strand = erv_annotation_primary$strand
)

# ERV情報をGRangesへ追加
mcols(erv_gr)$erv_fragment_id <-
  erv_annotation_primary$erv_fragment_id

mcols(erv_gr)$repName <-
  erv_annotation_primary$repName

mcols(erv_gr)$repFamily <-
  erv_annotation_primary$repFamily

mcols(erv_gr)$repeat_part <-
  erv_annotation_primary$repeat_part


## Step 11
# promoterとERV fragmentの重なりを調べる

# query  = promoter
# subject = ERV fragment
#
# ignore.strand = TRUE：
# 遺伝子とERVの向きが同じかは要求しない

promoter_hits <- findOverlaps(
  promoter_gr,
  erv_gr,
  ignore.strand = TRUE
)

# promoter–ERV overlapの総数
length(promoter_hits)

# ERVが重なったpromoterの数
length(
  unique(
    queryHits(promoter_hits)
  )
)


## Step 12
# 各geneのpromoterに重なるERV数を計算する

# promoterごとのERV fragment数
promoter_erv_counts <- tabulate(
  queryHits(promoter_hits),
  nbins = length(promoter_gr)
)

# 元の20,000 gene表へ追加する列を作る
#
# 解析対象外のchrM geneなどはNA
# 解析対象でoverlapがないgeneは0

geneModuleTable$promoter_ERV_n <- NA_integer_

geneModuleTable$promoter_ERV_n[
  gene_position_table$original_row
] <- promoter_erv_counts


# promoterに1つ以上ERVがあるか
geneModuleTable$promoter_ERV_any <- NA

geneModuleTable$promoter_ERV_any[
  gene_position_table$original_row
] <- promoter_erv_counts > 0


# 結果確認
table(
  geneModuleTable$promoter_ERV_any,
  useNA = "ifany"
)


## Step 13
# 各promoterに重なるERV familyをまとめる

# overlapしたERVのfamilyを、
# promoter番号ごとに分ける

family_by_promoter <- split(
  erv_annotation_primary$repFamily[
    subjectHits(promoter_hits)
  ],
  queryHits(promoter_hits)
)

# 各promoter用の結果列を作る
promoter_erv_families <- rep(
  NA_character_,
  length(promoter_gr)
)

# 同じfamilyが複数回出た場合は1つにまとめる
promoter_erv_families[
  as.integer(names(family_by_promoter))
] <- vapply(
  family_by_promoter,
  function(x) {
    paste(
      sort(unique(x)),
      collapse = ";"
    )
  },
  character(1)
)

# 元のgene表へ追加
geneModuleTable$promoter_ERV_families <-
  NA_character_

geneModuleTable$promoter_ERV_families[
  gene_position_table$original_row
] <- promoter_erv_families


## Step 14
# promoterとERVの詳細なpair表を作る

promoter_query_index <- queryHits(
  promoter_hits
)

erv_subject_index <- subjectHits(
  promoter_hits
)

promoter_ERV_pairs <- data.frame(
  
  gene = gene_position_table$gene[
    promoter_query_index
  ],
  
  symbol = gene_position_table$symbol[
    promoter_query_index
  ],
  
  module = gene_position_table$module[
    promoter_query_index
  ],
  
  chromosome = gene_position_table$chromosome[
    promoter_query_index
  ],
  
  gene_strand = gene_position_table$strand[
    promoter_query_index
  ],
  
  TSS = gene_position_table$TSS[
    promoter_query_index
  ],
  
  promoter_start = start(promoter_gr)[
    promoter_query_index
  ],
  
  promoter_end = end(promoter_gr)[
    promoter_query_index
  ],
  
  erv_fragment_id =
    erv_annotation_primary$erv_fragment_id[
      erv_subject_index
    ],
  
  repeat_start =
    erv_annotation_primary$repeat_start[
      erv_subject_index
    ],
  
  repeat_end =
    erv_annotation_primary$repeat_end[
      erv_subject_index
    ],
  
  erv_strand =
    erv_annotation_primary$strand[
      erv_subject_index
    ],
  
  repName =
    erv_annotation_primary$repName[
      erv_subject_index
    ],
  
  repFamily =
    erv_annotation_primary$repFamily[
      erv_subject_index
    ],
  
  repeat_part =
    erv_annotation_primary$repeat_part[
      erv_subject_index
    ],
  
  percent_divergence =
    erv_annotation_primary$percent_divergence[
      erv_subject_index
    ],
  
  stringsAsFactors = FALSE
)

# 内容確認
head(promoter_ERV_pairs)

# pair数
nrow(promoter_ERV_pairs)

# 完全に同じpairが重複していないか
sum(
  duplicated(
    promoter_ERV_pairs[
      , c(
        "gene",
        "erv_fragment_id"
      )
    ]
  )
)


## Step 15
# 各TSSに最も近いERV fragmentを調べる

# distanceToNearestは、
# 各TSSから最も近いERV fragmentを探す
#
# select = "arbitrary"：
# 同じ距離に複数ERVがある場合は、
# その中の1つを代表として返す

nearest_hits <- distanceToNearest(
  tss_gr,
  erv_gr,
  ignore.strand = TRUE,
  select = "arbitrary"
)

# nearest結果数
length(nearest_hits)


## Step 16
# nearest ERV情報を元のgene表へ追加する

nearest_gene_index <- queryHits(
  nearest_hits
)

nearest_erv_index <- subjectHits(
  nearest_hits
)

nearest_distance <- mcols(
  nearest_hits
)$distance


# 追加する列を作る
geneModuleTable$nearest_ERV_distance <-
  NA_integer_

geneModuleTable$nearest_ERV_fragment_id <-
  NA_character_

geneModuleTable$nearest_ERV_repName <-
  NA_character_

geneModuleTable$nearest_ERV_family <-
  NA_character_

geneModuleTable$nearest_ERV_part <-
  NA_character_

geneModuleTable$nearest_ERV_start <-
  NA_integer_

geneModuleTable$nearest_ERV_end <-
  NA_integer_

geneModuleTable$nearest_ERV_strand <-
  NA_character_


# 元のgeneModuleTableの行番号
nearest_original_rows <-
  gene_position_table$original_row[
    nearest_gene_index
  ]


# 距離
geneModuleTable$nearest_ERV_distance[
  nearest_original_rows
] <- nearest_distance


# ERV fragment ID
geneModuleTable$nearest_ERV_fragment_id[
  nearest_original_rows
] <- erv_annotation_primary$erv_fragment_id[
  nearest_erv_index
]


# repName
geneModuleTable$nearest_ERV_repName[
  nearest_original_rows
] <- erv_annotation_primary$repName[
  nearest_erv_index
]


# family
geneModuleTable$nearest_ERV_family[
  nearest_original_rows
] <- erv_annotation_primary$repFamily[
  nearest_erv_index
]


# internal / LTR_or_other
geneModuleTable$nearest_ERV_part[
  nearest_original_rows
] <- erv_annotation_primary$repeat_part[
  nearest_erv_index
]


# ERV start
geneModuleTable$nearest_ERV_start[
  nearest_original_rows
] <- erv_annotation_primary$repeat_start[
  nearest_erv_index
]


# ERV end
geneModuleTable$nearest_ERV_end[
  nearest_original_rows
] <- erv_annotation_primary$repeat_end[
  nearest_erv_index
]


# ERV strand
geneModuleTable$nearest_ERV_strand[
  nearest_original_rows
] <- erv_annotation_primary$strand[
  nearest_erv_index
]


## Step 17
# nearest ERVが上流か下流か判定する

# transcription方向を基準として、
#
# upstream：
#   TSSより転写方向の手前
#
# downstream：
#   TSSより転写方向の先
#
# overlaps_TSS：
#   ERV fragmentがTSSそのものに重なる

nearest_tss <- gene_position_table$TSS[
  nearest_gene_index
]

nearest_gene_strand <- gene_position_table$strand[
  nearest_gene_index
]

nearest_erv_start <-
  erv_annotation_primary$repeat_start[
    nearest_erv_index
  ]

nearest_erv_end <-
  erv_annotation_primary$repeat_end[
    nearest_erv_index
  ]


nearest_direction <- ifelse(
  
  # ERVがTSSを含む場合
  nearest_erv_start <= nearest_tss &
    nearest_erv_end >= nearest_tss,
  
  "overlaps_TSS",
  
  # ERVがTSSに重ならない場合
  ifelse(
    
    nearest_gene_strand == "+",
    
    # + strand
    ifelse(
      nearest_erv_end < nearest_tss,
      "upstream",
      "downstream"
    ),
    
    # - strand
    ifelse(
      nearest_erv_start > nearest_tss,
      "upstream",
      "downstream"
    )
  )
)


# 結果を追加
geneModuleTable$nearest_ERV_direction <-
  NA_character_

geneModuleTable$nearest_ERV_direction[
  nearest_original_rows
] <- nearest_direction


# 上流を負、下流を正とした距離も作る
#
# 例：
# upstream 500 bp   → -500
# downstream 500 bp →  500
# TSS overlap       →  0

nearest_signed_distance <- nearest_distance

nearest_signed_distance[
  nearest_direction == "upstream"
] <- -nearest_distance[
  nearest_direction == "upstream"
]

nearest_signed_distance[
  nearest_direction == "overlaps_TSS"
] <- 0


geneModuleTable$nearest_ERV_signed_distance <-
  NA_integer_

geneModuleTable$nearest_ERV_signed_distance[
  nearest_original_rows
] <- nearest_signed_distance


## Step 18
# nearest ERV結果を確認する

head(
  geneModuleTable[
    geneModuleTable$ERV_position_eligible,
    c(
      "gene",
      "symbol",
      "module",
      "chromosome",
      "TSS",
      "nearest_ERV_fragment_id",
      "nearest_ERV_repName",
      "nearest_ERV_family",
      "nearest_ERV_distance",
      "nearest_ERV_direction",
      "nearest_ERV_signed_distance"
    )
  ]
)

# family別件数
table(
  geneModuleTable$nearest_ERV_family,
  useNA = "ifany"
)

# 上流・下流・TSS overlap
table(
  geneModuleTable$nearest_ERV_direction,
  useNA = "ifany"
)

# 距離の要約
summary(
  geneModuleTable$nearest_ERV_distance
)


## Step 19
# moduleごとの記述的summaryを作る
#
# ここではまだ統計的enrichment検定は行わない
# moduleごとの件数と割合を確認するだけ

module_list <- sort(
  unique(
    geneModuleTable$module
  )
)

module_ERV_summary <- do.call(
  rbind,
  lapply(
    module_list,
    function(mod) {
      
      module_rows <- which(
        geneModuleTable$module == mod
      )
      
      eligible_module_rows <- module_rows[
        geneModuleTable$ERV_position_eligible[
          module_rows
        ]
      ]
      
      promoter_positive_n <- sum(
        geneModuleTable$promoter_ERV_any[
          eligible_module_rows
        ],
        na.rm = TRUE
      )
      
      data.frame(
        module = mod,
        
        module_gene_n =
          length(module_rows),
        
        ERV_position_eligible_n =
          length(eligible_module_rows),
        
        promoter_ERV_gene_n =
          promoter_positive_n,
        
        promoter_ERV_percent =
          promoter_positive_n /
          length(eligible_module_rows) *
          100,
        
        median_nearest_ERV_distance =
          median(
            geneModuleTable$nearest_ERV_distance[
              eligible_module_rows
            ],
            na.rm = TRUE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

# 内容確認
module_ERV_summary


## Step 20
# AGE関連候補moduleだけ確認する

target_modules <- c(
  "greenyellow",
  "red",
  "darkgreen",
  "turquoise",
  "midnightblue",
  "blue",
  "pink"
)

target_module_ERV_summary <- module_ERV_summary[
  module_ERV_summary$module %in%
    target_modules,
]

target_module_ERV_summary


## Step 21
# 座標・ERV情報付きWGCNA表を保存する

saveRDS(
  geneModuleTable,
  "clean/geneModuleTable_with_ERV_overlap_mad20k_power6.rds"
)

write.csv(
  geneModuleTable,
  "results/geneModuleTable_with_ERV_overlap_mad20k_power6.csv",
  row.names = FALSE
)


## Step 22
# promoter–ERV pair表を保存する

saveRDS(
  promoter_ERV_pairs,
  "clean/promoter_ERV_pairs_mad20k_power6.rds"
)

write.csv(
  promoter_ERV_pairs,
  "results/promoter_ERV_pairs_mad20k_power6.csv",
  row.names = FALSE
)


## Step 23
# module summaryを保存する

write.csv(
  module_ERV_summary,
  "results/module_ERV_overlap_descriptive_summary.csv",
  row.names = FALSE
)

write.csv(
  target_module_ERV_summary,
  "results/target_module_ERV_overlap_descriptive_summary.csv",
  row.names = FALSE
)


## Step 24
# 最終確認

cat(
  "WGCNA全gene数:",
  nrow(geneModuleTable),
  "\n"
)

cat(
  "ERV位置解析対象gene数:",
  sum(
    geneModuleTable$ERV_position_eligible
  ),
  "\n"
)

cat(
  "promoterにERVがあるgene数:",
  sum(
    geneModuleTable$promoter_ERV_any,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "promoter–ERV pair数:",
  nrow(promoter_ERV_pairs),
  "\n"
)

cat(
  "nearest ERVが付かなかった解析対象gene数:",
  sum(
    geneModuleTable$ERV_position_eligible &
      is.na(
        geneModuleTable$nearest_ERV_fragment_id
      )
  ),
  "\n"
)

target_module_ERV_summary





#確認
table(
  geneModuleTable$ERV_position_eligible
)

length(promoter_hits)

table(
  geneModuleTable$promoter_ERV_any,
  useNA = "ifany"
)

summary(
  geneModuleTable$nearest_ERV_distance
)

table(
  geneModuleTable$nearest_ERV_family,
  useNA = "ifany"
)

target_module_ERV_summary
