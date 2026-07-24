# ============================================================
# 08_ERV_annotation_preparation.R
#
# 目的：
# UCSC hg38のRepeatMasker annotationを取得し、
# GRCh38上のLTR・ERV座位情報を準備する
#
#
# 【RepeatMaskerとは】
#
# RepeatMaskerは、DNA配列に含まれる反復配列を検出・注釈する
# 代表的なbioinformatics toolである。
#
# RepeatMaskerが対象とするものには、以下が含まれる。
#
# ・LTR retrotransposons / endogenous retroviruses
# ・LINE
# ・SINE
# ・DNA transposons
# ・simple repeats
# ・low-complexity regions
#
# したがって、RepeatMaskerはERV専用ツールではない。
# RepeatMaskerの全annotationから、
#
#   repClass  = LTR
#   repFamily = ERV1、ERVK、ERVL、ERVL-MaLRなど
#
# を抽出することで、ヒトゲノム上のERV由来配列を取得する。
#
#
# 【今回利用するデータ】
#
# RepeatMaskerを自分で実行するのではなく、
# UCSC Genome Browserがhg38に対して作成した
# RepeatMasker annotation table（rmsk）を利用する。
#
# hg38とGRCh38は、今回の目的では同じassemblyを指す。
#
#
# 【重要な注意点1：1行が完全なprovirusとは限らない】
#
# RepeatMaskerの各行は、
# RepeatMaskerが認識したrepeat fragmentまたはrepeat instanceである。
#
# 例えば、1つのERV provirusが、
#
#   5' LTR
#   internal region
#   3' LTR
#
# の複数行に分かれてannotationされることがある。
#
# したがって、今回作成する表はまず
# 「RepeatMasker annotation fragment単位」の座標表である。
# 完全長ERV locusへの統合は後のStepで検討する。
#
#
# 【重要な注意点2：座標形式】
#
# UCSC rmsk table：
#   start = 0-based
#   end   = half-open
#
# R / GRanges：
#   start = 1-based
#   end   = inclusive
#
# そのため、UCSCのgenoStartには1を足す。
# genoEndはそのまま使用する。
#
#
# 入力：
# UCSC hg38 RepeatMasker rmsk table
#
# 出力：
# clean/RepeatMasker_hg38_LTR_GRCh38.rds
# clean/RepeatMasker_hg38_ERV_GRCh38.rds
# results/RepeatMasker_hg38_LTR_family_counts.csv
# results/RepeatMasker_hg38_ERV_family_counts.csv
# ============================================================


## Step 1
# 必要なフォルダを確認する
#すでにある

## Step 2
# UCSC RepeatMaskerファイルの保存場所を指定する

rmsk_file <- paste0(
  "/rds/homes/t/txk567/gtex_wgcna/raw/",
  "hg38_rmsk.txt.gz"
)

rmsk_url <- paste0(
  "https://hgdownload.soe.ucsc.edu/",
  "goldenPath/hg38/database/rmsk.txt.gz"
)


## Step 3
# UCSC hg38 RepeatMasker annotationをダウンロードする

# すでにファイルが存在する場合は再ダウンロードしない
if (!file.exists(rmsk_file)) {
  
  download.file(
    url = rmsk_url,
    destfile = rmsk_file,
    mode = "wb"
  )
}

# ファイルが存在するか確認
file.exists(rmsk_file)

# ファイルサイズを確認
file.info(rmsk_file)$size


## Step 4
# RepeatMasker rmsk tableを読み込む

# UCSCのdatabase dumpには列名がないため、
# header = FALSEで読み込む

rmsk_raw <- read.delim(
  rmsk_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

# 行数と列数を確認
dim(rmsk_raw)

# 最初の数行を確認
head(rmsk_raw)


## Step 5
# RepeatMasker tableに列名を付ける

# bin列を含む17列の場合
if (ncol(rmsk_raw) == 17) {
  
  colnames(rmsk_raw) <- c(
    "bin",
    "swScore",
    "milliDiv",
    "milliDel",
    "milliIns",
    "genoName",
    "genoStart",
    "genoEnd",
    "genoLeft",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "repStart",
    "repEnd",
    "repLeft",
    "repeat_id"
  )
  
  # bin列を含まない16列の場合
} else if (ncol(rmsk_raw) == 16) {
  
  colnames(rmsk_raw) <- c(
    "swScore",
    "milliDiv",
    "milliDel",
    "milliIns",
    "genoName",
    "genoStart",
    "genoEnd",
    "genoLeft",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "repStart",
    "repEnd",
    "repLeft",
    "repeat_id"
  )
  
} else {
  
  stop(
    paste0(
      "予想外の列数です：",
      ncol(rmsk_raw)
    )
  )
}

# 列名を確認
colnames(rmsk_raw)

# 内容を確認
head(rmsk_raw)


## Step 6
# UCSC座標をR用の1-based座標に変換する

rmsk_raw$chromosome <- rmsk_raw$genoName

# UCSC startは0-basedなので1を足す
rmsk_raw$repeat_start <- rmsk_raw$genoStart + 1

# UCSC endはそのまま使用する
rmsk_raw$repeat_end <- rmsk_raw$genoEnd

# RepeatMasker annotationの長さ
rmsk_raw$repeat_length <-
  rmsk_raw$repeat_end -
  rmsk_raw$repeat_start + 1


# divergenceなどはparts per thousandで記録されているため、
# 10で割ってpercentageにする

rmsk_raw$percent_divergence <-
  rmsk_raw$milliDiv / 10

rmsk_raw$percent_deletion <-
  rmsk_raw$milliDel / 10

rmsk_raw$percent_insertion <-
  rmsk_raw$milliIns / 10


# 内容確認
head(
  rmsk_raw[
    , c(
      "chromosome",
      "repeat_start",
      "repeat_end",
      "strand",
      "repName",
      "repClass",
      "repFamily",
      "percent_divergence"
    )
  ]
)



## Step 7
# 標準染色体だけを残す

# GENCODE v39のgene座標と比較しやすくするため、
# chr1-22、chrX、chrY、chrMだけを使用する
#
# alternative contigやrandom scaffoldは今回は除外する

standard_chromosomes <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY",
  "chrM"
)

rmsk_standard <- rmsk_raw[
  rmsk_raw$chromosome %in% standard_chromosomes,
]

# 全annotation数
nrow(rmsk_raw)

# 標準染色体上のannotation数
nrow(rmsk_standard)

# 染色体ごとのannotation数
table(rmsk_standard$chromosome)



## Step 8
# RepeatMaskerに含まれるrepeat classを確認する

repeat_class_counts <- sort(
  table(rmsk_standard$repClass),
  decreasing = TRUE
)

repeat_class_counts



## Step 9
# LTR classを抽出する

# ERV由来配列はRepeatMaskerでは主に
# repClass == "LTR"として分類される

ltr_annotation <- rmsk_standard[
  rmsk_standard$repClass == "LTR",
]

# LTR annotation数
nrow(ltr_annotation)

# LTR内のfamilyを確認
ltr_family_counts <- sort(
  table(ltr_annotation$repFamily),
  decreasing = TRUE
)

ltr_family_counts



## Step 10
# ERV familyを抽出する

# RepeatMaskerのヒトERVは主に次のfamilyに分類される
#
# ERV1
# ERVK
# ERVL
# ERVL-MaLR
#
# これらはすべてfamily名がERVから始まるため、
# "^ERV"で抽出する

erv_annotation <- ltr_annotation[
  grepl(
    "^ERV",
    ltr_annotation$repFamily
  ),
]

# ERV annotation数
nrow(erv_annotation)

# ERV familyごとのannotation数
erv_family_counts <- sort(
  table(erv_annotation$repFamily),
  decreasing = TRUE
)

erv_family_counts



## Step 11
# ERV familyに含めなかったLTR familyを確認する

non_erv_ltr <- ltr_annotation[
  !grepl(
    "^ERV",
    ltr_annotation$repFamily
  ),
]

non_erv_ltr_family_counts <- sort(
  table(non_erv_ltr$repFamily),
  decreasing = TRUE
)

non_erv_ltr_family_counts



## Step 12
# internal regionとそれ以外を暫定分類する

# repNameが「-int」で終わるものは、
# ERVのinternal regionとしてannotateされていることが多い

erv_annotation$repeat_part <- ifelse(
  grepl(
    "-int$",
    erv_annotation$repName,
    ignore.case = TRUE
  ),
  "internal",
  "LTR_or_other"
)

# internalとLTR_or_otherの数
table(erv_annotation$repeat_part)



## Step 13
# 各ERV annotationにlocus IDを付ける

# 染色体・座標・repeat名を組み合わせて、
# 各annotation fragmentのIDを作る

erv_annotation$erv_fragment_id <- paste(
  erv_annotation$chromosome,
  erv_annotation$repeat_start,
  erv_annotation$repeat_end,
  erv_annotation$repName,
  sep = "_"
)

# IDの重複を確認
sum(duplicated(erv_annotation$erv_fragment_id))



## Step 14
# 解析に使う列だけを整理する

erv_annotation_clean <- erv_annotation[
  , c(
    "erv_fragment_id",
    "chromosome",
    "repeat_start",
    "repeat_end",
    "repeat_length",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "repeat_part",
    "swScore",
    "percent_divergence",
    "percent_deletion",
    "percent_insertion"
  )
]

ltr_annotation_clean <- ltr_annotation[
  , c(
    "chromosome",
    "repeat_start",
    "repeat_end",
    "repeat_length",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "swScore",
    "percent_divergence",
    "percent_deletion",
    "percent_insertion"
  )
]

# 内容確認
head(erv_annotation_clean)



## Step 15
# LTR・ERV familyのsummary tableを作る

ltr_family_summary <- as.data.frame(
  ltr_family_counts,
  stringsAsFactors = FALSE
)

colnames(ltr_family_summary) <- c(
  "repeat_family",
  "annotation_n"
)


erv_family_summary <- as.data.frame(
  erv_family_counts,
  stringsAsFactors = FALSE
)

colnames(erv_family_summary) <- c(
  "ERV_family",
  "annotation_n"
)


non_erv_ltr_family_summary <- as.data.frame(
  non_erv_ltr_family_counts,
  stringsAsFactors = FALSE
)

colnames(non_erv_ltr_family_summary) <- c(
  "repeat_family",
  "annotation_n"
)


# 内容確認
ltr_family_summary
erv_family_summary
non_erv_ltr_family_summary


## Step 16
# WGCNA gene座標とのchromosome表記を確認する

geneModuleTable <- readRDS(
  "clean/geneModuleTable_with_GenAge_coordinates_mad20k_power6.rds"
)

# WGCNA側のchromosome
sort(
  unique(geneModuleTable$chromosome)
)

# ERV側のchromosome
sort(
  unique(erv_annotation_clean$chromosome)
)

# 両方に共通するchromosome
sort(
  intersect(
    unique(geneModuleTable$chromosome),
    unique(erv_annotation_clean$chromosome)
  )
)

# WGCNAにはあるがERV annotationにはないchromosome
setdiff(
  unique(geneModuleTable$chromosome),
  unique(erv_annotation_clean$chromosome)
)



## Step 17
# 整理したannotationを保存する

# 全LTR annotation
saveRDS(
  ltr_annotation_clean,
  "clean/RepeatMasker_hg38_LTR_GRCh38.rds"
)

# ERV familyに限定したannotation
saveRDS(
  erv_annotation_clean,
  "clean/RepeatMasker_hg38_ERV_GRCh38.rds"
)


# family別summary
write.csv(
  ltr_family_summary,
  "results/RepeatMasker_hg38_LTR_family_counts.csv",
  row.names = FALSE
)

write.csv(
  erv_family_summary,
  "results/RepeatMasker_hg38_ERV_family_counts.csv",
  row.names = FALSE
)

write.csv(
  non_erv_ltr_family_summary,
  "results/RepeatMasker_hg38_non_ERV_LTR_family_counts.csv",
  row.names = FALSE
)



## Step 18
# 最終確認

cat(
  "RepeatMasker全annotation数:",
  nrow(rmsk_raw),
  "\n"
)

cat(
  "標準染色体上のannotation数:",
  nrow(rmsk_standard),
  "\n"
)

cat(
  "LTR annotation数:",
  nrow(ltr_annotation_clean),
  "\n"
)

cat(
  "ERV annotation数:",
  nrow(erv_annotation_clean),
  "\n"
)

cat(
  "ERV fragment IDの重複数:",
  sum(duplicated(erv_annotation_clean$erv_fragment_id)),
  "\n"
)

# ERV family別件数
erv_family_summary

# 最初のERV annotation
head(erv_annotation_clean)


##Extra Step
## Step 19
# 主解析用のhigh-confidence ERV annotationを作る
#
# family名の末尾に「?」があるものは、
# RepeatMasker上で分類が不確実なannotationである
#
# 主解析では以下の4 familyだけを使用する
#
# ERV1
# ERVK
# ERVL
# ERVL-MaLR

primary_erv_families <- c(
  "ERV1",
  "ERVK",
  "ERVL",
  "ERVL-MaLR"
)

erv_annotation_primary <- erv_annotation_clean[
  erv_annotation_clean$repFamily %in%
    primary_erv_families,
]


# 主解析用ERV fragment数
nrow(erv_annotation_primary)

# familyごとの数
table(
  erv_annotation_primary$repFamily
)

# IDの重複確認
sum(
  duplicated(
    erv_annotation_primary$erv_fragment_id
  )
)


## Step 20
# 主解析用ERV annotationを保存する

saveRDS(
  erv_annotation_primary,
  "clean/RepeatMasker_hg38_ERV_primary_GRCh38.rds"
)

write.csv(
  as.data.frame(
    table(erv_annotation_primary$repFamily)
  ),
  "results/RepeatMasker_hg38_ERV_primary_family_counts.csv",
  row.names = FALSE
)


