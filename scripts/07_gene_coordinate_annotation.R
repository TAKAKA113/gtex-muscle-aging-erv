# ============================================================
# 07_gene_coordinate_annotation.R
#
# 目的：
# GTEx V10のWGCNA遺伝子に、
# GENCODE v39 / GRCh38のゲノム座標を追加する
#
# GTF：
# /rds/homes/t/txk567/gtex_wgcna/raw/gencode.v39.GRCh38.genes.gtf
#
# 入力：
# clean/geneModuleTable_with_GenAge_mad20k_power6.rds
#
# 出力：
# clean/geneModuleTable_with_GenAge_coordinates_mad20k_power6.rds
# results/geneModuleTable_with_GenAge_coordinates_mad20k_power6.csv
# ============================================================


## Step 1
# 必要なフォルダを作る
#ここはすでにあるので飛ばす

## Step 2
# 必要なpackageを準備する

# rtracklayerはGTFファイルをRで読み込むpackage
# packageがまだ入っていない場合のみインストールする

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  BiocManager::install("rtracklayer")
}

library(rtracklayer)


## Step 3
# GTFファイルの場所を指定する

gtf_file <- paste0(
  "/rds/homes/t/txk567/gtex_wgcna/raw/",
  "gencode.v39.GRCh38.genes.gtf"
)

# ファイルが存在するか確認
file.exists(gtf_file)

# ファイルが見つからなければ処理を停止
if (!file.exists(gtf_file)) {
  stop(
    paste0(
      "GTFファイルが見つかりません：\n",
      gtf_file
    )
  )
}


## Step 4
# GenAge annotation済みのWGCNA結果を読む

geneModuleTable <- readRDS(
  "clean/geneModuleTable_with_GenAge_mad20k_power6.rds"
)

# gene数を確認
nrow(geneModuleTable)

# 列名を確認
colnames(geneModuleTable)

# 内容を確認
head(
  geneModuleTable[
    , c(
      "gene",
      "module",
      "symbol",
      "is_GenAge"
    )
  ]
)

# 期待されるgene数
# 20000


## Step 5
# GENCODE v39 GTFを読み込む

gtf_v39 <- import(gtf_file)

# 読み込んだGTFを確認
gtf_v39

# featureの種類を確認
table(gtf_v39$type)

# GTFにはgene、transcript、exonなどが含まれる
# 今回必要なのはgene全体の座標なので、
# typeがgeneの行だけを残す

genes_v39 <- gtf_v39[
  gtf_v39$type == "gene"
]

# GTFに含まれるgene数を確認
length(genes_v39)


## Step 6
# GTFからgene座標表を作る

gene_coordinates <- data.frame(
  
  # version付きEnsembl gene ID
  gene = genes_v39$gene_id,
  
  # gene symbol
  symbol_gtf = genes_v39$gene_name,
  
  # protein_coding、lncRNAなど
  gene_type = genes_v39$gene_type,
  
  # 染色体
  chromosome = as.character(seqnames(genes_v39)),
  
  # geneの開始位置
  gene_start = start(genes_v39),
  
  # geneの終了位置
  gene_end = end(genes_v39),
  
  # DNA strand
  strand = as.character(strand(genes_v39)),
  
  stringsAsFactors = FALSE
)


## Step 7
# Ensembl IDのversion番号を外す

# WGCNA側
# ENSG00000168000.16
# ↓
# ENSG00000168000

geneModuleTable$gene_clean <- sub(
  "\\..*$",
  "",
  geneModuleTable$gene
)

# GTF側
gene_coordinates$gene_clean <- sub(
  "\\..*$",
  "",
  gene_coordinates$gene
)


## Step 8
# TSSとgene lengthを計算する

# TSS = Transcription Start Site
#
# + strand：
# gene_startがTSS
#
# - strand：
# gene_endがTSS

gene_coordinates$TSS <- ifelse(
  gene_coordinates$strand == "+",
  gene_coordinates$gene_start,
  gene_coordinates$gene_end
)

# gene全体の長さを計算
gene_coordinates$gene_length <-
  gene_coordinates$gene_end -
  gene_coordinates$gene_start + 1


## Step 9
# GTFから作った座標表を確認する

head(gene_coordinates)

# 同じgene_cleanが重複していないか確認
sum(duplicated(gene_coordinates$gene_clean))

# 0なら重複なし

## Step 10
# 完全なEnsembl gene IDに重複がないか確認する
#
# gene_cleanでは45個重複したが、
# version番号を含む完全なgene IDでは重複しないはず

sum(duplicated(gene_coordinates$gene))

# 期待される結果
# 0

## Step 11
# WGCNA遺伝子とGENCODE v39を照合する
#
# gene_cleanではなく、
# ENSG00000134184.13のような完全なgene IDを使う

match_index <- match(
  geneModuleTable$gene,
  gene_coordinates$gene
)

# 照合できたgene数
sum(!is.na(match_index))

# 照合できなかったgene数
sum(is.na(match_index))


## Step 12
# WGCNA表にゲノム座標を追加する

geneModuleTable$symbol_gtf <-
  gene_coordinates$symbol_gtf[match_index]

geneModuleTable$gene_type <-
  gene_coordinates$gene_type[match_index]

geneModuleTable$chromosome <-
  gene_coordinates$chromosome[match_index]

geneModuleTable$gene_start <-
  gene_coordinates$gene_start[match_index]

geneModuleTable$gene_end <-
  gene_coordinates$gene_end[match_index]

geneModuleTable$strand <-
  gene_coordinates$strand[match_index]

geneModuleTable$TSS <-
  gene_coordinates$TSS[match_index]

geneModuleTable$gene_length <-
  gene_coordinates$gene_length[match_index]


## Step 13
# 座標が正しく付いたか確認する

cat(
  "WGCNA gene数:",
  nrow(geneModuleTable),
  "\n"
)

cat(
  "座標が付いたgene数:",
  sum(!is.na(geneModuleTable$chromosome)),
  "\n"
)

cat(
  "座標が付かなかったgene数:",
  sum(is.na(geneModuleTable$chromosome)),
  "\n"
)

cat(
  "座標付与率:",
  mean(!is.na(geneModuleTable$chromosome)) * 100,
  "%\n"
)


## Step 14
# 座標付きWGCNA表を確認する

head(
  geneModuleTable[
    , c(
      "gene",
      "symbol",
      "symbol_gtf",
      "module",
      "chromosome",
      "gene_start",
      "gene_end",
      "strand",
      "TSS",
      "gene_length",
      "gene_type"
    )
  ]
)


## Step 15
# 以前付けたsymbolとGTFのsymbolが一致するか確認する

symbol_match <- geneModuleTable$symbol ==
  geneModuleTable$symbol_gtf

table(
  symbol_match,
  useNA = "ifany"
)


## Step 16
# 座標が付かなかった遺伝子を確認する

unmatched_genes <- geneModuleTable[
  is.na(geneModuleTable$chromosome),
  c(
    "gene",
    "gene_clean",
    "symbol",
    "module"
  )
]

nrow(unmatched_genes)

unmatched_genes


## Step 17
# 座標付きWGCNA表を保存する

saveRDS(
  geneModuleTable,
  "clean/geneModuleTable_with_GenAge_coordinates_mad20k_power6.rds"
)

write.csv(
  geneModuleTable,
  "results/geneModuleTable_with_GenAge_coordinates_mad20k_power6.csv",
  row.names = FALSE
)


## Step 18
# GENCODE v39のgene座標表を保存する

saveRDS(
  gene_coordinates,
  "clean/gencode_v39_GRCh38_gene_coordinates.rds"
)

write.csv(
  gene_coordinates,
  "results/gencode_v39_GRCh38_gene_coordinates.csv",
  row.names = FALSE
)

