# 02_prepare_Telescope_matrix.R
#
# 53 samplesのTelescope reportを読み込み、
# HERV locusごとのraw count matrixを作成する。
#
# TEcountがERV subfamily単位の発現量を表すのに対し、
# Telescopeは個々のHERV locus単位の発現量を推定する。
#
# 最終的に作るmatrix：
#   行 = HERV loci
#   列 = 53 samples
#
# このscriptでは低発現filterやnormalisationは行わない。


# 1. Project directory ----

# この解析projectの最上位directory。
# 01と同じproject rootを使用する。
project_root <- "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

setwd(project_root)

# 現在のworking directoryを確認する。
getwd()


# 2. Find Telescope report files ----

# Telescopeの結果が保存されているdirectory。
telescope_dir <- file.path(
  project_root,
  "05_counts",
  "telescope"
)

# list.files()は、指定したdirectory内にあるfileを探す関数。
#
# path：
#   探索を始めるdirectory
#
# pattern：
#   探したいfile名の特徴
#
# recursive = TRUE：
#   telescope/直下だけでなく、
#   その下にあるsampleごとのsubdirectoryも探索する
#
# full.names = TRUE：
#   file名だけではなく、完全なfile pathを取得する
telescope_report_files <- list.files(
  path = telescope_dir,
  pattern = "telescope_report\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


# 見つかったreport fileの数を確認する。
#
# length()は、objectに含まれる要素数を数える関数。
# 今回はfile pathが53個見つかれば、lengthは53になる。
n_telescope_reports <- length(telescope_report_files)

cat(
  "Telescope report files:",
  n_telescope_reports,
  "\n"
)


# 最初の6個のfile pathを確認する。
head(telescope_report_files)





# 3. Read one Telescope report --------------------------------------------


# telescope_report_filesには53個のfile pathが入っている。
#
# [1]は、その中の1番目の要素を取り出すという意味。
# Consoleに表示される[1]とは異なり、
# codeの中で使う[1]は実際の位置指定である。
first_report_file <- telescope_report_files[1]

cat("First report file:\n")
cat(first_report_file, "\n")


# 最初のTelescope reportをdata frameとして読み込む。
#
# Telescope reportはtab区切りのTSV fileなので、
# read.delim()を使用する。
first_report <- read.delim(
  first_report_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# reportの行数と列数を確認する。
dim(first_report)

# column名を確認する。
colnames(first_report)

# 最初の6行を確認する。
head(first_report)




# 4.Read Telescope report correctly ---------------------------------------


# Telescope reportの1行目には、
# versionやread数などのRunInfoが記録されている。
#
# 実際のtable headerは2行目にあるため、
# skip = 1を指定して最初の1行を読み飛ばす。
#
# skip = 1：
#   fileの先頭1行を無視してから読み込みを開始する。
first_report_clean <- read.delim(
  first_report_file,
  header = TRUE,
  sep = "\t",
  skip = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# 読み直したreportの構造を確認する。
dim(first_report_clean)

colnames(first_report_clean)

head(first_report_clean)




# 5.Extract locus counts from one report ----------------------------------


# Telescope reportには11種類のcolumnがあるが、
# raw count matrixに必要なのは次の2列。
#
#data_frame[行, 列]はdata frameから一部分を取り出す基本的なコマンド
# transcript：HERV locusの名前
#
# final_count：Telescopeが最終的にそのlocusへ割り当てたcount
#
# c()を使って、残したいcolumn名を2つまとめて指定する。
first_locus_counts <- first_report_clean[
  ,
  c("transcript", "final_count")
]


# __no_featureは実在するHERV locusではないため除外する。
#
# != は「等しくない」という比較演算子。
#
# transcriptが"__no_feature"ではない行を残す。
first_locus_counts <- first_locus_counts[
  first_locus_counts$transcript != "__no_feature",
]


# 1行削除した後の行番号を1から振り直す。
#
# row.names()はdata frameの左端に表示される行番号。
# NULLを代入すると、1, 2, 3...へ自動的に戻る。
row.names(first_locus_counts) <- NULL


#結果の確認
# 行数と列数を確認
dim(first_locus_counts)

# 最初の6行を確認
head(first_locus_counts)

# transcript IDの重複数を確認
n_duplicated_loci <- sum(
  duplicated(first_locus_counts$transcript)
)

cat("Duplicated locus IDs:", n_duplicated_loci, "\n")

# final_countに欠損値がないか確認
n_missing_locus_counts <- sum(
  is.na(first_locus_counts$final_count)
)

cat("Missing final counts:", n_missing_locus_counts, "\n")

# final_countが数値として読み込まれているか確認
class(first_locus_counts$final_count)




# 6.Extract sample ID -----------------------------------------------------


# 各Telescope reportは、
# sample IDと同じ名前のdirectory内に保存されている。
#
# 例：
# .../telescope/SRR13388732/SRR13388732-telescope_report.tsv
#
# dirname()：
#   file pathから、fileが入っているdirectoryのpathを取得する。
#
# basename()：
#   pathの最後の部分だけを取得する。
#
# したがって、
# dirname(first_report_file)でsample directoryを取得し、
# basename()で"SRR13388732"だけを取り出す。
sample_id <- basename(dirname(first_report_file))

sample_id




# 7. Rename count column with sample ID -----------------------------------


# 現在のcolumn名を確認する。
colnames(first_locus_counts)

# 2列目のcolumn名を、
# "final_count"からsample IDへ変更する。
#
# colnames(first_locus_counts)
#   すべてのcolumn名を表す
#
# [2]
#   その中の2番目のcolumn名を指定する
#
# <- sample_id
#   2番目のcolumn名に"SRR13388732"を代入する
colnames(first_locus_counts)[2] <- sample_id

# 変更後のcolumn名を確認する。
colnames(first_locus_counts)

# 最初の6行を確認する。
head(first_locus_counts)





# 8. Read all Telescope reports -------------------------------------------


# 各sampleのcount tableを保存するためのlistを作る。
#
# listは、複数のdata frameを1つのobject内に保存できる入れ物。
#
# 今回は53 sampleあるため、
# 53個の空き場所を持つlistを最初に用意する。
telescope_count_list <- vector(
  mode = "list",
  length = length(telescope_report_files)
)


# for loopは、同じ処理を繰り返すための構文。
#
# seq_along(telescope_report_files)は、
# file数に合わせて1, 2, 3, ... 53を作る。
#
# したがって、このloopは53回実行される。
for (i in seq_along(telescope_report_files)) {
  
  # i番目のreport fileを取り出す。
  report_file <- telescope_report_files[i]
  
  # report fileが入っているdirectory名からsample IDを取得する。
  sample_id <- basename(dirname(report_file))
  
  # Telescope reportの1行目はRunInfoなので、
  # skip = 1で読み飛ばす。
  report_data <- read.delim(
    report_file,
    header = TRUE,
    sep = "\t",
    skip = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # locus IDとfinal countの2列だけを残す。
  locus_counts <- report_data[, c("transcript", "final_count")]
  
  # __no_featureは実在するHERV locusではないため除外する。
  locus_counts <- locus_counts[
    locus_counts$transcript != "__no_feature",
  ]
  
  # 行を削除した後、行番号を1から振り直す。
  row.names(locus_counts) <- NULL
  
  # final_count列をsample IDへ変更する。
  colnames(locus_counts)[2] <- sample_id
  
  # 処理したdata frameをlistのi番目へ保存する。
  #
  # [[i]]はlistのi番目の中身を指定する書き方。
  telescope_count_list[[i]] <- locus_counts
}





# 9. Check imported report list -------------------------------------------


# list内に53個のdata frameが保存されたか確認する。
length(telescope_count_list)

# 1番目のdata frameを確認する。
head(telescope_count_list[[1]])

# 最後のdata frameも確認する。
head(telescope_count_list[[53]])




# 10. Check all imported samples ------------------------------------------


# 各sampleについて、
# sample IDとreport内のHERV locus数を保存するvectorを用意する。
#
# character()：
#   文字列を保存するための空のvectorを作る。
#
# integer()：
#   整数を保存するための空のvectorを作る。
sample_ids <- character(length(telescope_count_list))
n_loci_per_sample <- integer(length(telescope_count_list))


# 53個のdata frameを順番に確認する。
for (i in seq_along(telescope_count_list)) {
  
  # i番目のdata frameを取り出す。
  sample_data <- telescope_count_list[[i]]
  
  # 2列目のcolumn名がsample IDになっている。
  sample_ids[i] <- colnames(sample_data)[2]
  
  # そのsampleでreportに出力されたlocus数を記録する。
  n_loci_per_sample[i] <- nrow(sample_data)
}


# sample IDとlocus数を1つのtableにまとめる。
telescope_report_summary <- data.frame(
  sample_id = sample_ids,
  n_loci = n_loci_per_sample
)


# 最初の6 samplesを確認する。
head(telescope_report_summary)


# sample IDが全部で53個あるか確認する。
length(sample_ids)

# unique()は重複を除いた値を返す。
# sample IDがすべて異なれば、こちらも53になる。
length(unique(sample_ids))

# duplicated()は、同じsample IDが2回目以降に現れた場合にTRUEを返す。
# sum()で重複数を数える。
sum(duplicated(sample_ids))


# 各sampleのlocus数の分布を確認する。
#
# summary()は、
# 最小値、第一四分位、中央値、平均値、第三四分位、最大値
# をまとめて表示する。
summary(n_loci_per_sample)

# range()は最小値と最大値だけを表示する。
range(n_loci_per_sample)





# 11. Merge all Telescope samples -----------------------------------------


# 間違って作成したtelescope_rawがあれば削除する。
rm(telescope_raw)

# 1番目のsampleを結合の土台にする。
telescope_raw <- telescope_count_list[[1]]

# 2番目から53番目までをtranscript IDで結合する。
for (i in 2:length(telescope_count_list)) {
  
  telescope_raw <- merge(
    telescope_raw,
    telescope_count_list[[i]],
    by = "transcript",
    all = TRUE
  )
}


dim(telescope_raw)

head(telescope_raw[, 1:6])





# 12. Replace missing locus counts with zero ------------------------------


# merge()で全sampleを結合した際、
# あるsampleのreportに存在しなかったlocusにはNAが入った。
#
# このNAは「測定情報が失われた」という意味ではなく、
# そのsampleでは該当locusのfinal_countが報告されなかったことを示す。
# count matrixとして扱うため、これらを0へ置き換える。


# is.na()は、各値がNAかどうかをTRUE/FALSEで返す。
#
# TRUEは1、FALSEは0として数えられるため、
# sum()でmatrix内のNAの総数を確認できる。
n_missing_before <- sum(is.na(telescope_raw))

cat("Missing values before replacement:", n_missing_before, "\n")


# telescope_raw内でNAになっている場所を指定し、
# その値を0へ置き換える。
telescope_raw[is.na(telescope_raw)] <- 0


# 置き換え後はNAが0個になることを確認する。
n_missing_after <- sum(is.na(telescope_raw))

cat("Missing values after replacement:", n_missing_after, "\n")




# 13. Save output files ---------------------------------------------------


# Telescope raw matrixの保存先を指定する。
matrix_output_dir <- file.path(
  project_root,
  "results",
  "matrices"
)

object_output_dir <- file.path(
  project_root,
  "results",
  "objects"
)

table_output_dir <- file.path(
  project_root,
  "results",
  "tables"
)


# directoryが存在しない場合は作成する。
dir.create(
  matrix_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  object_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 保存するfile pathを指定する。
telescope_matrix_file <- file.path(
  matrix_output_dir,
  "Telescope_locus_raw.tsv"
)

telescope_summary_file <- file.path(
  table_output_dir,
  "Telescope_report_summary.tsv"
)


# TSV形式で保存する。
#
# telescope_raw：
#   14,177 loci × transcript列 + 53 sample列
#
# telescope_report_summary：
#   各sampleのreportに含まれていたlocus数
write.table(
  telescope_raw,
  file = telescope_matrix_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  telescope_report_summary,
  file = telescope_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# R objectとしても保存する。
#
# readRDS()で読み込むと、
# data frameの構造を保ったまま再利用できる。
saveRDS(
  telescope_raw,
  file.path(
    object_output_dir,
    "Telescope_locus_raw.rds"
  )
)




# 14.summary --------------------------------------------------------------


cat("\nStep 02 completed\n")
cat("Telescope reports:", length(telescope_report_files), "\n")
cat("HERV loci:", nrow(telescope_raw), "\n")
cat("Samples:", ncol(telescope_raw) - 1, "\n")
cat("Missing values:", sum(is.na(telescope_raw)), "\n")
cat("Raw matrix saved to:", telescope_matrix_file, "\n")
cat("Report summary saved to:", telescope_summary_file, "\n")
