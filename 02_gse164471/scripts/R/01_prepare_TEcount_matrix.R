# 01_prepare_TEcount_matrix.R
#
# TEcountで得られたraw count matrixを読み込み、
# protein-coding genesとERV subfamiliesを抽出するための準備
#
# TEcountの出力には通常のgeneとTE subfamilyが同じmatrix内に含まれる。
# そのため、feature annotationを使って各行が何を表すか確認し、
# 後の段階でprotein-coding genesとstrict ERVsを選択する。
#
# このscriptでは低発現filterやnormalisationはまだ行わない。
# まずraw countとannotationが正しく対応しているかを確認する。



# 1.Project directory -------------------------------------------------------
# project_rootは、この解析projectの最上位directoryを表す。
# 以降はこの場所を基準にしてinputやoutputのpathを指定する。
#
# working directoryをproject rootに固定すると、
# 「現在どのfolderにいるか」によってfile pathが変わる問題を避けられる。
project_root <- "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

setwd(project_root)

# getwd()は現在のworking directoryを返す。
# project_rootと同じpathが表示されればよい。
getwd()



# 2.Input files -----------------------------------------------------------

# TEcount raw count matrix
#
# 行：
#   geneまたはTE subfamily
#
# 列：
#   feature_idと53 samples
#
# 値：
#   TEcountが数えたraw read count
count_file <- file.path(
  project_root,
  "05_counts",
  "tecount",
  "gene_te_counts_raw.tsv"
)

# Feature annotation
#
# count matrixの各featureが、
# geneなのかTEなのかを判定するために使用する。
#
# geneの場合はgene biotype、
# TEの場合はclassやfamilyなどの情報が含まれている。
annotation_file <- file.path(
  project_root,
  "05_counts",
  "tecount",
  "gene_te_feature_annotation.tsv"
)

# ERV subfamily whitelist
#
# TEcount matrixにはERV以外のtransposable elementも含まれる。
# この一覧を使用して、ERVとして扱うsubfamilyだけを選択する。
erv_whitelist_file <- file.path(
  project_root,
  "03_reference",
  "te_annotation",
  "erv_subfamily_whitelist.tsv"
)


#
# file.exists()は、指定したfileが実際に存在するか確認する関数。
#
# TRUE：
#   fileが存在する
#
# FALSE：
#   pathまたはfile名が間違っている
file.exists(count_file)
file.exists(annotation_file)
file.exists(erv_whitelist_file)



# 3.Check input files -----------------------------------------------------

# 3つのinput fileを1つのcharacter vectorにまとめる。
# c()は複数の値を結合する関数。
input_files <- c(
  count_file,
  annotation_file,
  erv_whitelist_file
)

# file.exists(input_files)を実行すると、
# 3つのfileそれぞれについてTRUE/FALSEが返る。
#
# all()は、すべてがTRUEの場合だけTRUEを返す。
#
# ! はTRUEとFALSEを反転する記号なので、
# 1つでもfileが存在しない場合にstop()が実行される。
if (!all(file.exists(input_files))) {
  stop("Input fileが見つかりません。file pathを確認してください。")
}




# 4.Read input data -------------------------------------------------------


# read.delim()はtab区切りのtext fileをdata frameとして読み込む。
#
# TEcountのTSV fileはtab区切りなのでread.delim()を使用する。
#
# check.names = FALSE：
#   sample IDやcolumn名をRが自動的に変更しないようにする。
#
# stringsAsFactors = FALSE：
#   文字列をfactorへ自動変換せず、characterとして保持する。
tecount_raw <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

feature_annotation <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

erv_whitelist <- read.delim(
  erv_whitelist_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)




# 5.Check TEcount matrix --------------------------------------------------


# dim()はdata frameの「行数」と「列数」を返す。
#
# 今回はおよそ次の構造を想定している。
#
# 行数：
#   約59,608 features
#
# 列数：
#   feature ID 1列 + 53 sample列
dim(tecount_raw)


# colnames()で全column名を確認する。
# 最初のcolumnがfeature IDで、
# 残りがSRR sample IDになっているか確認する。
colnames(tecount_raw)

# head()は最初の6行を表示する。
#
# matrix全体を表示すると非常に大きいため、
# dataの内容や形式を確認するときはhead()を使用する。
head(tecount_raw)


# TEcount matrixの行数をfeature数として保存する。
n_features <- nrow(tecount_raw)

# 1列目がfeature IDであるため、
# 全column数から1を引いた値がsample数になる。
n_samples <- ncol(tecount_raw) - 1

cat("TEcount features:", n_features, "\n")
cat("Samples:", n_samples, "\n")

# 今回のdatasetは53 samplesで構成される。
# 53でなければ、sampleが欠けているか、
# 想定外のcolumnが追加されている可能性がある。
if (n_samples != 53) {
  warning("TEcount matrixのsample数が53ではありません。")
}




# 6.Check feature annotation ----------------------------------------------


# annotationの行数、column名、最初の6行を確認する。
#
# この段階ではcolumn名を決めつけない。
# 実際のcolumn構造を確認してから、
# protein-coding geneとERVを抽出するcodeを書く。
dim(feature_annotation)
colnames(feature_annotation)
head(feature_annotation)

cat("Annotation rows:", nrow(feature_annotation), "\n")

# count matrixとannotationの行数が同じでも、
# featureの並び順まで同じとは限らない。
#
# 最終的にはfeature IDを使って対応付ける必要がある。
# ここでは最初の基本確認として行数を比較する。
if (nrow(tecount_raw) != nrow(feature_annotation)) {
  warning("TEcount matrixとannotationの行数が一致していません。")
}




# 7.Check ERV whitelist ---------------------------------------------------


# whitelistの行数とcolumn構造を確認する。
#
# このfileにはstrict ERVとbroad ERVの情報が
# どのような形式で記録されているかを確認する必要がある。
dim(erv_whitelist)
colnames(erv_whitelist)
head(erv_whitelist)

cat("ERV whitelist rows:", nrow(erv_whitelist), "\n")
#ERV whitelist rows: 564




# 8. Check feature IDs ----------------------------------------------------


# TEcount matrixの最初のcolumn名を取得する。
#
# 現時点ではcolumn名を直接"feature_id"と決めつけず、
# 実際の最初のcolumn名を使用する。
feature_id_column <- colnames(tecount_raw)[1]

cat("Feature ID column:", feature_id_column, "\n")

# duplicated()は、以前に同じ値が出現した要素をTRUEにする。
#
# sum()でTRUEの数を数えることで、
# 重複したfeature IDの数を確認できる。
n_duplicated_features <- sum(
  duplicated(tecount_raw[[feature_id_column]])
)

cat("Duplicated feature IDs:", n_duplicated_features, "\n")

# feature IDが重複している場合、
# 同じfeatureが複数行に存在することになる。
#
# そのままWGCNAへ進めるとcolumn名が一意でなくなるため、
# 重複があれば原因を調べる必要がある。
if (n_duplicated_features > 0) {
  warning("重複したfeature IDがあります。")
}


# 9. Check count columns ----

# 1列目以外をsample columnとして取得する。
sample_columns <- colnames(tecount_raw)[-1]

# sapply()を使い、各sample columnにclass()を適用する。
#
# raw countは計数値なので、
# sample columnはintegerまたはnumericである必要がある。
sample_column_classes <- sapply(
  tecount_raw[, sample_columns, drop = FALSE],
  class
)

# 各classが何列あるか集計する。
table(sample_column_classes)

# numericまたはinteger以外のcolumnがないか確認する。
valid_count_classes <- sample_column_classes %in% c(
  "integer",
  "numeric"
)

if (!all(valid_count_classes)) {
  warning("数値として読み込まれていないsample columnがあります。")
}


# 10. Check missing values ----

# is.na()は欠損値をTRUEとして返す。
# sum()でmatrix全体の欠損値数を数える。
#
# raw count matrixにNAがあると、
# normalizationやWGCNAで問題になる。
n_missing_counts <- sum(
  is.na(tecount_raw[, sample_columns, drop = FALSE])
)

cat("Missing count values:", n_missing_counts, "\n")

if (n_missing_counts > 0) {
  warning("TEcount matrixにmissing valueがあります。")
}


# 11. Initial summary ----

# このsummaryは解析結果ではなく、
# input dataが想定した構造になっているかを確認するためのもの。
cat("\nTEcount input check completed\n")
cat("Features:", n_features, "\n")
cat("Samples:", n_samples, "\n")
cat("Annotation rows:", nrow(feature_annotation), "\n")
cat("Whitelist rows:", nrow(erv_whitelist), "\n")
cat("Duplicated feature IDs:", n_duplicated_features, "\n")
cat("Missing count values:", n_missing_counts, "\n")

#簡単なまとめ（自分用）
# TEcount matrixには59,608個のfeatureがある
# 予定どおり53 sampleが含まれている
# annotationも59,608行なので、count matrixと同じ数のfeatureが記録されている
# ERV whitelistには564 subfamiliesが含まれている
# feature IDの重複はない
# count matrixに欠損値NAはない


colnames(feature_annotation)
head(feature_annotation)

colnames(erv_whitelist)
head(erv_whitelist)




# 12.Check feature order --------------------------------------------------


# TEcount matrixとfeature annotationは、どちらも59,608行だった。
# ただし、行数が同じでもfeatureの並び順が同じとは限らない。
#
# 例えば、
#
# count matrix 1行目：Gene A
# annotation   1行目：Gene B
#
# となっていると、Gene AにGene Bのannotationを付けることになる。
# そのため、feature_idの内容と順番が完全に一致しているか確認する。


# identical()は、2つのobjectが完全に同じか確認する関数。
#
# 内容、順番、長さがすべて同じならTRUE、
# どこか1か所でも違えばFALSEを返す。
same_feature_order <- identical(
  tecount_raw$feature_id,
  feature_annotation$feature_id
)

cat(
  "Count matrix and annotation have the same feature order:",
  same_feature_order,
  "\n"
)






# 13. Identify protein-coding genes ---------------------------------------


# 今回のWGCNAでは、gene側はprotein-coding geneを使用する。
#
# protein-coding geneを選ぶ条件は次の2つ。
#
# 1. feature_type が "GENE"
# 2. gene_biotype が "protein_coding"
#
# 両方を満たす行だけをTRUEにする。


# まず、各行がGENEか確認する。
#
# == は「左右が同じか」を確認する比較演算子。
#
# feature_typeがGENEならTRUE、
# GENEでなければFALSEになる。
is_gene <- feature_annotation$feature_type == "GENE"

# TRUEとFALSEの数を確認する。
table(is_gene)


# 次に、各行がprotein-coding biotypeか確認する。
is_protein_coding_biotype <- (
  feature_annotation$gene_biotype == "protein_coding"
)

table(is_protein_coding_biotype)


# & は「両方の条件を満たす」という意味。
#
# GENEであり、
# かつ
# gene_biotypeがprotein_codingである行だけTRUEになる。
is_protein_coding_gene <- (
  is_gene &
    is_protein_coding_biotype
)

# TRUEになった行数を数える。
#
# Rでは論理値を数値として扱う場合、
# TRUEは1、FALSEは0として計算される。
#
# そのためsum()を使うとTRUEの数を数えられる。
n_protein_coding_genes <- sum(is_protein_coding_gene)

cat(
  "Protein-coding genes:",
  n_protein_coding_genes,
  "\n"
)
#Protein-coding genes: 19846となった





# 14.Identify strict ERV subfamilies -------------------------------------


# 今回のprimary analysisでは、次の条件を満たすTE subfamilyを
# strict ERVとして使用する。
#
# 1. te_class が "LTR"
# 2. te_family が "ERV1", "ERVK", "ERVL" のいずれか
#
# LTR配列には転写調節配列が含まれ、
# promoterやenhancerとして機能する可能性がある。
# そのため、ERVと近傍geneの発現制御を調べる今回の研究に関連性が高い。
#
# ただし、このfilterを通過したすべてのLTRが
# 実際にpromoterまたはenhancerとして機能している、
# という意味ではない。
#
# ERVL-MaLRはより広いERV関連LTRの定義には含められるが、
# primary analysisでは除外し、後のbroad analysis用に残している。


# include_primaryがTRUEの行だけを残す。
#
# data_frame[行の条件, columnの条件]
#
# コンマの左側：
#   残す行
#
# コンマの右側：
#   残すcolumn
#
# 今回は右側を空欄にしているため、
# すべてのcolumnを残す。
strict_erv_whitelist <- erv_whitelist[
  erv_whitelist$include_primary == TRUE,
]


# strict ERVとして選ばれた行数を確認する。
nrow(strict_erv_whitelist)


# ERV familyごとの数も確認する。
#table()は、各値が何回出てくるかを数える関数
table(strict_erv_whitelist$te_family)


# strict ERVのfeature IDだけを取り出す。
strict_erv_ids <- strict_erv_whitelist$feature_id


# whitelist内のstrict ERVが、
# feature annotationにも存在するか確認する。
#
# %in% は「左側の値が、右側の一覧に含まれているか」
# を確認する記号。
#
# 含まれていればTRUE、
# 含まれていなければFALSEになる。
strict_ids_found <- (
  strict_erv_ids %in% feature_annotation$feature_id
)


# whitelistに登録されたstrict ERVの総数
n_strict_erv_ids <- length(strict_erv_ids)

# feature annotation内で見つかったstrict ERVの数
n_strict_erv_found <- sum(strict_ids_found)

cat("Strict ERVs in whitelist:", n_strict_erv_ids, "\n")
cat("Strict ERVs found in annotation:", n_strict_erv_found, "\n")


# 次に、feature annotationの59,608行それぞれについて、
# strict ERVに該当するかをTRUE/FALSEで記録する。
is_strict_erv <- (
  feature_annotation$feature_id %in% strict_erv_ids
)

# TRUEになった行数を数える。
n_strict_erv <- sum(is_strict_erv)

cat("Strict ERV subfamilies:", n_strict_erv, "\n")






# 15.Build Gene + strict ERV raw matrix -----------------------------------


# これまでに、59,608 featuresの各行について、
#
# is_protein_coding_gene
#   protein-coding geneならTRUE
#
# is_strict_erv
#   strict ERVならTRUE
#
# という2つのTRUE/FALSE vectorを作成


# | は「または」を意味する。
#
# protein-coding geneである
# または
# strict ERVである
#
# どちらか一方でもTRUEなら、selected_featureはTRUEになる。
is_selected_feature <- (
  is_protein_coding_gene |
    is_strict_erv
)


# 選択されるfeature数を確認する。
#
# TRUEは1、FALSEは0として数えられるため、
# sum()でTRUEの総数を取得できる。
n_selected_features <- sum(is_selected_feature)

cat(
  "Selected Gene + strict ERV features:",
  n_selected_features,
  "\n"
)

# is_selected_featureがTRUEの行だけをTEcount matrixから取り出す。
#
# data_frame[行, 列]
#
# コンマの左側：
#   残す行
#
# コンマの右側：
#   残す列
#
# 右側が空欄なので、feature_idと53 sample列をすべて残す。
tecount_gene_erv_raw <- tecount_raw[
  is_selected_feature,
]


# 同じ行をfeature annotationからも取り出す。
#
# count matrixとannotationは同じ順番であることを
# Step 12で確認済みなので、同じTRUE/FALSE vectorを使用できる。
gene_erv_annotation <- feature_annotation[
  is_selected_feature,
]





# 16.Check selected matrix ------------------------------------------------


# 行数と列数を確認する。
dim(tecount_gene_erv_raw)
dim(gene_erv_annotation)


# 抽出後もcount matrixとannotationのfeature順が一致するか確認する。
selected_features_match <- identical(
  tecount_gene_erv_raw$feature_id,
  gene_erv_annotation$feature_id
)

cat(
  "Selected matrix and annotation have the same feature order:",
  selected_features_match,
  "\n"
)


# 抽出後のfeature種類を数える。
table(gene_erv_annotation$feature_type)


# 最終確認
cat("\nGene + strict ERV matrix created\n")
cat("Protein-coding genes:", sum(is_protein_coding_gene), "\n")
cat("Strict ERV subfamilies:", sum(is_strict_erv), "\n")
cat("Total selected features:", nrow(tecount_gene_erv_raw), "\n")
cat("Samples:", ncol(tecount_gene_erv_raw) - 1, "\n")






# 17.Save output files ----------------------------------------------------


# resultsの保存先を指定する。
#
# matrices：
#   発現matrixなど、行列形式のdataを保存する。
#
# tables：
#   feature annotationなど、説明情報を含む表を保存する。
#
# objects：
#   Rで作成したobjectを、そのままRDS形式で保存する。
matrix_output_dir <- file.path(
  project_root,
  "results",
  "matrices"
)

table_output_dir <- file.path(
  project_root,
  "results",
  "tables"
)

object_output_dir <- file.path(
  project_root,
  "results",
  "objects"
)


# dir.create()はdirectoryを作成する関数。
#
# recursive = TRUE：
#   results/が存在しない場合も、
#   results/matrices/までまとめて作成する。
#
# showWarnings = FALSE：
#   directoryがすでに存在していてもwarningを表示しない。
dir.create(
  matrix_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  object_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# 保存するfile pathを指定する。
count_output_file <- file.path(
  matrix_output_dir,
  "TEcount_gene_strict_ERV_raw.tsv"
)

annotation_output_file <- file.path(
  table_output_dir,
  "TEcount_gene_strict_ERV_annotation.tsv"
)


# write.table()でdata frameをtext fileとして保存する。
#
# sep = "\t"：
#   tab区切りのTSV形式で保存する。
#
# quote = FALSE：
#   characterを引用符 " " で囲まない。
#
# row.names = FALSE：
#   Rが自動的に作る行番号をfileに保存しない。
#
# row.names = FALSEにしないと、
# feature_idとは別に不要な1, 2, 3...の列が追加される。
write.table(
  tecount_gene_erv_raw,
  file = count_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_erv_annotation,
  file = annotation_output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# saveRDS()は、R objectをそのまま保存する関数。
#
# TSV：
#   人が開いて内容を確認しやすい。
#
# RDS：
#   data typeやcolumn構造を保ったまま、
#   後のR scriptで素早く読み込める。
saveRDS(
  tecount_gene_erv_raw,
  file.path(
    object_output_dir,
    "TEcount_gene_strict_ERV_raw.rds"
  )
)

saveRDS(
  gene_erv_annotation,
  file.path(
    object_output_dir,
    "TEcount_gene_strict_ERV_annotation.rds"
  )
)




# 18.summary --------------------------------------------------------------


cat("\nStep 01 completed\n")
cat("Protein-coding genes:", n_protein_coding_genes, "\n")
cat("Strict ERV subfamilies:", n_strict_erv, "\n")
cat("Total features:", nrow(tecount_gene_erv_raw), "\n")
cat("Samples:", ncol(tecount_gene_erv_raw) - 1, "\n")
cat("Count matrix saved to:", count_output_file, "\n")
cat("Annotation saved to:", annotation_output_file, "\n")


# 表のまとめ
# Expression matrix
# 20,324 features × 53 samples
# 
# Feature annotation
# 20,324 features × 10 columns
# 
# Sample metadata
# 53 samples × 年齢・性別など



# 流れの確認(自分用)
# TEcount raw matrixを読み込む
# ↓
# 入力構造を確認
# 59,608 features × 53 samples
# ↓
# annotationとの対応を確認
# feature順序：完全一致
# ↓
# protein-coding genesを特定
# 19,846 genes
# ↓
# LTR classかつERV1 / ERVK / ERVLをstrict ERVとして特定
# 478 ERV subfamilies
# ↓
# Gene + strict ERV raw matrixを作成
# 20,324 features × 53 samples
# ↓
# matrix・annotationを保存
