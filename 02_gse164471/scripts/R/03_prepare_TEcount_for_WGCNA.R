# 03_prepare_TEcount_for_WGCNA.R
#
# TEcountのGene + strict ERV raw count matrixを、
# WGCNAに使用できる形式へ準備する。
#
# このscriptで行うこと：
#   1. Sample metadataの作成
#   2. Count matrixとmetadataの対応付け
#   3. 低発現featureのfilter
#   4. DESeq2によるvariance stabilizing transformation（VST）
#   5. PCAによるsample QC
#
# 作成したVST matrixは2つの解析で共通して使用する。
#
# WGCNA1：
#   filter後の全Gene + ERVを使う。
#
# WGCNA2：
#   後のDESeq2 age associationで選ばれたfeatureだけを
#   同じVST matrixから抽出して使う。
#
# DESeq2による年齢関連検定そのものは、
# 05_DESeq2_age_association.Rで実行する。



# １. Project directory and packages ---------------------------------------


project_root <- "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"

setwd(project_root)

# DESeq2：
#   raw countの管理、normalisation、VSTに使用する。
#
# GEOquery：
#   GSE164471のGEO sample metadataを取得する。
library(DESeq2)
library(GEOquery)





# 2. Input files ----------------------------------------------------------


# Step 01で保存したGene + strict ERV raw count matrix
#
# 行：
#   19,846 protein-coding genes
#   478 strict ERV subfamilies
#
# 列：
#   feature_id
#   53 SRR samples
count_file <- file.path(
  project_root,
  "results",
  "objects",
  "TEcount_gene_strict_ERV_raw.rds"
)

# Step 01で保存したfeature annotation。
#
# 各行がGeneかERVか、
# Gene symbolやERV familyなどを記録している。
annotation_file <- file.path(
  project_root,
  "results",
  "objects",
  "TEcount_gene_strict_ERV_annotation.rds"
)

#SRRは、SRAで管理されるsequencing runのID
#SRA Run IDとGEO Sample IDを対応させるために使用する。
#
#Run：
#   SRR13388732など
#
#SampleName：
#   GSM5011616など
runinfo_file <- file.path(
  project_root,
  "00_metadata",
  "all_runinfo.csv"
)





# 3. Read input data ------------------------------------------------------


# readRDS()は、saveRDS()で保存したR objectを読み込む関数。
#
# TSVを読み直す方法と違い、
# columnの型などをR objectのまま復元できる。
tecount_raw <- readRDS(count_file)

feature_annotation <- readRDS(annotation_file)

runinfo <- read.csv(
  runinfo_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# Step 01で作成したmatrixの構造を確認する。
dim(tecount_raw)

# 20,324 features × 54 columnsを想定する。
#
# 54列の内訳：
#   feature_id 1列
#   sample 53列






# 4. Download and inspect GEO sample metadata -----------------------------


# GEOには、各sampleの年齢、性別、組織などの
# 生物学的なmetadataが登録されている。
#
# getGEO()でGSE164471のSeries Matrixを取得する。
gse_list <- getGEO(
  "GSE164471",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# getGEO()の結果はlistとして返される。
# 今回使用するExpressionSetは1番目に入っている。
gse_object <- gse_list[[1]]

# pData()はExpressionSetから、
# sample metadataだけをdata frameとして取り出す。
geo_metadata_raw <- pData(gse_object)


# GEO metadataの構造を確認する。
#
# dim()：
#   sample数とmetadata column数
#
# colnames()：
#   どのようなmetadata項目があるか
dim(geo_metadata_raw)
colnames(geo_metadata_raw)


# 今回必要な列の実際の中身を確認する。
#
# title：
#   GEO上のsample名
#
# age:ch1：
#   年齢
#
# gender:ch1：
#   性別
head(
  geo_metadata_raw[, c(
    "geo_accession",
    "title",
    "age:ch1",
    "gender:ch1"
  )]
)





# 5.Extract age and sex metadata ----------------------------------------


# GEO metadataから、今後の解析に必要な列だけを残す。
#
# geo_accession：
#   GSM sample ID
#
# age:ch1：
#   年齢
#
# gender:ch1：
#   性別
geo_metadata <- geo_metadata_raw[, c(
  "geo_accession",
  "age:ch1",
  "gender:ch1"
)]

# 扱いやすい列名へ変更する。
colnames(geo_metadata) <- c(
  "geo_id",
  "age",
  "sex"
)

# GEO accessionが左端のrow nameにも入っているため、
# 通常の1, 2, 3...という行番号へ戻す。
row.names(geo_metadata) <- NULL

# 整理後のmetadataを確認する。
head(geo_metadata)

# 列ごとのdata typeも確認する。
str(geo_metadata)



# 6. Check inconsistent metadata ----

# titleとage/sexの内容が一致しなかったsampleを指定する。
samples_to_check <- c(
  "GSM5011620",
  "GSM5011621"
)

# 各sampleが確認対象かどうかをTRUE/FALSEで記録する。
is_sample_to_check <- (
  geo_metadata_raw$geo_accession %in% samples_to_check
)

# GEOに登録された元のcharacteristicsと、
# 自動整理されたage:ch1・gender:ch1を並べて確認する。
geo_metadata_raw[
  is_sample_to_check,
  c(
    "geo_accession",
    "title",
    "characteristics_ch1",
    "characteristics_ch1.1",
    "characteristics_ch1.2",
    "age:ch1",
    "gender:ch1"
  )
]



#Metadata内でTitleと実際の数値で２つのsampleで不一致
#元論文で男性33人なので男性の数をカウントして
#33人ならmetadataの情報を採用
# GSM5011620
# title      ：AGE37_F
# age:ch1    ：35
# gender:ch1 ：F
# 
# GSM5011621
# title      ：AGE38_M
# age:ch1    ：38
# gender:ch1 ：F

table(geo_metadata$sex)

# F  M 
# 20 33 






# 7. Connect GEO metadata to SRR samples ----------------------------------


# all_runinfo.csvから必要な2列だけを取り出す。
#
# Run：
#   SRR13388732などのsequencing run ID。
#   TEcount matrixのsample列名として使われている。
#
# SampleName：
#   GSM5011616などのGEO sample ID。
#   GEO metadataと接続するために使う。
run_map <- runinfo[, c(
  "Run",
  "SampleName"
)]

# 分かりやすい列名へ変更する。
colnames(run_map) <- c(
  "sample_id",
  "geo_id"
)


# merge()は、2つのdata frameを共通する列で結合する関数。
#
# run_mapとgeo_metadataには、どちらにもgeo_id列がある。
# by = "geo_id"によって、同じGSM IDの行同士を結び付ける。
sample_metadata <- merge(
  run_map,
  geo_metadata,
  by = "geo_id"
)


# 列を解析で使いやすい順番へ並べる。
sample_metadata <- sample_metadata[, c(
  "sample_id",
  "geo_id",
  "age",
  "sex"
)]


# ageは文字列として読み込まれていたので、
# DESeq2で数値変数として使えるinteger型へ変換する。
sample_metadata$age <- as.integer(sample_metadata$age)


# 完成したsample metadataを確認する。
head(sample_metadata)

dim(sample_metadata)




# 8. Align metadata with TEcount samples ----------------------------------


# TEcount matrixの2列目以降が53 samplesのraw count。
# そのcolumn名から、sample IDの並び順を取得する。
sample_ids_in_counts <- colnames(tecount_raw)[-1]

# match()は、左側の各sample IDが、
# sample_metadata$sample_idの何番目にあるかを返す。
#
# 例えば、
# sample_ids_in_countsの1番目がSRR13388732で、
# metadataでも1行目にあれば、結果は1になる。
sample_order <- match(
  sample_ids_in_counts,
  sample_metadata$sample_id
)

# metadataに見つからなかったsampleがある場合はNAになる。
# すべて対応していれば0になる。
n_unmatched_samples <- sum(is.na(sample_order))

cat(
  "Samples not found in metadata:",
  n_unmatched_samples,
  "\n"
)


# TEcount matrixのsample順にmetadataを並べ替える。
sample_metadata <- sample_metadata[sample_order, ]


# DESeq2では、
# count matrixのcolumn名とmetadataのrow nameを対応させる。
#
# そのため、sample_idをmetadataのrow nameに設定する。
row.names(sample_metadata) <- sample_metadata$sample_id


# 並び順が完全に一致したか確認する。
same_sample_order <- identical(
  sample_ids_in_counts,
  row.names(sample_metadata)
)

cat(
  "Count matrix and metadata have the same sample order:",
  same_sample_order,
  "\n"
)






# 9. Create DESeq2 input objects ------------------------------------------


# TEcount dataからsampleのraw count列だけを取り出す。
#
# sample_ids_in_countsには、
# TEcount matrixに含まれる53個のSRR IDが保存されている。
#
# as.matrix()はdata frameをmatrixへ変換する関数。
# DESeq2のcountDataには数値matrixを渡す。
count_matrix <- as.matrix(
  tecount_raw[, sample_ids_in_counts]
)


# 現在のcount_matrixにはfeature_id列が含まれていない。
# その代わり、各行の名前としてfeature IDを設定する。
#
# 例：
# ENSG00000000003.14
# ERV24B_Prim-int:ERV1:LTR
row.names(count_matrix) <- tecount_raw$feature_id


# DESeq2で使用するmetadata列だけを取り出す。
#
# sample IDはすでにrow nameとして設定済みなので、
# ここではageとsexだけを使用する。
col_data <- sample_metadata[, c(
  "age",
  "sex"
)]


# ageを数値として明示する。
#
# DESeq2ではageを連続変数として扱う。
# つまり20代・30代などのgroupではなく、
# 年齢が1歳増えることと発現量の関係をモデル化する。
col_data$age <- as.numeric(col_data$age)


# sexをfactorへ変換する。
#
# factorは、MとFのようなカテゴリー変数を表す型。
#
# levels = c("F", "M")とすると、
# Fが基準カテゴリーになり、Mとの差をモデル内で考慮する。
col_data$sex <- factor(
  col_data$sex,
  levels = c("F", "M")
)






# 10. Create DESeqDataSet -------------------------------------------------


# DESeqDataSetFromMatrix()は、
# raw count matrixとsample metadataをまとめて、
# DESeq2用のobjectを作る関数。
#
# countData：
#   Gene・ERVのraw count matrix
#
# colData：
#   sampleごとのage・sex
#
# design = ~ sex + age：
#   性別の違いを考慮しながら、
#   年齢と発現量の関連を解析するモデル。
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ sex + age
)


# DESeqDataSetの大きさを確認する。
dim(dds)





# 11. Filter low-expression feature ---------------------------------------


# WGCNAではfeature間の発現相関を計算する。
#
# ほとんどのsampleでcountが0に近いfeatureは、
# 発現変動ではなく偶然の0や小さなcountによって
# 相関が決まる可能性がある。
#
# 今回は、少なくとも6 samplesでraw countが10以上の
# GeneまたはERVを残す。
#
# 6 samplesとした理由：
#   全53 samplesの約10%に相当する。
#
# このfilterにはageやsexの情報を使用しないため、
# WGCNA1のglobal networkを年齢方向へ誘導する処理ではない。


# featureを残すために必要なsample数
minimum_samples <- 6


# counts(dds)は、dds内のraw count matrixを取り出す。
#
# counts(dds) >= 10
#   各countが10以上ならTRUE、
#   10未満ならFALSEになる。
#
# rowSums()は、各featureについてTRUEの数を数える。
#
# その数が6以上なら、そのfeatureを残す。
keep_feature <- rowSums(
  counts(dds) >= 10
) >= minimum_samples


# filter前後のfeature数を確認する。
n_features_before <- nrow(dds)
n_features_after <- sum(keep_feature)

cat("Features before filtering:", n_features_before, "\n")
cat("Features after filtering:", n_features_after, "\n")
cat(
  "Features removed:",
  n_features_before - n_features_after,
  "\n"
)


# filterを通過したfeatureだけを残す。
dds_filtered <- dds[keep_feature, ]


# feature annotationも同じ行だけ残す。
#
# ddsとfeature_annotationは同じfeature順であることを
# Step 01で確認済み。
feature_annotation_filtered <- feature_annotation[
  keep_feature,
]


# filter後にGeneとERVがそれぞれ何個残ったか確認する。
table(feature_annotation_filtered$feature_type)



# Filter前：20,324 features
# Filter後：18,291 features
# 除外　　：2,033 features
# 
# Protein-coding genes：17,826
# Strict ERV subfamilies：465




# 12. Apply variance stabilizing transformation ---------------------------


# vst()は、raw countに対して
#
# 1. sample間のlibrary sizeの違いを補正する
# 2. 発現量に依存するvarianceの大きさを安定化する
#
# という変換を行う。
#
# blind = TRUE：
#   ageやsexを使わずに変換する。
#
# 今回のglobal WGCNAでは、
# ageを使ってnetwork構造を事前に誘導しないため、
# blind = TRUEを使用する。
vst_object <- vst(
  dds_filtered,
  blind = TRUE
)


# assay()は、DESeq2 objectの中から
# 実際のVST expression matrixを取り出す関数。
vst_matrix <- assay(vst_object)


# matrixの大きさを確認する。
dim(vst_matrix)




# 13. Check VST matrix ----------------------------------------------------


# missing valueが含まれていないか確認する。
sum(is.na(vst_matrix))


# 各featureについて、53 samples間のMADを計算する。
feature_mad <- apply(
  vst_matrix,
  1,
  mad
)


#just in case, MADが0、つまりsample間で全く変動しないfeature数を確認する。
sum(feature_mad == 0)

#61個がMAD=0と判断されたので除外する
#MAD=0は中央の50%以上のsampleで発現量が同じ
#つまり共発現ネットワークに寄与しない




# 13.Remove zero-MAD features for WGCNA -----------------------------------


keep_mad <- feature_mad > 0

vst_matrix_wgcna <- vst_matrix[keep_mad, ]

feature_annotation_wgcna <- feature_annotation_filtered[keep_mad, ]

#結果
dim(vst_matrix_wgcna)
#18230    53


table(feature_annotation_wgcna$feature_type)
# GENE    TE 
# 17765   465


# vst_matrix
# → 低発現filter後の全18,291 features
# 
# vst_matrix_wgcna
# → zero-MAD 61 featuresを除いたWGCNA入力





# 14. Perform PCA for sample QC -------------------------------------------

library(ggplot2)
# prcomp()ではrows = samples、columns = featuresが必要。
# 現在のmatrixはrows = features、columns = samplesなので、
# t()で転置してからPCAを行う。
pca_result <- prcomp(
  t(vst_matrix_wgcna),
  center = TRUE,
  scale. = FALSE
)

# 各principal componentが説明するvarianceの割合を計算する。
pca_variance <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

# PCA座標とsample metadataをまとめる。
pca_data <- data.frame(
  sample_id = row.names(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  age = sample_metadata[row.names(pca_result$x), "age"],
  sex = sample_metadata[row.names(pca_result$x), "sex"]
)

head(pca_data)
round(pca_variance[1:5], 1)


#結果
# PC1：40.4%
# PC2：23.8%
# 合計：64.2%



#PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = age,
    shape = sex
  )
) +
  geom_point(size = 3) +
  labs(
    x = paste0(
      "PC1 (",
      round(pca_variance[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2], 1),
      "%)"
    ),
    color = "Age",
    shape = "Sex"
  ) +
  theme_classic()

print(pca_plot)


# 15. Save PCA plot ----

# 保存先フォルダが存在しなければ作成する。
dir.create(
  file.path(project_root, "results", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

# PCA plotの保存先を指定する。
pca_file <- file.path(
  project_root,
  "results",
  "figures",
  "TEcount_VST_PCA.pdf"
)

# PCA plotをPDFとして保存する。
ggsave(
  filename = pca_file,
  plot = pca_plot,
  width = 7,
  height = 5
)

file.exists(pca_file)

#4sampleほど離れていたが一旦保留



# 16. Prepare WGCNA input and cluster samples -----------------------------

# WGCNAでは、
# rows = samples
# columns = features
# のmatrixを使用する。
#
# 現在はrows = features、columns = samplesなので、
# t()で転置する。
datExpr <- t(vst_matrix_wgcna)

dim(datExpr)

#sample metadataとの順番も確認
identical(
  row.names(datExpr),
  row.names(sample_metadata)
)


#続いてsample dendrogramを作成
# sample間の発現パターンの距離を計算し、
# 似ているsample同士を階層的にまとめる。
sample_tree <- hclust(
  dist(datExpr),
  method = "average"
)


# sample dendrogramをPDFとして保存する。
sample_tree_file <- file.path(
  project_root,
  "results",
  "figures",
  "TEcount_WGCNA_sample_dendrogram.pdf"
)

pdf(
  sample_tree_file,
  width = 10,
  height = 6
)

plot(
  sample_tree,
  main = "Sample clustering based on VST expression",
  xlab = "",
  sub = "",
  cex = 0.7
)

dev.off()




# 17. Inspect potential outlier sample ------------------------------------


outlier_candidates <- c(
  "SRR13388769",
  "SRR13388751",
  "SRR13388754"
)

# 年齢・性別を確認する。
sample_metadata[
  outlier_candidates,
  c("geo_id", "age", "sex")
]

# raw countのlibrary sizeを確認する。
library_size <- colSums(count_matrix)

data.frame(
  sample_id = outlier_candidates,
  library_size = library_size[outlier_candidates],
  relative_to_median = round(
    library_size[outlier_candidates] / median(library_size),
    2
  )
)



# 17. Save processed TEcount data ----

# PCAとsample dendrogramで外れ候補を確認したが、
# technical QC failureを示す根拠がないため、
# 53 samplesすべてを以降の解析に使用する。


# 低発現filter後のDESeq2 objectを保存する。
# 後のfeature-wise age association解析で使用する。
saveRDS(
  dds_filtered,
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_dds_low_expression_filtered.rds"
  )
)


# WGCNA用expression matrixを保存する。
#
# rows    = 53 samples
# columns = 18,230 Gene・ERV features
saveRDS(
  datExpr,
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_WGCNA_input.rds"
  )
)


# sample metadataを保存する。
saveRDS(
  sample_metadata,
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_sample_metadata.rds"
  )
)


# WGCNAに含まれる18,230 featuresのannotationを保存する。
saveRDS(
  feature_annotation_wgcna,
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_WGCNA_feature_annotation.rds"
  )
)


# DESeq2で使用する18,291 featuresのannotationを保存する。
saveRDS(
  feature_annotation_filtered,
  file.path(
    project_root,
    "results",
    "objects",
    "TEcount_DESeq2_feature_annotation.rds"
  )
)
