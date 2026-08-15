# 07_prepare_Telescope_loci.R
#
# 目的：
#   Telescopeで定量したindividual ERV lociについて、
#   downstream age association解析に使用可能な
#   locus-level expression matrixを作成する。
#
# このstepでは：
#   1. Telescope locus matrixの構造を確認する
#   2. 53 samplesが正しく含まれていることを確認する
#   3. 重複sample column (.x / .y) がないことを確認する
#   4. locus annotationを確認する
#   5. 低発現locusを除外する
#   6. TEcount ERV subfamilyとTelescope locusを対応付ける


# 1. Set project ----

setwd(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
)


# 2. Read Telescope locus matrix ----

telescope_raw <- readRDS(
  "results/objects/Telescope_locus_raw.rds"
)



# 3. Check Telescope matrix structure ----

cat("\n===== Telescope raw matrix =====\n")

cat("Class:\n")
print(class(telescope_raw))

cat("\nDimensions:\n")
print(dim(telescope_raw))

cat("\nColumn names (first 10):\n")
print(head(colnames(telescope_raw), 10))

cat("\nRow names (first 10):\n")
print(head(rownames(telescope_raw), 10))



# Check duplicated sample suffixes created by previous merge.

xy_columns <- grep(
  "\\.[xy]$",
  colnames(telescope_raw),
  value = TRUE
)

cat(
  "\nNumber of .x / .y columns:",
  length(xy_columns),
  "\n"
)

if (length(xy_columns) > 0) {
  print(xy_columns)
}


#確認まとめ
dim(telescope_raw)

head(colnames(telescope_raw), 10)

head(rownames(telescope_raw), 10)

length(
  grep(
    "\\.[xy]$",
    colnames(telescope_raw)
  )
)


# 4. Separate locus IDs and count matrix ----

# "transcript" columnにはTelescope locus IDが入っている。
# 残り53 columnsが各sampleのlocus-level counts。

telescope_locus_id <- telescope_raw$transcript

telescope_counts <- telescope_raw[
  ,
  colnames(telescope_raw) != "transcript"
]


# Locus IDをrow nameとして設定する。
rownames(telescope_counts) <- telescope_locus_id




# 5. Check locus IDs ----

cat("\n===== Telescope locus IDs =====\n")

cat("Number of loci:",
    length(telescope_locus_id), "\n")

cat("\nFirst 10 locus IDs:\n")
print(head(telescope_locus_id, 10))

cat("\nDuplicated locus IDs:\n")
print(sum(duplicated(telescope_locus_id)))

cat("\nMissing locus IDs:\n")
print(sum(is.na(telescope_locus_id) |
            telescope_locus_id == ""))


# 6. Check Telescope counts ----

cat("\n===== Telescope count matrix =====\n")

cat("Dimensions:\n")
print(dim(telescope_counts))

cat("\nAll sample columns numeric:\n")
print(
  all(
    vapply(
      telescope_counts,
      is.numeric,
      logical(1)
    )
  )
)

cat("\nNumber of NA values:\n")
print(sum(is.na(telescope_counts)))

cat("\nNumber of negative values:\n")
print(sum(as.matrix(telescope_counts) < 0))

cat("\nCount summary:\n")
print(summary(as.vector(as.matrix(telescope_counts))))



countまとめ
head(telescope_locus_id, 10)
sum(duplicated(telescope_locus_id))

dim(telescope_counts)
sum(is.na(telescope_counts))
sum(as.matrix(telescope_counts) < 0)

summary(as.vector(as.matrix(telescope_counts)))

sum(
  abs(
    as.matrix(telescope_counts) -
      round(as.matrix(telescope_counts))
  ) > 1e-8
)

#小数値の確認
sum(
  abs(
    as.matrix(telescope_counts) -
      round(as.matrix(telescope_counts))
  ) > 1e-8
)


# 7. Check sample IDs against metadata ----

# これまでTEcount / WGCNA / DESeq2で使用したsample metadata。
sample_metadata <- readRDS(
  "results/objects/TEcount_sample_metadata.rds"
)


cat("\n===== Sample metadata =====\n")

cat("Dimensions:\n")
print(dim(sample_metadata))

cat("\nColumn names:\n")
print(colnames(sample_metadata))

cat("\nFirst rows:\n")
print(head(sample_metadata))



# 8. Check Telescope sample IDs against metadata ----
#Telescopeの53列とmetadataの53行が「同じsample集合」かつ「同じ順番」かを確認

telescope_sample_ids <- colnames(telescope_counts)
metadata_sample_ids  <- sample_metadata$sample_id


# sample数を確認する。
cat("\n===== Sample ID check =====\n")

cat(
  "Telescope samples:",
  length(telescope_sample_ids),
  "\n"
)

cat(
  "Metadata samples:",
  length(metadata_sample_ids),
  "\n"
)


# 同じsample集合か確認する。

#setequal()は２つの集合が完全に一致してルカの確認


cat(
  "\nSame sample set:",
  setequal(
    telescope_sample_ids,
    metadata_sample_ids
  ),
  "\n"
)


# 順番まで完全一致しているか確認する。
#identical()は２つのdataやobjectが完全に一致しているか確認


cat(
  "Same sample order:",
  identical(
    telescope_sample_ids,
    metadata_sample_ids
  ),
  "\n"
)

cat("\nRow names:\n")
print(head(rownames(sample_metadata), 10))




# 9. Low-expression filtering ----

# TEcount解析と同じ基準を使用する。
#
# 各ERV locusについて、
# raw count >= 10となるsampleが
# 少なくとも6 samplesあるlocusのみを保持する。

n_samples_count10 <- rowSums(
  telescope_counts >= 10
)

#6sample以上なら残す
keep_locus <- n_samples_count10 >= 6



#ここまででraw count ≥ 10 が少なくとも6 samplesで観測されるlocusを保持した
#Filetrしたので確認
cat("\n===== Telescope low-expression filter =====\n")

cat(
  "Loci before filtering:",
  nrow(telescope_counts),
  "\n"
)




#filter後のMatrixを作成
telescope_counts_filtered <- telescope_counts[
  keep_locus,
  ,
  drop = FALSE
]
cat(
  "Loci retained:",
  sum(keep_locus),
  "\n"
)

cat(
  "Loci removed:",
  sum(!keep_locus),
  "\n"
)

#Matrixの構造確認
cat(
  "\nFiltered Telescope matrix dimensions:\n"
)

print(
  dim(telescope_counts_filtered)
)


#全部確認してみる
sum(keep_locus)
#[1] 11503    

sum(!keep_locus)
#[1] 2674    

dim(telescope_counts_filtered)
#[1] 11503    53




# 9.5 Save filtered Telescope matrix ----


saveRDS(
  telescope_counts_filtered,
  "results/objects/Telescope_locus_low_expression_filtered.rds"
)



# 10. TEcount subfamily ↔ Telescope locus mapping ----
#10.1 Telescope GTFを読む
# Telescope GTFのintModelと、
# TEcountで使用したERV subfamilyが
# 同じ命名体系であることを確認する。

wgcna1_feature_annotation <- readRDS(
  "results/objects/TEcount_WGCNA_feature_annotation.rds"
)

te_subfamilies <- unique(
  wgcna1_feature_annotation$te_subfamily[
    wgcna1_feature_annotation$feature_type == "TE"
  ]
)

head(te_subfamilies, 20)

"ERV3-16A3_I-int" %in% te_subfamilies



# Telescopeで使用したHERV annotationを読み込む。

telescope_gtf <- read.delim(
  "03_reference/telescope_annotation/HERV_rmsk.hg38.v2/transcripts.gtf",
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  quote = ""
)

colnames(telescope_gtf) <- c(
  "chr",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attribute"
)


#10.2 gene行だけ使う
# 1つのTelescope locusを1行で扱うため、
# feature == "gene" の行のみ使用する。

telescope_gene_gtf <- telescope_gtf[
  telescope_gtf$feature == "gene",
]


#10.3 locus IDとintModelを取り出す
# Telescope locus IDをattributeから抽出する。

telescope_gene_gtf$transcript_id <- sub(
  '.*transcript_id "([^"]+)".*',
  '\\1',
  telescope_gene_gtf$attribute
)


# Parent ERV subfamilyに相当するintModelを抽出する。

telescope_gene_gtf$intModel <- sub(
  '.*intModel "([^"]+)".*',
  '\\1',
  telescope_gene_gtf$attribute
)

# 10.4 crosswalkを作る
telescope_crosswalk <- telescope_gene_gtf[
  ,
  c(
    "transcript_id",
    "intModel",
    "chr",
    "start",
    "end",
    "strand"
  )
]


# 10.5 まず全体をQCする
cat("\n===== Telescope locus-subfamily crosswalk =====\n")

cat(
  "Rows:",
  nrow(telescope_crosswalk),
  "\n"
)

cat(
  "Unique Telescope loci:",
  length(unique(telescope_crosswalk$transcript_id)),
  "\n"
)

cat(
  "Duplicated Telescope loci:",
  sum(duplicated(telescope_crosswalk$transcript_id)),
  "\n"
)

cat(
  "Unique intModels:",
  length(unique(telescope_crosswalk$intModel)),
  "\n"
)

print(head(telescope_crosswalk, 10))


#次に、filtered 11,503 lociが全部GTFに存在するか確認
cat(
  "\nFiltered loci found in GTF:",
  sum(
    rownames(telescope_counts_filtered) %in%
      telescope_crosswalk$transcript_id
  ),
  "/",
  nrow(telescope_counts_filtered),
  "\n"
)


#さらにTelescope側のintModelのうち、
#TEcountで使った465 strict ERV subfamiliesと一致するものを確認
cat(
  "Telescope intModels matching TEcount strict ERV subfamilies:",
  sum(
    unique(telescope_crosswalk$intModel) %in%
      te_subfamilies
  ),
  "/",
  length(unique(telescope_crosswalk$intModel)),
  "\n"
)

#ここで
# Telescope GTF loci
# 14,968 loci
# ↓
# unique intModel = 60
# 
# Filtered Telescope loci
# 11,503 / 11,503
# → 全てGTFに存在
# 
# Telescope intModel
# 60 / 60
# → 全てTEcount strict ERV subfamily名と一致
#が確認できた
#つまり、Telescopeで使われている60種類のparent ERV modelについては、
#TEcount側の命名と100%対応できる



# 10.6 Read WGCNA1 module assignment 

wgcna1_module_assignment <- readRDS(
  "results/objects/WGCNA1_module_assignment.rds"
)

cat("\n===== WGCNA1 module assignment =====\n")

print(class(wgcna1_module_assignment))

print(dim(wgcna1_module_assignment))

print(colnames(wgcna1_module_assignment))

print(head(wgcna1_module_assignment))


# 10.7 Build TE subfamily-module map 

te_module_map <- merge(
  wgcna1_feature_annotation[
    wgcna1_feature_annotation$feature_type == "TE",
    c(
      "feature_id",
      "te_subfamily",
      "te_family",
      "te_class"
    )
  ],
  wgcna1_module_assignment[
    wgcna1_module_assignment$feature_type == "TE",
    c(
      "feature_id",
      "module"
    )
  ],
  by = "feature_id",
  all.x = TRUE
)


#10.8 465 ERVが全部入っているか確認
cat("\n===== TE subfamily-module map =====\n")

cat(
  "Number of TE subfamilies:",
  nrow(te_module_map),
  "\n"
)

cat(
  "Missing module assignments:",
  sum(is.na(te_module_map$module)),
  "\n"
)

cat("\nModule distribution:\n")

print(
  table(te_module_map$module)
)



#10.9 filtered Telescope lociにannotation
# 今の11,503 lociだけをcrosswalkから取り出します。
# ここではmerge()ではなくmatch()を使います。
# 理由は、Telescope count matrixと同じlocus順序を維持したいから

# Filter後のTelescope lociだけを、
# count matrixと同じ順序で取得する。

filtered_crosswalk <- telescope_crosswalk[
  match(
    rownames(telescope_counts_filtered),
    telescope_crosswalk$transcript_id
  ),
]


#10.10 intModelからWGCNA1 moduleを付ける

filtered_crosswalk$module <- te_module_map$module[
  match(
    filtered_crosswalk$intModel,
    te_module_map$te_subfamily
  )
]

# これで、以下の流れが完成
# Telescope locus
# ↓
# intModel
# ↓
# TEcount subfamily
# ↓
# WGCNA1 module
#



# 10.11 mapping結果を確認

cat("\n===== Filtered Telescope loci + WGCNA1 module =====\n")

cat(
  "Filtered loci:",
  nrow(filtered_crosswalk),
  "\n"
)

cat(
  "Loci with module assignment:",
  sum(!is.na(filtered_crosswalk$module)),
  "/",
  nrow(filtered_crosswalk),
  "\n"
)

cat(
  "Loci without module assignment:",
  sum(is.na(filtered_crosswalk$module)),
  "\n"
)

cat("\nLocus distribution by module:\n")

print(
  table(
    filtered_crosswalk$module,
    useNA = "ifany"
  )
)

print(
  head(filtered_crosswalk, 10)
)



#確認
table(te_module_map$module)

table(
  filtered_crosswalk$module,
  useNA = "ifany"
)
#sonokextuka
# 11,503 individual ERV loci
# ├── turquoise-derived  8,975
# └── yellow-derived     2,528



# 10.12 Check represented parent ERV subfamilies 
#11,503 lociが、何種類のparent subfamilyに由来するかを確認!!!!!!!!!
cat("\n===== Parent ERV subfamilies in filtered Telescope loci =====\n")

cat(
  "Total represented subfamilies:",
  length(unique(filtered_crosswalk$intModel)),
  "\n"
)

cat("\nSubfamilies by module:\n")

print(
  tapply(
    filtered_crosswalk$intModel,
    filtered_crosswalk$module,
    function(x) length(unique(x))
  )
)


# 11. Save Telescope locus preparation results ----

write.table(
  filtered_crosswalk,
  "results/tables/07_Telescope_locus_subfamily_module_mapping.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


mapping_summary <- data.frame(
  module = c("turquoise", "yellow"),
  n_subfamilies = c(44, 16),
  n_loci = c(8975, 2528)
)

write.table(
  mapping_summary,
  "results/tables/07_Telescope_mapping_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
