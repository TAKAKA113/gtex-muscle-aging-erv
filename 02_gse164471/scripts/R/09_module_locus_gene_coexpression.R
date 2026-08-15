# 09_module_locus_gene_coexpression.R
#
# 目的：
#   WGCNA1でERVが濃縮されたturquoise / yellow moduleについて、
#   Telescopeでindividual locusまで分解したERVと
#   同じmoduleに属するprotein-coding genesとの
#   coexpression relationshipを調べる。
#
# 対象：
#   turquoise:
#     8,975 Telescope loci
#     5,599 genes
#
#   yellow:
#     2,528 Telescope loci
#     1,437 genes

# WGCNA1でERVが濃縮されたturquoise/yellow moduleを対象に、年齢でfeatureを選別せず、
# Telescopeで分解したERV lociと同module genesのVST発現を使って、
# 53サンプル間のpairwise bicorを計算し、module内部のlocus–gene共発現構造を詳しく見る



setwd(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
)

library(DESeq2)


options(stringsAsFactors = FALSE)


# 2. Read input objects ----

dds_joint <- readRDS(
  "results/objects/08_DESeq2_joint_fitted.rds"
)

locus_annotation <- read.delim(
  "results/tables/07_Telescope_locus_subfamily_module_mapping.tsv"
)

wgcna1_module_assignment <- readRDS(
  "results/objects/WGCNA1_module_assignment.rds"
)


# 3. Create VST expression matrix ----
# raw countではなくVST後のexpressionを相関解析に使用する。
#
# blind = FALSE:
# DESeq2で推定済みのdesign/dispersion情報を利用する。
# ただしageやsexをexpressionから除去する処理ではない。

vsd_joint <- vst(
  dds_joint,
  blind = FALSE
)

vst_matrix <- assay(
  vsd_joint
)

#check
cat("\n===== VST matrix =====\n")

cat(
  "Features:",
  nrow(vst_matrix),
  "\n"
)

cat(
  "Samples:",
  ncol(vst_matrix),
  "\n"
)

cat(
  "Missing values:",
  sum(is.na(vst_matrix)),
  "\n"
)


# 4. Extract WGCNA1 module genes ----

turquoise_genes <- wgcna1_module_assignment$feature_id[
  wgcna1_module_assignment$feature_type == "GENE" &
    wgcna1_module_assignment$module == "turquoise"
]

yellow_genes <- wgcna1_module_assignment$feature_id[
  wgcna1_module_assignment$feature_type == "GENE" &
    wgcna1_module_assignment$module == "yellow"
]


#Check
cat("\n===== WGCNA1 genes =====\n")

cat(
  "Turquoise genes:",
  length(turquoise_genes),
  "\n"
)

cat(
  "Yellow genes:",
  length(yellow_genes),
  "\n"
)



# 5. Extract Telescope loci by parent WGCNA1 module ----

turquoise_loci <- locus_annotation$transcript_id[
  locus_annotation$module == "turquoise"
]

yellow_loci <- locus_annotation$transcript_id[
  locus_annotation$module == "yellow"
]


#Check
cat("\n===== Telescope loci =====\n")

cat(
  "Turquoise-associated loci:",
  length(turquoise_loci),
  "\n"
)

cat(
  "Yellow-associated loci:",
  length(yellow_loci),
  "\n"
)


# 6. Check feature IDs against VST matrix ----
#全featureがVST Matrixにあるか最終確認
cat("\n===== Feature ID check =====\n")

cat(
  "Turquoise genes found:",
  sum(
    turquoise_genes %in%
      rownames(vst_matrix)
  ),
  "/",
  length(turquoise_genes),
  "\n"
)

cat(
  "Yellow genes found:",
  sum(
    yellow_genes %in%
      rownames(vst_matrix)
  ),
  "/",
  length(yellow_genes),
  "\n"
)

cat(
  "Turquoise loci found:",
  sum(
    turquoise_loci %in%
      rownames(vst_matrix)
  ),
  "/",
  length(turquoise_loci),
  "\n"
)

cat(
  "Yellow loci found:",
  sum(
    yellow_loci %in%
      rownames(vst_matrix)
  ),
  "/",
  length(yellow_loci),
  "\n"
)



# 7. Prepare expression matrices for bicor ----
#相関用Matrixを作成
#bicor()では、行 = samples、列 = features にしておくと分かりやすいのでtransposeする
# 行 = sample
# 列 = ERV locus / gene

turquoise_locus_expr <- t(
  vst_matrix[
    turquoise_loci,
    ,
    drop = FALSE
  ]
)

turquoise_gene_expr <- t(
  vst_matrix[
    turquoise_genes,
    ,
    drop = FALSE
  ]
)


yellow_locus_expr <- t(
  vst_matrix[
    yellow_loci,
    ,
    drop = FALSE
  ]
)

yellow_gene_expr <- t(
  vst_matrix[
    yellow_genes,
    ,
    drop = FALSE
  ]
)

#Check
cat("\n===== Correlation input matrices =====\n")

cat(
  "Turquoise locus matrix:",
  dim(turquoise_locus_expr),
  "\n"
)

cat(
  "Turquoise gene matrix:",
  dim(turquoise_gene_expr),
  "\n"
)

cat(
  "Yellow locus matrix:",
  dim(yellow_locus_expr),
  "\n"
)

cat(
  "Yellow gene matrix:",
  dim(yellow_gene_expr),
  "\n"
)

# VST matrix
# 29,329 features × 53 samples
# ↓
# すべて欠損なし
# 
# turquoise
# 8,975 loci × 5,599 genes
# 
# yellow
# 2,528 loci × 1,437 genes






# 8. Calculate yellow locus-gene bicor ----
# WGCNA1と同じbicorを使用する
#特にturquoiseは約5,000万ペアなので、correlation matrixだけでもかなり大きくなる
#そこでlocusを500個ずつblockに分けて計算

#bicor
# 外れ値に比較的強い相関指標
# 正確には biweight midcorrelation で、極端なサンプルの影響をPearsonより受けにくくする
# robust correlation

#ここでの目的はageでERVを選別せず、全53サンプルにおけるERV locus–gene共発現を見る
#ためWGCNAで使ったmatrixの各数値をbicorに変換している！！！！！



# 行 = sample
# 列 = feature
#
# x = Telescope loci
# y = protein-coding genes

yellow_bicor <- WGCNA::bicor(
  yellow_locus_expr,
  yellow_gene_expr,
  maxPOutliers = 0.1
)


cat("\n===== Yellow bicor matrix =====\n")

cat(
  "Dimensions:",
  dim(yellow_bicor),
  "\n"
)

cat(
  "Missing correlations:",
  sum(is.na(yellow_bicor)),
  "\n"
)

cat(
  "Minimum bicor:",
  min(yellow_bicor, na.rm = TRUE),
  "\n"
)

cat(
  "Maximum bicor:",
  max(yellow_bicor, na.rm = TRUE),
  "\n"
)


dim(yellow_bicor)

sum(is.na(yellow_bicor))

range(yellow_bicor, na.rm = TRUE)

# 結果は
# NA               : 0
# 最小 bicor       : -0.7149
# 最大 bicor       : +0.9366
# 
# つまりyellow moduleの中には、
# ERV locusとgeneがかなり強く正に共発現しているpair,かなり強く負に共発現しているpair
# の両方が存在する



#yelloの数値確認
summary(
  as.vector(yellow_bicor)
)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.7149  0.1379  0.3058  0.2955  0.4630  0.9366 


quantile(
  abs(yellow_bicor),
  probs = c(
    0.90,
    0.95,
    0.99,
    0.999
  )
)
# 90%       95%       99%     99.9% 
#   0.5918425 0.6594113 0.7619328 0.8349032 




# 9. Flag loci with zero MAD ----
# zero MAD lociではbicorではなくPearson correlationが使用されたため、
# downstreamで通常のbicor pairと区別できるようflagを付ける

yellow_locus_method <- data.frame(
  transcript_id = colnames(yellow_locus_expr),
  correlation_method = ifelse(
    colnames(yellow_locus_expr) %in% yellow_zero_mad_loci,
    "Pearson_fallback",
    "bicor"
  )
)

print(
  table(yellow_locus_method$correlation_method)
)


# 10. Extract strongest yellow locus-gene correlations ----
# yellow module内の全3,632,736 pairから、
# |bicor|が上位0.1%に入るpairを探索的候補として抽出する。
#
# 注意：
# これは統計的有意性のthresholdではない。
# module内部で特に強いcoexpression pairを
# manageableな数に絞るためのranking基準として使用する。


# bicor matrixの行・列が想定したfeature IDになっているか確認
cat("\n===== Yellow bicor ID check =====\n")

cat(
  "Locus IDs identical:",
  identical(
    rownames(yellow_bicor),
    colnames(yellow_locus_expr)
  ),
  "\n"
)

cat(
  "Gene IDs identical:",
  identical(
    colnames(yellow_bicor),
    colnames(yellow_gene_expr)
  ),
  "\n"
)

#上位0.1%のthresholdを計算
# Absolute bicorの99.9 percentile
yellow_threshold <- quantile(
  abs(yellow_bicor),
  probs = 0.999,
  na.rm = TRUE
)

cat(
  "Yellow |bicor| 99.9% threshold:",
  yellow_threshold,
  "\n"
)

#Yellow |bicor| 99.9% threshold: 0.8349032

#上位pairを取り出す
yellow_top_index <- which(
  abs(yellow_bicor) >= yellow_threshold,
  arr.ind = TRUE
)

yellow_top_pairs <- data.frame(
  transcript_id = rownames(yellow_bicor)[
    yellow_top_index[, 1]
  ],
  
  gene_id = colnames(yellow_bicor)[
    yellow_top_index[, 2]
  ],
  
  bicor = yellow_bicor[
    yellow_top_index
  ]
)

yellow_top_pairs$abs_bicor <- abs(
  yellow_top_pairs$bicor
)

yellow_top_pairs$direction <- ifelse(
  yellow_top_pairs$bicor > 0,
  "positive",
  "negative"
)


#Pearson fallbackのflagを付ける
yellow_top_pairs$correlation_method <-
  yellow_locus_method$correlation_method[
    match(
      yellow_top_pairs$transcript_id,
      yellow_locus_method$transcript_id
    )
  ]

#parent ERV subfamilyも付る
yellow_top_pairs$intModel <-
  locus_annotation$intModel[
    match(
      yellow_top_pairs$transcript_id,
      locus_annotation$transcript_id
    )
  ]


#Gene symbolをつける
# Gene annotationを読み込む
wgcna1_feature_annotation <- readRDS(
  "results/objects/TEcount_WGCNA_feature_annotation.rds"
)

yellow_top_pairs$gene_symbol <-
  wgcna1_feature_annotation$gene_symbol[
    match(
      yellow_top_pairs$gene_id,
      wgcna1_feature_annotation$feature_id
    )
  ]


#相関が強い順に並べる
yellow_top_pairs <- yellow_top_pairs[
  order(
    yellow_top_pairs$abs_bicor,
    decreasing = TRUE
  ),
]



#結果確認
cat("\n===== Top yellow locus-gene pairs =====\n")

cat(
  "Number of top 0.1% pairs:",
  nrow(yellow_top_pairs),
  "\n"
)

cat("\nDirection:\n")

print(
  table(
    yellow_top_pairs$direction
  )
)

cat("\nCorrelation method:\n")

print(
  table(
    yellow_top_pairs$correlation_method
  )
)

print(
  head(
    yellow_top_pairs,
    20
  )
)


# 全3,632,736 pairのうち 3,633 pair が上位0.1%
# 3,633 pairはすべて正相関
# 3,633 pairはすべて通常のbicorで、Pearson fallbackは0


# 11. Save yellow coexpression results ----
# Yellow moduleのstrong ERV locus-gene pairを保存する
#
# Selection criterion:
# absolute bicorの上位0.1%
#
# 注意：
# 上位0.1%は統計的有意性を意味しない。
# downstream genomic annotationへ進めるcandidate prioritisationに使用する。


write.table(
  yellow_top_pairs,
  "results/tables/09_yellow_top01percent_locus_gene_pairs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  yellow_locus_method,
  "results/tables/09_yellow_locus_correlation_method.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


yellow_summary <- data.frame(
  module = "yellow",
  n_loci = ncol(yellow_locus_expr),
  n_genes = ncol(yellow_gene_expr),
  n_pairs = length(yellow_bicor),
  selection = "top 0.1% absolute bicor",
  bicor_cutoff = as.numeric(yellow_threshold),
  n_selected_pairs = nrow(yellow_top_pairs),
  n_zero_MAD_loci = length(yellow_zero_mad_loci)
)


write.table(
  yellow_summary,
  "results/tables/09_yellow_coexpression_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 2,528 loci
# ×
# 1,437 genes
# ↓
# 3,632,736 pairs
# ↓
# top 0.1%
# 3,633 strong pairs
# ↓
# cutoff |bicor| = 0.8349032




#__________________________________
#ここでturquoise作るよ
# 12. Check zero-MAD turquoise loci ----
# MAD = 0 のlocusではbicorが計算できないため、
# WGCNA::bicor()はそのlocusのみPearson correlationへfallbackする。
#
# これらは削除せず、
# downstreamで区別できるようflagを付ける。


turquoise_locus_mad <- apply(
  turquoise_locus_expr,
  2,
  mad
)


turquoise_zero_mad_loci <- names(
  turquoise_locus_mad[
    turquoise_locus_mad == 0
  ]
)


cat("\n===== Turquoise zero-MAD loci =====\n")

cat(
  "Turquoise loci with zero MAD:",
  length(turquoise_zero_mad_loci),
  "/",
  ncol(turquoise_locus_expr),
  "\n"
)


head(
  turquoise_zero_mad_loci,
  20
)



# correlation method flag

turquoise_locus_method <- data.frame(
  transcript_id = colnames(turquoise_locus_expr),
  
  correlation_method = ifelse(
    colnames(turquoise_locus_expr) %in%
      turquoise_zero_mad_loci,
    "Pearson_fallback",
    "bicor"
  )
)


print(
  table(
    turquoise_locus_method$correlation_method
  )
)



# 13. Calculate turquoise locus-gene correlations blockwise ----
#
# Turquoise:
#   8,975 loci × 5,599 genes
#
# Yellowと同様に、
# absolute bicorが上位0.1%のpairをcandidateとして残す。
#
# Turquoiseはpair数が非常に多いため、
# 500 lociずつblockに分けて計算する。
#
# このblock処理はmemory managementのためだけであり、
# bicorの計算方法自体はYellowと同じ。


turquoise_total_pairs <-
  ncol(turquoise_locus_expr) *
  ncol(turquoise_gene_expr)


turquoise_n_keep <- ceiling(
  turquoise_total_pairs * 0.001
)


cat("\n===== Turquoise correlation analysis =====\n")

cat(
  "Total pairs:",
  turquoise_total_pairs,
  "\n"
)

cat(
  "Top 0.1% pairs to retain:",
  turquoise_n_keep,
  "\n"
)



# 14. Run blockwise bicor ----
#blockごとにbicorを計算
block_size <- 500


block_starts <- seq(
  1,
  ncol(turquoise_locus_expr),
  by = block_size
)


turquoise_block_candidates <- vector(
  "list",
  length(block_starts)
)


# QC用
turquoise_min_bicor <- Inf
turquoise_max_bicor <- -Inf
turquoise_na_bicor <- 0


for (b in seq_along(block_starts)) {
  
  start_col <- block_starts[b]
  
  end_col <- min(
    start_col + block_size - 1,
    ncol(turquoise_locus_expr)
  )
  
  
  cat(
    "\nBlock",
    b,
    "/",
    length(block_starts),
    ": loci",
    start_col,
    "-",
    end_col,
    "\n"
  )
  
  
  # 500 lociだけ取り出す
  locus_block <- turquoise_locus_expr[
    ,
    start_col:end_col,
    drop = FALSE
  ]
  
  
  # ERV locus × gene bicor
  bicor_block <- WGCNA::bicor(
    locus_block,
    turquoise_gene_expr,
    maxPOutliers = 0.1
  )
  
  
  # QC
  turquoise_min_bicor <- min(
    turquoise_min_bicor,
    bicor_block,
    na.rm = TRUE
  )
  
  turquoise_max_bicor <- max(
    turquoise_max_bicor,
    bicor_block,
    na.rm = TRUE
  )
  
  turquoise_na_bicor <-
    turquoise_na_bicor +
    sum(is.na(bicor_block))
  
  
  # absolute bicor
  abs_block <- abs(
    bicor_block
  )
  
  
  # 各blockから十分多くの上位pairを残す
  #
  # Global top 0.1%に入るpairは、
  # 各block内でも少なくともglobal keep数以内に入るため、
  # 各blockからturquoise_n_keep個まで残せば
  # global top 0.1%を失わない。
  
  n_local_keep <- min(
    turquoise_n_keep,
    length(abs_block)
  )
  
  
  top_index <- order(
    abs_block,
    decreasing = TRUE
  )[
    seq_len(n_local_keep)
  ]
  
  
  top_position <- arrayInd(
    top_index,
    dim(bicor_block)
  )
  
  
  turquoise_block_candidates[[b]] <- data.frame(
    
    transcript_id =
      rownames(bicor_block)[
        top_position[, 1]
      ],
    
    gene_id =
      colnames(bicor_block)[
        top_position[, 2]
      ],
    
    bicor =
      bicor_block[
        top_index
      ],
    
    stringsAsFactors = FALSE
  )
  
  
  rm(
    locus_block,
    bicor_block,
    abs_block,
    top_index,
    top_position
  )
  
  gc()
}





# 15. Select global top 0.1% turquoise pairs ----

turquoise_candidates <- do.call(
  rbind,
  turquoise_block_candidates
)


turquoise_candidates$abs_bicor <- abs(
  turquoise_candidates$bicor
)


turquoise_candidates <- turquoise_candidates[
  order(
    turquoise_candidates$abs_bicor,
    decreasing = TRUE
  ),
]


turquoise_top_pairs <- turquoise_candidates[
  seq_len(turquoise_n_keep),
]



#cutoffも記録
turquoise_threshold <- min(
  turquoise_top_pairs$abs_bicor
)


cat(
  "\nTurquoise top 0.1% |bicor| cutoff:",
  turquoise_threshold,
  "\n"
)


# 16. Add annotation to turquoise top pairs ----


# positive / negative
turquoise_top_pairs$direction <- ifelse(
  turquoise_top_pairs$bicor > 0,
  "positive",
  "negative"
)


# bicor or Pearson fallback
turquoise_top_pairs$correlation_method <-
  turquoise_locus_method$correlation_method[
    match(
      turquoise_top_pairs$transcript_id,
      turquoise_locus_method$transcript_id
    )
  ]


# Parent ERV subfamily
turquoise_top_pairs$intModel <-
  locus_annotation$intModel[
    match(
      turquoise_top_pairs$transcript_id,
      locus_annotation$transcript_id
    )
  ]


# Gene symbol
turquoise_top_pairs$gene_symbol <-
  wgcna1_feature_annotation$gene_symbol[
    match(
      turquoise_top_pairs$gene_id,
      wgcna1_feature_annotation$feature_id
    )
  ]


#列を読みやすい順に
turquoise_top_pairs <- turquoise_top_pairs[
  ,
  c(
    "transcript_id",
    "intModel",
    "gene_id",
    "gene_symbol",
    "bicor",
    "abs_bicor",
    "direction",
    "correlation_method"
  )
]



# 17. Check turquoise results ----

cat("\n===== Turquoise top locus-gene pairs =====\n")


cat(
  "Total tested pairs:",
  turquoise_total_pairs,
  "\n"
)


cat(
  "Selected top 0.1% pairs:",
  nrow(turquoise_top_pairs),
  "\n"
)


cat(
  "Minimum bicor:",
  turquoise_min_bicor,
  "\n"
)


cat(
  "Maximum bicor:",
  turquoise_max_bicor,
  "\n"
)


cat(
  "Missing bicor:",
  turquoise_na_bicor,
  "\n"
)


cat(
  "Top 0.1% cutoff:",
  turquoise_threshold,
  "\n"
)


cat("\nDirection:\n")

print(
  table(
    turquoise_top_pairs$direction
  )
)


cat("\nCorrelation method:\n")

print(
  table(
    turquoise_top_pairs$correlation_method
  )
)


print(
  head(
    turquoise_top_pairs,
    20
  )
)



# 18. Save turquoise coexpression results ----

write.table(
  turquoise_top_pairs,
  "results/tables/09_turquoise_top01percent_locus_gene_pairs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


write.table(
  turquoise_locus_method,
  "results/tables/09_turquoise_locus_correlation_method.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


turquoise_summary <- data.frame(
  module = "turquoise",
  n_loci = ncol(turquoise_locus_expr),
  n_genes = ncol(turquoise_gene_expr),
  n_pairs = turquoise_total_pairs,
  selection = "top 0.1% absolute bicor",
  bicor_cutoff = turquoise_threshold,
  n_selected_pairs = nrow(turquoise_top_pairs),
  n_zero_MAD_loci = length(turquoise_zero_mad_loci)
)


write.table(
  turquoise_summary,
  "results/tables/09_turquoise_coexpression_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

