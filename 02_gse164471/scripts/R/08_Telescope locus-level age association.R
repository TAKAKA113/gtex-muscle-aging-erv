# 08_Telescope_age_association.R
#
# 目的：
#   low-expression filtering後のindividual ERV lociについて、
#   性別を調整した年齢関連解析を行う。
#
# model：
#   locus raw count ~ sex + age
#
# 背景：
#   TEcountを用いたsubfamily-level解析では、
#   FDR < 0.05のage-associated ERV subfamilyは
#   検出されなかった。
#
#   そこでTelescopeを用いてindividual ERV lociへ
#   解像度を上げ、subfamily aggregationでは見えない
#   locus-specificなage associationが存在するかを調べる。
#
# 対象：
#   11,503 filtered Telescope loci
#   53 skeletal muscle samples
#
#   parent ERV subfamily:
#     turquoise = 44 subfamilies / 8,975 loci
#     yellow    = 16 subfamilies / 2,528 loci
#
# model：
#   raw count ~ sex + age
#
# age coefficient：
#   beta > 0 : expression increases with age
#   beta < 0 : expression decreases with age
#
# Multiple testing：
#   Benjamini-Hochberg FDR



# 1. Set project and load package ----

setwd(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
)

library(DESeq2)

options(stringsAsFactors = FALSE)


# 2. Read input objects ----

# Low-expression filtering後のTelescope locus counts。
# 11,503 loci × 53 samples。
telescope_counts <- readRDS(
  "results/objects/Telescope_locus_low_expression_filtered.rds"
)

# Telescope locus → parent ERV subfamily → WGCNA1 module対応表。
locus_annotation <- read.delim(
  "results/tables/07_Telescope_locus_subfamily_module_mapping.tsv"
)

# sample metadata。
sample_metadata <- readRDS(
  "results/objects/TEcount_sample_metadata.rds"
)

# 以前のDESeq2解析で使用した、
# low-expression filtering後のGene + ERV subfamily count object。
tecount_dds <- readRDS(
  "results/objects/TEcount_dds_low_expression_filtered.rds"
)

# feature annotation。
tecount_annotation <- readRDS(
  "results/objects/TEcount_DESeq2_feature_annotation.rds"
)


# 3. Extract gene raw counts ----

tecount_raw <- counts(
  tecount_dds,
  normalized = FALSE
)



#annotationとrow orderが一致しているか確認
cat("\n===== TEcount raw counts =====\n")

cat(
  "Annotation order identical:",
  identical(
    rownames(tecount_raw),
    tecount_annotation$feature_id
  ),
  "\n"
)


#Geneだけ抽出
gene_ids <- tecount_annotation$feature_id[
  tecount_annotation$feature_type == "GENE"
]

gene_counts <- tecount_raw[
  gene_ids,
  ,
  drop = FALSE
]

cat(
  "Gene count matrix dimensions:",
  dim(gene_counts),
  "\n"
)

#Gene count matrix dimensions: 17826 53 


# 4. Check sample order ----
#TelescopeとGeneでsample orderを確認
cat("\n===== Joint matrix sample check =====\n")

cat(
  "Gene vs Telescope:",
  identical(
    colnames(gene_counts),
    colnames(telescope_counts)
  ),
  "\n"
)

cat(
  "Gene vs metadata:",
  identical(
    colnames(gene_counts),
    rownames(sample_metadata)
  ),
  "\n"
)

cat(
  "Telescope vs metadata:",
  identical(
    colnames(telescope_counts),
    rownames(sample_metadata)
  ),
  "\n"
)



# 5. Gene + Telescope locus count matrixを作る -------------------------------
#まずTelescope側をmatrixにする
telescope_counts_matrix <- as.matrix(
  telescope_counts
)

#
storage.mode(telescope_counts_matrix) <- "integer"

#GeneとTelescope locusでrow nameが衝突していないことも確認
cat(
  "\nDuplicated IDs between genes and loci:",
  sum(
    rownames(gene_counts) %in%
      rownames(telescope_counts_matrix)
  ),
  "\n"
)



#ここで結合
joint_counts <- rbind(
  gene_counts,
  telescope_counts_matrix
)


#確認
cat("\n===== Joint count matrix =====\n")

cat(
  "Genes:",
  nrow(gene_counts),
  "\n"
)

cat(
  "Telescope loci:",
  nrow(telescope_counts_matrix),
  "\n"
)

cat(
  "Total features:",
  nrow(joint_counts),
  "\n"
)

cat(
  "Samples:",
  ncol(joint_counts),
  "\n"
)


# 29,329 features × 53 samples
# 
# 17,826 host genes
# +
#   11,503 individual ERV loci
#というDESeq2 inputを確認できた




# 6. metadata for DESeq2 ---------------------------------------------------
#raw count ~ sex + age

coldata <- sample_metadata[
  colnames(joint_counts),
  c(
    "sex",
    "age"
  )
]


#確認
identical(
  rownames(tecount_raw),
  tecount_annotation$feature_id
)

dim(gene_counts)

identical(
  colnames(gene_counts),
  colnames(telescope_counts)
)

sum(
  rownames(gene_counts) %in%
    rownames(telescope_counts_matrix)
)

dim(joint_counts)

identical(
  colnames(joint_counts),
  rownames(coldata)
)

# 7. Run DESeq2 age-association model ----
#DESeq2を実行
dds_joint <- DESeqDataSetFromMatrix(
  countData = joint_counts,
  colData = coldata,
  design = ~ sex + age
)

dds_joint <- DESeq(dds_joint)


# Check model coefficients
resultsNames(dds_joint)



# 8. Extract Telescope locus dataset ----
# Telescope lociだけ取り出す
# 29,329 features全部から、11,503 Telescope lociだけ抜きだす
#
# 修正メモ：
# results()を29,329 features全体に実行してからTelescope lociを抜き出すと、
# independent filteringとBH multiple-testing correctionも
# Gene + Telescope loci全体を対象として行われてしまう。
#
# 今回のprimary testは11,503 Telescope lociなので、
# DESeq2 model fitting後にまずTelescope lociだけをsubsetし、
# その後results()を実行する。
#
# dds_jointはすでにDESeq()済みなので、
# size factor・dispersion・GLMを最初から再計算しているわけではない。

dds_locus <- dds_joint[
  rownames(telescope_counts_matrix),
]



# 9. Extract age-association results ----
#age coefficientを取り出す
# 11,503 Telescope lociを対象にresults()を実行する
#
# 修正メモ：
# 以前のage_results_allは不要。
# 今後はlocus-specific FDRを得るため、
# dds_locusに対して直接results()を実行する。

locus_age_results <- results(
  dds_locus,
  name = "age",
  alpha = 0.05
)

locus_age_results <- as.data.frame(
  locus_age_results
)

locus_age_results$transcript_id <- rownames(
  locus_age_results
)


#check
cat("\n===== Telescope age-association results =====\n")

cat(
  "Number of loci:",
  nrow(locus_age_results),
  "\n"
)

cat(
  "Missing p-values:",
  sum(is.na(locus_age_results$pvalue)),
  "\n"
)

cat(
  "Missing adjusted p-values:",
  sum(is.na(locus_age_results$padj)),
  "\n"
)



# 10. Add locus annotation ----
#locus annotationを付ける
annotation_index <- match(
  locus_age_results$transcript_id,
  locus_annotation$transcript_id
)

locus_age_results$intModel <-
  locus_annotation$intModel[annotation_index]

locus_age_results$chr <-
  locus_annotation$chr[annotation_index]

locus_age_results$start <-
  locus_annotation$start[annotation_index]

locus_age_results$end <-
  locus_annotation$end[annotation_index]

locus_age_results$strand <-
  locus_annotation$strand[annotation_index]

locus_age_results$module <-
  locus_annotation$module[annotation_index]


#Check
cat(
  "Missing locus annotations:",
  sum(is.na(locus_age_results$intModel)),
  "\n"
)

print(
  table(
    locus_age_results$module,
    useNA = "ifany"
  )
)

# turquoise    yellow 
# 8975         2528 



# 11. Calculate 10-year age effect ----
#10年間のeffect sizeを追加
#DESeq2のage coefficientは1歳あたりなので、修論では10年間あたりも見る

locus_age_results$log2FC_per_10y <-
  locus_age_results$log2FoldChange * 10

locus_age_results$fold_change_per_10y <-
  2 ^ locus_age_results$log2FC_per_10y

locus_age_results$percent_change_per_10y <-
  (
    locus_age_results$fold_change_per_10y - 1
  ) * 100



# 12. Summarize age-associated Telescope loci ----
#結果全体を確認
cat("\n===== Age-associated Telescope loci =====\n")

cat(
  "FDR < 0.05:",
  sum(
    locus_age_results$padj < 0.05,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FDR < 0.05, beta > 0:",
  sum(
    locus_age_results$padj < 0.05 &
      locus_age_results$log2FoldChange > 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "FDR < 0.05, beta < 0:",
  sum(
    locus_age_results$padj < 0.05 &
      locus_age_results$log2FoldChange < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Nominal p < 0.05:",
  sum(
    locus_age_results$pvalue < 0.05,
    na.rm = TRUE
  ),
  "\n"
)



#module別に確認
sig_loci <- locus_age_results[
  !is.na(locus_age_results$padj) &
    locus_age_results$padj < 0.05,
]

cat("\nSignificant loci by module:\n")

print(
  table(sig_loci$module)
)


#上位の確認
sig_loci <- sig_loci[
  order(sig_loci$padj),
]

print(
  head(
    sig_loci[
      ,
      c(
        "transcript_id",
        "intModel",
        "module",
        "baseMean",
        "log2FoldChange",
        "log2FC_per_10y",
        "percent_change_per_10y",
        "pvalue",
        "padj"
      )
    ],
    20
  )
)

# 13. Check DESeq2 model convergence ----
# DESeq2でbeta coefficientが収束しなかったfeatureを確認する

nonconv_ids <- rownames(dds_joint)[
  which(mcols(dds_joint)$betaConv == FALSE)
]

cat("\n===== Non-converged features =====\n")

cat(
  "Number of non-converged features:",
  length(nonconv_ids),
  "\n"
)

print(nonconv_ids)


# Telescope lociか確認
cat(
  "\nNon-converged Telescope loci:",
  sum(
    nonconv_ids %in%
      rownames(telescope_counts_matrix)
  ),
  "\n"
)

print(
  nonconv_ids[
    nonconv_ids %in%
      rownames(telescope_counts_matrix)
  ]
)


# Significant MER41 locusが正常に収束しているか確認
cat(
  "\nMER41_5q14.3a beta convergence:",
  mcols(dds_joint)[
    "MER41_5q14.3a",
    "betaConv"
  ],
  "\n"
)


# 14. Save Telescope age-association results ----

# Full locus-level age-association results
saveRDS(
  locus_age_results,
  "results/objects/08_Telescope_locus_age_results.rds"
)

write.table(
  locus_age_results,
  "results/tables/08_Telescope_locus_age_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# FDR-significant loci only
write.table(
  sig_loci,
  "results/tables/08_Telescope_FDR05_loci.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Fitted DESeq2 object
# 後でVSTを作る時にDESeq()をやり直さなくて済むよう保存
saveRDS(
  dds_joint,
  "results/objects/08_DESeq2_joint_fitted.rds"
)



# ============================================================
# Additional analysis: MER41_5q14.3a follow-up
# ============================================================
#
# MER41_5q14.3a was the only Telescope locus showing
# significant age association after FDR correction.
#
# This section characterises this locus separately from the
# 9 ERV-gene candidates identified by coexpression and
# genomic-overlap filtering.
# ============================================================


# 15. Extract MER41_5q14.3a age-association result ----

mer41_age <- locus_age_results[
  locus_age_results$transcript_id == "MER41_5q14.3a",
]


cat("\n===== MER41_5q14.3a age association =====\n")

print(mer41_age)


cat(
  "\nNumber of matching loci:",
  nrow(mer41_age),
  "\n"
)

cat(
  "FDR:",
  mer41_age$padj,
  "\n"
)

cat(
  "log2FC per year:",
  mer41_age$log2FoldChange,
  "\n"
)

cat(
  "log2FC per 10 years:",
  mer41_age$log2FC_per_10y,
  "\n"
)

cat(
  "Percent change per 10 years:",
  mer41_age$percent_change_per_10y,
  "%\n"
)


# 16. Annotate genomic context of MER41_5q14.3a ----

library(GenomicRanges)
library(rtracklayer)

# Same GENCODE annotation used in the main pipeline
gencode_gtf <- import(
  "03_reference/gencode_v26/gencode.v26.primary_assembly.annotation.gtf"
)

gene_gtf <- gencode_gtf[gencode_gtf$type == "gene"]
exon_gtf <- gencode_gtf[gencode_gtf$type == "exon"]


# MER41 genomic interval
mer41_gr <- GRanges(
  seqnames = mer41_age$chr,
  ranges = IRanges(
    start = mer41_age$start,
    end   = mer41_age$end
  ),
  strand = mer41_age$strand
)


# 16.1 Check direct gene-body overlap ----

gene_hits <- findOverlaps(
  mer41_gr,
  gene_gtf,
  ignore.strand = TRUE
)

cat("\n===== MER41_5q14.3a genomic context =====\n")
cat("Gene-body overlaps:", length(gene_hits), "\n")


if (length(gene_hits) > 0) {
  
  overlapping_genes <- gene_gtf[subjectHits(gene_hits)]
  
  overlap_gene_ids <- mcols(overlapping_genes)$gene_id
  
  # Check whether MER41 overlaps any exon of those genes
  candidate_exons <- exon_gtf[
    mcols(exon_gtf)$gene_id %in% overlap_gene_ids
  ]
  
  exon_hits <- findOverlaps(
    mer41_gr,
    candidate_exons,
    ignore.strand = TRUE
  )
  
  genomic_context <- ifelse(
    length(exon_hits) > 0,
    "exonic",
    "intronic"
  )
  
  mer41_genomic_context <- data.frame(
    transcript_id = "MER41_5q14.3a",
    erv_chr = as.character(seqnames(mer41_gr)),
    erv_start = start(mer41_gr),
    erv_end = end(mer41_gr),
    erv_strand = as.character(strand(mer41_gr)),
    gene_id = mcols(overlapping_genes)$gene_id,
    gene_symbol = mcols(overlapping_genes)$gene_name,
    genomic_context = genomic_context,
    stringsAsFactors = FALSE
  )
  
} else {
  
  # 16.2 If no gene-body overlap, identify nearest gene TSS ----
  
  gene_tss <- resize(
    gene_gtf,
    width = 1,
    fix = "start"
  )
  
  nearest_hit <- distanceToNearest(
    mer41_gr,
    gene_tss,
    ignore.strand = TRUE
  )
  
  nearest_gene <- gene_gtf[subjectHits(nearest_hit)]
  
  mer41_genomic_context <- data.frame(
    transcript_id = "MER41_5q14.3a",
    erv_chr = as.character(seqnames(mer41_gr)),
    erv_start = start(mer41_gr),
    erv_end = end(mer41_gr),
    erv_strand = as.character(strand(mer41_gr)),
    gene_id = mcols(nearest_gene)$gene_id,
    gene_symbol = mcols(nearest_gene)$gene_name,
    genomic_context = "non-overlapping",
    tss_distance = mcols(nearest_hit)$distance,
    stringsAsFactors = FALSE
  )
}


print(mer41_genomic_context)


# transcript_id erv_chr erv_start  erv_end erv_strand           gene_id gene_symbol
# 1 MER41_5q14.3a    chr5  89917678 89919657          - ENSG00000214942.4  AC113167.1
# genomic_context tss_distance
# 1 non-overlapping        17013


# 16.3 Save MER41 genomic annotation ----

write.table(
  mer41_genomic_context,
  "results/tables/08_MER41_5q14.3a_genomic_context.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

