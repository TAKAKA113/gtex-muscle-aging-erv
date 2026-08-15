#10_genomic_annotation.R
#
# 目的：
#   Step 09で抽出したstrong ERV locus-gene coexpression pairsについて、
#   ERV locusとpaired geneのgenomic relationshipを調べる。
#
# 評価項目：
#   1. Same chromosome
#   2. Gene-body overlap
#   3. Exon overlap
#   4. Intron overlap
#   5. Distance to paired-gene TSS
#   6. Nearest gene
#
# Gene annotation：
#   GENCODE v26 primary assembly
#
#   03_reference/gencode_v26/
#   gencode.v26.primary_assembly.annotation.gtf
#
# このGTFはSTAR / featureCounts / TEcountで使用したものと同じ。
# expression analysisとgenomic annotationで同一annotationを使用する。


# 1. Set project ----

setwd(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
)


# 2. Load packages ----

library(GenomicRanges)
library(rtracklayer)
# 3. Read strong ERV-gene coexpression pairs ----

yellow_pairs <- read.delim(
  "results/tables/09_yellow_top01percent_locus_gene_pairs.tsv"
)

turquoise_pairs <- read.delim(
  "results/tables/09_turquoise_top01percent_locus_gene_pairs.tsv"
)


cat("\n===== Strong coexpression pairs =====\n")

cat(
  "Yellow:",
  nrow(yellow_pairs),
  "\n"
)

cat(
  "Turquoise:",
  nrow(turquoise_pairs),
  "\n"
)


# Yellow:    3633
# Turquoise: 50252



# 4. Import GENCODE v26 gene annotation ----

gencode_gtf <- rtracklayer::import(
  "03_reference/gencode_v26/gencode.v26.primary_assembly.annotation.gtf"
)


#v26 annotationの構造の確認
cat("\n===== GENCODE v26 annotation =====\n")

print(
  class(gencode_gtf)
)

print(
  length(gencode_gtf)
)

print(
  names(
    mcols(gencode_gtf)
  )
)

print(
  head(gencode_gtf)
)



# 5. Combine strong coexpression pairs ----
# Yellow / Turquoiseを同じ形式でまとめる
# まだfilterはしない

yellow_pairs$module <- "yellow"
turquoise_pairs$module <- "turquoise"

strong_pairs <- rbind(
  yellow_pairs,
  turquoise_pairs
)

strong_pairs$pair_id <- seq_len(
  nrow(strong_pairs)
)


cat("\n===== Strong ERV-gene pairs =====\n")

cat(
  "Total:",
  nrow(strong_pairs),
  "\n"
)

print(
  table(strong_pairs$module)
)


# 6. Extract protein-coding gene annotation ----
#GENCODEからprotein-cordingのみ抽出
gene_gr <- gencode_gtf[
  gencode_gtf$type == "gene" &
    gencode_gtf$gene_type == "protein_coding"
]


cat("\n===== Protein-coding gene annotation =====\n")

cat(
  "Genes:",
  length(gene_gr),
  "\n"
)

cat(
  "Duplicated gene IDs:",
  sum(
    duplicated(gene_gr$gene_id)
  ),
  "\n"
)


#geneがGENCODEに全部あるか確認
cat(
  "Strong-pair genes found in GENCODE:",
  sum(
    unique(strong_pairs$gene_id) %in%
      gene_gr$gene_id
  ),
  "/",
  length(unique(strong_pairs$gene_id)),
  "\n"
)


# 7. Extract protein-coding exon annotation ----

exon_gr <- gencode_gtf[
  gencode_gtf$type == "exon" &
    gencode_gtf$gene_type == "protein_coding"
]


cat(
  "\nProtein-coding exon records:",
  length(exon_gr),
  "\n"
)


# 8. Create GRanges for ERV loci ----
#Strong pairに含まれるERV lociをGRanges化
#ERV座標はTelescope annotationから既に作ったlocus_annotationを使う

strong_loci <- unique(
  strong_pairs$transcript_id
)


strong_locus_annotation <- locus_annotation[
  match(
    strong_loci,
    locus_annotation$transcript_id
  ),
]


erv_gr <- GRanges(
  seqnames = strong_locus_annotation$chr,
  
  ranges = IRanges(
    start = strong_locus_annotation$start,
    end = strong_locus_annotation$end
  ),
  
  strand = strong_locus_annotation$strand
)


erv_gr$transcript_id <-
  strong_locus_annotation$transcript_id


cat("\n===== Strong ERV loci =====\n")

cat(
  "Unique ERV loci:",
  length(erv_gr),
  "\n"
)

cat(
  "Missing locus IDs:",
  sum(is.na(erv_gr$transcript_id)),
  "\n"
)



# 9. Find ERV-gene body overlaps ----
# ERV × paired gene のgene-body overlap

# ここ本命
 
# findOverlaps()は2つのgenomic range集合のoverlapを検索するGenomicRangesの基本操作
# 今回はorientationではなく物理的位置の重なりを見たいのでignore.strand = TRUE



gene_hits <- findOverlaps(
  erv_gr,
  gene_gr,
  ignore.strand = TRUE
)


gene_overlap_keys <- paste(
  erv_gr$transcript_id[
    queryHits(gene_hits)
  ],
  
  gene_gr$gene_id[
    subjectHits(gene_hits)
  ],
  
  sep = "|||"
)


strong_pairs$pair_key <- paste(
  strong_pairs$transcript_id,
  strong_pairs$gene_id,
  sep = "|||"
)


strong_pairs$gene_body_overlap <-
  strong_pairs$pair_key %in%
  gene_overlap_keys


# 10. Find ERV-exon overlaps ----
# paired geneのannotated exonのいずれかと重なるかを確認する

exon_hits <- findOverlaps(
  erv_gr,
  exon_gr,
  ignore.strand = TRUE
)


exon_overlap_keys <- paste(
  erv_gr$transcript_id[
    queryHits(exon_hits)
  ],
  
  exon_gr$gene_id[
    subjectHits(exon_hits)
  ],
  
  sep = "|||"
)


strong_pairs$exon_overlap <-
  strong_pairs$pair_key %in%
  exon_overlap_keys



# 11. Define intronic ERV-gene pairs ----
#
# 今回の定義：
# gene body overlap = TRUE
# exon overlap      = FALSE
#
# つまりpaired geneのgene body内部だが、
# GENCODE v26でannotatedされたexonとは重ならない領域。

strong_pairs$intron_overlap <-
  strong_pairs$gene_body_overlap &
  !strong_pairs$exon_overlap



# 12. Add paired-gene coordinates ----
#TSS distanceを見る前にpaired geneと同じchromosomeかを付ける
gene_index <- match(
  strong_pairs$gene_id,
  gene_gr$gene_id
)


strong_pairs$gene_chr <-
  as.character(
    seqnames(gene_gr)[gene_index]
  )

strong_pairs$gene_start <-
  start(gene_gr)[gene_index]

strong_pairs$gene_end <-
  end(gene_gr)[gene_index]

strong_pairs$gene_strand <-
  as.character(
    strand(gene_gr)[gene_index]
  )


locus_index <- match(
  strong_pairs$transcript_id,
  locus_annotation$transcript_id
)


strong_pairs$erv_chr <-
  locus_annotation$chr[locus_index]

strong_pairs$erv_start <-
  locus_annotation$start[locus_index]

strong_pairs$erv_end <-
  locus_annotation$end[locus_index]


strong_pairs$same_chr <-
  strong_pairs$erv_chr ==
  strong_pairs$gene_chr


# 13. Build transcript TSS annotation ----
#
# 
# GENCODEでは1 geneに複数transcriptがあります。そのため今回、
# paired geneのannotated transcript TSSのうち、ERVに最も近いTSSまでの距離を使う
# + strand:
#   transcript start = TSS
#
# - strand:
#   transcript end = TSS

transcript_gr <- gencode_gtf[
  gencode_gtf$type == "transcript" &
    gencode_gtf$gene_type == "protein_coding"
]


transcript_tss <- data.frame(
  gene_id = transcript_gr$gene_id,
  
  chr = as.character(
    seqnames(transcript_gr)
  ),
  
  tss = ifelse(
    as.character(strand(transcript_gr)) == "+",
    start(transcript_gr),
    end(transcript_gr)
  ),
  
  stringsAsFactors = FALSE
)


# 同じgene・同じTSSが複数recordにある場合を除く
transcript_tss <- unique(
  transcript_tss
)



# 14. Calculate distance to nearest annotated TSS
#     of the paired gene ----
#TSSまでの最短距離を計算

pair_tss <- merge(
  strong_pairs[
    ,
    c(
      "pair_id",
      "gene_id",
      "erv_chr",
      "erv_start",
      "erv_end"
    )
  ],
  
  transcript_tss,
  
  by = "gene_id",
  all.x = TRUE
)


# ERVとTSSが別chromosomeならNA
pair_tss$tss_distance <- NA_real_


same_chr_tss <-
  pair_tss$erv_chr ==
  pair_tss$chr


# TSSがERVより左側
left_tss <-
  same_chr_tss &
  pair_tss$tss <
  pair_tss$erv_start


pair_tss$tss_distance[left_tss] <-
  pair_tss$erv_start[left_tss] -
  pair_tss$tss[left_tss]


# TSSがERVより右側
right_tss <-
  same_chr_tss &
  pair_tss$tss >
  pair_tss$erv_end


pair_tss$tss_distance[right_tss] <-
  pair_tss$tss[right_tss] -
  pair_tss$erv_end[right_tss]


# TSSがERV range内
inside_tss <-
  same_chr_tss &
  pair_tss$tss >=
  pair_tss$erv_start &
  pair_tss$tss <=
  pair_tss$erv_end


pair_tss$tss_distance[inside_tss] <- 0


#かくpairについて最小値
min_tss_distance <- tapply(
  pair_tss$tss_distance,
  pair_tss$pair_id,
  function(x) {
    
    if (all(is.na(x))) {
      return(NA_real_)
    }
    
    min(
      x,
      na.rm = TRUE
    )
  }
)


strong_pairs$tss_distance <- as.numeric(
  min_tss_distance[
    as.character(strong_pairs$pair_id)
  ]
)



# 15. Check whether paired gene is also nearest gene ----

nearest_hits <- GenomicRanges::nearest(
  erv_gr,
  gene_gr,
  select = "all",
  ignore.strand = TRUE
)


nearest_keys <- paste(
  erv_gr$transcript_id[
    queryHits(nearest_hits)
  ],
  
  gene_gr$gene_id[
    subjectHits(nearest_hits)
  ],
  
  sep = "|||"
)


strong_pairs$paired_gene_is_nearest <-
  strong_pairs$pair_key %in%
  nearest_keys



# 16. Summarize genomic relationships ----


cat("\n====================================\n")
cat("GENOMIC ANNOTATION SUMMARY\n")
cat("====================================\n")


for (mod in c("yellow", "turquoise")) {
  
  x <- strong_pairs[
    strong_pairs$module == mod,
  ]
  
  
  cat(
    "\n---",
    toupper(mod),
    "---\n"
  )
  
  
  cat(
    "Strong pairs:",
    nrow(x),
    "\n"
  )
  
  
  cat(
    "Same chromosome:",
    sum(x$same_chr),
    "\n"
  )
  
  
  cat(
    "Gene-body overlap:",
    sum(x$gene_body_overlap),
    "\n"
  )
  
  
  cat(
    "Exon overlap:",
    sum(x$exon_overlap),
    "\n"
  )
  
  
  cat(
    "Intron overlap:",
    sum(x$intron_overlap),
    "\n"
  )
  
  
  cat(
    "Non-overlap + TSS <= 500 bp:",
    sum(
      !x$gene_body_overlap &
        x$tss_distance <= 500,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "Non-overlap + TSS <= 2 kb:",
    sum(
      !x$gene_body_overlap &
        x$tss_distance <= 2000,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "Non-overlap + TSS <= 10 kb:",
    sum(
      !x$gene_body_overlap &
        x$tss_distance <= 10000,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "Paired gene is nearest gene:",
    sum(
      x$paired_gene_is_nearest
    ),
    "\n"
  )
}




# --- YELLOW ---
#   Strong pairs: 3633 
# Same chromosome: 198 
# Gene-body overlap: 1 
# Exon overlap: 0 
# Intron overlap: 1 
# Non-overlap + TSS <= 500 bp: 0 
# Non-overlap + TSS <= 2 kb: 0 
# Non-overlap + TSS <= 10 kb: 0 
# Paired gene is nearest gene: 3 
# 
# --- TURQUOISE ---
#   Strong pairs: 50252 
# Same chromosome: 3109 
# Gene-body overlap: 8 
# Exon overlap: 4 
# Intron overlap: 4 
# Non-overlap + TSS <= 500 bp: 0 
# Non-overlap + TSS <= 2 kb: 0 
# Non-overlap + TSS <= 10 kb: 0 
# Paired gene is nearest gene: 28 


# 17. Count unique ERV loci by genomic relationship ----

for (mod in c("yellow", "turquoise")) {
  
  x <- strong_pairs[
    strong_pairs$module == mod,
  ]
  
  
  cat(
    "\n---",
    toupper(mod),
    "UNIQUE ERV LOCI ---\n"
  )
  
  
  cat(
    "All strong-pair loci:",
    length(
      unique(x$transcript_id)
    ),
    "\n"
  )
  
  
  cat(
    "Gene-body overlap:",
    length(
      unique(
        x$transcript_id[
          x$gene_body_overlap
        ]
      )
    ),
    "\n"
  )
  
  
  cat(
    "Exonic:",
    length(
      unique(
        x$transcript_id[
          x$exon_overlap
        ]
      )
    ),
    "\n"
  )
  
  
  cat(
    "Intronic:",
    length(
      unique(
        x$transcript_id[
          x$intron_overlap
        ]
      )
    ),
    "\n"
  )
  
  
  cat(
    "Non-overlap + TSS <= 10 kb:",
    length(
      unique(
        x$transcript_id[
          !x$gene_body_overlap &
            !is.na(x$tss_distance) &
            x$tss_distance <= 10000
        ]
      )
    ),
    "\n"
  )
}

# --- YELLOW UNIQUE ERV LOCI ---
#   All strong-pair loci: 524 
# Gene-body overlap: 1 
# Exonic: 0 
# Intronic: 1 
# Non-overlap + TSS <= 10 kb: 0 
# 
# --- TURQUOISE UNIQUE ERV LOCI ---
#   All strong-pair loci: 2668 
# Gene-body overlap: 8 
# Exonic: 4 
# Intronic: 4 
# Non-overlap + TSS <= 10 kb: 0 


# 18. Check local candidate overlap ----

cat("\n===== Local candidate relationships =====\n")

cat(
  "Gene-body overlap AND paired gene is nearest:",
  sum(
    strong_pairs$gene_body_overlap &
      strong_pairs$paired_gene_is_nearest
  ),
  "\n"
)

cat(
  "Nearest gene but no gene-body overlap:",
  sum(
    !strong_pairs$gene_body_overlap &
      strong_pairs$paired_gene_is_nearest
  ),
  "\n"
)


#distanceも見る
nearest_nonoverlap <- strong_pairs[
  !strong_pairs$gene_body_overlap &
    strong_pairs$paired_gene_is_nearest,
]

summary(
  nearest_nonoverlap$tss_distance
)

sort(
  nearest_nonoverlap$tss_distance
)



# 19. Extract primary overlapping ERV-gene candidates ----
# Step 09でstrong coexpressionを示したERV locus-gene pairのうち、
# Step 10のGenomicRanges annotationで
# ERV locusがそのpaired gene bodyと直接overlapするpairを抽出する。
#
# ここでは新しいfilterは追加しない。
#
# 目的：
#   9個のprimary candidateについて、
#   ・どのERV locusか
#   ・どのERV subfamilyか
#   ・Yellow / Turquoiseのどちらか
#   ・どのgeneとstrong coexpressionを示したか
#   ・bicorはいくつか
#   ・exonic / intronicのどちらか
#   を確認する。
#
# 注意：
#   strong coexpression + genomic overlapは
#   ERVによるgene regulationの証明ではない。
#   high-priority ERV-gene candidate relationshipとして扱う。


overlap_candidates <- strong_pairs[
  strong_pairs$gene_body_overlap,
]


# Genomic contextを分かりやすい1列にまとめる
overlap_candidates$genomic_context <- ifelse(
  overlap_candidates$exon_overlap,
  "exonic",
  "intronic"
)


# bicorが強い順に並べる
overlap_candidates <- overlap_candidates[
  order(
    overlap_candidates$abs_bicor,
    decreasing = TRUE
  ),
]


# 結果確認
cat("\n===== Primary overlapping ERV-gene candidates =====\n")

cat(
  "Number of candidate pairs:",
  nrow(overlap_candidates),
  "\n"
)

cat(
  "Unique ERV loci:",
  length(
    unique(overlap_candidates$transcript_id)
  ),
  "\n"
)

cat(
  "Unique genes:",
  length(
    unique(overlap_candidates$gene_id)
  ),
  "\n"
)

cat("\nModule:\n")

print(
  table(
    overlap_candidates$module
  )
)

cat("\nGenomic context:\n")

print(
  table(
    overlap_candidates$genomic_context
  )
)


# 9 candidateの実体を見る
print(
  overlap_candidates[
    ,
    c(
      "module",
      "transcript_id",
      "intModel",
      "gene_id",
      "gene_symbol",
      "bicor",
      "abs_bicor",
      "genomic_context",
      "paired_gene_is_nearest"
    )
  ]
)



# 20. Save genomic annotation results ----
# Step 09で抽出したstrong ERV-gene coexpression pairsに
# GenomicRangesによるgenomic annotationを付加した結果を保存する。
#
# 保存するもの：
#   1. 全53,885 strong pairs + genomic annotation
#   2. gene-body overlapを示したprimary 9 candidates
#   3. nearest geneだがgene-body overlapしない22 pairs
#
# 注意：
#   primary candidateは
#   strong coexpression + direct gene-body overlap
#   を満たすpairとして保存する。
#
#   ただしERVによるgene regulationを証明するものではない。


# 1. Full strong-pair genomic annotation ----

saveRDS(
  strong_pairs,
  "results/objects/10_strong_pairs_genomic_annotation.rds"
)

write.table(
  strong_pairs,
  "results/tables/10_strong_pairs_genomic_annotation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 2. Primary overlapping candidates ----
# 9 ERV loci × 9 genes
# 4 exonic + 5 intronic

saveRDS(
  overlap_candidates,
  "results/objects/10_primary_overlapping_ERV_gene_candidates.rds"
)

write.table(
  overlap_candidates,
  "results/tables/10_primary_overlapping_ERV_gene_candidates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 3. Secondary nearest-gene pairs ----
# paired geneはnearest geneだがgene-body overlapしない22 pairs
# 最短TSS distance = 約24 kbなので、
# primary local candidateには含めないが記録として保存する。

saveRDS(
  nearest_nonoverlap,
  "results/objects/10_nearest_nonoverlap_ERV_gene_pairs.rds"
)

write.table(
  nearest_nonoverlap,
  "results/tables/10_nearest_nonoverlap_ERV_gene_pairs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# 4. Save genomic annotation summary ----

genomic_summary <- data.frame(
  metric = c(
    "Strong coexpression pairs",
    "Same chromosome",
    "Gene-body overlap",
    "Exon overlap",
    "Intron overlap",
    "Non-overlap + TSS <= 10 kb",
    "Paired gene is nearest",
    "Nearest gene but non-overlap"
  ),
  
  yellow = c(
    3633,
    198,
    1,
    0,
    1,
    0,
    3,
    2
  ),
  
  turquoise = c(
    50252,
    3109,
    8,
    4,
    4,
    0,
    28,
    20
  ),
  
  total = c(
    53885,
    3307,
    9,
    4,
    5,
    0,
    31,
    22
  )
)


write.table(
  genomic_summary,
  "results/tables/10_genomic_annotation_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(genomic_summary)

# metric yellow turquoise total
# 1    Strong coexpression pairs   3633     50252 53885
# 2              Same chromosome    198      3109  3307
# 3            Gene-body overlap      1         8     9
# 4                 Exon overlap      0         4     4
# 5               Intron overlap      1         4     5
# 6   Non-overlap + TSS <= 10 kb      0         0     0
# 7       Paired gene is nearest      3        28    31
# 8 Nearest gene but non-overlap      2        20    22