# ============================================================
# 11_Telescope_figures.R
# ============================================================
#
# Purpose:
# Telescope locus-level analysisの主要結果を
# dissertation用Figureとして可視化する。
#
# このscriptでは新しいfilterや統計解析は行わない。
# Step 08–10で確定・保存した結果のみを使用する。
#
# Planned figures:
#
# Figure 1: Telescope locus-level analysis overview
#   A. Telescope / analysis workflow
#   B. Locus-level age-association volcano plot
#   C. Candidate prioritisation flow
#
# Figure 2: Primary ERV-gene candidates
#   A. 9 candidate pairs ranked by bicor
#   B-E. Representative expression/genomic examples
# ============================================================


# 1. Setup ----

setwd(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
)

library(ggplot2)


# 2. Load saved results ----

# Telescope locus-level age association
locus_age_results <- readRDS(
  "results/objects/08_Telescope_locus_age_results.rds"
)

# All strong ERV-gene pairs with genomic annotation
strong_pairs <- readRDS(
  "results/objects/10_strong_pairs_genomic_annotation.rds"
)

# Primary 9 overlapping ERV-gene candidates
overlap_candidates <- readRDS(
  "results/objects/10_primary_overlapping_ERV_gene_candidates.rds"
)


# 3. Input QC ----

cat("\n===== 11 Telescope figures: input QC =====\n")

cat(
  "Telescope age results:",
  nrow(locus_age_results),
  "rows\n"
)

cat(
  "Strong ERV-gene pairs:",
  nrow(strong_pairs),
  "rows\n"
)

cat(
  "Primary overlap candidates:",
  nrow(overlap_candidates),
  "rows\n"
)


cat("\nTelescope age-result columns:\n")
print(
  names(locus_age_results)
)

cat("\nPrimary candidate columns:\n")
print(
  names(overlap_candidates)
)



names(locus_age_results)
names(overlap_candidates)


# 4. Figure 1B: Telescope locus-level age association ----
#
# x-axis:
#   age 10年あたりのlog2 fold change
#
# y-axis:
#   -log10(p-value)
#
# FDR < 0.05のlocusのみラベルする。
#
# 注意：
# padj < 0.05で有意性を判定するため、
# p = 0.05のhorizontal lineは描かない。
# p = 0.05はFDR thresholdではないため。


volcano_df <- locus_age_results


# Plottingできないp-value NAを除外
volcano_df <- volcano_df[
  !is.na(volcano_df$pvalue),
]


# FDR significance
volcano_df$FDR_significant <-
  !is.na(volcano_df$padj) &
  volcano_df$padj < 0.05


# y-axis
volcano_df$minus_log10_p <-
  -log10(volcano_df$pvalue)


cat("\n===== Telescope age-association volcano =====\n")

cat(
  "Plotted loci:",
  nrow(volcano_df),
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(volcano_df$FDR_significant),
  "\n"
)

print(
  volcano_df[
    volcano_df$FDR_significant,
    c(
      "transcript_id",
      "intModel",
      "module",
      "log2FC_per_10y",
      "pvalue",
      "padj"
    )
  ]
)


#plot
p_telescope_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FC_per_10y,
    y = minus_log10_p
  )
) +
  
  geom_point(
    aes(
      colour = FDR_significant
    ),
    alpha = 0.6,
    size = 1.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  geom_text(
    data = volcano_df[
      volcano_df$FDR_significant,
    ],
    aes(
      label = transcript_id
    ),
    vjust = -0.7,
    size = 3.5,
    show.legend = FALSE
  ) +
  
  labs(
    title = "Age association of individual ERV loci",
    subtitle = paste0(
      "Telescope loci; FDR < 0.05: ",
      sum(volcano_df$FDR_significant)
    ),
    x = "log2 fold change per 10 years of age",
    y = expression(-log[10](italic(p))),
    colour = "FDR < 0.05"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(
      face = "bold"
    )
  )


print(p_telescope_volcano)

# 5. Save Telescope age volcano ----

ggsave(
  "results/figures/11_Telescope_age_association_volcano.pdf",
  p_telescope_volcano,
  width = 7,
  height = 5
)

ggsave(
  "results/figures/11_Telescope_age_association_volcano.png",
  p_telescope_volcano,
  width = 7,
  height = 5,
  dpi = 300
)



# 6. Figure 1: Telescope candidate prioritisation flow ----
#
# Telescopeによるindividual ERV locus解析から、
# 最終9 ERV-gene candidate pairsまでの流れを示す。
#
# 注意：
# 途中で単位が
# "ERV loci" → "ERV-gene pairs"
# に変わるため、Figure内に明記する。


# Raw Telescope loci
telescope_raw <- readRDS(
  "results/objects/Telescope_locus_raw.rds"
)

# Low-expression filtered Telescope loci
telescope_filtered <- readRDS(
  "results/objects/Telescope_locus_low_expression_filtered.rds"
)


# 数字をobjectから取得
n_raw_loci <- nrow(telescope_raw)

n_filtered_loci <- nrow(telescope_filtered)

module_locus_counts <- table(
  locus_age_results$module
)

n_turquoise_loci <- unname(
  module_locus_counts["turquoise"]
)

n_yellow_loci <- unname(
  module_locus_counts["yellow"]
)

n_strong_pairs <- nrow(
  strong_pairs
)

n_primary_pairs <- nrow(
  overlap_candidates
)

n_exonic <- sum(
  overlap_candidates$genomic_context == "exonic"
)

n_intronic <- sum(
  overlap_candidates$genomic_context == "intronic"
)


cat("\n===== Telescope candidate flow =====\n")

cat("Raw Telescope loci:", n_raw_loci, "\n")
cat("Filtered loci:", n_filtered_loci, "\n")
cat("Turquoise-associated loci:", n_turquoise_loci, "\n")
cat("Yellow-associated loci:", n_yellow_loci, "\n")
cat("Strong ERV-gene pairs:", n_strong_pairs, "\n")
cat("Primary overlapping pairs:", n_primary_pairs, "\n")
cat("Exonic:", n_exonic, "\n")
cat("Intronic:", n_intronic, "\n")


# 7. Draw Telescope analysis flow ----

flow_df <- data.frame(
  x = c(
    0,
    0,
    -1.4,
    1.4,
    0,
    0,
    -0.8,
    0.8
  ),
  
  y = c(
    7,
    5.8,
    4.5,
    4.5,
    3.0,
    1.7,
    0.4,
    0.4
  ),
  
  label = c(
    paste0(
      "Telescope quantification\n",
      n_raw_loci,
      " individual ERV loci"
    ),
    
    paste0(
      "Low-expression filtering\n",
      n_filtered_loci,
      " expressed loci"
    ),
    
    paste0(
      "Turquoise\n",
      n_turquoise_loci,
      " loci"
    ),
    
    paste0(
      "Yellow\n",
      n_yellow_loci,
      " loci"
    ),
    
    paste0(
      "Module-constrained locus–gene coexpression\n",
      n_strong_pairs,
      " strong ERV–gene pairs"
    ),
    
    paste0(
      "Genomic annotation\n",
      n_primary_pairs,
      " overlapping ERV–gene pairs"
    ),
    
    paste0(
      n_exonic,
      " exonic"
    ),
    
    paste0(
      n_intronic,
      " intronic"
    )
  )
)


p_telescope_flow <- ggplot() +
  
  # boxes
  geom_label(
    data = flow_df,
    aes(
      x = x,
      y = y,
      label = label
    ),
    size = 3.6,
    label.size = 0.3,
    lineheight = 1.05
  ) +
  
  # raw → filtered
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 6.65,
    yend = 6.15,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  # filtered → turquoise
  annotate(
    "segment",
    x = 0,
    xend = -1.4,
    y = 5.45,
    yend = 4.85,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  # filtered → yellow
  annotate(
    "segment",
    x = 0,
    xend = 1.4,
    y = 5.45,
    yend = 4.85,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  # modules → coexpression
  annotate(
    "segment",
    x = -1.4,
    xend = 0,
    y = 4.15,
    yend = 3.35,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  annotate(
    "segment",
    x = 1.4,
    xend = 0,
    y = 4.15,
    yend = 3.35,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  # coexpression → genomic annotation
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 2.65,
    yend = 2.05,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  # genomic annotation → exon
  annotate(
    "segment",
    x = 0,
    xend = -0.8,
    y = 1.35,
    yend = 0.75,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  # genomic annotation → intron
  annotate(
    "segment",
    x = 0,
    xend = 0.8,
    y = 1.35,
    yend = 0.75,
    arrow = arrow(length = unit(0.15, "cm"))
  ) +
  
  labs(
    title = "Locus-specific ERV candidate prioritisation using Telescope"
  ) +
  
  xlim(-3, 3) +
  ylim(-0.2, 7.7) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0.5
    )
  )


print(p_telescope_flow)


# 8. Save Telescope flow figure ----

ggsave(
  "results/figures/11_Telescope_candidate_flow.pdf",
  p_telescope_flow,
  width = 7,
  height = 7
)

ggsave(
  "results/figures/11_Telescope_candidate_flow.png",
  p_telescope_flow,
  width = 7,
  height = 7,
  dpi = 300
)


# 9. Figure: Primary ERV-gene candidates ranked by bicor ----

candidate_plot_df <- overlap_candidates

candidate_plot_df$candidate_label <- paste0(
  candidate_plot_df$transcript_id,
  "  –  ",
  candidate_plot_df$gene_symbol
)

candidate_plot_df$candidate_label <- factor(
  candidate_plot_df$candidate_label,
  levels = candidate_plot_df$candidate_label[
    order(candidate_plot_df$bicor)
  ]
)


p_candidate_rank <- ggplot(
  candidate_plot_df,
  aes(
    x = bicor,
    y = candidate_label,
    shape = genomic_context,
    colour = module
  )
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text(
    aes(
      label = sprintf("%.3f", bicor)
    ),
    nudge_x = 0.006,
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  ) +
  
  labs(
    title = "Primary ERV–gene candidates ranked by coexpression",
    x = "Biweight midcorrelation (bicor)",
    y = NULL,
    shape = "Genomic context",
    colour = "WGCNA module"
  ) +
  
  coord_cartesian(
    xlim = c(0.85, 0.965),
    clip = "off"
  ) +
  
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 2)
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.position = "bottom",
    
    legend.box = "horizontal",
    
    axis.text.y = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      t = 15,
      r = 35,
      b = 15,
      l = 15
    )
  )


print(p_candidate_rank)


# 10. Save candidate ranked plot ----

ggsave(
  "results/figures/11_Telescope_primary_candidates_bicor.pdf",
  p_candidate_rank,
  width = 9.5,
  height = 6
)

ggsave(
  "results/figures/11_Telescope_primary_candidates_bicor.png",
  p_candidate_rank,
  width = 9.5,
  height = 6,
  dpi = 300
)



# 11. Representative candidate expression scatter ----
#
# Primary 9 candidatesのうち、
# exonic / intronicそれぞれbicor最大のpairを代表例として可視化する。
#
# 新しいcandidate selectionではなく、
# 既に確定した9 candidatesの代表的なexpression patternを示す。


library(DESeq2)


# Step 08で保存したjoint fitted DESeq2 object
dds_joint <- readRDS(
  "results/objects/08_DESeq2_joint_fitted.rds"
)


# Step 09と同じVSTを再現
vsd_joint <- vst(
  dds_joint,
  blind = FALSE
)

vst_matrix <- assay(vsd_joint)


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

# Highest-bicor exonic candidate
top_exonic <- overlap_candidates[
  overlap_candidates$genomic_context == "exonic",
]

top_exonic <- top_exonic[
  which.max(top_exonic$bicor),
]


# Highest-bicor intronic candidate
top_intronic <- overlap_candidates[
  overlap_candidates$genomic_context == "intronic",
]

top_intronic <- top_intronic[
  which.max(top_intronic$bicor),
]


cat("\nTop exonic candidate:\n")
print(
  top_exonic[
    ,
    c(
      "transcript_id",
      "gene_symbol",
      "gene_id",
      "bicor"
    )
  ]
)


cat("\nTop intronic candidate:\n")
print(
  top_intronic[
    ,
    c(
      "transcript_id",
      "gene_symbol",
      "gene_id",
      "bicor"
    )
  ]
)

# 12. Exonic candidate scatter ----

exonic_df <- data.frame(
  erv_expression =
    vst_matrix[
      top_exonic$transcript_id,
    ],
  
  gene_expression =
    vst_matrix[
      top_exonic$gene_id,
    ]
)


p_exonic_scatter <- ggplot(
  exonic_df,
  aes(
    x = erv_expression,
    y = gene_expression
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.75
  ) +
  
  labs(
    title = paste0(
      top_exonic$transcript_id,
      " – ",
      top_exonic$gene_symbol
    ),
    
    subtitle = paste0(
      "Exonic overlap; bicor = ",
      sprintf("%.3f", top_exonic$bicor)
    ),
    
    x = paste0(
      top_exonic$transcript_id,
      " VST expression"
    ),
    
    y = paste0(
      top_exonic$gene_symbol,
      " VST expression"
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


print(p_exonic_scatter)

# 13. Intronic candidate scatter ----

intronic_df <- data.frame(
  erv_expression =
    vst_matrix[
      top_intronic$transcript_id,
    ],
  
  gene_expression =
    vst_matrix[
      top_intronic$gene_id,
    ]
)


p_intronic_scatter <- ggplot(
  intronic_df,
  aes(
    x = erv_expression,
    y = gene_expression
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.75
  ) +
  
  labs(
    title = paste0(
      top_intronic$transcript_id,
      " – ",
      top_intronic$gene_symbol
    ),
    
    subtitle = paste0(
      "Intronic overlap; bicor = ",
      sprintf("%.3f", top_intronic$bicor)
    ),
    
    x = paste0(
      top_intronic$transcript_id,
      " VST expression"
    ),
    
    y = paste0(
      top_intronic$gene_symbol,
      " VST expression"
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    )
  )


print(p_intronic_scatter)

# 14. Save representative candidate scatter plots ----

ggsave(
  "results/figures/11_HERVH_5q32b_SCGB3A2_scatter.pdf",
  p_exonic_scatter,
  width = 5,
  height = 4.5
)

ggsave(
  "results/figures/11_HERVL_18p11.21a_PIEZO2_scatter.pdf",
  p_intronic_scatter,
  width = 5,
  height = 4.5
)



# 15. Scatter plots for all 9 primary ERV-gene candidates ----
#
# 9 primary candidatesすべてについて、
# 53 samplesにおけるERV locus expressionと
# paired-gene expressionの共発現patternを可視化する。
#
# x = ERV locus VST expression
# y = paired gene VST expression
#
# 各panelにはbicorとgenomic contextを表示する。


# Step 08で保存したjoint fitted DESeq2 object
dds_joint <- readRDS(
  "results/objects/08_DESeq2_joint_fitted.rds"
)

# Step 09と同じVST
vsd_joint <- DESeq2::vst(
  dds_joint,
  blind = FALSE
)

vst_matrix <- assay(vsd_joint)

# 16. Build plotting data ----

scatter_list <- lapply(
  seq_len(nrow(overlap_candidates)),
  function(i) {
    
    x <- overlap_candidates[i, ]
    
    data.frame(
      sample_id = colnames(vst_matrix),
      
      erv_expression =
        vst_matrix[
          x$transcript_id,
        ],
      
      gene_expression =
        vst_matrix[
          x$gene_id,
        ],
      
      transcript_id = x$transcript_id,
      gene_symbol = x$gene_symbol,
      bicor = x$bicor,
      genomic_context = x$genomic_context,
      module = x$module,
      
      stringsAsFactors = FALSE
    )
  }
)


scatter_df <- do.call(
  rbind,
  scatter_list
)


cat("\n===== Primary candidate scatter data =====\n")

cat(
  "Rows:",
  nrow(scatter_df),
  "\n"
)

cat(
  "Candidates:",
  length(unique(scatter_df$transcript_id)),
  "\n"
)


# 17. Candidate labels ----

scatter_df$panel_label <- paste0(
  scatter_df$transcript_id,
  " – ",
  scatter_df$gene_symbol,
  "\n",
  scatter_df$genomic_context,
  "; bicor = ",
  sprintf("%.3f", scatter_df$bicor)
)


# bicorが高い順に並べる
candidate_order <- overlap_candidates$transcript_id[
  order(
    overlap_candidates$bicor,
    decreasing = TRUE
  )
]


panel_order <- sapply(
  candidate_order,
  function(id) {
    unique(
      scatter_df$panel_label[
        scatter_df$transcript_id == id
      ]
    )
  }
)


scatter_df$panel_label <- factor(
  scatter_df$panel_label,
  levels = panel_order
)



# 18. 3 × 3 scatter plot --------------------------------------------------

p_candidate_scatter <- ggplot(
  scatter_df,
  aes(
    x = erv_expression,
    y = gene_expression
  )
) +
  
  geom_point(
    alpha = 0.7,
    size = 1.7
  ) +
  
  facet_wrap(
    ~ panel_label,
    ncol = 3,
    scales = "free"
  ) +
  
  labs(
    title = "Expression relationships of primary ERV–gene candidates",
    x = "ERV locus expression (VST)",
    y = "Gene expression (VST)"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    
    panel.spacing = unit(
      1,
      "lines"
    )
  )


print(p_candidate_scatter)

# 19. Save 9-candidate scatter figure ----

ggsave(
  "results/figures/11_Telescope_primary_candidates_scatter.pdf",
  p_candidate_scatter,
  width = 10,
  height = 9
)

ggsave(
  "results/figures/11_Telescope_primary_candidates_scatter.png",
  p_candidate_scatter,
  width = 10,
  height = 9,
  dpi = 300
)



# 20. Telescope locus composition by parent ERV model ----
#
# 11,503 expressed Telescope lociが、
# どのparent ERV model (intModel) に属するかを集計する。
#
# Turquoise / Yellowは、
# そのlocus自体をWGCNAした結果ではなく、
# parent ERV subfamilyが所属したWGCNA1 moduleを示す。


# intModel × module のlocus数
intmodel_counts <- as.data.frame(
  table(
    locus_age_results$intModel,
    locus_age_results$module
  )
)

colnames(intmodel_counts) <- c(
  "intModel",
  "module",
  "n_loci"
)


# 0 locusの組み合わせを除く
intmodel_counts <- intmodel_counts[
  intmodel_counts$n_loci > 0,
]


# intModelごとのtotal locus数
intmodel_totals <- aggregate(
  n_loci ~ intModel,
  data = intmodel_counts,
  FUN = sum
)

colnames(intmodel_totals)[2] <- "total_loci"


# total locus数を付与
intmodel_counts <- merge(
  intmodel_counts,
  intmodel_totals,
  by = "intModel"
)


# QC
cat("\n===== Telescope intModel composition =====\n")

cat(
  "Total loci:",
  sum(intmodel_counts$n_loci),
  "\n"
)

cat(
  "Unique intModels:",
  length(unique(intmodel_counts$intModel)),
  "\n"
)

cat("\nModule totals:\n")
print(
  aggregate(
    n_loci ~ module,
    data = intmodel_counts,
    FUN = sum
  )
)



# 21. Plot all 60 parent ERV models ----

# total locus数が多い順
intmodel_order <- intmodel_totals$intModel[
  order(
    intmodel_totals$total_loci,
    decreasing = FALSE
  )
]

intmodel_counts$intModel <- factor(
  intmodel_counts$intModel,
  levels = intmodel_order
)


p_intmodel_composition <- ggplot(
  intmodel_counts,
  aes(
    x = n_loci,
    y = intModel,
    fill = module
  )
) +
  
  geom_col(
    width = 0.8
  ) +
  
  scale_fill_manual(
    values = c(
      turquoise = "turquoise3",
      yellow = "gold"
    )
  ) +
  
  labs(
    title = "Telescope-resolved ERV loci by subfamily",
    subtitle = "11,503 expressed loci across 60 ERV subfamilies",
    x = "Number of ERV loci",
    y = "ERV subfamily",
    fill = "WGCNA module"
  )
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    
    legend.position = "top",
    
    axis.text.y = element_text(
      size = 7
    )
  )


print(p_intmodel_composition)


# 22. Save intModel composition figure ----

ggsave(
  "results/figures/11_Telescope_intModel_composition.pdf",
  p_intmodel_composition,
  width = 8,
  height = 12
)

ggsave(
  "results/figures/11_Telescope_intModel_composition.png",
  p_intmodel_composition,
  width = 8,
  height = 12,
  dpi = 300
)





# 23. Representative genomic context of primary ERV-gene candidates ----
#
# exonic / intronic candidateのうちbicor最大のpairを1つずつ選択し、
# GENCODE v26 gene structureとTelescope ERV locusを重ねて表示する。
#
# Exon:
#   paired geneにannotateされた全exonをgene-levelでunionして表示。
#
# 注意：
#   genomic overlap + coexpressionはregulationの証明ではない。


library(GenomicRanges)
library(rtracklayer)
library(ggplot2)


# 1. Select representative candidates 

top_exonic <- overlap_candidates[
  overlap_candidates$genomic_context == "exonic",
]

top_exonic <- top_exonic[
  which.max(top_exonic$bicor),
]


top_intronic <- overlap_candidates[
  overlap_candidates$genomic_context == "intronic",
]

top_intronic <- top_intronic[
  which.max(top_intronic$bicor),
]


representative_candidates <- rbind(
  top_exonic,
  top_intronic
)


print(
  representative_candidates[
    ,
    c(
      "transcript_id",
      "gene_symbol",
      "gene_id",
      "bicor",
      "genomic_context",
      "erv_chr",
      "erv_start",
      "erv_end"
    )
  ]
)

# 2. Load GENCODE v26 annotation 

gencode_gtf <- rtracklayer::import(
  "03_reference/gencode_v26/gencode.v26.primary_assembly.annotation.gtf"
)


# gene rows
gene_gr <- gencode_gtf[
  gencode_gtf$type == "gene"
]


# exon rows
exon_gr <- gencode_gtf[
  gencode_gtf$type == "exon"
]

# 3. Build gene structure data 

gene_plot_list <- list()
exon_plot_list <- list()
erv_plot_list <- list()


for (i in seq_len(nrow(representative_candidates))) {
  
  candidate <- representative_candidates[i, ]
  
  candidate_gene_id <- candidate$gene_id
  
  
  # paired gene
  this_gene <- gene_gr[
    gene_gr$gene_id == candidate_gene_id
  ]
  
  
  # paired gene exons
  this_exons <- exon_gr[
    exon_gr$gene_id == candidate_gene_id
  ]
  
  
  if (length(this_gene) != 1) {
    stop(
      paste(
        "Expected exactly one gene row for",
        candidate_gene_id,
        "but found",
        length(this_gene)
      )
    )
  }
  
  
  if (length(this_exons) == 0) {
    stop(
      paste(
        "No exons found for",
        candidate_gene_id
      )
    )
  }
  
  
  # transcript間で重複するexonをunion
  this_exons_reduced <- reduce(
    this_exons,
    ignore.strand = TRUE
  )
  
  
  panel_label <- paste0(
    candidate$transcript_id,
    " – ",
    candidate$gene_symbol,
    "\n",
    candidate$genomic_context,
    "; bicor = ",
    sprintf("%.3f", candidate$bicor)
  )
  
  
  # gene body
  gene_plot_list[[i]] <- data.frame(
    panel = panel_label,
    chr = as.character(seqnames(this_gene)),
    start = start(this_gene),
    end = end(this_gene),
    strand = as.character(strand(this_gene)),
    gene_symbol = candidate$gene_symbol
  )
  
  
  # exons
  exon_plot_list[[i]] <- data.frame(
    panel = panel_label,
    start = start(this_exons_reduced),
    end = end(this_exons_reduced)
  )
  
  
  # ERV
  erv_plot_list[[i]] <- data.frame(
    panel = panel_label,
    start = candidate$erv_start,
    end = candidate$erv_end,
    midpoint = (
      candidate$erv_start +
        candidate$erv_end
    ) / 2,
    transcript_id = candidate$transcript_id
  )
}


gene_plot_df <- do.call(
  rbind,
  gene_plot_list
)

exon_plot_df <- do.call(
  rbind,
  exon_plot_list
)

erv_plot_df <- do.call(
  rbind,
  erv_plot_list
)

# 4. Plot genomic context 

p_genomic_context <- ggplot() +
  
  # gene body
  geom_segment(
    data = gene_plot_df,
    aes(
      x = start,
      xend = end,
      y = 1,
      yend = 1
    ),
    linewidth = 0.7
  ) +
  
  # exons
  geom_rect(
    data = exon_plot_df,
    aes(
      xmin = start,
      xmax = end,
      ymin = 0.85,
      ymax = 1.15
    )
  ) +
  
  # actual ERV interval
  geom_rect(
    data = erv_plot_df,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1.45,
      ymax = 1.75
    ),
    fill = "firebrick",
    alpha = 0.8
  ) +
  
  # ERV midpoint
  geom_vline(
    data = erv_plot_df,
    aes(
      xintercept = midpoint
    ),
    colour = "firebrick",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # gene label
  geom_text(
    data = gene_plot_df,
    aes(
      x = (start + end) / 2,
      y = 0.55,
      label = gene_symbol
    ),
    fontface = "bold",
    size = 4
  ) +
  
  # ERV label
  geom_text(
    data = erv_plot_df,
    aes(
      x = midpoint,
      y = 1.95,
      label = transcript_id
    ),
    colour = "firebrick",
    fontface = "bold",
    size = 3.5
  ) +
  
  facet_wrap(
    ~ panel,
    ncol = 1,
    scales = "free_x"
  ) +
  
  scale_y_continuous(
    breaks = c(
      1,
      1.6
    ),
    labels = c(
      "Gene / exons",
      "ERV locus"
    ),
    limits = c(
      0.3,
      2.15
    )
  ) +
  
  scale_x_continuous(
    labels = scales::label_number(
      big.mark = ",",
      accuracy = 1
    )
  ) +
  
  labs(
    title = "Genomic context of representative ERV–gene candidates",
    x = "Genomic coordinate (bp)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    panel.spacing = unit(
      1.3,
      "lines"
    )
  )


print(p_genomic_context)


# 5. Save genomic context figure 

ggsave(
  "results/figures/11_Telescope_representative_genomic_context.pdf",
  p_genomic_context,
  width = 10,
  height = 6
)

ggsave(
  "results/figures/11_Telescope_representative_genomic_context.png",
  p_genomic_context,
  width = 10,
  height = 6,
  dpi = 300
)

