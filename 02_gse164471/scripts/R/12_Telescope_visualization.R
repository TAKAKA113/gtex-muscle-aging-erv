# 12_Telescope_visualization.R
#
# Purpose:
# Visualize all 9 primary Telescope-resolved ERV–gene candidates.
#
# Each panel shows:
# - paired gene body
# - annotated exons from GENCODE v26
# - Telescope-resolved ERV locus
# - exonic / intronic genomic context
# - bicor coexpression value
#
# Important:
# Genomic overlap + coexpression supports a candidate relationship,
# but does not demonstrate causal regulation.


library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

setwd(
  "/rds/projects/z/zhoujz-gnn-chem-mixture/gse164471_erv_aging"
)


# 1. Load primary ERV-gene candidates ----

overlap_candidates <- read.delim(
  "results/tables/10_primary_overlapping_ERV_gene_candidates.tsv"
)


cat("\n===== Primary ERV-gene candidates =====\n")

cat(
  "Candidate pairs:",
  nrow(overlap_candidates),
  "\n"
)

cat(
  "Exonic:",
  sum(overlap_candidates$genomic_context == "exonic"),
  "\n"
)

cat(
  "Intronic:",
  sum(overlap_candidates$genomic_context == "intronic"),
  "\n"
)



# 2. Load GENCODE v26 ----

gencode_gtf <- rtracklayer::import(
  "03_reference/gencode_v26/gencode.v26.primary_assembly.annotation.gtf"
)


gene_gr <- gencode_gtf[
  gencode_gtf$type == "gene"
]


exon_gr <- gencode_gtf[
  gencode_gtf$type == "exon"
]


# 3. Build genomic plotting data ----

gene_plot_list <- list()
exon_plot_list <- list()
erv_plot_list <- list()


for (i in seq_len(nrow(overlap_candidates))) {
  
  candidate <- overlap_candidates[i, ]
  
  candidate_gene_id <- candidate$gene_id
  
  
  # paired gene
  this_gene <- gene_gr[
    gene_gr$gene_id == candidate_gene_id
  ]
  
  
  # exons belonging to the paired gene
  this_exons <- exon_gr[
    exon_gr$gene_id == candidate_gene_id
  ]
  
  
  if (length(this_gene) != 1) {
    
    stop(
      paste(
        "Expected one gene row for",
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
  
  
  # Merge overlapping exons across annotated transcripts
  this_exons_reduced <- reduce(
    this_exons,
    ignore.strand = TRUE
  )
  
  
  # Panel title
  panel_label <- paste0(
    candidate$transcript_id,
    " – ",
    candidate$gene_symbol,
    "\n",
    candidate$genomic_context,
    "; bicor = ",
    sprintf("%.3f", candidate$bicor)
  )
  
  
  # Gene body
  gene_plot_list[[i]] <- data.frame(
    panel = panel_label,
    chr = as.character(seqnames(this_gene)),
    start = start(this_gene),
    end = end(this_gene),
    strand = as.character(strand(this_gene)),
    gene_symbol = candidate$gene_symbol,
    genomic_context = candidate$genomic_context,
    bicor = candidate$bicor
  )
  
  
  # Exons
  exon_plot_list[[i]] <- data.frame(
    panel = panel_label,
    start = start(this_exons_reduced),
    end = end(this_exons_reduced)
  )
  
  
  # ERV locus
  erv_plot_list[[i]] <- data.frame(
    panel = panel_label,
    start = candidate$erv_start,
    end = candidate$erv_end,
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


# 4. Order panels ----

candidate_order_df <- overlap_candidates[
  order(
    overlap_candidates$genomic_context,
    -overlap_candidates$bicor
  ),
]


candidate_order_df$panel <- paste0(
  candidate_order_df$transcript_id,
  " – ",
  candidate_order_df$gene_symbol,
  "\n",
  candidate_order_df$genomic_context,
  "; bicor = ",
  sprintf("%.3f", candidate_order_df$bicor)
)


panel_levels <- candidate_order_df$panel


gene_plot_df$panel <- factor(
  gene_plot_df$panel,
  levels = panel_levels
)

exon_plot_df$panel <- factor(
  exon_plot_df$panel,
  levels = panel_levels
)

erv_plot_df$panel <- factor(
  erv_plot_df$panel,
  levels = panel_levels
)


# 5. Plot all 9 genomic contexts ----

p_all_genomic_context <- ggplot() +
  
  # Gene body
  geom_segment(
    data = gene_plot_df,
    aes(
      x = start,
      xend = end,
      y = 1,
      yend = 1
    ),
    linewidth = 0.6,
    colour = "black"
  ) +
  
  # Exons
  geom_rect(
    data = exon_plot_df,
    aes(
      xmin = start,
      xmax = end,
      ymin = 0.88,
      ymax = 1.12
    ),
    fill = "grey35"
  ) +
  
  # ERV genomic interval
  geom_rect(
    data = erv_plot_df,
    aes(
      xmin = start,
      xmax = end,
      ymin = 1.45,
      ymax = 1.72
    ),
    fill = "firebrick",
    alpha = 0.85
  ) +
  
  # Gene symbol
  geom_text(
    data = gene_plot_df,
    aes(
      x = (start + end) / 2,
      y = 0.58,
      label = gene_symbol
    ),
    fontface = "bold",
    size = 3.2
  ) +
  
  # ERV label
  geom_text(
    data = erv_plot_df,
    aes(
      x = (start + end) / 2,
      y = 1.92,
      label = transcript_id
    ),
    colour = "firebrick",
    fontface = "bold",
    size = 2.8
  ) +
  
  facet_wrap(
    ~ panel,
    ncol = 3,
    scales = "free_x"
  ) +
  
  scale_y_continuous(
    breaks = c(
      1,
      1.58
    ),
    labels = c(
      "Gene / exons",
      "ERV locus"
    ),
    limits = c(
      0.35,
      2.1
    )
  ) +
  
  scale_x_continuous(
    labels = scales::label_number(
      big.mark = ",",
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0.03,
        0.03
      )
    )
  ) +
  
  labs(
    title = "Genomic context of primary Telescope ERV–gene candidates",
    subtitle = "Nine coexpressed ERV–gene pairs with direct gene-body overlap",
    x = "Genomic coordinate (bp)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 8
    ),
    
    axis.text.x = element_text(
      size = 7
    ),
    
    axis.text.y = element_text(
      size = 7
    ),
    
    panel.spacing = unit(
      1.1,
      "lines"
    )
  )


print(p_all_genomic_context)


# 6. Save figure ----

ggsave(
  "results/figures/12_Telescope_all9_genomic_context.pdf",
  p_all_genomic_context,
  width = 13,
  height = 10
)

ggsave(
  "results/figures/12_Telescope_all9_genomic_context.png",
  p_all_genomic_context,
  width = 13,
  height = 10,
  dpi = 300
)


# 7. Build final candidate summary table ----

candidate_summary <- overlap_candidates[
  order(-overlap_candidates$bicor),
  c(
    "transcript_id",
    "intModel",
    "gene_symbol",
    "module",
    "bicor",
    "genomic_context"
  )
]


# Rename columns for manuscript presentation
names(candidate_summary) <- c(
  "ERV_locus",
  "ERV_subfamily",
  "Host_gene",
  "WGCNA_module",
  "bicor",
  "Genomic_context"
)


# Presentation formatting
candidate_summary$Genomic_context <- ifelse(
  candidate_summary$Genomic_context == "exonic",
  "Exonic",
  "Intronic"
)

candidate_summary$bicor <- round(
  candidate_summary$bicor,
  3
)


# Confirm candidate number
stopifnot(nrow(candidate_summary) == 9)


cat("\n===== Final 9 ERV-gene candidate associations =====\n")

print(candidate_summary)


# ERV_locus   ERV_subfamily Host_gene WGCNA_module bicor Genomic_context
# 4000        HERVH_5q32b       HERVH-int   SCGB3A2    turquoise 0.944          Exonic
# 8642    HERVL_18p11.21a       HERVL-int    PIEZO2    turquoise 0.922        Intronic
# 13904     HERVH_5q35.1a       HERVH-int     GABRP    turquoise 0.914          Exonic
# 15966    HERVH_19p13.3b       HERVH-int      FUT3    turquoise 0.912          Exonic
# 28468 ERV316A3_8p11.22b ERV3-16A3_I-int    ADAM18    turquoise 0.904        Intronic
# 28579     HERVL_3q13.31       HERVL-int     LSAMP    turquoise 0.904        Intronic
# 30472    HERVH_15q26.3b       HERVH-int     LRRK1    turquoise 0.903          Exonic
# 46345  HERVL40_12q21.31     HERVL40-int    LRRIQ1    turquoise 0.896        Intronic
# 403       PABLA_7q21.13      PABL_A-int    CFAP69       yellow 0.873        Intronic

# 8. Save final candidate summary table ----

write.table(
  candidate_summary,
  "results/tables/12_Telescope_primary_ERV_gene_candidates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

