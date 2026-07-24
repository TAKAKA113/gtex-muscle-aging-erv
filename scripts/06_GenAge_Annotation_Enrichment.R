#AGE関連候補module
#→ GenAge Humanとの重なりを確認
#→ module全体のenrichment
#→ hub gene候補内の既知老化geneを確認

##Step1
#GenAge Human list を取得して読む

# dataフォルダを作る
dir.create("data", showWarnings = FALSE)

# GenAge Humanをダウンロードする
download.file(
  "https://genomics.senescence.info/genes/human_genes.zip",
  "data/human_genes.zip",
  mode = "wb"
)

# ZIP内のCSVを直接読み込む
genage <- read.csv(
  unz("data/human_genes.zip", "genage_human.csv"),
  stringsAsFactors = FALSE
)

# 読み込み確認
dim(genage)
head(genage)


# GenAge Humanのgene symbolを取り出す
genage_symbols <- genage$symbol

# 内容を確認
length(genage_symbols)
head(genage_symbols)


##Step2
#WGCNAのresultとGenAgeの照合
# WGCNAの各geneがGenAge geneか判定する

#Fileを読む
gene_info <- readRDS(
  "clean/gene_info_v10_muscle_skeletal.rds"
)

geneModuleTable <- readRDS(
  "clean/geneModuleTable_mad20k_power6.rds"
)

hub_summary_all <- readRDS(
  "clean/top30_hub_candidates_summary_mad20k_power6.rds"
)


#WGCNA遺伝子へのgene symbol付与
# WGCNA側のgene列はEnsembl IDのまま
# 例：ENSG00000134184.13
# 末尾のversion番号を外す
# ENSG00000134184.13 → ENSG00000134184
geneModuleTable$gene_clean <- sub(
  "\\..*$", "", geneModuleTable$gene
)

hub_summary_all$gene_clean <- sub(
  "\\..*$", "", hub_summary_all$gene
)

gene_info$gene_clean <- sub(
  "\\..*$", "", gene_info$Name
)


# gene_infoを使ってEnsembl IDにgene symbolを追加！！！
# ENSG00000134184のようなEnsenble IDではGenAgeとmatchできない
#GSTM1のようなSymbol ID列を対応表により追加
geneModuleTable$symbol <- gene_info$Description[
  match(geneModuleTable$gene_clean, gene_info$gene_clean)
]

hub_summary_all$symbol <- gene_info$Description[
  match(hub_summary_all$gene_clean, gene_info$gene_clean)
]


# symbolが付かなかったgene数を確認
# 両方0なら成功
sum(is.na(geneModuleTable$symbol))
sum(is.na(hub_summary_all$symbol))


# 内容確認
head(geneModuleTable[, c("gene", "module", "symbol")])

head(
  hub_summary_all[
    , c("gene", "module", "symbol", "MM", "GS_AGE")
  ]
)



# WGCNA遺伝子とGenAge Humanの照合
# 各gene symbolがGenAgeの307 genesに含まれるか判定
#
# TRUE  = GenAgeに登録されている
# FALSE = GenAgeに登録されていない
geneModuleTable$is_GenAge <-
  geneModuleTable$symbol %in% genage_symbols

hub_summary_all$is_GenAge <-
  hub_summary_all$symbol %in% genage_symbols


# 全20,000 genesのうち、GenAge geneが何個あるか
table(geneModuleTable$is_GenAge)

# hub候補210 genesのうち、GenAge geneが何個あるか
table(hub_summary_all$is_GenAge)


#ModuleごとのGenAge遺伝子の割合
# target modules
target_modules <- c(
  "greenyellow", "red", "darkgreen", "turquoise",
  "midnightblue", "blue", "pink"
)

# moduleごとのGenAge overlap summary
module_overlap_summary <- do.call(
  rbind,
  lapply(target_modules, function(mod) {
    
    # module内の重複しないgene symbol
    module_symbols <- unique(
      geneModuleTable$symbol[
        geneModuleTable$module == mod
      ]
    )
    
    # GenAgeに含まれるsymbol
    genage_hits <- module_symbols[
      module_symbols %in% genage_symbols
    ]
    
    data.frame(
      module = mod,
      module_size = length(module_symbols),
      GenAge_overlap_n = length(genage_hits),
      GenAge_overlap_symbols = paste(genage_hits, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
)

module_overlap_summary

#Step3
#moduleごとのGenAge enrichment解析

# AGE関連候補module
target_modules <- c(
  "greenyellow", "red", "darkgreen", "turquoise",
  "midnightblue", "blue", "pink"
)


# WGCNAに入った全geneを背景とする
# 同じsymbolが重複していた場合は1つとして数える
universe_symbols <- unique(
  geneModuleTable$symbol[
    !is.na(geneModuleTable$symbol) &
      geneModuleTable$symbol != ""
  ]
)

# N = 背景の全gene数
# K = 背景に含まれるGenAge gene数
N <- length(universe_symbols)
K <- sum(universe_symbols %in% genage_symbols)

# 背景を確認
N
K
K / N


# 各moduleについてGenAge enrichmentを検定
module_genage_results <- do.call(
  rbind,
  lapply(target_modules, function(mod) {
    
    # 現在のmoduleに含まれるgene symbol
    module_symbols <- unique(
      geneModuleTable$symbol[
        geneModuleTable$module == mod &
          !is.na(geneModuleTable$symbol) &
          geneModuleTable$symbol != ""
      ]
    )
    
    # module内に含まれるGenAge gene
    genage_hits <- module_symbols[
      module_symbols %in% genage_symbols
    ]
    
    # n = moduleの全gene数
    # k = module内のGenAge gene数
    n <- length(module_symbols)
    k <- length(genage_hits)
    
    # Fisher's exact test用の2×2表
    fisher_table <- matrix(
      c(
        k,                  # module内・GenAge
        n - k,              # module内・非GenAge
        K - k,              # module外・GenAge
        N - K - (n - k)     # module外・非GenAge
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    # 背景よりGenAge geneが多いか検定
    ft <- fisher.test(
      fisher_table,
      alternative = "greater"
    )
    
    # moduleごとの結果を1行にまとめる
    data.frame(
      module = mod,
      module_size = n,
      GenAge_overlap_n = k,
      expected_overlap = n * K / N,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      GenAge_overlap_symbols = paste(
        genage_hits,
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  })
)


# 7 modulesのp値を多重検定補正
module_genage_results$FDR <- p.adjust(
  module_genage_results$p_value,
  method = "BH"
)

# FDRが小さい順に並べる
module_genage_results <- module_genage_results[
  order(
    module_genage_results$FDR,
    module_genage_results$p_value
  ),
]

# 結果を確認
module_genage_results


# 結果を保存
write.csv(
  module_genage_results,
  "results/module_GenAge_enrichment_mad20k_power6.csv",
  row.names = FALSE
)

saveRDS(
  module_genage_results,
  "clean/module_GenAge_enrichment_mad20k_power6.rds"
)


#Figure
library(ggplot2)

module_genage_results <- readRDS(
  "clean/module_GenAge_enrichment_mad20k_power6.rds"
)

# odds ratio順に並べる
module_genage_results$module <- factor(
  module_genage_results$module,
  levels = module_genage_results$module[
    order(module_genage_results$odds_ratio, decreasing = TRUE)
  ]
)

module_genage_results$neglog10p <- -log10(module_genage_results$p_value)

p_or <- ggplot(module_genage_results,
               aes(x = module, y = odds_ratio, fill = neglog10p)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_gradient(
    low = "grey85", high = "firebrick",
    name = "-log10(p)"
  ) +
  labs(
    title = "GenAge enrichment across candidate modules",
    subtitle = "No module reached FDR < 0.05; pink shows a modest trend",
    x = "Module",
    y = "Odds ratio (Fisher's exact test)"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(
  "results/06_GenAge_enrichment_light.pdf",
  p_or, width = 8, height = 5
)


## Step4
# hub gene候補内の既知GenAge geneを確認
# 210個のhub候補から、
# GenAgeに登録されているgeneだけを抽出
hub_genage_only <- hub_summary_all[
  hub_summary_all$is_GenAge,
  c(
    "gene",
    "module",
    "symbol",
    "MM",
    "GS_AGE",
    "rank_in_module"
  )
]

# 結果を確認
hub_genage_only

# gene数を確認
nrow(hub_genage_only)

#保存
write.csv(
  hub_genage_only,
  "results/top30_hub_GenAge_only_mad20k_power6.csv",
  row.names = FALSE
)

saveRDS(
  hub_genage_only,
  "clean/top30_hub_GenAge_only_mad20k_power6.rds"
)


## Step5
# GenAge annotation済みの全データを保存

# WGCNA全20,000 genes
write.csv(
  geneModuleTable,
  "results/geneModuleTable_with_GenAge_mad20k_power6.csv",
  row.names = FALSE
)

saveRDS(
  geneModuleTable,
  "clean/geneModuleTable_with_GenAge_mad20k_power6.rds"
)


# hub候補210 genes
write.csv(
  hub_summary_all,
  "results/hub_summary_with_GenAge_mad20k_power6.csv",
  row.names = FALSE
)

saveRDS(
  hub_summary_all,
  "clean/hub_summary_with_GenAge_mad20k_power6.rds"
)


## Step6
# 06_GenAge解析の最終確認

cat("WGCNA gene数:", nrow(geneModuleTable), "\n")
cat("WGCNA内のGenAge gene数:",
    sum(geneModuleTable$is_GenAge), "\n")

cat("hub候補数:", nrow(hub_summary_all), "\n")
cat("hub候補内のGenAge gene数:",
    sum(hub_summary_all$is_GenAge), "\n")

cat("FDR < 0.05のmodule数:",
    sum(module_genage_results$FDR < 0.05), "\n")

# enrichment結果
module_genage_results

# GenAgeに登録されていたhub候補
hub_genage_only
