#この段階は基本コピペでいい
##Step1
#goodSamplesGenes() でQC
#欠損が多い遺伝子や、異常なサンプルがないかを機械的に確認
library(WGCNA)

options(stringsAsFactors = FALSE)

# WGCNA用に作った発現データを読み込む
# datExpr は「行 = サンプル、列 = 遺伝子」の形になっている
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")

# データの大きさを確認
# 期待される形：818 samples × 20000 genes
dim(datExpr)


# WGCNAに入れる前の最低限チェック
# goodSamplesGenes() は、WGCNAに使えないサンプルや遺伝子がないか確認する関数
gsg <- goodSamplesGenes(datExpr, verbose = 3)


# TRUEなら、全サンプル・全遺伝子がWGCNAに使える状態
# FALSEなら、何か問題のあるサンプルまたは遺伝子がある
gsg$allOK


# 問題のある遺伝子が何個あるか確認
sum(!gsg$goodGenes)

# 問題のあるサンプルが何個あるか確認
sum(!gsg$goodSamples)


# もし問題があった場合、どの遺伝子・サンプルが問題か表示する
if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) {
    cat("問題のある遺伝子:\n")
    print(colnames(datExpr)[!gsg$goodGenes])
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    cat("問題のあるサンプル:\n")
    print(rownames(datExpr)[!gsg$goodSamples])
  }
}


##Step2
#Sample crustering

# datExprを読み込む
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")

# サンプル同士の発現パターンの違いをもとにクラスタリング
sampleTree <- hclust(dist(datExpr), method = "average")


# PDFとして保存
pdf("results/01_sample_clustering_mad20k_no_labels.pdf", width = 10, height = 8)

plot(sampleTree,
     main = "Sample clustering for WGCNA input QC",
     xlab = "",
     sub = "",
     labels = FALSE)

dev.off()

##Step3
#AGE / SEX との関係を見る

# 読み込み
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")
datTraits <- readRDS("clean/datTraits_v11_muscle_skeletal.rds")

# サンプル順を datExpr に合わせる
datTraits <- datTraits[match(rownames(datExpr), datTraits$SAMPID), ]

# sample IDの順番確認
all(rownames(datExpr) == datTraits$SAMPID)

# sample clustering
sampleTree <- hclust(dist(datExpr), method = "average")

# trait を数値化
sex_num <- ifelse(datTraits$SEX_LABEL == "Male", 0, 1)

age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
age_num <- match(datTraits$AGE, age_levels)

# Ageは白〜赤のグラデーション
age_color <- numbers2colors(age_num, signed = FALSE)

# Sexは手動で色指定
# Male = blue, Female = red
sex_color <- ifelse(datTraits$SEX_LABEL == "Male", "blue", "red")

traitColors <- data.frame(
  Age = age_color,
  Sex = sex_color
)

# 保存先フォルダ
dir.create("results", showWarnings = FALSE)

# 保存
pdf("results/02_sample_clustering_traits_mad20k_no_labels.pdf", width = 12, height = 8)

plotDendroAndColors(
  sampleTree,
  traitColors,
  groupLabels = c("Age", "Sex"),
  main = "Sample clustering with trait heatmap",
  dendroLabels = FALSE
)

dev.off()


##Step4
#Pick β based on SoftThreshold
#Soft-thresholding powerはWGCNAでgene同士の相関をネットワークの接続強度に変換するためのパラメータ
#power = 相関をどれくらい厳しくネットワーク接続に変換するかを決める数字
#Scale independence=その power を使って作るネットワークが、どれくらい “scale-free network っぽい形” になっているか
#WGCNAでは「生物学的ネットワークはscale-freeっぽい構造を持つことが多い」という前提があり
#「少数の遺伝子がたくさんつながり、多くの遺伝子は少数しかつながらない」ほうが良い
#R^2=「そのネットワークが scale-free topology にどれだけうまく当てはまっているか」を数値化
options(stringsAsFactors = FALSE)

# datExpr を読み込む
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")

# power 候補
powers <- c(1:10, 12, 14, 16, 18, 20)

# soft-threshold を評価
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

# 結果確認
sft$fitIndices

# 図を保存
pdf("results/02_soft_threshold_mad20k.pdf", width = 12, height = 6)

par(mfrow = c(1, 2))

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers,
     col = "red")
abline(h = 0.8, col = "blue", lty = 2)

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     col = "red")

dev.off()　　　　#power = 6

# 保存
saveRDS(sft, file = "clean/sft_mad20k.rds")

##Step5
#blockwiseModulesでmoduleと作成
#power = 6
#発現パターンが似ている遺伝子同士をまとめて、moduleとして分類する


options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# datExpr を読み込む
# datExpr は WGCNA 用の発現行列
# 行 = sample、列 = gene
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")

dim(datExpr)


# module detection
# blockwiseModules() は WGCNA の本体
# 発現パターンが似ている gene をまとめて module に分類する
#
# power = 6:
# pickSoftThreshold() で選んだ soft-threshold power
#
# networkType = "signed":
# 正の相関を重視する設定
# 一緒に上がる / 一緒に下がる gene を同じ module にまとめやすい
#
# TOMType = "signed":
# signed network に合わせた TOM を使う
# TOM は単純な相関だけでなく、network上での近さも考慮する
#
# minModuleSize = 30:
# 30 genes 未満の小さすぎる module は作らない
#
# mergeCutHeight = 0.25:
# 似ている module を統合する基準
# module eigengene の相関が高い module 同士は merge される
#
# numericLabels = FALSE:
# module名を数字ではなく色で表示する
# 例：blue, brown, turquoise
#
# saveTOMs = FALSE:
# TOM matrix は非常に大きいので保存しない

net <- blockwiseModules(
  datExpr,
  power = 6,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)


# module の数とサイズを確認
# 各moduleに何個のgeneが入ったかを見る
table(net$colors)


# 重要な結果を取り出す
# moduleColors: 各geneがどのmoduleに入ったか
# MEs: module eigengene。各moduleの代表的な発現パターン
moduleColors <- net$colors
MEs <- net$MEs
MEs <- orderMEs(MEs)


# gene と module の対応表を作る
# 後で hub gene を見るときに使う
geneModuleTable <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors
)

head(geneModuleTable)


# 保存
# net はWGCNA全体の結果
# moduleColors はgeneごとのmodule情報
# MEs はmodule eigengene
# geneModuleTable はgene-module対応表

saveRDS(net, file = "clean/WGCNA_net_mad20k_power6.rds")
saveRDS(moduleColors, file = "clean/moduleColors_mad20k_power6.rds")
saveRDS(MEs, file = "clean/moduleEigengenes_mad20k_power6.rds")
saveRDS(geneModuleTable, file = "clean/geneModuleTable_mad20k_power6.rds")

write.csv(geneModuleTable,
          file = "results/gene_module_assignment_mad20k_power6.csv",
          row.names = FALSE)


# check
dim(datExpr)
table(moduleColors)
dim(MEs)
head(geneModuleTable)
list.files("clean", pattern = "mad20k_power6")
list.files("results", pattern = "mad20k_power6")


##Step6
#module eigengene と trait の相関

options(stringsAsFactors = FALSE)

# 読み込み
# MEs は module eigengene
# 各moduleの代表的な発現パターン
MEs <- readRDS("clean/moduleEigengenes_mad20k_power6.rds")

# datTraits は sample metadata
datTraits <- readRDS("clean/datTraits_v11_muscle_skeletal.rds")

# datExpr は WGCNAに使った発現行列
# 行 = samples, 列 = genes
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")


# サンプル順を datExpr に揃える
datTraits <- datTraits[match(rownames(datExpr), datTraits$SAMPID), ]

# sample ID の順番確認
all(rownames(datExpr) == datTraits$SAMPID)

# MEs のsample順も確認
all(rownames(MEs) == rownames(datExpr))


# trait を数値化
# SEX_LABEL: Male = 0, Female = 1
sex_num <- ifelse(datTraits$SEX_LABEL == "Male", 0, 1)

# AGE は "20-29", "30-39" のような年齢階級なので、
# 若い順に 1,2,3,4,5,6 として数値化する
age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
age_num <- match(datTraits$AGE, age_levels)

# module eigengene と相関させる trait table を作る
traitData <- data.frame(
  AGE = age_num,
  SEX = sex_num
)

rownames(traitData) <- datTraits$SAMPID


# module-trait correlation
# 各module eigengene と AGE / SEX の相関を計算する
moduleTraitCor <- cor(MEs, traitData, use = "p")

# 相関のp-valueを計算する
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))


# 確認
moduleTraitCor
moduleTraitPvalue


# heatmapの中に表示する文字を作る
# 上段 = correlation
# 下段 = p-value
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(moduleTraitCor)



# heatmapをPDFとして保存
# module名を短くする
# 例：MEblue → blue
moduleLabels <- gsub("^ME", "", names(MEs))

pdf("results/04_module_trait_relationships_mad20k_power6_fixed.pdf",
    width = 10,
    height = 12)

# 余白を調整
# c(下, 左, 上, 右)
par(mar = c(8, 12, 4, 6))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),
  yLabels = moduleLabels,
  ySymbols = moduleLabels,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.65,
  cex.lab = 0.8,
  xLabelsAngle = 45,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)

dev.off()

# 保存
saveRDS(moduleTraitCor,
        file = "clean/moduleTraitCor_AGE_SEX_mad20k_power6.rds")

saveRDS(moduleTraitPvalue,
        file = "clean/moduleTraitPvalue_AGE_SEX_mad20k_power6.rds")

write.csv(moduleTraitCor,
          file = "results/moduleTraitCor_AGE_SEX_mad20k_power6.csv")

write.csv(moduleTraitPvalue,
          file = "results/moduleTraitPvalue_AGE_SEX_mad20k_power6.csv")


#Step7
#Hub Gene


# 読み込み
# datExpr: WGCNAに使った発現行列
# 行 = sample, 列 = gene
datExpr <- readRDS("clean/datExpr_v11_muscle_skeletal_MAD20k.rds")

# MEs: module eigengene
# 各moduleの代表的な発現パターン
MEs <- readRDS("clean/moduleEigengenes_mad20k_power6.rds")
MEs <- orderMEs(MEs)

# datTraits: sample metadata
datTraits <- readRDS("clean/datTraits_v11_muscle_skeletal.rds")

# moduleColors: 各geneがどのmoduleに属するか
moduleColors <- readRDS("clean/moduleColors_mad20k_power6.rds")


# datTraits の順番を datExpr に合わせる
datTraits <- datTraits[match(rownames(datExpr), datTraits$SAMPID), ]

# sample順が正しいか確認
all(rownames(datExpr) == datTraits$SAMPID)
all(rownames(MEs) == rownames(datExpr))


# AGE を数値化
# GTExのAGEは "20-29", "30-39" のような年齢階級なので、
# 若い順に 1,2,3,4,5,6 として数値化する
age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79")
age_num <- match(datTraits$AGE, age_levels)


# MM: module membership
# 各geneの発現パターンが、各module eigengeneとどれくらい似ているかを見る
# MMが高いgeneは、そのmoduleの中心的なgene、つまりhub gene候補になりやすい
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(geneModuleMembership),
    nSamples = nrow(datExpr)
  )
)

# 列名を分かりやすくする
# 例：MM_MEblue, MM_MEred
names(geneModuleMembership) <- paste0("MM_", names(MEs))
names(MMPvalue) <- paste0("p.MM_", names(MEs))


# GS: gene significance for AGE
# 各geneの発現量がAGEとどれくらい相関するかを見る
# GS_AGEが高いgeneは、年齢と関連するgeneと考えられる
geneTraitSignificance <- as.data.frame(cor(datExpr, age_num, use = "p"))

GSPvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(geneTraitSignificance),
    nSamples = nrow(datExpr)
  )
)

names(geneTraitSignificance) <- "GS_AGE"
names(GSPvalue) <- "p.GS_AGE"


# geneごとの統合表を作る
# 各geneについて、
# gene ID, module, MM, MM p-value, GS_AGE, GS p-value
# をまとめる
allGeneInfo <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  geneModuleMembership,
  MMPvalue,
  geneTraitSignificance,
  GSPvalue,
  stringsAsFactors = FALSE
)

head(allGeneInfo)


# 保存先フォルダ
dir.create("clean", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)


# 全gene情報を保存
# 後でhub gene抽出やannotationに使う
saveRDS(allGeneInfo, file = "clean/allGeneInfo_mad20k_power6.rds")

write.csv(
  allGeneInfo,
  file = "results/allGeneInfo_mad20k_power6.csv",
  row.names = FALSE
)


# AGE関連module候補を指定
# AGEとのmodule-trait correlationで比較的目立っていたmodule
# greyは通常解釈対象外なので含めない
target_modules <- c("greenyellow", "red", "darkgreen", "turquoise",
                    "midnightblue", "blue", "pink")


# 各moduleごとにhub gene候補を抽出する
# 基本方針：
# 1. そのmoduleに属するgeneだけを取り出す
# 2. MMの絶対値が高い順に並べる
# 3. さらにGS_AGEの絶対値も参考にする
# 4. 上位30 genesをhub gene候補として保存する

all_top_hubs <- list()

for (mod in target_modules) {
  
  me_col <- paste0("MM_ME", mod)
  pme_col <- paste0("p.MM_ME", mod)
  
  # 指定moduleに属するgeneだけを取り出す
  module_df <- allGeneInfo[allGeneInfo$module == mod, ]
  
  # 念のため、moduleが存在するか確認
  if (nrow(module_df) == 0) {
    cat("\nModule", mod, "has no genes. Skipped.\n")
    next
  }
  
  # MMが高く、かつGS_AGEも高いgeneを上に並べる
  module_df <- module_df[
    order(
      abs(module_df[[me_col]]),
      abs(module_df$GS_AGE),
      decreasing = TRUE
    ),
  ]
  
  # 上位30 geneをhub gene候補として抽出
  top_hub <- head(module_df, 30)
  
  # module名を追加して、後でまとめやすくする
  top_hub$target_module <- mod
  
  all_top_hubs[[mod]] <- top_hub
  
  # moduleごとの全候補を保存
  write.csv(
    module_df,
    file = paste0("results/hub_candidates_", mod, "_full_mad20k_power6.csv"),
    row.names = FALSE
  )
  
  # moduleごとのtop30を保存
  write.csv(
    top_hub,
    file = paste0("results/hub_candidates_", mod, "_top30_mad20k_power6.csv"),
    row.names = FALSE
  )
  
  # R consoleに上位30を表示
  cat("\n====================\n")
  cat("Module:", mod, "\n")
  cat("Number of genes in module:", nrow(module_df), "\n")
  cat("====================\n")
  
  print(
    top_hub[, c("gene", "module", me_col, pme_col, "GS_AGE", "p.GS_AGE")]
  )
}


# 全target moduleのtop30 hub genesを1つの表にまとめる
all_top_hubs_df <- do.call(rbind, all_top_hubs)

# 保存
saveRDS(
  all_top_hubs_df,
  file = "clean/top30_hub_candidates_target_modules_mad20k_power6.rds"
)

write.csv(
  all_top_hubs_df,
  file = "results/top30_hub_candidates_target_modules_mad20k_power6.csv",
  row.names = FALSE
)


# check
dim(allGeneInfo)
table(allGeneInfo$module)
target_modules
dim(all_top_hubs_df)
head(all_top_hubs_df)
list.files("results", pattern = "hub_candidates")



#tableを綺麗に整形
# hub gene候補を見やすい形に整理する

hub_summary_list <- list()

for (mod in target_modules) {
  
  me_col <- paste0("MM_ME", mod)
  pme_col <- paste0("p.MM_ME", mod)
  
  module_df <- all_top_hubs_df[all_top_hubs_df$target_module == mod, ]
  
  hub_summary <- data.frame(
    gene = module_df$gene,
    module = module_df$module,
    MM = module_df[[me_col]],
    p_MM = module_df[[pme_col]],
    GS_AGE = module_df$GS_AGE,
    p_GS_AGE = module_df$p.GS_AGE,
    stringsAsFactors = FALSE
  )
  
  hub_summary$abs_MM <- abs(hub_summary$MM)
  hub_summary$abs_GS_AGE <- abs(hub_summary$GS_AGE)
  hub_summary$rank_in_module <- 1:nrow(hub_summary)
  
  hub_summary_list[[mod]] <- hub_summary
}

hub_summary_all <- do.call(rbind, hub_summary_list)

# 確認
head(hub_summary_all)
dim(hub_summary_all)

# 保存
write.csv(
  hub_summary_all,
  file = "results/top30_hub_candidates_summary_mad20k_power6.csv",
  row.names = FALSE
)

saveRDS(
  hub_summary_all,
  file = "clean/top30_hub_candidates_summary_mad20k_power6.rds"
)

