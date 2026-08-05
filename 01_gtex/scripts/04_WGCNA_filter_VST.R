#DESeq2 パッケージを読む（なければ install）
#低発現遺伝子を捨てる
#vst() で変換する
#行と列を入れ替えて WGCNA 用の向きにする
#サンプル外れ値の確認（hclust）
#保存する

library(DESeq2)

#Edit matrix
expr_mat <- round(as.matrix(expr_mat))　　　#make expression matrix
mode(expr_mat) <- "integer"         #make values into integer just in case

expr_mat[1:5, 1:3]


#Remove low read counts gene
keep <- rowSums(expr_mat) >= 50  #filter out gene with fewer than 50 total reads
sum(keep)

expr_mat <- expr_mat[keep, ]
dim(expr_mat)     #59033 => 43186


#build DEseq2 object
coldata <- data.frame(
  row.names = muscle_info$SAMPID,
  sex = muscle_info$SEX_LABEL,
  age = muscle_info$AGE
)

dds <- DESeqDataSetFromMatrix(
  countData = expr_mat,
  colData   = coldata,
  design    = ~ 1   # DESeq is not needed, just for VST only
)

dds


#Variance Stabilized Data　VST後のData
#分散は平均からの差だから、同じスケール（相対変化）で比較しないと、真の姿を捉えられない。
#発現レベルごとに異なるスケールを補正し、相対的な生物学的変動を捉える
vsd <- vst(dds, blind = TRUE)    #VST(Variance Stabilizing Transformation:分散安定化変換)
vst_mat <- assay(vsd)      #Extract Expressed matrix from vsd
                           #vstには発現量、サンプル情報、遺伝情報、その他が入ってる
dim(vst_mat)
vst_mat[1:5, 1:3]


#Transpose the expression matrix
#WGCNA expects samples as rows and genes as columns
datExpr <- t(vst_mat)    #t() :transpose row and columns
dim(datExpr)             #WGCNA用の発現行列（行＝サンプル,列＝遺伝子）
datExpr[1:3, 1:5]


#Save
saveRDS(datExpr,    file = "clean/datExpr_v10_muscle_skeletal.rds")
saveRDS(muscle_info, file = "clean/datTraits_v10_muscle_skeletal.rds")

list.files("clean")


#Check
sum(keep)       #43186
dim(expr_mat)   #43186   818
dim(vst_mat)    #43186   818
dim(datExpr)    #818 43186


