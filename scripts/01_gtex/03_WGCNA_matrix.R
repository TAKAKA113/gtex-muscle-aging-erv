#発現の表（行: 遺伝子, 列: サンプル）を R に読み込み、metadata の 818 sample と完全に対応させる

#Read the expression data
#GCT file contain two metadata lines at beginning, so skip the first two lines when reading the file

expr <- read.delim(
  gzfile("raw/gene_reads_v10_muscle_skeletal.gct.gz"),
  skip = 2,             # 上の2行（#1.2 と 行数 列数）をとばす
  header = TRUE,        # 3行目を列名にする
  sep = "\t",
  check.names = FALSE
)


dim(expr)     #59033   820  so i got 59033 genes
expr[1:5, 1:5]



#Separate express data into gene_info and expr_matrix
#Since WGCNA requiers a numeric expression matrix
#Gean ID and Name are moved to another file for downstream analysis

gene_info <- expr[, 1:2]           # 1〜2列目: Name と Description
expr_mat  <- expr[, -c(1, 2)]      # それ以外: サンプルごとの発現

dim(gene_info)
dim(expr_mat)  #59033   818



#Set gene ID as row names
#WGCNAやDEseqは発現量のみの行列、Nameという文字列が邪魔、解析ソフトが文字があるとerror
#だから遺伝子IDを行名にして、ラベルとして扱う
rownames(expr_mat) <- gene_info$Name    #At this point,I got pure expression matrix
expr_mat[1:5, 1:3]

#cheack
nrow(muscle_info)
head(colnames(expr_mat))
head(muscle_info$SAMPID)

#Reorder the Sample metadata to match the column order of the expression matrix
#発現データとサンプル情報を一致させる。列・行の一致
# metadata 側を expression と同じ並びにする
#これをしないとサンプルに別のサンプルの発現量がついてしまう
#メタデータはsampleIDを基準に発現行列の順番に合わせる
muscle_info <- muscle_info[match(colnames(expr_mat), muscle_info$SAMPID), ]

# 完全一致を確認
all(colnames(expr_mat) == muscle_info$SAMPID)　　#TRUE

expr_mat[1:5, 1:5]
muscle_info[1:5, ]


#Save
saveRDS(expr_mat,    file = "clean/expr_mat_v10_muscle_skeletal.rds")
saveRDS(gene_info,   file = "clean/gene_info_v10_muscle_skeletal.rds")
saveRDS(muscle_info, file = "clean/muscle_info_v10_muscle_skeletal.rds")

list.files("clean")



dim(expr)
dim(expr_mat)
all(colnames(expr_mat) == muscle_info$SAMPID)
