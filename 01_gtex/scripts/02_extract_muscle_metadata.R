# 発現ファイルの上3行だけ読む
top3 <- readLines(gzfile("raw/gene_reads_v10_muscle_skeletal.gct.gz"), n = 3)

# 3行目がサンプルIDの行
header_line <- top3[3]

# タブで分けて、最初の2つ（Name, Description）を捨てる
expr_samples <- strsplit(header_line, "\t")[[1]][-c(1, 2)]

# 確認
length(expr_samples)
head(expr_samples)

#make DonorID from SampleID
# 例: "GTEX-1117F-0526-SM-5GZZ7" → "GTEX-1117F"
subj_id <- sub("^(GTEX-[^-]+).*", "\\1", expr_samples)

head(subj_id)

#Make the matrix only in Sample and DonorID
muscle_info <- data.frame(
  SAMPID = expr_samples,
  SUBJID = subj_id
)

head(muscle_info)
nrow(muscle_info)



#Check the infomation of Donor
subj <- read.delim(
  "meta/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

head(subj)
colnames(subj)

#
#I got "SUBJID"  "SEX"     "AGE"     "DTHHRDY"
#But SUBJID、SEX、AGE are only needed
subj_small <- subj[, c("SUBJID", "SEX", "AGE")]
head(subj_small)


#Conbine two matrix (SAMPID, SUBJID）and（SUBJID, SEX, AGE）
muscle_info <- merge(muscle_info, subj_small, by = "SUBJID")
head(muscle_info)
nrow(muscle_info). #nrow(muscle_info) = 818



#Change the rownames, 1 -> Male, 2 -> Female
muscle_info$SEX_LABEL <- ifelse(muscle_info$SEX == 1, "Male", "Female")

table(muscle_info$SEX_LABEL)
table(muscle_info$AGE)

#Save
write.table(
  muscle_info,
  file = "clean/muscle_skeletal_metadata_v10.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

list.files("clean")
