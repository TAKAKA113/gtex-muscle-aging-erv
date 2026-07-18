dir.create("raw")
dir.create("meta")
dir.create("clean")
dir.create("results")
list.files()

#
url1 <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/counts-by-tissue/gene_reads_v10_muscle_skeletal.gct.gz"
file1 <- "raw/gene_reads_v10_muscle_skeletal.gct.gz"
download.file(url1, file1, mode = "wb")

file.exists(file1)
file.size(file1)

#TPM file
url2 <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_muscle_skeletal.gct.gz"
file2 <- "raw/gene_tpm_v10_muscle_skeletal.gct.gz"
download.file(url2, file2, mode = "wb")

file.exists(file2)
file.size(file2)
list.files("raw")


#Metadata
url3 <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
file3 <- "meta/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
download.file(url3, file3, mode = "wb")

file.exists(file3)
file.size(file3)
list.files("meta")


#subject_phenotype
url4 <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
file4 <- "meta/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
download.file(url4, file4, mode = "wb")

file.exists(file4)
file.size(file4)
list.files("meta")

