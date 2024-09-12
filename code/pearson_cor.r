library(WGCNA)
library(dplyr)
library(HiClimR)
library(reshape2)
library(DESeq2)
library(Matrix)
library(data.table)
library(spam)
library(spam64)
options(spam.force64=TRUE)
vst <- as.matrix(read.table("../results/matrices/sorghum_RUVr_k1_raw_counts.tsv", sep = "\t", header = T, row.names = 1))
#vst <- as.matrix(read.table("/home/dmpachon/jorge/comparative_cane/results/matrices/sorghum_RUVr_k1_vst_counts.tsv", sep = "\t", header = T, row.names = 1))

row_means <- rowMeans(vst)
row_sds <- rowSds(vst)

# Calculate coefficient of variation (CV) for each gene
# 100 = 3rd quartil for sorghum
# 31 = 3rd quartil for cane
#cv <- (row_sds / row_means) * 100
#summary(cv)
# Set a threshold for the minimum coefficient of variation you want to keep
#min_cv_threshold <- 31

# Filter rows (genes) based on the coefficient of variation threshold
#filtered_genes <- vst[cv >= min_cv_threshold, ]
#dim(filtered_genes)

pcor <- round(as.data.table(fastCor(t(vst), nSplit = 50, upperTri = TRUE, optBLAS = TRUE, verbose = TRUE)), digits = 3)
rownames(pcor) <- colnames(pcor)
write.table(pcor, file = "../results/networks/matrix_full_sorghum.txt", quote = F, row.names = T, col.names = T, sep="\t")
#sparse_pcor <- as(as.matrix(pcor), "sparseMatrix")
#dim(pcor)
#melt <- na.exclude(melt(sparse_pcor))

#person_fil <- melt %>% filter(abs(value) > 0.8 & Var1 != Var2) %>% mutate(value = abs(value*1000))
#dim(person_fil)

#write.table(person_fil, file = "../results/networks/cane/triples_0.8_cv_0.txt", quote = F, row.names = F, col.names = F, sep="\t")
