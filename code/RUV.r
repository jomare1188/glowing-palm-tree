#BiocManager::install("RUVSeq")
library("RUVSeq")
library("DESeq2")
library("tximport")
library("tidyverse")
library("wesanderson")

# load data direct from deseq
#txgen2="../data/sugarcane_tx2gen.txt"
txgen2 = "../data/sorghum_tx2gen.txt"
#metadata="../data/sugarcane_metadata.csv"
metadata = "../data/sorghum_metadata.csv"


# matrices
#sorghum_raw <- read.table("../results/sorghum_raw_counts.tsv", header = T, row.names = 1, sep = "\t")
#cane_raw <- t(read.table("../results/sugarcane_raw_counts.tsv", header = T, row.names = 1, sep = "\t"))
# metadata
#cane_metadata <- read.table("../data/sugarcane_metadata.csv", header = T,  sep = ",")
sorghum_metadata <- read.table("../data/sorghum_metadata.csv", header = T,  sep = ",")

sweet_genotypes <- c("BTx623", "BTx642", "RTx430") # for sorghum
sorghum_metadata$group <- ifelse((sorghum_metadata$Cultivar %in% sweet_genotypes), "Sweet", "Grain")

#sugar_genotypes <- c("SRA5", "SRA3", "Q238") # for cane
#cane_metadata$group  <-  ifelse((cane_metadata$genotype %in% sugar_genotypes), "sugar", "fiber")


sorghum_metadata_files = paste0("/home/dmpachon/jorge/comparative_cane/data/", pull(sorghum_metadata, "Cultivar"),"/3_salmon_quant_longest/" ,pull(sorghum_metadata , "Run"), "/quant.sf")
names(sorghum_metadata_files) = pull(sorghum_metadata, "Run")

tx2gene = read.table(txgen2, sep = ",", col.names =c("transid","geneid"))

# import count data to tximport
count_data = tximport( files = sorghum_metadata_files,
	type = "salmon",
        tx2gene =  tx2gene,
        txOut = T,
        ignoreTxVersion = F)

raw <- DESeqDataSetFromTximport(txi = count_data,
        colData = sorghum_metadata,
        design = ~ group)

# filter not expressed genes
dim(counts(raw))
keep<- rowSums(counts(raw))>0
table(keep)
raw <- raw[keep,]
print("datos raw filtrados")
dim(raw)

# RUVs
cv_row <- function(row) {
  sd(row) / mean(row) * 100
}

cv <- apply(counts(raw), 1, cv_row)
negative_control <- names(head(cv[order(cv)],1000))
differences <- makeGroups(sorghum_metadata$Cultivar)

# RUVr
# some magic
design <- model.matrix(~sorghum_metadata$group)
y <- DGEList(counts=counts(raw), group=sorghum_metadata$group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

# try for different k values
set_RUVr <- RUVr(y$counts,rownames(y), k=1, res)
dataRUVr <- DESeqDataSetFromMatrix(countData = set_RUVr$normalizedCounts,
                                   colData = sorghum_metadata,
                                   design = ~ group)

RUVr_vst <- varianceStabilizingTransformation(dataRUVr)

# plot PCA for RUVr(normalized) and vst(transformed) data
colors = wes_palette("Darjeeling1",2 , type = "discrete")
pca_data <- plotPCA(RUVr_vst, intgroup = c("group", "Cultivar"), returnData = TRUE)

percentVar <- round(100 * attr(pca_data, "percentVar"))
# plot R color condition
colors = wes_palette("Darjeeling1",2 , type = "discrete")
# Plooot 
p <- ggplot(pca_data, aes(x = PC1, y = PC2, shape = as.factor(Cultivar),  color = as.factor(group.1)))
p <- p + geom_point(size=3)
p <- p + labs(title = "sorghum", color="Group", shape = "Genotype" )
p <- p + xlab((paste0("PC1 : ", percentVar[1],"%")))
p <- p + ylab((paste0("PC2 : ", percentVar[2],"%")))
#p <- p + geom_text_repel(aes(label=name), size=2.5, max.overlaps = 20 )
p <- p + scale_colour_manual(values = colors)
#p <- p + scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8))
p <- p + theme_bw()
p <- p + theme(text = element_text(family = "Times New Roman", size=22),
		panel.border = element_blank(),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		axis.line = element_line(colour = "black"))
ggsave(p, filename = "../results/RUVr/sorghum_k1_RUVr_groups.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

#
#	theme_bw() +
#	geom_point(size=3.5, aes(shape = RUVr_vst$genotype)) +
#	labs(title = "RUVr cane", col="Group", shape="Genotype") +
#	scale_colour_manual(values = colors) +
#	scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 8))+
#  	theme(text = element_text(size=22),
#		panel.border = element_blank(),
#		panel.grid.major = element_blank(),
#		panel.grid.minor = element_blank(),
#		plot.title = element_text(hjust = 0.5),
#		axis.line = element_line(colour = "black"))
#ggsave(ruvr_condition, filename = "../results/RUVr/cane/k2_RUVr_groups.png",units = "cm",width = 15*1.3, height = 15,dpi = 320)

write.table(assay(dataRUVr), file = "../results/matrices/sorghum_RUVr_k1_raw_counts.tsv", sep = "\t" ,quote = F, col.names = T)
write.table(assay(RUVr_vst), file = "../results/matrices/sorghum_RUVr_k1_vst_counts.tsv", sep = "\t" ,quote = F, col.names = T)

