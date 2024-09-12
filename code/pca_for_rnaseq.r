## variables
cores = 1
# set working directory

## install packages
# install.packages(c("BiocManager","backports","tidyverse"),
# lib = "./libraries",
# repos = 'http://cran.us.r-project.org')
#BiocManager::install("tximport", lib = "./../libraries")
#BiocManager::install("DESeq2", lib = "./../libraries")
#BiocManager::install("BiocParallel", lib = "./../libraries")

## load libraries 
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggrepel)

#txgen2="../data/sugarcane_tx2gen.txt"
txgen2 = "../data/sorghum_tx2gen.txt"
#metadata="../data/sugarcane_metadata.csv"
metadata = "../data/sorghum_metadata.csv"

get_full_trans_count_matrixes <- function(metadata, txgen2) {
        sample_table <-read.table(metadata, sep = ",", header = T)
	sugar_genotypes <- c("SRA5", "SRA3", "Q238")
	fiber_genotypes <- c("KQB09-20432", "Q157", "Q151", "Q135")

	grain_genotypes <- c("Rio", "Keller", "Della")
	sweet_genotypes <- c("BTx623", "BTx642", "RTx430")

	sample_table$group <- ifelse((sample_table$Cultivar %in% sweet_genotypes), "Sweet", "Grain")
#	sample_table$group <- ifelse((sample_table$genotype %in% sugar_genotypes), "sugar", "fiber")
	
	separated_wider <- separate_wider_delim(sample_table , cols = ConditionDescription, delim = "_", names = c("time", "genotype2", "tissue"))
	separated_wider$ConditionDescription <- sample_table$ConditionDescription
	sample_table <- separated_wider

        # load files paths
#       sample_files = paste0("/home/dmpachon/jorge/comparative_cane/data/", pull(sample_table, "genotype"),"/3_salmon_quant_longest/" ,pull(sample_table , "SampleID"), "/quant.sf")
	sample_files = paste0("/home/dmpachon/jorge/comparative_cane/data/", pull(sample_table, "Cultivar"),"/3_salmon_quant_longest/" ,pull(sample_table , "Run"), "/quant.sf")
        # name table columns
#       names(sample_files) = pull(sample_table, "SampleID")
	names(sample_files) = pull(sample_table, "Run")
        # relate genes to transcripts
        tx2gene = read.table(txgen2, sep = ",", col.names =c("transid","geneid"))
        # import count data to tximport
        count_data = tximport( files = sample_files,
                               type = "salmon",
                               tx2gene =  tx2gene,
                               txOut = T,
                               ignoreTxVersion = F)

        raw <- DESeqDataSetFromTximport(txi = count_data,
                                        colData = sample_table,
                                        design = ~ group)
        # full raw and normalized count matrixes for all experiments
        data <- estimateSizeFactors(raw)
        vst <- varianceStabilizingTransformation(data)
        df_data <- assay(raw)
        df_vst <- assay(vst)
        write.table(df_data, file = paste0("../results/", "cane_raw_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
        write.table(df_vst, file = paste0("../results/", "cane_vst_counts.tsv", sep = ""), sep = "\t" ,quote = F, col.names = T)
	# Make PCA with DESEQ2 function
#	pca_data <- plotPCA(vst, intgroup = c("group", "genotype", "tissue", "ConditionDescription", "time"), returnData = TRUE)
	pca_data <- plotPCA(vst, intgroup = c("group", "Cultivar", "tissue", "Age", "BioProject", "dev_stage"), returnData = TRUE)
	# Plot customized PCA 
	percentVar <- round(100 * attr(pca_data, "percentVar"))
	# plot R color condition
	colors = wes_palette("Darjeeling1",2 , type = "discrete")
	# Plooot 
	p <- ggplot(pca_data, aes(x = PC1, y = PC2, shape = as.factor(genotype),  color = as.factor(group.1)))
#	p <- p + geom_jitter(aes(color = group), size = 3, show.legend = FALSE) +
	p <- p + facet_wrap(~ ConditionDescription)
	p <- p + geom_point(size=3)
	p <- p + scale_shape_manual(values = c(15,16,17,18,3,4,6))
	p <- p + labs(title = "sugarcane", color="Group", shape = "Genotype" )
	p <- p + xlab((paste0("PC1 : ", percentVar[1],"%")))
	p <- p + ylab((paste0("PC2 : ", percentVar[2],"%")))
#	p <- p + geom_text_repel(aes(label=name), size=2.5, max.overlaps = 20 )
	p <- p + scale_colour_manual(values = colors)
	p <- p + theme_bw()
	p <- p + theme( text = element_text(family = "Times New Roman", size=22),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"))

ggsave(p, filename = "../results/pca_sugarcane_facet.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

}
#get_full_gene_count_matrixes("./../data/metadata_complete.csv",  "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")
#get_full_trans_count_matrixes("./../data/metadata_complete.csv",  "/Storage/data1/jorge.munoz/good_t2gene/results/tx2gen_DESEQ.csv")

