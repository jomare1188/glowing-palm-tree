library(cogeqc)
library(tidyverse)

# read OG tsv file
orthogroups <- read_orthogroups("/home/dmpachon/jorge/comparative_cane/results/orthofinder/Results_Aug30/Orthogroups/Orthogroups.tsv")
# format sorghum proteins (delete final .p* from name)
orthogroups$Gene <- sub("\\.p.*$", "", orthogroups$Gene)
colnames(orthogroups)[colnames(orthogroups) == "Gene"] <- "gene"

# Sorghum
# read cor genes from sorghum
sorghum_table <- read.table("../results/correlation_analysis/abs/abs_cor_sorghum_sweet_rho80_padj001.csv", header = T, sep = ",")
# read modules table for sorghum
modules_sorghum <- read.table("../results/networks/abs_mcl_p0.6_i2.7.mcl.formated.csv", col.names = c("gene", "module"))
# get the modules correlated 
sorghum_cor_genes_and_modules <- modules_sorghum %>% filter(module %in% sorghum_table$Module)
# merge correlated modules with orthogroups table
sorghum_cor_OG_table <- merge(sorghum_cor_genes_and_modules, orthogroups, by = "gene")

# Sugarcane
# read cor genes from sugarcane
sugarcane_table <- read.table("../results/correlation_analysis/abs/abs_cor_sugarcane_sweet_rho80_padj001.csv", header = T, sep = "," )
modules_sugarcane <- read.table("../results/networks/abs_mcl_cane_p0.6_i2.7.mcl.formated.csv", col.names = c("gene", "module"))
# get the modules correlated 
sugarcane_cor_genes_and_modules <- modules_sugarcane %>% filter(module %in% sugarcane_table$Module)
# merge correlated modules with orthogroups table
sugarcane_cor_OG_table <- merge(sugarcane_cor_genes_and_modules, orthogroups, by = "gene")
# ??
# making full table with orthologues correlated genes 
test <- sorghum_cor_OG_table %>% filter(Orthogroup %in% sugarcane_cor_OG_table$Orthogroup)
clusters_expattern_sorghum <- read.table("../results/correlation_analysis/abs/abs_clusters_correlated_modules_sorghum_rho08.csv", header = T, sep = ",")
clusters_expattern_sorghum$module <- sub("module_", "", rownames(clusters_expattern_sorghum))
tmp <- merge(test, clusters_expattern_sorghum, by = "module")
sorghum_rich_table <- tmp %>% select(gene, Orthogroup, Species, module, cluster)

test2 <- sugarcane_cor_OG_table %>% filter(Orthogroup %in% sorghum_cor_OG_table$Orthogroup)
clusters_expattern_sugarcane <- read.table("../results/correlation_analysis/abs/abs_clusters_correlated_modules_sugarcane_rho08.csv", header = T, sep = ",")
clusters_expattern_sugarcane$module <- sub("module_", "", rownames(clusters_expattern_sugarcane))
tmp <- merge(test2, clusters_expattern_sugarcane, by = "module")
sugarcane_rich_table <- tmp %>% select(gene, Orthogroup, Species, module, cluster)

full_cor_table <- rbind(sugarcane_rich_table, sorghum_rich_table)

result_comparative <- full_cor_table %>%
  group_by(Orthogroup) %>%
  summarize(different_clusters = sum(cluster != cluster[1]))
## See in which expression pattern cluster are the orthologues in sugarcane and sorghum
table(result_comparative$different_clusters)/sum(table(result_comparative$different_clusters))
# 12 OG have the genes from sorghum and sugcarcane in the same cluster
# 22 OG have one gene from another cluster
# 2 OG have two genes from differents cluster

## See in which groups are orthologues tf
tf_sugarcane <- read.table("../results/tf_annot/sugarcane/tf_formatted_sugarcane.txt", header = T)
colnames(tf_sugarcane)[1] <- "gene"

tf_sorghum <- read.table("../results/tf_annot/sorghum/tf_formatted_sorghum.txt", header = T)
colnames(tf_sorghum)[1] <- "gene"

ttt <- merge(tf_sorghum, sorghum_rich_table, by = "gene")
zzz <- merge(tf_sugarcane, sugarcane_rich_table, by = "gene")

full_cor_tf_OG_table <- rbind(ttt, zzz)

result_comparative_tf <- full_cor_tf_OG_table %>%
  group_by(Orthogroup) %>%
  summarize(different_clusters = sum(cluster != cluster[1]))
## ??
## see in which cluster are the orthologues in sugarcane and sorghum  i.e. the orthologues preserve its expression pattern in sorghum and sugarcane?
## read matrices from expression patterns of eigengene of modules correlated to sweet - fiber in sugarcane and sorghum

eigen_sorghum_mods <- read.table("../results/matrices/eigen_sorghum_cor_rho08.csv", row.names = 1, header = T, sep = ",") 
eigen_cane_mods <- read.table("../results/matrices/eigen_sugarcane_cor_rho08.csv", row.names = 1, header = T, sep = ",")
library(factoextra)
df <- eigen_sorghum_mods

fviz_nbclust(df, kmeans, method = "silhouette")
fviz_nbclust(df, kmeans, method = "wss") +
geom_vline(xintercept = 3, linetype = 2)
library(cluster)
gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
 K.max = 10, B = 10)
fviz_gap_stat(gap_stat)

##


# make heatmaps
library(scico)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(wesanderson)
library(parallel)

vst_path <- "../results/matrices/sugarcane_RUVr_k2_vst_counts.tsv"
# "../results/matrices/sorghum_RUVr_k1_vst_counts.tsv"
metadata_path <-"../data/sugarcane_metadata.csv"
# "../data/sorghum_metadata.csv"

vst <- read.table(vst_path)
sample_table <- read.table(metadata_path, sep = ",", header = T)
# only for sugarcae
colnames(sample_table)[2] <- "Cultivar"
colnames(sample_table)[1] <- "Run"

# define groups
grain_genotypes <- c("Rio", "Keller", "Della")
sweet_genotypes <- c("BTx623", "BTx642", "RTx430")
sample_table$Group <- ifelse((sample_table$Cultivar %in% sweet_genotypes), "Sweet", "Grain")
# Cane
sweet_genotypes <- c("SRA5", "SRA3", "Q238")
fiber_genotypes <- c("KQB09-20432", "Q157", "Q151", "Q135")
sample_table$Group <- ifelse((sample_table$Cultivar %in% sweet_genotypes), "Sugar", "Fiber")

sweet <- ifelse((sample_table$Cultivar %in% sweet_genotypes), 1, 0)

# preprate anotation for heatmap
annotation_col <- sample_table
anot <- select(annotation_col, "Cultivar", "Group")
anot$Cultivar <- as.factor(anot$Cultivar)
anot$Group <- as.factor(anot$Group)
rownames(anot) <- annotation_col$Run

colorGroup=wes_palette("Darjeeling1", 2, type = "discrete")
#colorCultivar= c("#E69F00", "#56B4E9", "#9D9C85", "#F0E442", "#D55E00", "#CC79A7")
colorCultivar = brewer.pal(n = 7, name = "Set2")

# CANE
ann_colors = list(
  Cultivar = c(SRA5=colorCultivar[1], SRA3=colorCultivar[2], Q238=colorCultivar[3], "KQB09-20432"=colorCultivar[4], Q157=colorCultivar[5], Q151=colorCultivar[6], Q135=colorCultivar[7]),
  Group = c(Sugar = colorGroup[1], Fiber = colorGroup[2])
)

# SORGHUM
ann_colors = list(
  Cultivar = c(BTx623=colorCultivar[1], BTx642=colorCultivar[2], Della=colorCultivar[3], Keller=colorCultivar[4], Rio=colorCultivar[5], RTx430=colorCultivar[6]),
  Group = c(Sweet = colorGroup[1], Grain = colorGroup[2])
)


names <- test2
"test"
df <- vst[names$gene,]
colnames(df) <- colnames(vst)
rownames(anot) <- colnames(df)
# heat map with ean values for column i.e media por modulo
png("../results/correlation_analysis/heatmap_sorghum_orthologues.png", res = 300, width = 2*1200, height = 0.65*2850)
pheatmap(df,
         main =paste0("Sorghum Orthologues"),
         scale = "row",
         annotation_col = anot,
         show_rownames = F,
         show_colnames = F,
         annotation_colors = ann_colors,
         col = rev(scico(100, palette = 'roma')),
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = 3,
         cellheight = 2)
dev.off()


## make GO enrichment

library(topGO)
library(ggplot2)
library(tidyverse)
library(parallel)
library(clusterProfiler)
#results_path_cluster <-"/home/dmpachon/jorge/comparative_cane/results/correlation_analysis/"
dir.create(results_path_cluster)

GO <- read.table("/home/dmpachon/jorge/comparative_cane/data/GO_sugarcane/GO_cane_formated.txt", header=FALSE, stringsAsFactors=FALSE)
#GO <- read.table("/home/dmpachon/jorge/comparative_cane/results/panzzer/GO_formated.txt", header=FALSE, stringsAsFactors=FALSE)


# Rename columns if necessary (assuming your data has two columns)
colnames(GO) <- c("Gene", "GO_term")

# For cane
ids_longest <- read.table("/home/dmpachon/jorge/comparative_cane/data/GO_sugarcane/ids_longest_protein_per_OG.txt",  col.names = "Gene")
tmp <- GO %>% filter(Gene %in% ids_longest$Gene)
GO <- tmp
###

# Group by Gene and aggregate GO terms into a single string separated by spaces
formatted_GO <- aggregate(GO_term ~ Gene, GO, function(x) paste(x, collapse=" "))
gene2GO <- strsplit(formatted_GO$GO_term , " ")
names(gene2GO) <- formatted_GO$Gene

#
corr_mods <- read.table("../results/correlation_analysis/abs/abs_clusters_correlated_modules_sugarcane_rho08.csv", header = T, , sep = ",")
corr_mods$module <- sub("module_", "", rownames(corr_mods))

corr_mods <- corr_mods %>% filter(cluster == 2)

dynamicMods <- read.table(file = "../results/networks/abs_mcl_cane_p0.6_i2.7.mcl.formated.csv", header = F, row.names=1)
colnames(dynamicMods) <- "module_No"
geneNames <- rownames(dynamicMods)


myInterestingGenes <- geneNames[dynamicMods$module_No %in% corr_mods$module]
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)
allGO=usedGO(GOdata)
# MIRAR SI DA VALOR Pcorregido el TOP GO
# classic ingnora la topologia del go
# revisar algoritmos 

Classic <<- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#resultsWeight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
# Make results  table
#table <- GenTable(GOdata, Classic = resultClassic, Weight01 = resultsWeight01, topNodes = length(allGO), orderBy = 'Classic')
table <- GenTable(GOdata, Classic = Classic, topNodes = length(allGO), orderBy = 'Classic')
# Filter not significant values for classic algorithm
table1 <- filter(table, Classic < 0.01 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"), digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <<- cbind(table1,p.adj)
all_res_final <<- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.01),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.01),]
# Save first top 50 ontolgies sorted by adjusted pvalues
write.table(results.table.bh, file = "../results/GO_enrichment/abs/correlation/group2_sugarcane_GO.csv", quote=FALSE, row.names=FALSE, sep = ",")

ntop <- 22

ggdata <- results.table.bh[1:ntop,]
ggdata <- ggdata[complete.cases(ggdata), ]
#ggdata <- all_res_final
aux <- go2term(ggdata$GO.ID)
colnames(aux) <- c("GO.ID", "Lterm")

ggdata <- merge(ggdata, aux, by = "GO.ID")
ggdata$p.adj <- as.numeric(ggdata$p.adj)

ggdata <- ggdata[order(ggdata$p.adj),]
ggdata$Lterm <- factor(ggdata$Lterm, levels = rev(ggdata$Lterm)) # fixes order
gg1 <- ggplot(ggdata, aes(x = Lterm, y = -log10(p.adj) ))+
  geom_point(size = 6, colour = "black") +
  scale_size(range = c(2.5,12.5)) +
  xlab('GO Term') +
  ylab('-log(p)') +
  labs(title = 'GO Biological processes')+
  theme_bw(base_size = 24) +
  coord_flip()
ggsave("../results/GO_enrichment/abs/correlation/group2_sugarcane_GO.png", device = "png", width = 40, height = 30, dpi = 300, units = "cm")

