
# files
modules_path <- "../results/networks/abs_mcl_sorghum_p0.6_i2.7.mcl.formated.csv"
vst_path <- "../results/matrices/sorghum_RUVr_k1_vst_counts.tsv"
metadata_path <-"../data/sorghum_metadata.csv"


#libraries
library(scico)
library(RColorBrewer)
library(tidyverse)
library(pheatmap) 
library(wesanderson)
library(parallel)

modules <- read.table(modules_path, col.names = c("gene", "module_No"))
rownames(modules) <- modules$gene

#colnames(modules) <- c("module_No")
# read vst matrix 
vst <- read.table(vst_path)
sample_table <- read.table(metadata_path, sep = ",", header = T)
# only for sugarcae
colnames(sample_table)[2] <- "Cultivar"
colnames(sample_table)[1] <- "Run"
#


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

# prepare colors for heatmap
#colorCultivar=wes_palette("AsteroidCity2", 6, type = "discrete")
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


# make heatmap
make_heatmap_sorghum <- function(i){
df <- vst[modules[modules$module_No == i,]$gene,]
names <- modules[modules$module_No == i,]
df <- vst[names$gene,]
colnames(df) <- colnames(vst)
rownames(anot) <- colnames(df)
# heat map with ean values for column i.e media por modulo
png(paste0("../results/heatmaps/sorghum/abs/", "heatmap_samples_module_", i, ".png", sep = ""), res = 300, width = 2*1200, height = 2*2850)
pheatmap(df,
         main =paste0("Sorghum Module ",i , sep = ""),
         scale = "row",
         annotation_col = anot,
         show_rownames = F,
	 show_colnames = F,
	 annotation_colors = ann_colors,
         col = rev(scico(100, palette = 'roma')),
         cluster_cols = T,
         cluster_rows = T,
         cellwidth = 3,
         cellheight = 0.2) 
dev.off()
}


numCores <- 100
mclapply(modules$module_No, make_heatmap_sorghum , mc.cores = numCores)

## make heatmaps from eigen genes

module_list <- read.table("../results/correlation_analysis/abs/abs_cor_sorghum_sweet_rho80_padj001.csv", header = T, sep = ",")
make_heatmaps_eigengen <- function(module){
df <- vst[modules[modules$module_No == module,]$gene,]
names <- modules[modules$module_No == module,]
df <- vst[names$gene,]
colnames(df) <- colnames(vst)
pca <- prcomp(t(df))
eigen <-pca$x[,1]

raw <- data.frame(eigen=t(eigen), row.names = paste0("module_",module))

return(raw)
}

results_list <- list()

for (module in module_list$Module) {
  result <- make_heatmaps_eigengen(module)
  results_list[[module]] <- result
}


## anot row
joined_results <- do.call(rbind, results_list)
tmp <- pheatmap(joined_results, cluster_rows = T, cluster_cols = T, scale = "row")
anot_row <- as.data.frame(cutree(tmp$tree_row, k = 2))
colnames(anot_row) <- "cluster"

colnames(joined_results) <- colnames(vst)
write.table(joined_results, "../results/matrices/abs_eigen_sugarcane_cor_rho08.csv", col.names = T, row.names = T, sep = ",", quote = F)
write.table(anot_row, "../results/correlation_analysis/abs/abs_clusters_correlated_modules_sugarcane_rho08.csv", col.names = T, row.names = T, sep = ",", quote = F)

#png(paste0("../results/correlation_analysis/heatmap_eigen_cor_mod_sweet.png"), res = 300, width = 2.4*1200, height = 4200)
png(paste0("../results/correlation_analysis/abs/abs_heatmap_sugarcane_eigen_cor_mod_sweet.png"), width = 2.2*1200, height = 1600, res = 300)
pheatmap(joined_results,
         main = "",
         scale = "row",
         annotation_col = anot,
         show_rownames = T,
         show_colnames = F,
         annotation_colors = ann_colors,
         col = rev(scico(100, palette = 'roma')),
         cluster_cols = T,
	 cutree_rows = 2,
         cluster_rows = T,
         cellwidth = 3,
         cellheight = 5)
#8
dev.off()

