# files
modules_path <- "../results/networks/abs_mcl_p0.6_i2.7.mcl.formated.csv"
# "../results/networks/abs_mcl_cane_p0.6_i2.7.mcl.formated.csv"
vst_path <- "../results/matrices/sugarcane_RUVr_k2_vst_counts.tsv"
# "../results/matrices/sugarcane_RUVr_k2_vst_counts.tsv"
metadata_path <-"../data/sugarcane_metadata.csv"

#libraries
library(tidyverse)

modules <- read.table(modules_path, col.names = c("gene", "module_No"))
rownames(modules) <- modules$gene

#colnames(modules) <- c("module_No")
# read vst matrix 
vst <- read.table(vst_path)
sample_table <- read.table(metadata_path, sep = ",", header = T)

# only for sugarcae
colnames(sample_table)[2] <- "Cultivar"
#

# define groups
# SORGHUM
#grain_genotypes <- c("Rio", "Keller", "Della")
#sweet_genotypes <- c("BTx623", "BTx642", "RTx430")
#sample_table$Group <- ifelse((sample_table$Cultivar %in% sweet_genotypes), "Sweet", "Grain")
# CANE 
sweet_genotypes <- c("SRA5", "SRA3", "Q238")
fiber_genotypes <- c("KQB09-20432", "Q157", "Q151", "Q135")

sweet <- ifelse((sample_table$Cultivar %in% sweet_genotypes), 1, 0)

# make correlation analysis
make_correlation_analysis <- function(module, trait){
df <- vst[modules[modules$module_No == module,]$gene,]
names <- modules[modules$module_No == module,]
df <- vst[names$gene,]
colnames(df) <- colnames(vst)
pca <- prcomp(t(df))
eigen <-pca$x[,1]
cor <- cor.test(eigen, trait, method = "spearman", exact=FALSE)
pvalue <- cor$p.value
rho <- cor$estimate
raw <- data.frame(rho=rho, p = pvalue, Module = module, row.names = NULL, trait = "sweet")

return(raw)
}

results_list <- list()

for (module in 1:max(modules$module_No)) {
  result <- make_correlation_analysis(module, sweet)
  results_list[[module]] <- result
}

joined_results <- do.call(rbind, results_list)
filtered_results <- joined_results %>% mutate(bonferroni = p.adjust(p, method = "bonferroni")) %>% select (Module, trait, rho, p, bonferroni) %>% filter( abs(rho) > 0.80 & bonferroni < 0.01)

write.table(filtered_results, "../results/correlation_analysis/abs_cor_sugarcane_sweet_rho80_padj001.csv", sep =  ",",  quote = F, row.names = F)


