library(topGO)
library(ggplot2)
library(dplyr)
library(parallel)
library(clusterProfiler)
results_path_cluster <-"/home/dmpachon/jorge/comparative_cane/results/GO_enrichment/sugarcane/modules/"
dir.create(results_path_cluster)
numCores <- 30
## read clusters size
#number.clusters <- read.table("../results/networks/sorghum/mcl_p90_i2.7.mcl.formated_numbers.csv", header = F)

GO <- read.table("/home/dmpachon/jorge/comparative_cane/data/GO_sugarcane/GO_cane_formated.txt", header=FALSE, stringsAsFactors=FALSE)

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
dynamicMods <- read.table(file = "../results/networks/mcl_cane_p0.6_i2.7.mcl.formated.csv", header = F, row.names=1)
colnames(dynamicMods) <- "module_No"
geneNames <- rownames(dynamicMods)

# make annot fun
anot_modules <- function(Module_no, results_path){
#paste0("MODULO:", Module_no,"      ","NODOS:",number.clusters[Module_no,1])
myInterestingGenes <- geneNames[dynamicMods$module_No==Module_no]
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
table1 <- filter(table, Classic < 0.05 )
# Performing BH correction on our p values FDR
p.adj <- round(p.adjust(table1$Classic,method="BH"),digits = 4)
# Create the file with all the statistics from GO analysis
all_res_final <<- cbind(table1,p.adj)
all_res_final <<- all_res_final[order(all_res_final$p.adj),]
# Get list of significant GO before multiple testing correction
results.table.p = all_res_final[which(all_res_final$Classic <=0.05),]
# Get list of significant GO after multiple testing correction
results.table.bh = all_res_final[which(all_res_final$p.adj<=0.05),]
# Save first top 50 ontolgies sorted by adjusted pvalues
write.table(results.table.bh, file = paste0(results_path_cluster, "module_", Module_no, ".csv"), quote=FALSE, row.names=FALSE, sep = ",")

ntop <- 20

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
ggsave(paste0(results_path_cluster, "module_", Module_no, ".png", sep = ""), device = "png", width = 40, height = 30, dpi = 300, units = "cm")



}

mclapply(1:dim(table(dynamicMods)), anot_modules, mc.cores = numCores)

#for (i in 1:dim(table(dynamicMods))){ 
#	if(number.clusters[i,1]>=4){
		print(i)
		anot_modules(i,results_path_cluster)
		
#}
#}


