library(wesanderson)
library(tidyverse)
library(WGCNA)
library(cogeqc)
library(ggplot2)
# read vst expression matrix
vst_path_sugarcane <- "../results/matrices/sugarcane_RUVr_k2_vst_counts.tsv"
vst_path_sorghum <- "../results/matrices/sorghum_RUVr_k1_vst_counts.tsv"

vst_sorghum <- read.table(vst_path_sorghum)
vst_sugarcane <- read.table(vst_path_sugarcane)
# read metadata 
metadata_path_sugarcane <- "../data/sugarcane_metadata.csv"
metadata_path_sorghum <- "../data/sorghum_metadata.csv"

sample_table_sorghum <- read.table(metadata_path_sorghum, sep = ",", header = T)
sample_table_sugarcane <- read.table(metadata_path_sugarcane, sep = ",", header = T)

colnames(sample_table_sugarcane)[2] <- "Cultivar"
# make groups

# CANE 
sweet_genotypes_sugarcane <- c("SRA5", "SRA3", "Q238")
fiber_genotypes_sorghum <- c("KQB09-20432", "Q157", "Q151", "Q135")
# SORGHUM
grain_genotypes_sorghum <- c("Rio", "Keller", "Della")
sweet_genotypes_sorghum <- c("BTx623", "BTx642", "RTx430")

sample_table_sorghum$Group <- ifelse((sample_table_sorghum$Cultivar %in% sweet_genotypes_sorghum), "Sweet", "Grain")

sweet_sorghum <- ifelse((sample_table_sorghum$Cultivar %in% sweet_genotypes_sorghum), 1, 0)
sweet_sugarcane <- ifelse((sample_table_sugarcane$Cultivar %in% sweet_genotypes_sugarcane), 1, 0)

#
gene_correlation_trait_sorghum <- cor(t(vst_sorghum), sweet_sorghum, method = "spearman")
gene_correlation_trait_sugarcane <- cor(t(vst_sugarcane), sweet_sugarcane, method = "spearman")

gene_correlation_trait_sorghum.pvalue <- corPvalueStudent(gene_correlation_trait_sorghum, ncol(vst_sorghum))
gene_correlation_trait_sugarcane.pvalue <- corPvalueStudent(gene_correlation_trait_sugarcane, ncol(vst_sugarcane))

raw_sorghum <- data.frame(rho=gene_correlation_trait_sorghum, p = gene_correlation_trait_sorghum.pvalue, gene = rownames(gene_correlation_trait_sorghum))
raw_sugarcane <- data.frame(rho = gene_correlation_trait_sugarcane, p= gene_correlation_trait_sugarcane.pvalue, gene = rownames(gene_correlation_trait_sugarcane))

filtered_gene_correlation_trait_sorghum <- raw_sorghum %>% filter (p < 0.01 & abs(rho) > 0.8)
filtered_gene_correlation_trait_sugarcane <- raw_sugarcane %>% filter (p < 0.01 & abs(rho) > 0.8)


#


# read OG tsv file
orthogroups <- read_orthogroups("/home/dmpachon/jorge/comparative_cane/results/orthofinder/Results_Aug30/Orthogroups/Orthogroups.tsv")
# format sorghum proteins (delete final .p* from name)
orthogroups$Gene <- sub("\\.p.*$", "", orthogroups$Gene)
colnames(orthogroups)[colnames(orthogroups) == "Gene"] <- "gene"

# add correlation values to orthofinder results
gene_correlation_trait_sorghum_OG <- merge(raw_sorghum, orthogroups, by = "gene")
gene_correlation_trait_sugarcane_OG <- merge(raw_sugarcane, orthogroups, by = "gene")

# read tf
tf_sugarcane <- read.table("../results/tf_annot/sugarcane/tf_formatted_sugarcane.txt", header = T)
tf_sorghum <- read.table("../results/tf_annot/sorghum/tf_formatted_sorghum.txt", header = T)
colnames(tf_sugarcane)[1] <- "gene"
colnames(tf_sorghum)[1] <- "gene"

tf_correlation_trait_sorghum_OG <- merge(gene_correlation_trait_sorghum_OG, tf_sorghum, by = "gene")
tf_correlation_trait_sugarcane_OG <- merge(gene_correlation_trait_sugarcane_OG, tf_sugarcane, by = "gene")

#tf_correlation_trait_OG_sugarcane <- tf_correlation_trait_OG
#tf_correlation_trait_OG_sorghum <- tf_correlation_trait_OG

# we filter sorghum orthologues with sugarcane and abs(rho) > 0.8
tf_correlation_trait_OG_sorghum_in_sugarcane <- tf_correlation_trait_sorghum_OG  %>% filter(Orthogroup %in% tf_correlation_trait_sugarcane_OG$Orthogroup) %>% filter(abs(rho) > 0.5 & p < 0.01)
# we filter sugarcane orthologues with sorghum and abs(rho) > 0.8
tf_correlation_trait_OG_sugarcane_in_sorghum <- tf_correlation_trait_sugarcane_OG %>% filter(Orthogroup %in% tf_correlation_trait_sorghum_OG$Orthogroup) %>% filter(abs(rho) > 0.5 & p < 0.01)

tmp <- tf_correlation_trait_OG_sorghum_in_sugarcane %>% filter(Orthogroup %in% tf_correlation_trait_OG_sugarcane_in_sorghum$Orthogroup)
tmp2 <- tf_correlation_trait_OG_sugarcane_in_sorghum %>% filter(Orthogroup %in% tmp$Orthogroup)

full <- rbind(tmp, tmp2)

full$Species <- as.character(full$Species)
full$Species[full$Species == "sorghum_protein_of_longest_cds_per_orthogroup"] <- "Sorghum"
full$Species[full$Species == "sugarcane_proteins_of_longest_cds_per_OG"] <- "Sugarcane"


full <- full %>%
  mutate(corr_sign = ifelse(rho > 0, "Positive", "Negative"))
full$Species <- as.factor(full$Species)

# Summarize the data by orthogroup and species
df <- full %>%
  group_by(Orthogroup, Species, corr_sign) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = corr_sign, values_from = count, values_fill = 0)


tmp3 <- full %>%
  group_by(Orthogroup) %>%
  summarize(different_clusters = sum(rho < 0 )) %>%
  mutate(group = ifelse(different_clusters == 0, "Conserved", "Not Conserved"))

full_df <- merge(df, tmp3, by = "Orthogroup")
full_df <- full_df[order(full_df$group), ]

colors = wes_palette("Darjeeling1",2 , type = "discrete")

# Create the stacked bar chart
ggplot(full_df, aes(x = factor(Orthogroup, levels=unique(Orthogroup)), y = Positive, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_bar(aes(y = -Negative, fill = Species), stat = "identity", position = "stack") +
  geom_hline(yintercept = 0, "black") +
# geom_vline(xintercept = "OG0050036", "black")+
  labs(x = "Orthogroup", y = "Correlation sign/Number of genes", fill = "Species") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(text = element_text(family = "Times", size=16),
	axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
	panel.border = element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	axis.line = element_line(colour = "black"))

ggsave(filename = "/home/dmpachon/jorge/comparative_cane/results/correlation_analysis/stacked_bar_comparative_tf.svg" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)

full_all <- merge(full_df, full, by = "Orthogroup") %>% select (-Species.x)
# this gives you the frequency of the tf families in the genes of the "conserved" and "not conserved orthologues tf"
# Conserved are the genes in an OG which share the same sign (+/-) of the correlation coefficient in both species 
not_conserved_tf_families <- data.frame(full_all %>% filter(group != "Conserved") %>% select(Family) %>% table())
colnames(not_conserved_tf_families) <- c("Family", "Genes")


conserved_tf_families <- as.data.frame(full_all %>% filter(group == "Conserved") %>% select(Family) %>% table())
colnames(conserved_tf_families) <- c("Family", "Genes")

write.table(full_all, "/home/dmpachon/jorge/comparative_cane/results/correlation_analysis/classification_TAP_families.csv", sep = ",", header = T, quote = F, rownames = F)

ggplot(not_conserved_tf_families, aes(x = reorder(Family, -Genes) , y = Genes)) +
  geom_bar(stat = "identity", position = "stack", fill = "black") +
  labs(x = "TAP Family", y = "Genes") +
  theme_bw() +
  theme(text = element_text(family = "Times", size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave(filename = "/home/dmpachon/jorge/comparative_cane/results/correlation_analysis/not_conserved_tf_families.png" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)


ggplot(conserved_tf_families, aes(x = reorder(Family, -Genes) , y = Genes)) +
  geom_bar(stat = "identity", position = "stack", fill = "black") +
  labs(x = "TAP Family", y = "Genes") +
  theme_bw() +
  theme(text = element_text(family = "Times", size=16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave(filename = "/home/dmpachon/jorge/comparative_cane/results/correlation_analysis/conserved_tf_families.svg" ,units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)




#ggplot(full) +
#  geom_point(aes(x = Orthogroup, y = rho, color = Species)) +
#  labs(x = "Orthogroup", y = "Correlation value") +
#  theme_bw() + 
#  theme(text = element_text( size=22),
#                panel.border = element_blank(),
#                panel.grid.major = element_blank(),
#                panel.grid.minor = element_blank(),
#                axis.line = element_line(colour = "black"))

#ggsave(filename = "test.png" ,device = "png", units = "cm", width = 15*1.3, height = 15, dpi = 320, limitsize = F)




#test <- full %>%
#  group_by(Orthogroup) %>%
#  summarize(different_clusters = sum(rho < 0 )) %>%
#  mutate(group = ifelse(different_clusters == 0, "Conserved", "Not Conserved"))


library(igraph)
library(readr)

# Read the TSV file
adjacency_data <- read_tsv("../results/networks/abs_triplets_p60_sorghum_fancy.tsv", 
                           col_names = c("from", "to", "weight"))

# Create an igraph object
graph <- graph_from_data_frame(adjacency_data, directed = FALSE)

# Step 1: Create a set of your selected genes
selected_genes <- full_all %>% filter(Species.y == "Sorghum") %>% select(gene)
# Replace with your actual list of genes

# Step 1: Get the indices of the selected genes in the graph
selected_vertices <- V(graph)[name %in% selected_genes$gene]

# Step 2: Find all neighbors of the selected genes at distance 1
# Including the selected genes themselves
neighbors_list <- neighborhood(graph, order = 1, nodes = selected_vertices)

# Step 3: Flatten the list of neighbors into a single vector of vertex IDs
neighbors_vertices <- unique(unlist(neighbors_list))

# Step 4: Create the subgraph induced by the selected genes and their neighbors
subgraph <- induced_subgraph(graph, vids = selected_vertices)

# Step 5: View basic information about the subgraph
print(subgraph)


write_graph(subgraph, file = "subgraph_sorghum_orthogonal_correlated_tf_gephi.graphml", format = "graphml")















# Step 2: Find all neighbors at distance 1
selected_indices <- which(V(graph)$name %in% selected_genes)

# Step 2: Find all neighbors at distance 1
neighbor_indices <- unique(unlist(lapply(selected_indices, function(i) neighbors(graph, i))))

# Combine selected genes and their neighbors
all_indices <- unique(c(selected_indices, neighbor_indices))

# Step 3: Create the subgraph
subgraph <- induced_subgraph(graph, all_indices)

