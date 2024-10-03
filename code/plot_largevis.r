library(gridExtra)
library(igraph)
library(tidyverse)
library(readr)
library(wesanderson)
library(viridis)

# Read correlations file 
# Read the TSV file
adjacency_data <- read_tsv("../results/networks/abs_triplets_p60_cane_fancy.tsv", 
                           col_names = c("from", "to", "weight"))
# Create an igraph object
graph <- graph_from_data_frame(adjacency_data, directed = FALSE)
graph <- simplify(graph)

# Step 1: Extract edge list and weights
edge_list <- as.data.frame(get.edgelist(graph, names = FALSE))
vertex_name <- as.data.frame(V(graph))
anot_df <- data.frame(gene = row.names(vertex_name), index = vertex_name$x)

edge_weights <- ifelse(is.null(E(graph)$weight), rep(1, ecount(graph)), E(graph)$weight)

# Step 2: Add weights to the edge list
edge_list$weight <- edge_weights

# Step 3: Adjust node indexing to start from 0
edge_list$V1 <- edge_list$V1 - 1
edge_list$V2 <- edge_list$V2 - 1

# Step 4: Aggregate nodes with their neighbors and weights for LargeVis
knn_format <- aggregate(edge_list$V2 ~ edge_list$V1, data = edge_list, FUN = function(x) paste(x, collapse=" "))
knn_weights <- aggregate(edge_list$weight ~ edge_list$V1, data = edge_list, FUN = function(x) paste(x, collapse=" "))
knn_format$V2_with_weights <- paste(knn_format$`edge_list$V2`, knn_weights$weight)

# Step 5: Write the output in LargeVis format to a file
writeLines(paste(knn_format$`edge_list$V1`, knn_format$V2_with_weights), "graph_for_largevis_sorghum.txt")

## run largevis outside R
# load correlation file
cor <- read.table("../results/correlation_analysis/all_genes_cor.csv", sep = ",", header = T)
# load embeddings
embedings <- read.table("../../LargeVis/Linux/weight_sorghum_largevis.txt", skip = 1, header = F, col.names = c("index", "dim1", "dim2"))
# fix index to start from 1 and not from 0 like largevis needed
embedings$index <- embedings$index + 1
# merge data from correlations and embedings using index
temp <- merge(anot_df, cor, by = "gene")
full <- merge(temp, embedings, by = "index")
## plot embedings from largevis
## choose species: sorghum, sugarcane
species <- "sorghum"
df <- full %>% filter(species == species)
df$correlation <- ifelse(abs(df$rho) < 0.8, 0.1, 1)

plot2 <- ggplot(data = df, aes(x = dim1, y = dim2, color = abs(correlation))) +
	geom_point(size = abs(df$correlation)*3, alpha = abs(df$correlation)) +
	scale_color_gradientn(colours = viridis(250)) +
	theme_void() +
        theme(legend.position="none",
             panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_blank(),
	     panel.background = element_rect(fill = 'white', color = 'white')  	
)
ggsave("../results/correlation_analysis/largevis_cor_combined.png", arrangeGrob(plot, plot2, ncol = 2),  device = "svg", width = 1600*1.5*2, height = 1200*1.5, units = "px")


ggsave("../results/correlation_analysis/largevis_cor_sugarcane.png", device = "png", dpi = 300 , width = 1600*1.5, height = 1200*1.5, units = "px")





 

## 
