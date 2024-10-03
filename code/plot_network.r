library(igraph)
library(readr)
library(grDevices)  # For the adjustcolor function

graph_plot <- function(adjacency,out){
# Read the TSV file
adjacency_data <- read_tsv(adjacency,
                           col_names = c("from", "to", "weight"))
graph <- graph_from_data_frame(adjacency_data, directed = FALSE)
graph <- simplify(graph)

# Assuming your graph is named 'graph'
# Calculate connected components
comps <- components(graph)

# Select only the largest connected component (or iterate over components as needed)
# For example, we take the largest component for simplicity
largest_comp <- which.max(comps$csize)
subgraph <- induced_subgraph(graph, which(comps$membership == largest_comp))

library(viridis)

# Calculate node degrees
node_degrees <- degree(subgraph)

# Normalize the degree values for coloring (between 0 and 1)
degree_norm <- (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees))

# Use a continuous color palette (viridis in this case)
node_colors <- viridis(length(node_degrees))[rank(degree_norm)]  # Map degrees to colors

# Set the transparency of the edges (alpha = 0.4)
E(subgraph)$color <- adjustcolor("gray", alpha.f = 0.4)

# Save the plot to a file (e.g., PNG or PDF)
svg(filename = out, width = 1600, height = 1200)

# Plot the graph with colored nodes by degree
plot(subgraph, 
     vertex.color = node_colors,   # Apply the continuous color palette based on degree
     vertex.size = 0.5,              # Adjust vertex size if needed
     vertex.label = NA,            # Remove vertex labels for clarity
#     edge.arrow.size = 0.5,        # Adjust edge arrow size for undirected graphs
     layout = layout_nicely)      # Use Fruchterman-Reingold layout (or other)

# Close the file
dev.off()
}
graph_plot("../results/networks/abs_triplets_p60_sorghum_fancy.tsv", "sorghum_network_simplify_degree_colored.svg")
graph_plot("../results/networks/abs_triplets_p60_cane_fancy.tsv", "sugarcane_network_simplify_degree_colored.svg")



