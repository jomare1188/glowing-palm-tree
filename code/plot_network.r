library(igraph)
library(readr)

# Read the TSV file
adjacency_data <- read_tsv("../results/networks/abs_triplets_p60_sorghum_fancy.tsv",
                           col_names = c("from", "to", "weight"))
graph <- graph_from_data_frame(adjacency_data, directed = FALSE)

library(viridis)

# Calculate node degrees
node_degrees <- degree(graph)

# Normalize the degree values for coloring (between 0 and 1)
degree_norm <- (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees))

# Use a continuous color palette (viridis in this case)
node_colors <- viridis(length(node_degrees))[rank(degree_norm)]  # Map degrees to colors

# Save the plot to a file (e.g., PNG or PDF)
png(filename = "sorghum_network_degree_colored.png", width = 1600, height = 1200)

# Plot the graph with colored nodes by degree
plot(graph, 
     vertex.color = node_colors,   # Apply the continuous color palette based on degree
     vertex.size = 1,              # Adjust vertex size if needed
     vertex.label = NA,            # Remove vertex labels for clarity
#     edge.arrow.size = 0.5,        # Adjust edge arrow size for undirected graphs
     layout = layout_nicely)      # Use Fruchterman-Reingold layout (or other)

# Close the file
dev.off()
