library(igraph)
library(dplyr)

# Function to generate a graph with preferential attachment rule
generate_graph <- function(num_nodes, m) {
  g <- erdos.renyi.game(num_nodes, p = m / num_nodes, directed = TRUE)
  return(g)
}

# Generate a graph with preferential attachment (100 nodes, m=2 for preferential attachment)
graph <- generate_graph(100, m = 2)

# Function to assign node conditions (Susceptible or Infectious)
assign_node_conditions <- function(graph, prob_infectious = 0.1) {
  nodes <- V(graph)$name
  conditions <- sample(c("Susceptible", "Infectious"), length(nodes), replace = TRUE, prob = c(1 - prob_infectious, prob_infectious))
  node_conditions <- data.frame(node_name = nodes, condition = conditions)
  return(node_conditions)
}

# Assign node conditions to the graph
node_conditions <- assign_node_conditions(graph, prob_infectious = 0.1)

# Function to create communities using Louvain community detection
create_communities <- function(graph) {
  communities <- cluster_louvain(graph)
  nodes <- V(graph)$name
  community_data <- data.frame(node_name = nodes, community = communities$membership)
  return(community_data)
}

# Create communities based on the graph
community_data <- create_communities(graph)

# Function to generate edge data with edge weights
generate_edge_data <- function(graph) {
  edge_data <- as_data_frame(as_edgelist(graph, names = TRUE, weights = TRUE))
  colnames(edge_data) <- c("Start", "End", "Weight")
  return(edge_data)
}

# Generate edge data with edge weights
edge_data <- generate_edge_data(graph)

# Write data to files
write.csv(edge_data, "edge_data.csv", row.names = FALSE)
write.csv(node_conditions, "node_conditions.csv", row.names = FALSE)
write.csv(community_data, "communities.csv", row.names = FALSE)

