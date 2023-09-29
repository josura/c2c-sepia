library(igraph)
library(dplyr)

# Function to generate a graph with preferential attachment rule
generate_graph <- function(num_nodes, m) {
  g <- erdos.renyi.game(num_nodes, p = m / num_nodes, directed = FALSE)
  return(g)
}

# Generate a graph with preferential attachment (100 nodes, m=2 for preferential attachment)
graph <- generate_graph(100, m = 2)

# Function to assign node conditions (Susceptible or Infectious)
assign_node_conditions <- function(graph, prob_infectious = 0.1) {
  nodes <- V(graph)
  conditions <- sample(c("Susceptible", "Infectious"), length(nodes), replace = TRUE, prob = c(1 - prob_infectious, prob_infectious))
  node_conditions <- data.frame(node_name = nodes, condition = conditions)
  return(node_conditions)
}

# Assign node conditions to the graph
node_conditions <- assign_node_conditions(graph, prob_infectious = 0.1)

# Function to create communities using Louvain community detection
create_communities <- function(graph) {
  communities <- cluster_louvain(graph)
  nodes <- V(graph)
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

# Plot the graph
plot_graph <- function(graph, node_conditions) {
  # Assign node colors based on conditions
  node_colors <- ifelse(node_conditions$condition == "Susceptible", "blue", "red")
  
  # Assign edge widths based on edge weights
  edge_weights <- E(graph)$weight * 5  # Multiply by a factor for better visualization
  
  plot(
    graph,
    layout = layout_with_fr(graph),
    vertex.color = node_colors,
    vertex.size = 6,
    edge.width = edge_weights,
    vertex.label.cex = 0.5
  )
}

# Plot the graph with node colors and edge widths
plot_graph(graph, node_conditions)

# Create the input graphs for single types(communities) and create file for interactions between communities
save.graph.communities.and.interactions <- function(graph, node_conditions){
  # Create a graph for each community
  communities <- unique(node_conditions$community)
  for (community in communities) {
    # Get the nodes in the community
    nodes <- node_conditions %>% filter(community == community) %>% select(node_name)
    nodes <- nodes$node_name
    
    # Create a subgraph for the community
    subgraph <- induced_subgraph(graph, nodes)
    
    # Save the subgraph as a tsv file
    write.graph(subgraph, file = paste0("community_", community, ".graphml"), format = "graphml")
  }
  
  # Create a graph for interactions between communities
  # Get the edges between communities
  edges <- edge_data %>% filter(Start %in% communities & End %in% communities)
  
  # Create a graph with the edges
  interactions <- graph_from_data_frame(edges, directed = FALSE)
  
  # Save the graph as a tsv file
  write.graph(interactions, file = "interactions.graphml", format = "graphml")
}

# Write data to tsv files
write.csv(edge_data, "edge_data.csv", sep = "\t" , row.names = FALSE)
write.csv(node_conditions, "node_conditions.csv", sep = "\t", row.names = FALSE)
write.csv(community_data, "communities.csv", sep = "\t", row.names = FALSE)

