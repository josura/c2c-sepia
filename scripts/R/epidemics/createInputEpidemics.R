library(igraph)
library(dplyr)

# Function to generate a graph with preferential attachment rule, and random edge weights from 0 to 1
generate_graph <- function(num_nodes, m) {
  g <- barabasi.game(num_nodes, m = m, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
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
  #edge_data <- as_data_frame(as_edgelist(graph, names = TRUE, weights = TRUE))
  edge_data <- get.data.frame(graph)
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
    vertex.size = 3,
    edge.width = edge_weights,
    vertex.label.cex = 0.5
  )
}

# Plot the graph with node colors and edge widths
plot_graph(graph, node_conditions)

# Create the input graphs for single types(communities) and create file for interactions between communities
# The interaction file should have the following format:
# startType	startNodeName	endType	endNodeName	weight
create_input_graphs <- function(graph, node_conditions, community_data){
  # create a subgraph for each community
  communities <- unique(community_data$community)
  subgraphs <- list()
  for (community in communities) {
    subgraph <- induced_subgraph(graph, V(graph)[community_data$community == community])
    subgraphs[[community]] <- subgraph
  }

  # create a file for each community
  #the files should be:
  # a file for each community graph, has the format o a graph with Start End and Weight
  # a file for the interactions between communities, has the format of the interaction file, that is startType	startNodeName	endType	endNodeName	weight
  # a file for each community,  with the conditions of each node , has the format of the node conditions file, that is name condition
  for (community in communities) {
    print(paste0("saving community: ",community))
    subgraph <- subgraphs[[community]]
    subgraph_edge_data <- generate_edge_data(subgraph)
    subgraph_node_conditions <- node_conditions[as_ids(node_conditions$node_name) %in% as_ids(V(subgraph)), ]
    subgraph_community_data <- community_data[as_ids(community_data$node_name) %in% as_ids(V(subgraph)), ]
    write.csv(subgraph_edge_data, paste0("edge_data_community_", community, ".tsv"), sep = "\t", row.names = FALSE)
    write.csv(subgraph_node_conditions, paste0("node_conditions_community_", community, ".tsv"), sep = "\t", row.names = FALSE)
    write.csv(subgraph_community_data, paste0("communities_community_", community, ".tsv"), sep = "\t", row.names = FALSE)
  }

  # create a file for the interactions between communities
  # the file should be:
  # startType	startNodeName	endType	endNodeName	weight
  # where startType and endType are the community names
  # startNodeName and endNodeName are the node names
  # weight is the edge weight

  # get the edges between communities
  edges_between_communities <- get.data.frame(graph, what = "edges")
  # filter out edges that are not between communities, that is the edge is not between nodes that belong to different communities, so the edge is between nodes that belong to the same community
  edges_between_communities <- edges_between_communities[community_data[edges_between_communities$from, ]$community != community_data[edges_between_communities$to, ]$community, ]

  community_data.asids <- community_data
  community_data.asids$node_name <- as_ids(community_data.asids$node_name)
  
  # get the community names for each node
  edges_between_communities <- edges_between_communities %>%
    #mutate_at(c('from', 'to'), as.character) %>%
    left_join(community_data.asids, by = c("from" = "node_name")) %>%
    rename(startType = community) %>%
    left_join(community_data.asids, by = c("to" = "node_name")) %>%
    rename(endType = community)

  # get the edge weights
  edges_between_communities <- edges_between_communities %>%
    rename(startNodeName = from, endNodeName = to)
  
  # write the file
  write.csv(edges_between_communities, "interactions.tsv", sep = "\t", row.names = FALSE)
}

create_input_graphs(graph, node_conditions, community_data)
# Write data to tsv files
write.csv(edge_data, "edge_data.csv", sep = "\t" , row.names = FALSE)
write.csv(node_conditions, "node_conditions.csv", sep = "\t", row.names = FALSE)
write.csv(community_data, "communities.csv", sep = "\t", row.names = FALSE)

