library(igraph)
library(dplyr)
library(readr)

# Function to generate a graph with preferential attachment rule, and random edge weights from 0 to 1
generate_graph <- function(num_nodes, m) {
  g <- barabasi.game(num_nodes, m = m, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
  #assign names to nodes since they are not assigned by default and the ids are used as names
  V(g)$name <- V(g)
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
  # use the name of the node in the graph instead of the id
  # edge_data <- get.data.frame(graph)
  # colnames(edge_data) <- c("Start", "End", "Weight")
  edge_data <- as_data_frame(as_edgelist(graph, names = TRUE)) %>%
    rename(Start = V1, End = V2)
  # add weights to the edge data manually
  edge_data$Weight <- E(graph)$weight
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
create_input_graphs <- function(graph, node_conditions, community_data, output_dir){
  # create a subgraph for each community
  # the subgraphs should have the same names as the original graph
  # the subgraphs should have the same node conditions as the original graph
  communities <- unique(community_data$community)
  subgraphs <- list()
  for (community in communities) {
    subgraph <- induced_subgraph(graph, V(graph)[community_data$community == community])
    subgraphs[[community]] <- subgraph
  }

  # create the output directories
  dir.create(paste0(output_dir,"graphs"), showWarnings = FALSE)
  dir.create(paste0(output_dir,"node_conditions"), showWarnings = FALSE)
  dir.create(paste0(output_dir,"node_conditions_discr"), showWarnings = FALSE)
  dir.create(paste0(output_dir,"communities"), showWarnings = FALSE)
  dir.create(paste0(output_dir,"interactions"), showWarnings = FALSE)

  # create a file for each community
  #the files should be:
  # a file for each community graph, has the format of a graph with Start End and Weight
  # a file for the interactions between communities, has the format of the interaction file, that is startType	startNodeName	endType	endNodeName	weight
  # a file for each community,  with the conditions of each node , has the format of the node conditions file, that is name condition
  for (community in communities) {
    print(paste0("saving community: ",community))
    subgraph <- subgraphs[[community]]
    subgraph_edge_data <- generate_edge_data(subgraph)
    subgraph_node_conditions <- node_conditions[as_ids(node_conditions$node_name) %in% as_ids(V(subgraph) ), ]
    # filter the community data to only have the nodes in the subgraph
    subgraph_community_data <- community_data[as_ids(community_data$node_name) %in% as_ids(V(subgraph)), ]
    # write the files, no quotes
    write_tsv(subgraph_edge_data, paste0(output_dir,"graphs/", community, ".tsv"), quote = FALSE)
    write_tsv(subgraph_node_conditions, paste0(output_dir,"node_conditions/", community, ".tsv"), quote = FALSE)
    write_tsv(subgraph_community_data, paste0(output_dir,"communities/", community, ".tsv"), quote = FALSE)

    # also save the node conditions for the community as 0 for Susceptible and 1 for Infectious
    subgraph_node_conditions$conditiondiscr <- ifelse(subgraph_node_conditions$condition == "Susceptible", 0, 1)
    subgraph_node_conditions_discr <- subgraph_node_conditions %>%
      select(node_name, conditiondiscr) %>%
      rename(name = node_name) %>%
      rename(value = conditiondiscr)

    write_tsv(subgraph_node_conditions_discr, paste0(output_dir,"node_conditions_discr/", community, ".tsv"), quote = FALSE)

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
  write_tsv(edges_between_communities, paste0(output_dir,"interactions/interactions.tsv"))
}

create_input_graphs(graph, node_conditions, community_data, "syntheticGraphs/100Nodes/")
# Write data to tsv files
write_tsv(edge_data, "edge_data.tsv")
write_tsv(node_conditions, "node_conditions.tsv")
write_tsv(community_data, "communities.tsv")

