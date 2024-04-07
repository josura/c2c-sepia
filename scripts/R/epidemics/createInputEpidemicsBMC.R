library(igraph)
library(dplyr)
library(readr)

max_number_of_iteration <- 10

# Function to generate a graph with preferential attachment rule, and random edge weights from 0 to 1
generate_graph_barabasi <- function(num_nodes, m) {
  #g <- barabasi.game(num_nodes, m = m, directed = FALSE)
  g <- sample_pa(num_nodes, m = m, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
  #assign names to nodes since they are not assigned by default and the ids are used as names
  V(g)$name <- V(g)
  return(g)
}

# Function to generate a graph with preferential attachment and aging rule, and random edge weights from 0 to 1
generate_graph_barabasi_aging <- function(num_nodes, m, aging) {
  g <- sample_pa_age(n = num_nodes, pa.exp = 1, aging.exp = aging, m = m, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
  #assign names to nodes since they are not assigned by default and the ids are used as names
  V(g)$name <- V(g)
  return(g)
}

# Function to generate a random graph with erdos renyi model, and random edge weights from 0 to 1
generate_graph_erdos <- function(num_nodes, prob) {
  #g <- erdos.renyi.game(num_nodes, prob, directed = FALSE)
  g <- sample_gnp(num_nodes, prob, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
  #assign names to nodes since they are not assigned by default and the ids are used as names
  V(g)$name <- V(g)
  return(g)
}

# Function to generate a deterministic graph with a fixed number of nodes and edges, in a lattice structure
generate_graph_lattice <- function(num_nodes, num_edges) {
  g <- make_lattice(dim = c(num_nodes/6,6), nei = 1, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
  #assign names to nodes since they are not assigned by default and the ids are used as names
  V(g)$name <- V(g)
  return(g)
}

# Function to generate a deterministic regular graph
generate_graph_regular <- function(num_nodes, degree) {
  #g <- make_full_graph(num_nodes) %du% make_full_graph(degree)
  g <- sample_k_regular(num_nodes, degree, directed = FALSE)
  E(g)$weight <- runif(ecount(g), min = 0, max = 1)
  #assign names to nodes since they are not assigned by default and the ids are used as names
  V(g)$name <- V(g)
  return(g)
}

# Function to assign node conditions (Susceptible or Infectious)
assign_node_conditions <- function(graph, prob_infectious = 0.1) {
  nodes <- V(graph)
  conditions <- sample(c("Susceptible", "Infectious"), length(nodes), replace = TRUE, prob = c(1 - prob_infectious, prob_infectious))
  node_conditions <- data.frame(node_name = nodes, condition = conditions)
  return(node_conditions)
}


# Function to create communities using Louvain community detection
create_communities <- function(graph) {
  communities <- cluster_louvain(graph)
  nodes <- V(graph)
  community_data <- data.frame(node_name = nodes, community = communities$membership)
  return(community_data)
}

# Function to generate edge data with edge weights
generate_edge_data <- function(graph) {
  #edge_data <- as_data_frame(as_edgelist(graph, names = TRUE, weights = TRUE))
  # use the name of the node in the graph instead of the id
  # edge_data <- get.data.frame(graph)
  # colnames(edge_data) <- c("Start", "End", "Weight")
  edge_data <- data.frame(as_edgelist(graph, names = TRUE)) %>%
    rename(Start = X1, End = X2)
  # add weights to the edge data manually
  edge_data$Weight <- E(graph)$weight
  return(edge_data)
}

# Plot the graph
plot_graph <- function(graph, node_conditions) {
  # Assign node colors based on conditions
  node_colors <- ifelse(node_conditions$condition == "Susceptible", "white", "blue")
  
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
    # write the files, no quotes, overwrite the files if they exist
    write_tsv(subgraph_edge_data, paste0(output_dir,"graphs/", community, ".tsv"), quote = "none",append = FALSE)
    write_tsv(subgraph_node_conditions, paste0(output_dir,"node_conditions/", community, ".tsv"), quote = "none",append = FALSE)
    write_tsv(subgraph_community_data, paste0(output_dir,"communities/", community, ".tsv"), quote = "none",append = FALSE)
    
    # also save the node conditions for the community as 0 for Susceptible and 1 for Infectious
    subgraph_node_conditions$conditiondiscr <- ifelse(subgraph_node_conditions$condition == "Susceptible", 0, 1)
    subgraph_node_conditions_discr <- subgraph_node_conditions %>%
      select(node_name, conditiondiscr) %>%
      rename(name = node_name) %>%
      rename(value = conditiondiscr)
    
    write_tsv(subgraph_node_conditions_discr, paste0(output_dir,"node_conditions_discr/", community, ".tsv"), quote = "none", append = FALSE)
    
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

  # add random contact times (unique integers for every edge, increasing order,) from 0 to max_number_of_iteration to every edge between communities
  # the contact times for every edge are strings separated by a comma
  contact_times_column <- list()
  for (i in 1:nrow(edges_between_communities)){
    number_of_contacts <- sample(1:max_number_of_iteration, 1)
    contact_times_column[[i]] <- paste(sample(0:(max_number_of_iteration-1), number_of_contacts), collapse = ",")
  }

  edges_between_communities$contactTimes <- unlist(contact_times_column)

  
  # write the file
  write_tsv(edges_between_communities, paste0(output_dir,"interactions/interactions.tsv"), quote = "none",append = FALSE)

}

#create the directory to contain the significance analysis input data
dir.create("significanceAnalysys", showWarnings = FALSE)

# create the directory for the graphs generated with preferential attachment (barabasi model)
dir.create("significanceAnalysys/preferentialAttachment", showWarnings = FALSE)

# create the directory for the graphs generated with erdos renyi model
dir.create("significanceAnalysys/erdosRenyi", showWarnings = FALSE)

# create the directory for the graphs generated with lattice model
dir.create("significanceAnalysys/regular", showWarnings = FALSE)

# create the directory for the graphs generated with preferential attachment and aging model
dir.create("significanceAnalysys/preferentialAttachmentAging", showWarnings = FALSE)

nodes.list <- c(100, 1000, 10000)

for(numNodes in nodes.list){
  # create the directory for the defined nodes in the list nodes.list, with barabasi model, generate 30 graphs
  dir.create(paste0("significanceAnalysys/preferentialAttachment/",numNodes,"nodes"), showWarnings = FALSE)

  # create the directory for the defined nodes in the list nodes.list, with erdos renyi model, generate 30 graphs
  dir.create(paste0("significanceAnalysys/erdosRenyi/",numNodes,"nodes"), showWarnings = FALSE)

  # create the directory for the defined nodes in the list nodes.list, with lattice model, generate 30 graphs
  dir.create(paste0("significanceAnalysys/regular/",numNodes,"nodes"), showWarnings = FALSE)

  # create the directory for the defined nodes in the list nodes.list, with preferential attachment and aging model, generate 30 graphs
  dir.create(paste0("significanceAnalysys/preferentialAttachmentAging/",numNodes,"nodes"), showWarnings = FALSE)

  for (i in 1:30){
    ## PREFERENTIAL ATTACHMENT
    current_dir <- paste0("significanceAnalysys/preferentialAttachment/",numNodes,"nodes/",i,"/")
    dir.create(current_dir, showWarnings = FALSE)
    # Generate a graph with preferential attachment (100 nodes, m=2 for preferential attachment)
    graph <- generate_graph_barabasi(numNodes, m = 2)
    
    # Assign node conditions to the graph
    node_conditions <- assign_node_conditions(graph, prob_infectious = 0.1)
    
    # Create communities based on the graph
    community_data <- create_communities(graph)
    
    # Generate edge data with edge weights
    edge_data <- generate_edge_data(graph)
    
    # Plot the graph with node colors and edge widths
    #plot_graph(graph, node_conditions)
    
    create_input_graphs(graph, node_conditions, community_data, current_dir)
    # Write data to tsv files
    write_tsv(edge_data, paste0(current_dir,"edge_data.tsv"))
    write_tsv(node_conditions, paste0(current_dir,"node_conditions.tsv"))
    write_tsv(community_data, paste0(current_dir,"communities.tsv"))

    print(paste0("saved graph ",i," for ",numNodes," nodes with preferential attachment model."))

    ## PREFERENTIAL ATTACHMENT AGING
    current_dir <- paste0("significanceAnalysys/preferentialAttachmentAging/",numNodes,"nodes/",i,"/")
    dir.create(current_dir, showWarnings = FALSE)
    # Generate a graph with preferential attachment and aging (100 nodes, m=2 for preferential attachment, aging=1 for aging)
    graph <- generate_graph_barabasi_aging(numNodes, m = 2, aging = 1)
    
    # Assign node conditions to the graph
    node_conditions <- assign_node_conditions(graph, prob_infectious = 0.1)

    # Create communities based on the graph
    community_data <- create_communities(graph)

    # Generate edge data with edge weights
    edge_data <- generate_edge_data(graph)

    # Plot the graph with node colors and edge widths
    #plot_graph(graph, node_conditions)
    
    create_input_graphs(graph, node_conditions, community_data, current_dir)
    # Write data to tsv files
    write_tsv(edge_data, paste0(current_dir,"edge_data.tsv"))
    write_tsv(node_conditions, paste0(current_dir,"node_conditions.tsv"))
    write_tsv(community_data, paste0(current_dir,"communities.tsv"))

    print(paste0("saved graph ",i," for ",numNodes," nodes with preferential attachment and aging model."))
    

    ## ERDOS RENYI
    current_dir <- paste0("significanceAnalysys/erdosRenyi/",numNodes,"nodes/",i,"/")
    dir.create(current_dir, showWarnings = FALSE)
    # Generate a graph with erdos renyi model (100 nodes, prob=0.1 for erdos renyi)
    graph <- generate_graph_erdos(numNodes, prob = 0.1)

    # Assign node conditions to the graph
    node_conditions <- assign_node_conditions(graph, prob_infectious = 0.1)

    # Create communities based on the graph
    community_data <- create_communities(graph)

    # Generate edge data with edge weights
    edge_data <- generate_edge_data(graph)

    # Plot the graph with node colors and edge widths
    #plot_graph(graph, node_conditions)

    create_input_graphs(graph, node_conditions, community_data, current_dir)
    # Write data to tsv files
    write_tsv(edge_data, paste0(current_dir,"edge_data.tsv"))
    write_tsv(node_conditions, paste0(current_dir,"node_conditions.tsv"))
    write_tsv(community_data, paste0(current_dir,"communities.tsv"))
  }

  ## Regular graph
  current_dir <- paste0("significanceAnalysys/regular/",numNodes,"nodes/")
  #dir.create(current_dir, showWarnings = FALSE)

  # Generate a graph with regular model (100 nodes, degree=2 for regular)
  graph <- generate_graph_regular(numNodes, degree = 2)

  # Assign node conditions to the graph
  node_conditions <- assign_node_conditions(graph, prob_infectious = 0.1)

  # Create communities based on the graph
  community_data <- create_communities(graph)

  # Generate edge data with edge weights
  edge_data <- generate_edge_data(graph)

  # Plot the graph with node colors and edge widths
  #plot_graph(graph, node_conditions)
  
  create_input_graphs(graph, node_conditions, community_data, current_dir)
  # Write data to tsv files
  write_tsv(edge_data, paste0(current_dir,"edge_data.tsv"))
  write_tsv(node_conditions, paste0(current_dir,"node_conditions.tsv"))
  write_tsv(community_data, paste0(current_dir,"communities.tsv"))
}

