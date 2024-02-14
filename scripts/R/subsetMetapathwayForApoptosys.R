# Read edges.tsv file
edges <- read.table("../../data/testdata/edges.tsv", header = TRUE, sep = "\t")

# Read nodes.tsv file
nodes <- read.table("../../data/testdata/nodes.tsv", header = TRUE, sep = "\t")

# Define the string to match in the Aliases feature
search_string <- "your_search_string"

# Get the subset of nodes that contain the search string in the Aliases feature
subgraph_nodes <- subset(nodes, grepl(search_string, Aliases))

# Get the subset of edges that connect the subgraph nodes
subgraph_edges <- subset(edges, from %in% subgraph_nodes$ID | to %in% subgraph_nodes$ID)

# Print the subgraph nodes
print(subgraph_nodes)

# Print the subgraph edges
print(subgraph_edges)
