if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("KEGGgraph")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

library(KEGGgraph)

# libraries to convert entrez to ensemble
library(AnnotationDbi)
library(org.Hs.eg.db)

# Get the KEGG pathway data
tmp <- tempfile()
retrieveKGML("04210", organism="hsa", destfile=tmp, method="auto", quiet=TRUE)

# Read the pathways data from the file into a KGML file 
pathways <- parseKGML2Graph(tmp)

# Expand the embedded pathways data to a single graph object
mapkGembed <- parseKGMLexpandMaps(tmp)

# Get the nodes and edges of the graph
nodes.apoptosys <- nodes(mapkGembed)
edgesList.apoptosys <- edges(mapkGembed)

#remove the hsa prefix from the nodes
nodes.apoptosys.noorg <- gsub("hsa:", "", nodes.apoptosys)

# Read edges.tsv file
edges <- read.table("../../resources/graphs/metapathwayNew/edges.tsv", header = TRUE, sep = "\t")

library(readr)

nodes <- read_tsv("../../resources/graphs/metapathwayNew/nodes.tsv")

library(dplyr)

apoptosys.nodes.filtered <- nodes %>% filter(Id %in% nodes.apoptosys.noorg)


# Print the subgraph nodes
print(subgraph_nodes)

# Print the subgraph edges
print(subgraph_edges)
