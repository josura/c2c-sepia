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

mapkpathway <- parseKGML(tmp)
mapkpathway.expanded <- expandKEGGPathway(mapkpathway)

# Read the pathways data from the file into a KGML file 
pathways <- parseKGML2Graph(tmp)

# Expand the embedded pathways data to a single graph object
mapkGembed <- parseKGMLexpandMaps(tmp)

# Get the nodes and edges of the graph
nodes <- nodes(mapkGembed)
edgesList <- edges(mapkGembed)

# convert nodes names from KEGG to ensemble gene names
nodes.ensemble <- lapply(nodes, function(x) {
    if (startsWith(x, "hsa:")) {
        x <- sub("hsa:", "", x)
        x <- unlist(mapIds(org.Hs.eg.db, keys=x, column="SYMBOL", keytype="ENTREZID"))
    }
    return(x)
})

# extract the entrez names from the nodes.ensemble
nodes.names.entrez <- list()
for (i in 1:length(nodes.ensemble)) {
    nodes.names.entrez <- c(nodes.names.entrez, names(nodes.ensemble[[i]]))
}
# extract the values from the nodes.ensemble
nodes.names.ensemble <- list()
for (i in 1:length(nodes.ensemble)) {
  nodes.names.ensemble <- c(nodes.names.ensemble, unname(nodes.ensemble[[i]]))
}

# convert edges list names from KEGG to ensemble gene names
edgesList.ensemble <- lapply(edgesList, function(x) {
    x <- lapply(x, function(y) {
        if (startsWith(y, "hsa:")) {
            y <- sub("hsa:", "", y)
            y <- nodes.names.ensemble[nodes.names.entrez == y][[1]]
        }
        return(y)
    })
    return(x)
})
names(edgesList.ensemble) <- nodes.names.ensemble


create_apoptosys_graph_inputs <- function(nodes, edgesList){
    # Create the graph
    g <- graph_from_edgelist(edgesList, directed=TRUE)
    # Add the nodes
    V(g)$label <- nodes
    return(g)
}