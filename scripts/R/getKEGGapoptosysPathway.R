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

subtypes <- lapply(mapkpathway.expanded@edges, function(x) {
    return(x@subtype[[1]]@name)
})

subtypes <- unique(unlist(subtypes))

edges.dataframe <- data.frame(start=character(), end=character(), type=character(), subtype=character(), weight=numeric(), stringsAsFactors=FALSE)

for(i in 1:length(mapkpathway.expanded@edges)){
    edge <- mapkpathway.expanded@edges[[i]]
    start <- sub("hsa:", "", edge@entry1ID)
    end <- sub("hsa:", "", edge@entry2ID)
    subtype <- edge@subtype[[1]]@name
    type <- ""
    weight <- 0.0
    if (subtype == "activation") {
        type <- "PPREL"
        weight <- 1.0
    } else if (subtype == "inhibition") {
        type <- "PPREL"
        weight <- -1.0
    } else if (subtype == "expression"){
        type <- "GEREL"
        weight <- 1.0
    } else if (subtype == "binding/association"){
        type <- "PPREL"
        weight <- 0.0
    } else if (subtype == "dissociation"){
        type <- "PPREL"
        weight <- 0.0
    }
    edges.dataframe <- rbind(edges.dataframe, data.frame(start=start, end=end, type=type, subtype=subtype, weight=weight))
}

create_apoptosys_graph_inputs <- function(nodes, edgesList){
    # Create the graph
    g <- graph_from_edgelist(edgesList, directed=TRUE)
    # Add the nodes
    V(g)$label <- nodes
    return(g)
}