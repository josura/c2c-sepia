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
# pathways <- parseKGML2Graph(tmp)

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

# convert edges list names from KEGG to ensemble gene names
edgesList.ensemble <- lapply(edgesList, function(x) {
    x <- lapply(x, function(y) {
        if (startsWith(y, "hsa:")) {
            y <- sub("hsa:", "", y)
            y <- unlist(mapIds(org.Hs.eg.db, keys=y, column="SYMBOL", keytype="ENTREZID"))
        }
        return(y)
    })
    return(x)
})

# Plot the graph
plot(mapkGembed)
