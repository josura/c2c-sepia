library(dplyr)


edges.all <- read.csv("testdata/edges.tsv",sep = "\t")
nodes <- read.csv("testdata/nodes.tsv",sep = "\t")


edges.filtered <- edges.all[edges.all$X.Start %in% nodes$X.Id & edges.all$End %in% nodes$X.Id,] 
write.table(edges.filtered,"testdata/edges-filtered.tsv",sep = "\t",quote = FALSE)
