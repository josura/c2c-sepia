install.packages('Seurat')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

tisch.CCI <- read.csv("tisch/ESCA_GSE160269_CCI_LR.tsv",sep = "\t")
tisch.logFoldChange <- read.csv("tisch/ESCA_GSE160269_AllDiffGenes_table.tsv",sep = "\t")
tisch.metainfo <- read.csv("tisch/ESCA_GSE160269_CellMetainfo_table.tsv",sep = "\t")

dim(table(tisch.logFoldChange$Gene))

library(dplyr)

table(tisch.logFoldChange$Celltype..minor.lineage.)
tisch.dataframe <- tisch.logFoldChange
cell.list <- unique(tisch.dataframe$Celltype..minor.lineage.)
cell.logfold.list <- list()

cell.name <- cell.list[2]
sub.tisch.dataframe <- tisch.dataframe %>%
  filter(Celltype..minor.lineage. == cell.name)


gene.names.noversion <- lapply(sub.tisch.dataframe$Gene,  sub, pattern = "\\.\\d+$", replacement = "")
rownames(sub.tisch.dataframe) <- gene.names.noversion
logfold.embed <- data.frame(t(unlist(sub.tisch.dataframe$log2FC)),row.names = cell.name)
colnames(logfold.embed) <- gene.names.noversion


create.logfold.input <- function(tisch.dataframe){
  #TODO control input
  cell.list <- unique(tisch.dataframe$Celltype..minor.lineage.)
  cell.logfold.list <- list()
  for (cell.name in cell.list) {
    sub.tisch.dataframe <- tisch.dataframe %>%
      filter(Celltype..minor.lineage. == cell.name)
    #gene.names.noversion <- lapply(sub.tisch.dataframe$Gene,  gsub, pattern = "\"", replacement = "")
    gene.names.noversion <- lapply(sub.tisch.dataframe$Gene,  sub, pattern = "\\.\\d+$", replacement = "")
    sub.tisch.dataframe$Genefiltered <- gene.names.noversion
    sub.tisch.dataframe.filtered <- sub.tisch.dataframe %>%
      group_by(Genefiltered) %>%
      #arrange(log2FC) %>%
      slice(1) %>% ungroup
    
    gene.names.noversion <- sub.tisch.dataframe.filtered$Genefiltered
    logfold.embed <- data.frame(t(unlist(sub.tisch.dataframe.filtered$log2FC)),row.names = cell.name)
    
    colnames(logfold.embed) <- gene.names.noversion
    #logfold.embed$gene.name <- cell.name
    cell.logfold.list[[length(cell.logfold.list)+1]] <- logfold.embed
    print(paste("logFolds done for cell",cell.name,sep = ": "))
  }
  merged.embeddings <- cell.logfold.list[[1]]
  cell.logfold.list <- cell.logfold.list[-1]
  for (path.embed in cell.logfold.list) {
    #merged.embeddings <- bind_rows(merged.embeddings,path.embed)
    merged.embeddings[setdiff(names(path.embed), names(merged.embeddings))] <- NA
    path.embed[setdiff(names(merged.embeddings), names(path.embed))] <- NA
    merged.embeddings <- bind_rows(merged.embeddings, path.embed)
    print(paste("merge done for cell",rownames(path.embed),sep = ": "))
  }
  merged.embeddings[is.na(merged.embeddings)] <- 0
  merged.embeddings
}

test <- create.logfold.input(tisch.logFoldChange)
write.table(t(test),"logfoldchangeMatrix/tisch_ESCA_GSE160269-logfold-noversion.tsv",sep = "\t",quote = FALSE)

devtools::install_github("sqjin/CellChat")

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)