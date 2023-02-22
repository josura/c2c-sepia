install.packages('Seurat')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

tisch.CCI <- read.csv("tisch/ESCA_GSE160269_CCI_LR.tsv",sep = "\t")
tisch.logFoldChange <- read.csv("tisch/ESCA_GSE160269_AllDiffGenes_table.tsv",sep = "\t")
tisch.metainfo <- read.csv("tisch/ESCA_GSE160269_CellMetainfo_table.tsv",sep = "\t")

dim(table(tisch.logFoldChange$Gene))