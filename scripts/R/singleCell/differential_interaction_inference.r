library(CellChat)
library(Seurat)
library(readr)

# load seurat object
seurat_file.7d = "/home/josura/Projects/ccc/datiIdo/seurat_obj.7d_MGI.rds" #PN 7d
seurat_file.10h = "/home/josura/Projects/ccc/datiIdo/seurat_obj.10h_MGI.rds" #PN 10h

seurat_obj.7d = readRDS(seurat_file.7d)
seurat_obj.10h = readRDS(seurat_file.10h)

# merge the Seurat objects, adding a new column to the metadata to distinguish the two datasets
seurat_obj.7d$batchTime <- "7d"
seurat_obj.10h$batchTime <- "10h"

seurat_obj <- merge(seurat_obj.7d, seurat_obj.10h)
seurat_obj <- JoinLayers(seurat_obj)

seurat_obj <- NormalizeData(seurat_obj)
data.input <- seurat_obj[["RNA"]]$data # normalized data matrix
Idents(seurat_obj) <- "cell_type"
labels <- Idents(seurat_obj)
meta <- data.frame(labels = labels, row.names = names(labels), samples = seurat_obj@meta.data$Amp_batch_ID) # create a dataframe of the cell labels

# differential analysis between the two time points and the common cell types
# filter the uncommon cell types
intersected_cell_types <- intersect(levels(as.factor(meta$labels[meta$samples == "7d"])), levels(as.factor(meta$labels[meta$samples == "10h"])))

# all genes
all_genes <- rownames(data.input)

#load the metapathway files
#the edges file has the following columns separated by \t: Source	Target	Type	Subtype	Weight
metapathwayEdgesFile <- "~/Projects/ccc/datiIdo/inputs/mmu_metapathway/edges_mmu.txt"
#the nodes file has the following columns separated by \t: Id	Name	Type	Aliases
nodesInfoFile <- "~/Projects/ccc/datiIdo/inputs/mmu_metapathway/nodes_mmu.txt"


metaPathwayEdges <- read_tsv(metapathwayEdgesFile)
nodesinfo <- read_tsv(nodesInfoFile)

# create map that maps The Id in nodesInfoFile to the Name in nodesInfoFile, this map will be used to change the source and target in the metapathwayEdgesFile
nodesInfoMap <- nodesinfo %>% select(Id, Name) %>% as.data.frame()
nodesInfoMap <- nodesInfoMap[!duplicated(nodesInfoMap$Id),]
rownames(nodesInfoMap) <- nodesInfoMap$Id
nodesInfoMap$Id <- NULL

# change the source and target in the metapathwayEdgesFile
metaPathwayEdges$Source <- nodesInfoMap[as.character(metaPathwayEdges$Source),]
metaPathwayEdges$Target <- nodesInfoMap[as.character(metaPathwayEdges$Target),]
# remove duplicate rows on source,target and weight values
metaPathwayEdges <- metaPathwayEdges[!duplicated(metaPathwayEdges),]
# remove edges with weight equal to 0
metaPathwayEdges <- metaPathwayEdges[metaPathwayEdges$Weight != 0,]
# filter nodes that are not in the edges
nodesInfo.filtered <- nodesinfo[nodesinfo$Name %in% c(metaPathwayEdges$Source, metaPathwayEdges$Target),] 

# save the metapathway and the nodes info
write_tsv(metaPathwayEdges,"differentialInputs/metapathwayEdges.tsv")
write_tsv(nodesInfo.filtered,"differentialInputs/nodesInfo.tsv")


# # DC differential analysis
# seurat_obj.filtered.DC <- seurat_obj[,seurat_obj@meta.data$cell_type %in% intersected_cell_types]
# seurat_obj.filtered.DC <- seurat_obj.filtered.DC[,seurat_obj.filtered.DC@meta.data$cell_type == "DC"]
# Idents(seurat_obj.filtered.DC) <- "batchTime"
# diff_res.DC <- FindMarkers(seurat_obj.filtered.DC, ident.1 = "7d", ident.2 = "10h", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.1)
# diff_res.DC$gene_names <- rownames(diff_res.DC) 
# # column names are "p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj"  "gene_names"
# 
# 
# # join all genes with the differential analysis results, p_val, p_val_adj, pct.1 and pct.2 will be -1 for the genes that are not differentially expressed, while avg_log2FC will be 0
# diff_res.DC_all_genes <- data.frame(p_val = -1, avg_log2FC = 0, pct.1 = -1, pct.2 = -1, p_val_adj = -1, gene_names = all_genes)
# # replace the values of the differentially expressed genes
# diff_res.DC_all_genes[match(diff_res.DC$gene_names, diff_res.DC_all_genes$gene_names),] <- diff_res.DC

# filter the genes in seurat_obj, only the genes in nodesInfo.filtered will be kept
seurat_obj.clean = seurat_obj
Idents(seurat_obj.clean) <- "cell_type"
genesToKeep <- rownames(seurat_obj.clean)[rownames(seurat_obj.clean) %in% nodesInfo.filtered$Name]
seurat_obj.filtered <- seurat_obj.clean[genesToKeep,]

createInputsFilesFromSeuratobject <- function(seurat_obj, outputFolder){
    cell_types <- unique( Idents(seurat_obj) )
    # all genes
    all_genes <- rownames(seurat_obj)
    for (cell_type in cell_types){
        seurat_obj.filtered <- seurat_obj[ , Idents(seurat_obj) == cell_type ]
        Idents(seurat_obj.filtered) <- "batchTime"
        diff_res <- FindMarkers(seurat_obj.filtered, ident.1 = "7d", ident.2 = "10h", test.use = "MAST", logfc.threshold = 0.25, min.pct = 0.1, min.cells.group = 1)
        diff_res$gene_names <- rownames(diff_res)
        # join all genes
        diff_res_all_genes <- data.frame(p_val = -1, avg_log2FC = 0, pct.1 = -1, pct.2 = -1, p_val_adj = -1, gene_names = all_genes)
        # replace the values of the differentially expressed genes
        diff_res_all_genes[match(diff_res$gene_names, diff_res_all_genes$gene_names),] <- diff_res
        # save the results in a file, they will not be used as inputs but they will be used after the simulation to select the significant differentially expressed genes
        write.table(diff_res_all_genes, file = paste0(outputFolder, "/diffExprTables/", cell_type, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        # save the inputs, only names and values(log2FC) are needed
        input.table <- data.frame(name = diff_res_all_genes$gene_names, value = diff_res_all_genes$avg_log2FC)
        write.table(input.table, file = paste0(outputFolder, "/inputValues/", cell_type, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
}

createInputsFilesFromSeuratobject(seurat_obj.filtered, "/home/josura/Projects/ccc/datiIdo/differentialInputs")

# inference of differential interactions
# the interactions will be inferred for the whole seurat object

# filter the genes in seurat_obj, only the genes in nodesInfo.filtered will be kept
seurat_obj.clean = seurat_obj
Idents(seurat_obj.clean) <- "cell_type"
genesToKeep <- rownames(seurat_obj.clean)[rownames(seurat_obj.clean) %in% nodesInfo.filtered$Name]
seurat_obj.filtered <- seurat_obj.clean[genesToKeep,] 
seurat_obj.filtered <- NormalizeData(seurat_obj.filtered)
data.input <- seurat_obj.filtered[["RNA"]]$data # normalized data matrix
labels <- Idents(seurat_obj.filtered)
meta <- data.frame(labels = labels, row.names = names(labels), samples = seurat_obj.filtered@meta.data$Amp_batch_ID) # create a dataframe of the cell labels

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)


# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")
## Users can filter out the cell-cell communication if there are only few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 2)

df.net <- subsetCommunication(cellchat)

# save the inferred cellular communication network as an object
saveRDS(df.net, file = "differentialInputs/cellchat_inferred_network_filtered.rds")

# load the inferred cellular communication network
df.net <- readRDS("differentialInputs/cellchat_inferred_network_filtered.rds")



# create the interaction input file
createInteractionsFile <- function(interactionNetwork,outputDir){
    #create a file with the interactions
    # the file has five columns
    # startType, startNodeName , endType, endNodeName, weight
    # weight is the communication probability
    # the file is saved in the outputDir
    startTypes <- interactionNetwork$source
    # convert factors to characters
    startTypes <- as.character(startTypes)
    endTypes <- interactionNetwork$target
    # convert factors to characters
    endTypes <- as.character(endTypes)
    startNodeNames <- interactionNetwork$ligand
    endNodeNames <- interactionNetwork$receptor
    weights <- interactionNetwork$prob
    interactionDF <- data.frame(startType = startTypes, startNodeName = startNodeNames, endType = endTypes, endNodeName = endNodeNames, weight = weights)
    # change node names since some of them are in all caps, while official names have only the first letter in caps
    interactionDF$startNodeName <- tolower(interactionDF$startNodeName)
    interactionDF$endNodeName <- tolower(interactionDF$endNodeName)
    # make the first letter of each word in the node name uppercase
    interactionDF$startNodeName <- sapply(interactionDF$startNodeName, function(x) {
        paste0(toupper(substring(x, 1, 1)), substring(x, 2))
    })
    interactionDF$endNodeName <- sapply(interactionDF$endNodeName, function(x) {
        paste0(toupper(substring(x, 1, 1)), substring(x, 2))
    })
    # in the case of double interactions between one ligand and more than one receptor(identified by the two or more receptors contatenated with _ in the receptor feature(endNodeName in the dataframe) in the interaction network), the interaction is split in two or more rows, with the same ligand and other features
    splittedReceptors <- strsplit(interactionDF$endNodeName,"_")
    interactionDF.splitted <- transform(interactionDF[rep(seq_len(nrow(interactionDF)), lengths(splittedReceptors)),],c = unlist(splittedReceptors))
    interactionDF.splitted$fullname <- interactionDF.splitted$endNodeName
    interactionDF.splitted$endNodeName <-  interactionDF.splitted$c
    interactionDF.splitted$c <- NULL
    # endNodeName should also have the first letter of each word in uppercase
    interactionDF.splitted$endNodeName <- sapply(interactionDF.splitted$endNodeName, function(x) {
        paste0(toupper(substring(x, 1, 1)), substring(x, 2))
    })

    interactionFile <- paste0(outputDir,"interactions.tsv")
    write_tsv(interactionDF.splitted,interactionFile, quote = "none", append = FALSE)
  
}

createInteractionsFile(df.net, "/home/josura/Projects/ccc/datiIdo/differentialInputs/interactions/")
