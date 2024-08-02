require(dplyr)

.read.scdbtab=function(amp_batches.7d, scdb_path=".")
{
  files=paste(scdb_path,"/",amp_batches.7d,".txt",sep="")
  l=lapply(files,FUN=function(x){read.delim(x,header=T,check.names=F,stringsAsFactors=F,row.names=1)})
  a=do.call("cbind", l)
  return(a)
}

#batches.7d = c(paste("AB", c(4005:4006,4082,4083), sep=""))
batches.7d = c(paste("AB", c(1185,1186), sep=""))  #PN 7d
batches.10h = c(paste("AB", c(1245,1246), sep=""))  #PN 10h

umidir = "/home/josura/Projects/ccc/datiIdo/AB_UMI"

umis_all.7d  = .read.scdbtab(batches.7d, scdb_path = umidir)
umis_all.10h  = .read.scdbtab(batches.10h, scdb_path = umidir)

umis_all.7d.filtered = umis_all.7d
umis_all.10h.filtered = umis_all.10h

#  load metadata from tsv
metadata = read.delim("/home/josura/Projects/ccc/datiIdo/metadata_2.tsv", header=T, check.names=F, stringsAsFactors=F)
well_metadata = read.delim("/home/josura/Projects/ccc/datiIdo/lung_map.txt", header=T, check.names=F, stringsAsFactors=F)

# filter metadata to the well metadata
metadata.filtered = metadata[metadata$well %in% well_metadata$well,]
# remove columns with NA values
metadata.filtered = metadata.filtered[, colSums(is.na(metadata.filtered)) == 0]
# join metadata.filtered with well_metadata
metadata.filtered = merge(metadata.filtered, well_metadata, by="well")

# 7 days PN

# change column names in umis_all.7d.filtered to match metadata.filtered by well, changing well to cell_type
# first select the columns in umis_all.7d.filtered that are in metadata.filtered$well
umis_all.7d.filtered.columns = umis_all.7d.filtered[, colnames(umis_all.7d.filtered) %in% metadata.filtered$well]
# change the column names to match metadata.filtered$cell_type
colnames(umis_all.7d.filtered.columns) = metadata.filtered$cell_type[match(colnames(umis_all.7d.filtered.columns), metadata.filtered$well)]

# filter metadata.filtered rows to the columns in umis_all.7d.filtered via the well column (colnames in umis_all.7d.filtered)
metadata.filtered.onlyPresentRowsinBatch7d = metadata.filtered[metadata.filtered$well %in% colnames(umis_all.7d.filtered),]

# filter the umis_all.7d.filtered to the metadata.filtered.onlyPresentRowsinBatch7d, selecting the columns with names in metadata.filtered.onlyPresentRowsinBatch7d$well
umis_all.7d.filtered.onlyPresentRowsinBatch7d = umis_all.7d.filtered[, colnames(umis_all.7d.filtered) %in% metadata.filtered.onlyPresentRowsinBatch7d$well]

# load the umis_all.7d.filtered.columns into a seurat object
library(Seurat)
seurat_obj.7d.nocelltypes = CreateSeuratObject(counts = umis_all.7d.filtered.onlyPresentRowsinBatch7d, project = "ABlung")
# add cell types to the seurat object
seurat_obj.7d = AddMetaData(object = seurat_obj.7d.nocelltypes, metadata = metadata.filtered.onlyPresentRowsinBatch7d, col.name = "cell_type")
# add sample id to seurat object
seurat_obj.7d = AddMetaData(object = seurat_obj.7d, metadata = metadata.filtered.onlyPresentRowsinBatch7d, col.name = "Amp_batch_ID")

# save the seurat object
saveRDS(seurat_obj.7d, file="/home/josura/Projects/ccc/datiIdo/seurat_obj.7d_MGI.rds")


# 8 weeks PN

# change column names in umis_all.10h.filtered to match metadata.filtered by well, changing well to cell_type
# first select the columns in umis_all.10h.filtered that are in metadata.filtered$well
umis_all.10h.filtered.columns = umis_all.10h.filtered[, colnames(umis_all.10h.filtered) %in% metadata.filtered$well]
# change the column names to match metadata.filtered$cell_type
colnames(umis_all.10h.filtered.columns) = metadata.filtered$cell_type[match(colnames(umis_all.10h.filtered.columns), metadata.filtered$well)]

# filter metadata.filtered rows to the columns in umis_all.10h.filtered via the well column (colnames in umis_all.10h.filtered)
metadata.filtered.onlyPresentRowsinBatch10h = metadata.filtered[metadata.filtered$well %in% colnames(umis_all.10h.filtered),]

# filter the umis_all.10h.filtered to the metadata.filtered.onlyPresentRowsinBatch10h, selecting the columns with names in metadata.filtered.onlyPresentRowsinBatch10h$well
umis_all.10h.filtered.onlyPresentRowsinBatch10h = umis_all.10h.filtered[, colnames(umis_all.10h.filtered) %in% metadata.filtered.onlyPresentRowsinBatch10h$well]

# load the umis_all.10h.filtered.columns into a seurat object
seurat_obj.10h.nocelltypes = CreateSeuratObject(counts = umis_all.10h.filtered.onlyPresentRowsinBatch10h, project = "ABlung")

# add cell types to the seurat object
seurat_obj.10h = AddMetaData(object = seurat_obj.10h.nocelltypes, metadata = metadata.filtered.onlyPresentRowsinBatch10h, col.name = "cell_type")
# add sample id to seurat object
seurat_obj.10h = AddMetaData(object = seurat_obj.10h, metadata = metadata.filtered.onlyPresentRowsinBatch10h, col.name = "Amp_batch_ID")

# save the seurat object
saveRDS(seurat_obj.10h, file="/home/josura/Projects/ccc/datiIdo/seurat_obj.10h_MGI.rds")

# show the seurat object, with the cell types and UMAP
seurat_obj.7d <- readRDS("/home/josura/Projects/ccc/datiIdo/seurat_obj.7d_MGI.rds")
    
seurat_obj.7d.umap <- NormalizeData(seurat_obj.7d)
seurat_obj.7d.umap <- FindVariableFeatures(seurat_obj.7d.umap, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(seurat_obj.7d.umap)
seurat_obj.7d.umap <- ScaleData(seurat_obj.7d.umap, features = all.genes)
seurat_obj.7d.umap <- RunPCA(seurat_obj.7d.umap, features = VariableFeatures(object = seurat_obj.7d.umap))
seurat_obj.7d.umap <- RunUMAP(seurat_obj.7d.umap, dims = 1:17)
Idents(seurat_obj.7d.umap) <- "cell_type"
DimPlot(seurat_obj.7d.umap, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 10h

seurat_obj.10h <- readRDS("/home/josura/Projects/ccc/datiIdo/seurat_obj.10h_MGI.rds")

seurat_obj.10h.umap <- NormalizeData(seurat_obj.10h)
seurat_obj.10h.umap <- FindVariableFeatures(seurat_obj.10h.umap, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(seurat_obj.10h.umap)
seurat_obj.10h.umap <- ScaleData(seurat_obj.10h.umap, features = all.genes)
seurat_obj.10h.umap <- RunPCA(seurat_obj.10h.umap, features = VariableFeatures(object = seurat_obj.10h.umap))
seurat_obj.10h.umap <- RunUMAP(seurat_obj.10h.umap, dims = 1:17)
Idents(seurat_obj.10h.umap) <- "cell_type"
DimPlot(seurat_obj.10h.umap, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
