# testing the code at https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html
import numpy as np
import liana as li
import mudata as mu
import scanpy as sc
import pandas as pd

from matplotlib import pyplot as plt
from plotly import express as px

# Load the data
# prot = sc.read('citeseq_prot.h5ad', backup_url='https://figshare.com/ndownloader/files/47625196')
# rna = sc.read('citeseq_rna.h5ad', backup_url='https://figshare.com/ndownloader/files/47625193')
rna_matriceIPF_filepath = "/home/josura/Projects/ccc/datiGrete/matrice_IPF_agr.txt"
rna_matriceNormal_filepath = "/home/josura/Projects/ccc/datiGrete/matrice_Normal_agr.txt"

#IPF
rna_IPF_pd = pd.read_csv(rna_matriceIPF_filepath, sep="\t")
## first column is the gene names
rna_gene_names_IPF = rna_IPF_pd.iloc[:, 0].values
## all the other columns are the cellnames 
rna_cell_names_IPF = rna_IPF_pd.columns[1:].values
## remove the subscriptions from the cell names
rna_cell_names_IPF = [name.split(".")[0] for name in rna_cell_names_IPF]
## transpose the dataframe to have the cells as rows and genes as columns, since the column names are repeated, use some incremental index
rna_IPF_pd = rna_IPF_pd.T
## deleting the non transposed dataframe
#del rna_IPF_pd
## remove the first row which contains the gene names
rna_IPF_pd.columns = rna_gene_names_IPF
rna_IPF_pd = rna_IPF_pd.iloc[1:, :]
## set the index to an incremental index
rna_IPF_pd.index = pd.RangeIndex(start=0, stop=rna_IPF_pd.shape[0], step=1)
## cast the dataframe to float
rna_IPF_pd = rna_IPF_pd.astype(float)
## create the AnnData object
rna_IPF = sc.AnnData(rna_IPF_pd)
rna_IPF.obs["celltype"] = rna_cell_names_IPF
rna_IPF.obs["condition"] = "IPF"
## create the umap embedding
sc.pp.neighbors(rna_IPF, n_neighbors=10)
sc.tl.umap(rna_IPF)
## show the data
X = list(map(lambda x: x[0], rna_IPF.obsm["X_umap"]))
Y = list(map(lambda x: x[1], rna_IPF.obsm["X_umap"]))
celltypes_IPF = list(rna_IPF.obs["celltype"].values)

fig = px.scatter(x=X, y=Y, color=celltypes_IPF)
fig.show()

# Normal
rna_Normal_pd = pd.read_csv(rna_matriceNormal_filepath, sep="\t")
## first column is the gene names
rna_gene_names_Normal = rna_Normal_pd.iloc[:, 0].values
## all the other columns are the cellnames
rna_cell_names_Normal = rna_Normal_pd.columns[1:].values
## remove the subscriptions from the cell names
rna_cell_names_Normal = [name.split(".")[0] for name in rna_cell_names_Normal]
## transpose the dataframe to have the cells as rows and genes as columns, since the column names
# are repeated, use some incremental index
rna_Normal_pd = rna_Normal_pd.T
## deleting the non transposed dataframe
# del rna_Normal_pd
## remove the first row which contains the gene names
rna_Normal_pd.columns = rna_gene_names_Normal
rna_Normal_pd = rna_Normal_pd.iloc[1:, :]
## set the index to an incremental index
rna_Normal_pd.index = pd.RangeIndex(start=0, stop=rna_Normal_pd.shape[0], step=1)
## cast the dataframe to float
rna_Normal_pd = rna_Normal_pd.astype(float)


# TODO create the mdata object with the rna and the inferred metabolites
