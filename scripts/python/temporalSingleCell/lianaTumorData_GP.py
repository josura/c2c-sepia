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
rna_gene_names = rna_IPF_pd.iloc[:, 0].values
## all the other columns are the cellnames 
rna_cell_names = rna_IPF_pd.columns[1:].values
## remove the subscriptions from the cell names
rna_cell_names = [name.split(".")[0] for name in rna_cell_names]
# transpose the dataframe to have the cells as rows and genes as columns, since the column names are repeated, use some incremental index
rna_IPF_pd = rna_IPF_pd.T
## deleting the non transposed dataframe
#del rna_IPF_pd
# remove the first row which contains the gene names
rna_IPF_pd.columns = rna_gene_names
rna_IPF_pd = rna_IPF_pd.iloc[1:, :]
# set the index to an incremental index
rna_IPF_pd.index = pd.RangeIndex(start=0, stop=rna_IPF_pd.shape[0], step=1)
# cast the dataframe to float
rna_IPF_pd = rna_IPF_pd.astype(float)
# create the AnnData object
rna_IPF = sc.AnnData(rna_IPF_pd)



# TODO create the mdata object with the rna and the inferred metabolites
