import numpy as np
import liana as li
import mudata as mu
import scanpy as sc

from matplotlib import pyplot as plt
from plotly import express as px

# Load the data
# prot = sc.read('temporalSingleCell/citeseq_prot.h5ad', backup_url='https://figshare.com/ndownloader/files/47625196')
# rna = sc.read('temporalSingleCell/citeseq_rna.h5ad', backup_url='https://figshare.com/ndownloader/files/47625193')
prot = sc.read('temporalSingleCell/citeseq_prot.h5ad')
rna = sc.read('temporalSingleCell/citeseq_rna.h5ad')

mdata = mu.MuData({'rna': rna, 'prot': prot})
# make sure that cell type is accessible
mdata.obs['celltype'] = mdata.mod['rna'].obs['celltype'].astype('category')
# inspect the object
mdata

# show the data
X = list(map(lambda x: x[0], rna.obsm["X_umap"]))
Y = list(map(lambda x: x[1], rna.obsm["X_umap"]))
celltypes = list(rna.obs["celltype"].values)

fig = px.scatter(x=X, y=Y, color=celltypes)
fig.show()