import numpy as np
import liana as li
import mudata as mu
import scanpy as sc

from matplotlib import pyplot as plt
from adjustText import adjust_text


kwargs = {'frameon':False, 'size':1.5, 'img_key':'lowres'}

# Load the data
prot = sc.read('temporalSingleCell/citeseq_prot.h5ad', backup_url='https://figshare.com/ndownloader/files/47625196')
rna = sc.read('temporalSingleCell/citeseq_rna.h5ad', backup_url='https://figshare.com/ndownloader/files/47625193')

mdata = mu.MuData({'rna': rna, 'prot': prot})
# make sure that cell type is accessible
mdata.obs['celltype'] = mdata.mod['rna'].obs['celltype'].astype('category')
# inspect the object
mdata

# show the data
fig = plt.figure(figsize=(10, 10))
sc.pl.embedding(rna, basis='umap', color='celltype', ax=fig.gca(), **kwargs)