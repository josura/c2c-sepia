# testing the code at https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html
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

# Obtain a ligand-receptor resource of interest
resource = li.rs.select_resource(resource_name='consensus')
# Append AB: to the receptor names
resource['receptor'] = 'AB:' + resource['receptor']

# Append AB: to the protein modality
mdata.mod['prot'].var_names = 'AB:' + mdata.mod['prot'].var['gene_ids']

li.mt.rank_aggregate(adata=mdata,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource,
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     # Ligand-Receptor pairs are directed so we need to correctly pass
                     # `RNA` with ligands as `x_mod` and receptors as `y_mod`
                     'x_mod': 'rna',
                     'y_mod': 'prot',
                     # We use .X from the x_mod
                     'x_use_raw':False,
                     # We use .X from the y_mod
                     'y_use_raw':False,
                     # NOTE: we need to ensure that the modalities are correctly transformed
                     'x_transform':li.ut.zi_minmax,
                     'y_transform':li.ut.zi_minmax,
                    },
                  verbose=True
                  )

mdata.uns['liana_res'].head()

li.pl.dotplot(adata = mdata,
              colour='lr_means',
              size='specificity_rank',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['CD4+ naïve T', 'NK', 'Treg', 'CD8+ memory T'],
              target_labels=['CD14 mono', 'mature B', 'CD8+ memory T', 'CD16 mono'],
              figure_size=(9, 5),
              # finally, since cpdbv2 suggests using a filter to FPs
              # we filter the pvals column to <= 0.05
              filter_fun=lambda x: x['cellphone_pvals'] <= 0.05,
              cmap='plasma'
             )


# Create the object to use for the matbolite-liana analysis
adata = mdata.mod['rna']
