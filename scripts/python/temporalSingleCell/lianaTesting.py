# testing the code at https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html
import numpy as np
import liana as li
import mudata as mu
import scanpy as sc

from matplotlib import pyplot as plt
from plotly import express as px

# Load the data
# prot = sc.read('citeseq_prot.h5ad', backup_url='https://figshare.com/ndownloader/files/47625196')
# rna = sc.read('citeseq_rna.h5ad', backup_url='https://figshare.com/ndownloader/files/47625193')
prot = sc.read('/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/citeseq_prot.h5ad')
rna = sc.read('/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/citeseq_rna.h5ad')

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

plot = li.pl.dotplot(adata = mdata,
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
plot.show()


# Create the object to use for the matbolite-liana analysis
## Focus on Transcriptomics Data
adata = mdata.mod['rna']


## convert to murine symbols from human to mouse, vignette available at: https://liana-py.readthedocs.io/en/latest/notebooks/sma.html
map_df = li.rs.get_hcop_orthologs(columns=['human_symbol', 'mouse_symbol'],
                                  min_evidence=3
                                  ).rename(columns={'human_symbol':'source',
                                                   'mouse_symbol':'target'})

## Obtain MetalinksDB Prior Knowledge
metalinks = li.resource.get_metalinks(biospecimen_location='Blood',
                                      source=['CellPhoneDB', 'Cellinker', 'scConnect', # Ligand-Receptor resources
                                              'recon', 'hmr', 'rhea', 'hmdb' # Production-Degradation resources
                                              ],
                                      types=['pd', 'lr'], # NOTE: we obtain both ligand-receptor and production-degradation sets
                                     )

metalinks_translated = li.rs.translate_column(resource=metalinks,
                                   map_df=map_df,
                                   column='gene_symbol',
                                   one_to_many=1)
metalinks_translated.head()



## preparing the metabolite-receptor pairs
resource = metalinks_translated[metalinks_translated['type']=='lr'].copy()
resource = resource[['metabolite', 'gene_symbol']]\
    .rename(columns={'gene_symbol':'receptor'}).drop_duplicates()
resource.head()

## Prepare the Production-Degradation Network
pd_net = metalinks_translated[metalinks_translated['type'] == 'pd']
# we need to aggregate the production-degradation values
pd_net = pd_net[['metabolite', 'gene_symbol', 'mor']].groupby(['metabolite', 'gene_symbol']).agg('mean').reset_index()
pd_net.head()


## Prepare the Transport Network
t_net = metalinks_translated[metalinks_translated['type'] == 'pd']
t_net = t_net[['metabolite', 'gene_symbol', 'transport_direction']].dropna()
# Note that we treat export as positive and import as negative
t_net['mor'] = t_net['transport_direction'].apply(lambda x: 1 if x == 'out' else -1 if x == 'in' else None)
t_net = t_net[['metabolite', 'gene_symbol', 'mor']].dropna().groupby(['metabolite', 'gene_symbol']).agg('mean').reset_index()
t_net = t_net[t_net['mor']!=0]

