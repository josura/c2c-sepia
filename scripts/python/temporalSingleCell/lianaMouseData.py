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
rna_1hFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-1h.tsv"
rna_1h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-1h-metadata.tsv"
metabolitesFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_1h_metabolites_module168_cell1646_20241014-125146.csv" 

rna_pd = pd.read_csv(rna_1hFile, sep="\t", index_col=0)
rna_metadata_pd = pd.read_csv(rna_1h_metadataFile, sep="\t", index_col=0)
metabolites_pd = pd.read_csv(metabolitesFile, sep=",", index_col=0)
#drop first row in metabolites_pd
metabolites_pd = metabolites_pd.iloc[1:]

rna = sc.AnnData(rna_pd)
rna.obs = rna_metadata_pd
metabolites = sc.AnnData(metabolites_pd)

mdata = mu.MuData({'rna': rna, 'metabolites': metabolites})
# make sure that cell type is accessible
mdata.obs['celltype'] = mdata.mod['rna'].obs['cell_type'].astype('category')
# inspect the object
mdata

# add UMAP coordinates to the RNA data
sc.pp.neighbors(rna, n_neighbors=10)
sc.tl.umap(rna)

# show the data
X = list(map(lambda x: x[0], rna.obsm["X_umap"]))
Y = list(map(lambda x: x[1], rna.obsm["X_umap"]))
celltypes = list(rna.obs["cell_type"].values)

fig = px.scatter(x=X, y=Y, color=celltypes)
fig.show()

# Obtain a ligand-receptor resource of interest
resource = li.rs.select_resource(resource_name='consensus')

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

resource_translated = li.rs.translate_column(resource=resource,
                                      map_df=map_df,
                                      column='ligand',
                                      one_to_many=1)

resource_translated = li.rs.translate_column(resource=resource_translated,
                                      map_df=map_df,
                                      column='receptor',
                                      one_to_many=1)


# ligand receptors analysis
li.mt.rank_aggregate(adata=mdata,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource_translated,
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     # Ligand-Receptor pairs are directed so we need to correctly pass
                     # `RNA` with ligands as `x_mod` and receptors as `y_mod`
                     'x_mod': 'rna',
                     'y_mod': 'rna',
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

# metabolite receptors analysis TODO

## control the intersection between pd.unique(metalinks_translated["metabolite"]) and metabolites.to_df().columns
list(set(metalinks_translated["metabolite"]) & set(metabolites.to_df().columns))
## the list is very small, ['Ornithine', 'Glycine', 'Glutathione', 'Deoxyadenosine', 'Hypoxanthine', 'Cholesterol', 'Choline']
## we need to use other IDs to match the metabolites

name_map_70 = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_70_mouse.csv", sep=",")
name_map_metalinks = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_metalinks.csv", sep=",")

# filter the rows that have NaN value in the HMDB column
name_map_metalinks = name_map_metalinks[name_map_metalinks['HMDB'].notna()]
name_map_70 = name_map_70[name_map_70['HMDB'].notna()]

# filter the column that are not in the name_map_70 for the metabolites data (match column in the name_map_70)
## last column is not in the name_map_70( controlled by observing the data)
metabolites_filtered_df = metabolites.to_df()[metabolites.to_df().columns[:-1]]
## changing the column names to the HMDB ids
metabolites_filtered_df.columns = name_map_70['HMDB'].values

metabolites_filtered = sc.AnnData(metabolites_filtered_df)

mdata = mu.MuData({'rna': rna, 'metabolites': metabolites_filtered})
# make sure that cell type is accessible
mdata.obs['celltype'] = mdata.mod['rna'].obs['cell_type'].astype('category')
# inspect the object
mdata

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

#adata translated, changing row names with the mouse symbols
adata_translated = adata.copy()
adata_translated.var.index = 



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

meta = li.mt.fun.estimate_metalinks(adata,
                                    resource,
                                    pd_net=pd_net,
                                    t_net=t_net, # (Optional)
                                    use_raw=False,
                                    # keyword arguments passed to decoupler-py
                                    source='metabolite', target='gene_symbol',
                                    weight='mor', min_n=3)
# pass cell type information
meta.obs['celltype'] = adata.obs['celltype']

with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": (100)}):
    sc.pl.umap(meta.mod['metabolite'], color=['Prostaglandin J2', 'Metanephrine', 'celltype'], cmap='coolwarm')


# to get the dataframe representing the metabolite/protein data
mdata.mod["prot"].to_df()