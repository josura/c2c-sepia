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
rna_6hFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-6h.tsv"
rna_7hFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-7h.tsv"
rna_10hFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-10h.tsv"

rna_1h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-1h-metadata.tsv"
rna_6h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-6h-metadata.tsv"
rna_7h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-7h-metadata.tsv"
rna_10h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-10h-metadata.tsv"

metabolites_1hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_1h_metabolites_module168_cell1646_20241014-125146.csv" 
metabolites_6hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_6h_metabolites_module168_cell1037_20241014-125629.csv"
metaoblites_7hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_7h_metabolites_module168_cell1722_20241014-123944.csv"
metabolites_10hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_10h_metabolites_module168_cell1240_20241014-171445.csv"

#1h
rna_1h_pd = pd.read_csv(rna_1hFile, sep="\t", index_col=0)
rna_1h_metadata_pd = pd.read_csv(rna_1h_metadataFile, sep="\t", index_col=0)
metabolites_1h_pd = pd.read_csv(metabolites_1hFile, sep=",", index_col=0)

rna_1h = sc.AnnData(rna_1h_pd)
rna_1h.obs = rna_1h_metadata_pd
metabolites_1h = sc.AnnData(metabolites_1h_pd)

mdata_1h = mu.MuData({'rna': rna_1h, 'metabolites': metabolites_1h})
# make sure that cell type is accessible
mdata_1h.obs['celltype'] = mdata_1h.mod['rna'].obs['cell_type'].astype('category')
# inspect the object
mdata_1h

# add UMAP coordinates to the RNA data
sc.pp.neighbors(rna_1h, n_neighbors=10)
sc.tl.umap(rna_1h)

# show the data
X = list(map(lambda x: x[0], rna_1h.obsm["X_umap"]))
Y = list(map(lambda x: x[1], rna_1h.obsm["X_umap"]))
celltypes = list(rna_1h.obs["cell_type"].values)

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
li.mt.rank_aggregate(adata=mdata_1h,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource_translated,
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_1h_kwargs={
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

mdata_1h.uns['liana_res'].head()

# metabolite receptors analysis TODO

## control the intersection between pd.unique(metalinks_translated["metabolite"]) and metabolites_1h.to_df().columns
list(set(metalinks_translated["metabolite"]) & set(metabolites_1h.to_df().columns))
## the list is very small, ['Ornithine', 'Glycine', 'Glutathione', 'Deoxyadenosine', 'Hypoxanthine', 'Cholesterol', 'Choline']
## we need to use other IDs to match the metabolites_1h

name_map_70 = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_70_mouse_fixed.csv", sep=",")  #using the fixed name map for changing the glucose
name_map_metalinks = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_metalinks.csv", sep=",")

# filter the rows that have NaN value in the HMDB column
name_map_metalinks = name_map_metalinks[name_map_metalinks['HMDB'].notna()]
name_map_70 = name_map_70[name_map_70['HMDB'].notna()]

# filter the column that are not in the name_map_70 for the metabolites_1h data (match column in the name_map_70)
## last column is not in the name_map_70( controlled by observing the data)
metabolites_1h_filtered_df = metabolites_1h.to_df()[metabolites_1h.to_df().columns[:-1]]
## changing the column names to the HMDB ids
metabolites_1h_filtered_df.columns = name_map_70['HMDB'].values

# control the intersection between pd.unique(metalinks_translated["HMDB"]) and metabolites_filtered_df.columns
list(set(metalinks_translated["hmdb"]) & set(metabolites_1h_filtered_df.columns))

metabolites_1h_filtered = sc.AnnData(metabolites_1h_filtered_df)

mdata_1h = mu.MuData({'rna': rna_1h, 'metabolites': metabolites_1h_filtered})
# make sure that cell type is accessible
mdata_1h.obs['celltype'] = mdata_1h.mod['rna'].obs['cell_type'].astype('category')
# inspect the object
mdata_1h

# Create the object to use for the matbolite-liana analysis
## Focus on Transcriptomics Data
adata = mdata_1h.mod['rna']


## preparing the metabolite-receptor pairs
resource = metalinks_translated[metalinks_translated['type']=='lr'].copy()
resource = resource[['hmdb', 'gene_symbol']]\
    .rename(columns={'gene_symbol':'receptor'}).drop_duplicates()
resource.head()

## Prepare the Production-Degradation Network
pd_net = metalinks_translated[metalinks_translated['type'] == 'pd']
# we need to aggregate the production-degradation values
pd_net = pd_net[['hmdb', 'gene_symbol', 'mor']].groupby(['hmdb', 'gene_symbol']).agg('mean').reset_index()
pd_net.head()


## Prepare the Transport Network
t_net = metalinks_translated[metalinks_translated['type'] == 'pd']
t_net = t_net[['hmdb', 'gene_symbol', 'transport_direction']].dropna()
# Note that we treat export as positive and import as negative
t_net['mor'] = t_net['transport_direction'].apply(lambda x: 1 if x == 'out' else -1 if x == 'in' else None)
t_net = t_net[['hmdb', 'gene_symbol', 'mor']].dropna().groupby(['hmdb', 'gene_symbol']).agg('mean').reset_index()
t_net = t_net[t_net['mor']!=0]

meta = li.mt.fun.estimate_metalinks(adata,
                                    resource,
                                    pd_net=pd_net,
                                    t_net=t_net, # (Optional)
                                    use_raw=False,
                                    # keyword arguments passed to decoupler-py
                                    x_name='hmdb',
                                    source='hmdb', target='gene_symbol',
                                    weight='mor', min_n=3)
# pass cell type information
meta.obs['celltype'] = adata.obs['cell_type']


with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": (100)}):
    sc.pl.umap(meta.mod['hmdb'], color=['HMDB0000122', 'HMDB0000123', 'cell_type'], cmap='coolwarm')
# plotting D-glucose and glycine for the inferred metabolites 


mdata_1h.mod["metabolites"].obs["cell_type"] = adata.obs['cell_type']
# add UMAP coordinates to the metabolites data
mdata_1h.mod["metabolites"].obsm["X_umap"] = mdata_1h.mod["rna"].obsm["X_umap"]

with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": (100)}):
    sc.pl.umap(mdata_1h.mod['metabolites'], color=['HMDB0000122', 'HMDB0000123', 'cell_type'], cmap='coolwarm')
# plotting D-glucose and glycine for the original metabolites

# add scFEA metabolites to meta 
meta.mod["scFEA"] = mdata_1h.mod["metabolites"]

# filter the resources to only include the metabolites that are in the scFEA data (column names)
resource_filtered = resource[resource['hmdb'].isin(mdata_1h.mod["metabolites"].to_df().columns)]


#Infer Metabolite-Receptor Interactions
#We will next infer the putative ligand-receptor interactions between these two modalities.
li.mt.rank_aggregate(adata=meta,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource.rename(columns={'hmdb':'ligand'}),
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     'x_mod': 'scFEA',
                     'y_mod': 'receptor',
                     'x_use_raw':False,
                     'y_use_raw':False,
                     'x_transform':li.ut.zi_minmax,
                     'y_transform':li.ut.zi_minmax,
                    },
                  verbose=True
                  )
# li.mt.rank_aggregate(adata=mdata_1h,
#                      groupby='celltype',
#                      # pass our modified resource
#                      resource=resource.rename(columns={'hmdb':'ligand'}),
#                      # NOTE: Essential arguments when handling multimodal data
#                      mdata_kwargs={
#                      'x_mod': 'metabolites',
#                      'y_mod': 'rna',
#                      'x_use_raw':False,
#                      'y_use_raw':False,
#                      'x_transform':li.ut.zi_minmax,
#                      'y_transform':li.ut.zi_minmax,
#                     },
#                   verbose=True
#                   )

meta.uns['liana_res'].head()
# celltypes available for source are 'DC', 'AT1', 'T', 'Neut', 'MacIII', 'Mon', 'MacII', 'Fibro', 'B', 'Endothel', 'Baso', 'NK', 'Smooth', 'Clara', 'Matrix', 'AT2'
# celltypes available for target are 'Clara', 'NK', 'Smooth', 'AT2', 'Matrix', 'Endothel', 'DC', 'Fibro', 'Baso', 'AT1'
interactionPlot = li.pl.dotplot(adata = meta,
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['MacII', 'Neut', 'Baso', 'MacIII'],
              target_labels=['AT1' ,'AT2', 'DC', 'NK'],
              figure_size=(12, 6),
              # Filter to top 10 acc to magnitude rank
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              cmap='plasma'
             )

interactionPlot.show()

# save the interactions with the following format
# startType	startNodeName	endType	endNodeName	weight contactTimes
results_1h = meta.uns['liana_res']
## change the name of columns source target ligand_complex receptor_complex and scaled_weight
results_1h = results_1h.rename(columns={'source':'startType', 'target':'endType', 'ligand_complex':'startNodeName', 'receptor_complex':'endNodeName', 'scaled_weight':'weight'})
results_1h['contactTimes'] = 1
results_1h.to_csv("/home/josura/Projects/ccc/datiIdo/inputGraphs/1h/interactions/results_metabolite_1h.tsv", sep="\t", index=False)

# 6h
rna_6h_pd = pd.read_csv(rna_6hFile, sep="\t", index_col=0)
rna_6h_metadata_pd = pd.read_csv(rna_6h_metadataFile, sep="\t", index_col=0)
metabolites_6h_pd = pd.read_csv(metabolites_6hFile, sep=",", index_col=0)

## last column is not in the name_map_70( controlled by observing the data)
metabolites_6h_filtered_df = metabolites_6h_pd[metabolites_6h_pd.columns[:-1]]
## changing the column names to the HMDB ids
metabolites_6h_filtered_df.columns = name_map_70['HMDB'].values

rna_6h = sc.AnnData(rna_6h_pd)
rna_6h.obs = rna_6h_metadata_pd
metabolites_6h = sc.AnnData(metabolites_6h_filtered_df)

# add UMAP coordinates to the RNA data
sc.pp.neighbors(rna_6h, n_neighbors=10)
sc.tl.umap(rna_6h)

# show the data
# X = list(map(lambda x: x[0], rna_6h.obsm["X_umap"]))
# Y = list(map(lambda x: x[1], rna_6h.obsm["X_umap"]))
# celltypes = list(rna_6h.obs["cell_type"].values)
# fig = px.scatter(x=X, y=Y, color=celltypes)
# fig.show()

# No need to obtain the ligand-receptor resource of interest, since it is already obtained in the previous step

# No need to obtain MetalinksDB Prior Knowledge, since it is already obtained in the previous step

# No need to translate the resource, since it is already translated in the previous step

mdata_6h = mu.MuData({'rna': rna_6h, 'metabolites': metabolites_6h})
# make sure that cell type is accessible
mdata_6h.obs['celltype'] = mdata_6h.mod['rna'].obs['cell_type'].astype('category')
# inspect the object
mdata_6h


# ligand receptors analysis
li.mt.rank_aggregate(adata=mdata_6h,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource.rename(columns={'hmdb':'ligand'}),
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     'x_mod': 'metabolites',
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

mdata_6h.uns['liana_res'].head()

interactionPlot = li.pl.dotplot(adata = meta,
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['MacII', 'Neut', 'Baso', 'MacIII'],
              target_labels=['AT1' ,'AT2', 'DC', 'NK'],
              figure_size=(12, 6),
              # Filter to top 10 acc to magnitude rank
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              cmap='plasma'
             )

interactionPlot.show()

# save the interactions with the following format
# startType	startNodeName	endType	endNodeName	weight contactTimes
results_6h = mdata_6h.uns['liana_res']
## change the name of columns source target ligand_complex receptor_complex and scaled_weight
results_6h = results_6h.rename(columns={'source':'startType', 'target':'endType', 'ligand_complex':'startNodeName', 'receptor_complex':'endNodeName', 'scaled_weight':'weight'})
results_6h['contactTimes'] = 6
results_6h.to_csv("/home/josura/Projects/ccc/datiIdo/inputGraphs/6h/interactions/results_metabolite_6h.tsv", sep="\t", index=False)

# 7h
rna_7h_pd = pd.read_csv(rna_7hFile, sep="\t", index_col=0)
rna_7h_metadata_pd = pd.read_csv(rna_7h_metadataFile, sep="\t", index_col=0)
metabolites_7h_pd = pd.read_csv(metaoblites_7hFile, sep=",", index_col=0)

## last column is not in the name_map_70( controlled by observing the data)
metabolites_7h_filtered_df = metabolites_7h_pd[metabolites_7h_pd.columns[:-1]]
## changing the column names to the HMDB ids
metabolites_7h_filtered_df.columns = name_map_70['HMDB'].values

rna_7h = sc.AnnData(rna_7h_pd)
rna_7h.obs = rna_7h_metadata_pd
metabolites_7h = sc.AnnData(metabolites_7h_filtered_df)

# add UMAP coordinates to the RNA data
sc.pp.neighbors(rna_7h, n_neighbors=10)
sc.tl.umap(rna_7h)

# show the data
# X = list(map(lambda x: x[0], rna_7h.obsm["X_umap"]))
# Y = list(map(lambda x: x[1], rna_7h.obsm["X_umap"]))
# celltypes = list(rna_7h.obs["cell_type"].values)
# fig = px.scatter(x=X, y=Y, color=celltypes)
# fig.show()

# No need to obtain the ligand-receptor resource of interest, since it is already obtained in the previous step

# No need to obtain MetalinksDB Prior Knowledge, since it is already obtained in the previous step

# No need to translate the resource, since it is already translated in the previous step

mdata_7h = mu.MuData({'rna': rna_7h, 'metabolites': metabolites_7h})
# make sure that cell type is accessible
mdata_7h.obs['celltype'] = mdata_7h.mod['rna'].obs['cell_type'].astype('category')
# inspect the object
mdata_7h


# ligand receptors analysis
li.mt.rank_aggregate(adata=mdata_7h,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource.rename(columns={'hmdb':'ligand'}),
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     'x_mod': 'metabolites',
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

mdata_7h.uns['liana_res'].head()

interactionPlot = li.pl.dotplot(adata = meta,
                colour='lr_means',
                size='cellphone_pvals',
                inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                source_labels=['MacII', 'Neut', 'Baso', 'MacIII'],
                target_labels=['AT1' ,'AT2', 'DC', 'NK'],
                figure_size=(12, 6),
                # Filter to top 10 acc to magnitude rank
                top_n=10,
                orderby='magnitude_rank',
                orderby_ascending=True,
                cmap='plasma'
                 )

interactionPlot.show()

# save the interactions with the following format
# startType	startNodeName	endType	endNodeName	weight contactTimes
results_7h = mdata_7h.uns['liana_res']
## change the name of columns the needed names
results_7h = results_7h.rename(columns={'source':'startType', 'target':'endType', 'ligand_complex':'startNodeName', 'receptor_complex':'endNodeName', 'scaled_weight':'weight'})
results_7h['contactTimes'] = 7
results_7h.to_csv("/home/josura/Projects/ccc/datiIdo/inputGraphs/7h/interactions/results_metabolite_7h.tsv", sep="\t", index=False)

# 10h
rna_10h_pd = pd.read_csv(rna_10hFile, sep="\t", index_col=0)
rna_10h_metadata_pd = pd.read_csv(rna_10h_metadataFile, sep="\t", index_col=0)
metabolites_10h_pd = pd.read_csv(metabolites_10hFile, sep=",", index_col=0)

## last column is not in the name_map_70( controlled by observing the data)
metabolites_10h_filtered_df = metabolites_10h_pd[metabolites_10h_pd.columns[:-1]]
## changing the column names to the HMDB ids
metabolites_10h_filtered_df.columns = name_map_70['HMDB'].values

rna_10h = sc.AnnData(rna_10h_pd)
rna_10h.obs = rna_10h_metadata_pd
metabolites_10h = sc.AnnData(metabolites_10h_filtered_df)

# add UMAP coordinates to the RNA data
sc.pp.neighbors(rna_10h, n_neighbors=10)
sc.tl.umap(rna_10h)

# show the data
# X = list(map(lambda x: x[0], rna_10h.obsm["X_umap"]))
# Y = list(map(lambda x: x[1], rna_10h.obsm["X_umap"]))
# celltypes = list(rna_10h.obs["cell_type"].values)
# fig = px.scatter(x=X, y=Y, color=celltypes)


# No need to obtain the ligand-receptor resource of interest, since it is already obtained in the previous step

# No need to obtain MetalinksDB Prior Knowledge, since it is already obtained in the previous step

# No need to translate the resource, since it is already translated in the previous step

mdata_10h = mu.MuData({'rna': rna_10h, 'metabolites': metabolites_10h})

# make sure that cell type is accessible
mdata_10h.obs['celltype'] = mdata_10h.mod['rna'].obs['cell_type'].astype('category')

# inspect the object
mdata_10h

# ligand receptors analysis
li.mt.rank_aggregate(adata=mdata_10h,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource.rename(columns={'hmdb':'ligand'}),
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     'x_mod': 'metabolites',
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

mdata_10h.uns['liana_res'].head()

interactionPlot = li.pl.dotplot(adata = meta,
                colour='lr_means',
                size='cellphone_pvals',
                inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                source_labels=['MacII', 'Neut', 'Baso', 'MacIII'],
                target_labels=['AT1' ,'AT2', 'DC', 'NK'],
                figure_size=(12, 6),
                # Filter to top 10 acc to magnitude rank
                top_n=10,
                orderby='magnitude_rank',
                orderby_ascending=True,
                cmap='plasma'
                 )

interactionPlot.show()

# intra-cellular communication creation