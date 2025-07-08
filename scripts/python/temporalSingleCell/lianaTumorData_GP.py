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
# rna_matriceIPF_filepath = "/home/josura/Projects/ccc/datiGrete/matrice_IPF_agr.txt"
rna_matriceIPF_filepath = "/home/josura/Projects/ccc/datiGrete/matrice_IPF_specifica_agr.txt"
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
## create the AnnData object
rna_Normal = sc.AnnData(rna_Normal_pd)
rna_Normal.obs["celltype"] = rna_cell_names_Normal
rna_Normal.obs["condition"] = "Normal"
## create the umap embedding
sc.pp.neighbors(rna_Normal, n_neighbors=10)
sc.tl.umap(rna_Normal)
## show the data
X = list(map(lambda x: x[0], rna_Normal.obsm["X_umap"]))
Y = list(map(lambda x: x[1], rna_Normal.obsm["X_umap"]))
celltypes_Normal = list(rna_Normal.obs["celltype"].values)

fig = px.scatter(x=X, y=Y, color=celltypes_Normal)
fig.show()

# TODO combine the two AnnData into one single AnnData object, although maybe it is better to estimate the metabolites with the two separately and then combine them

# TODO create the mdata object with the rna and the inferred metabolites


#---- METABOLITES ESTIMATION ----#

# Obtain a ligand-receptor resource of interest
## Obtain MetalinksDB Prior Knowledge
metalinks = li.resource.get_metalinks(biospecimen_location='Blood',
                                      source=['CellPhoneDB', 'Cellinker', 'scConnect', # Ligand-Receptor resources
                                              'recon', 'hmr', 'rhea', 'hmdb' # Production-Degradation resources
                                              ],
                                      types=['pd', 'lr'], # NOTE: we obtain both ligand-receptor and production-degradation sets
                                     )

# Filter the resource to keep only ligand-receptor interactions
resource = metalinks[metalinks['type']=='lr'].copy()
resource = resource[['metabolite', 'gene_symbol']]\
    .rename(columns={'gene_symbol':'receptor'}).drop_duplicates()
resource.head()

# Preparing the production-degradation resource
prod_degr_net = metalinks[metalinks['type'] == 'pd']
# we need to aggregate the production-degradation values
prod_degr_net = prod_degr_net[['metabolite', 'gene_symbol', 'mor']].groupby(['metabolite', 'gene_symbol']).agg('mean').reset_index()
prod_degr_net.head()

# Preparing the transporter resource
transporter_net = metalinks[metalinks['type'] == 'pd']
transporter_net = transporter_net[['metabolite', 'gene_symbol', 'transport_direction']].dropna()
# Note that we treat export as positive and import as negative
transporter_net['mor'] = transporter_net['transport_direction'].apply(lambda x: 1 if x == 'out' else -1 if x == 'in' else None)
transporter_net = transporter_net[['metabolite', 'gene_symbol', 'mor']].dropna().groupby(['metabolite', 'gene_symbol']).agg('mean').reset_index()
transporter_net = transporter_net[transporter_net['mor']!=0]                           


## IPF
# Estimating the metabolites
meta_IPF = li.mt.fun.estimate_metalinks(rna_IPF,
                                    resource,
                                    pd_net=prod_degr_net,
                                    t_net=transporter_net, # (Optional)
                                    use_raw=False,
                                    # keyword arguments passed to decoupler-py
                                    source='metabolite', target='gene_symbol',
                                    weight='mor', min_n=3)
# pass cell type information
meta_IPF.obs['celltype'] = rna_IPF.obs['celltype']

#dataset with two modalities, one for RNA and one for Metabolites. The metabolites are estimated as t-values.
# visualization 
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": (100)}):
    sc.pl.umap(meta_IPF.mod['metabolite'], color=['Prostaglandin J2', 'Metanephrine', 'celltype'], cmap='coolwarm')

# Save the inferred metabolites in a table
meta_IPF_tobesave = meta_IPF.mod['metabolite'].to_df()
meta_IPF_tobesave['celltype'] = meta_IPF.obs['celltype']
meta_IPF_tobesave.to_csv('inferred_metabolites_IPF.csv', sep='\t', index=False)

# infer the metabolites interactions
li.mt.rank_aggregate(adata=meta_IPF,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource.rename(columns={'metabolite':'ligand'}),
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     'x_mod': 'metabolite',
                     'y_mod': 'receptor',
                     'x_use_raw':False,
                     'y_use_raw':False,
                     'x_transform':li.ut.zi_minmax,
                     'y_transform':li.ut.zi_minmax,
                    },
                  verbose=True
                  )

#visualization of the inferred metabolites interactions

meta_IPF.uns['liana_res'].head()

# Save the inferred interactions in a table
meta_IPF.uns['liana_res'].to_csv('inferred_interactions_IPF.csv', sep='\t', index=False)

# sources and target for IPF are:
# sources: 'Epithelial_cells', 'none', 'T_cells', 'Pre-B_cell_CD34-', 'NK_cell', 'Tissue_stem_cells', 'B_cell', 'CMP', 'Chondrocytes', 'Endothelial_cells', 'Fibroblasts', 'Smooth_muscle_cells', 'Monocyte', 'DC', 'Macrophage', 'GMP'
# targets: ['GMP', 'DC', 'CMP', 'Tissue_stem_cells', 'Pre-B_cell_CD34-', 'Monocyte', 'Chondrocytes', 'NK_cell', 'T_cells', 'Smooth_muscle_cells', 'Endothelial_cells', 'Fibroblasts', 'none', 'Macrophage', 'B_cell', 'Epithelial_cells']
interactionPlot = li.pl.dotplot(adata = meta_IPF,
                colour='lr_means',
                size='cellphone_pvals',
                inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                source_labels=list(pd.unique(meta_IPF.uns['liana_res']['source'])),
                target_labels=list(pd.unique(meta_IPF.uns['liana_res']['target'])),
                figure_size=(24, 12),
                # Filter to top 10 acc to magnitude rank
                top_n=10,
                orderby='magnitude_rank',
                orderby_ascending=True,
                cmap='plasma'
                 )

interactionPlot.show()


## Normal
# Estimating the metabolites
meta_Normal = li.mt.fun.estimate_metalinks(rna_Normal,
                                    resource,
                                    pd_net=prod_degr_net,
                                    t_net=transporter_net, # (Optional)
                                    use_raw=False,
                                    # keyword arguments passed to decoupler-py
                                    source='metabolite', target='gene_symbol',
                                    weight='mor', min_n=3)

# pass cell type information
meta_Normal.obs['celltype'] = rna_Normal.obs['celltype']

#dataset with two modalities, one for RNA and one for Metabolites. The metabolites are estimated as t-values.
# visualization
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": (100)}):
    sc.pl.umap(meta_Normal.mod['metabolite'], color=['Prostaglandin J2', 'Metanephrine', 'celltype'], cmap='coolwarm')

# Save the inferred metabolites in a table
meta_Normal_tobesave = meta_Normal.mod['metabolite'].to_df()
meta_Normal_tobesave['celltype'] = meta_Normal.obs['celltype']
meta_Normal_tobesave.to_csv('inferred_metabolites_Normal.csv', sep='\t')

# infer the metabolites interactions
li.mt.rank_aggregate(adata=meta_Normal,
                     groupby='celltype',
                     # pass our modified resource
                     resource=resource.rename(columns={'metabolite':'ligand'}),
                     # NOTE: Essential arguments when handling multimodal data
                     mdata_kwargs={
                     'x_mod': 'metabolite',
                     'y_mod': 'receptor',
                     'x_use_raw':False,
                     'y_use_raw':False,
                     'x_transform':li.ut.zi_minmax,
                     'y_transform':li.ut.zi_minmax,
                    },
                  verbose=True
                  )

#visualization of the inferred metabolites interactions
meta_Normal.uns['liana_res'].head()

# Save the inferred interactions in a table
meta_Normal.uns['liana_res'].to_csv('inferred_interactions_Normal.csv', sep='\t',index=False)

# sources and target for Normal are:
# sources:
# targets:
interactionPlot = li.pl.dotplot(adata = meta_Normal,
                colour='lr_means',
                size='cellphone_pvals',
                inverse_size=True, # we inverse sign since we want small p-values to have large sizes
                source_labels=list(pd.unique(meta_Normal.uns['liana_res']['source'])),
                target_labels=list(pd.unique(meta_Normal.uns['liana_res']['target'])),
                figure_size=(24, 12),
                # Filter to top 10 acc to magnitude rank
                top_n=10,
                orderby='magnitude_rank',
                orderby_ascending=True,
                cmap='plasma'
                 )

interactionPlot.show()
