import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import plotly.express as px

# get the outputs(iteration matrices) from the folder
iterationMatrixFolder = "/home/josura/Projects/ccc/datiIdo/inputGraphs/1h/outputWithLR/iterationMatrices"
iterationMatrices = {}


# get the files in the folder, every single file is a type(the name is specified by the first part of the filename), and the content is the iteration matrix
for file in os.listdir(iterationMatrixFolder):
    if file.endswith(".tsv"):
        temp_iterationMatrix = pd.read_csv(iterationMatrixFolder + '/' + file, sep='\t')
        ## drop the last column (useless)
        temp_iterationMatrix = temp_iterationMatrix.drop(temp_iterationMatrix.columns[-1], axis=1)
        ## first row of the matrix is the node names (nodeNames), the rest of the columns are the timepoints (the names go from 0 to the end of the simulation timepoint)
        nodeNames = temp_iterationMatrix['nodeNames']
        ## drop the nodeNames column
        temp_iterationMatrix = temp_iterationMatrix.drop('nodeNames', axis=1)
        ## get the timepoints
        timepoints = temp_iterationMatrix.columns
        ## we need to transpose the matrix so that the columns are the node names and the rows are the timepoints iteration results
        temp_iterationMatrix = temp_iterationMatrix.transpose()
        ## set the column names to be the node names
        temp_iterationMatrix.columns = nodeNames
        ## set the index to be the timepoints
        temp_iterationMatrix.index = timepoints
        ## add a column to be the timepoints
        # temp_iterationMatrix['time'] = timepoints
        ## change index name to be 'time'
        temp_iterationMatrix.index.name = 'time'
        ## add the iteration matrix to the dictionary
        iterationMatrices[file.split('.')[0]] = temp_iterationMatrix

# get the namemap for the nodes, translating the metabolite names to the original names
namemap = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_full_mouse_fixed.csv", sep=',')
## drop rows with NaN values on HMDB column
namemap = namemap.dropna(subset=['HMDB'])


# get the iteration matrix for the AT1-metabolites
AT1_metabolites_iterationMatrix = iterationMatrices['AT1-metabolites']

# read the metabolites data with the different time-points
metabolites_1hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_1h_metabolites_module168_cell1646_20241014-125146.csv" 
metabolites_6hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_6h_metabolites_module168_cell1037_20241014-125629.csv"
metaoblites_7hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_7h_metabolites_module168_cell1722_20241014-123944.csv"
metabolites_10hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_10h_metabolites_module168_cell1240_20241014-171445.csv"

# read the information of the transcriptomics about the cells at different timepoints
rna_1h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-1h-metadata.tsv"
rna_6h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-6h-metadata.tsv"
rna_7h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-7h-metadata.tsv"
rna_10h_metadataFile = "/home/josura/Projects/ccc/datiIdo/lianaInputs/rna-10h-metadata.tsv"

#1h
rna_1h_metadata_pd = pd.read_csv(rna_1h_metadataFile, sep="\t", index_col=0)
## read well data from the seurat obj
well_celltype_df_1h = pd.DataFrame()
well_celltype_df_1h["celltype"] = rna_1h_metadata_pd["cell_type"]
well_celltype_df_1h["well"] = rna_1h_metadata_pd.index
## reading the metabolites and join them with the well information
metabolites_1h = pd.read_csv(metabolites_1hFile, sep=",", index_col=0)
metabolites_1h["well"] = metabolites_1h.index
metabolites_1h = metabolites_1h.merge(well_celltype_df_1h, on="well", how="left")
#6h
rna_6h_metadata_pd = pd.read_csv(rna_6h_metadataFile, sep="\t", index_col=0)
## read well data from the seurat obj
well_celltype_df_6h = pd.DataFrame()
well_celltype_df_6h["celltype"] = rna_6h_metadata_pd["cell_type"]
well_celltype_df_6h["well"] = rna_6h_metadata_pd.index
## reading the metabolites and join them with the well information
metabolites_6h = pd.read_csv(metabolites_6hFile, sep=",", index_col=0)
metabolites_6h["well"] = metabolites_6h.index
metabolites_6h = metabolites_6h.merge(well_celltype_df_6h, on="well", how="left")
#7h
rna_7h_metadata_pd = pd.read_csv(rna_7h_metadataFile, sep="\t", index_col=0)
## read well data from the seurat obj
well_celltype_df_7h = pd.DataFrame()
well_celltype_df_7h["celltype"] = rna_7h_metadata_pd["cell_type"]
well_celltype_df_7h["well"] = rna_7h_metadata_pd.index
## reading the metabolites and join them with the well information
metabolites_7h = pd.read_csv(metabolites_7hFile, sep=",", index_col=0)
metabolites_7h["well"] = metabolites_7h.index
metabolites_7h = metabolites_7h.merge(well_celltype_df_7h, on="well", how="left")
#10h
rna_10h_metadata_pd = pd.read_csv(rna_10h_metadataFile, sep="\t", index_col=0)
## read well data from the seurat obj
well_celltype_df_10h = pd.DataFrame()
well_celltype_df_10h["celltype"] = rna_10h_metadata_pd["cell_type"]
well_celltype_df_10h["well"] = rna_10h_metadata_pd.index
## reading the metabolites and join them with the well information
metabolites_10h = pd.read_csv(metabolites_10hFile, sep=",", index_col=0)
metabolites_10h["well"] = metabolites_10h.index
metabolites_10h = metabolites_10h.merge(well_celltype_df_10h, on="well", how="left")



# example plotting for the AT1-metabolites iteration matrix
# plot the iteration matrix
## plot it in different subplots of 3 rows and 3 columns to show all the nodes
current_completed_plots = 0
while current_completed_plots < len(AT1_metabolites_iterationMatrix.columns):
    fig, axs = plt.subplots(3, 3)
    fig.suptitle('AT1-metabolites')
    for i in range(3):
        for j in range(3):
            if current_completed_plots < len(AT1_metabolites_iterationMatrix.columns):
                axs[i, j].plot(AT1_metabolites_iterationMatrix.index, AT1_metabolites_iterationMatrix[AT1_metabolites_iterationMatrix.columns[current_completed_plots]])
                HMDB_id = AT1_metabolites_iterationMatrix.columns[current_completed_plots]
                ## get the original name of the metabolite
                original_name_Series = namemap[namemap['HMDB'] == HMDB_id]['Match']
                if len(original_name_Series) > 0:
                    original_name = original_name_Series.values[0]
                else:
                    original_name = HMDB_id
                axs[i, j].set_title(original_name)
                # change the orientation of the x labels to be slighly rotated
                axs[i, j].tick_params(axis='x', rotation=90)
                current_completed_plots += 1
            else:
                break
    plt.show()

# example plotting for the AT1-metabolites iteration matrix with plotly
# plot the iteration matrix
## plot it in different subplots of 3 rows and 3 columns to show all the nodes
current_completed_plots = 0
while current_completed_plots < len(AT1_metabolites_iterationMatrix.columns):
    fig = px.line()
    fig.update_layout(title='AT1-metabolites')
    for i in range(3):
        for j in range(3):
            if current_completed_plots < len(AT1_metabolites_iterationMatrix.columns):
                HMDB_id = AT1_metabolites_iterationMatrix.columns[current_completed_plots]
                ## get the original name of the metabolite
                original_name_Series = namemap[namemap['HMDB'] == HMDB_id]['Match']
                if len(original_name_Series) > 0:
                    original_name = original_name_Series.values[0]
                else:
                    original_name = HMDB_id
                fig.add_scatter(x=AT1_metabolites_iterationMatrix.index, y=AT1_metabolites_iterationMatrix[AT1_metabolites_iterationMatrix.columns[current_completed_plots]], name=original_name)
                current_completed_plots += 1
            else:
                break
    fig.show()