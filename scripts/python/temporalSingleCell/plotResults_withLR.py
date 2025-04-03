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

# read the metabolites data with the different time-points
metabolites_1hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_1h_metabolites_module168_cell1646_20241014-125146.csv" 
metabolites_6hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_6h_metabolites_module168_cell1037_20241014-125629.csv"
metabolites_7hFile = "/home/josura/Projects/ccc/fluxes/scFEA/output/scRNA_7h_metabolites_module168_cell1722_20241014-123944.csv"
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
metabolites_1h.index = metabolites_1h["well"]
## drop the well column
metabolites_1h = metabolites_1h.drop("well", axis=1)
## averaging the metabolites for the same cell type
metabolites_1h_averaged = metabolites_1h.groupby(["celltype"]).mean()
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
metabolites_6h.index = metabolites_6h["well"]
## drop the well column
metabolites_6h = metabolites_6h.drop("well", axis=1)
## averaging the metabolites for the same cell type
metabolites_6h_averaged = metabolites_6h.groupby(["celltype"]).mean()
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
metabolites_7h.index = metabolites_7h["well"]
## drop the well column
metabolites_7h = metabolites_7h.drop("well", axis=1)
## averaging the metabolites for the same cell type
metabolites_7h_averaged = metabolites_7h.groupby(["celltype"]).mean()
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
metabolites_10h.index = metabolites_10h["well"]
## drop the well column
metabolites_10h = metabolites_10h.drop("well", axis=1)
## averaging the metabolites for the same cell type
metabolites_10h_averaged = metabolites_10h.groupby(["celltype"]).mean()

# reading the names for the metabolites as ids from this file, since the full names don't match the original names sometimes
metabolites_names = pd.read_csv("/home/josura/Projects/ccc/fluxes/scFEA/data/cName_complete_mouse_c70_m168.csv", sep=",")
## there is only one row in the dataframe with the CNAMES for every entry in the row, create a map from the metabolite id to the name
### first let's change the column names to be the metabolite ids in the average metabolites dataframes
### let's control if they have the same column names
if len(metabolites_1h_averaged.columns) == len(metabolites_names.columns):
    columns_match = True
    for i in range(len(metabolites_1h_averaged.columns)):
        if metabolites_1h_averaged.columns[i] != metabolites_names.columns[i]:
            print("The columns don't match: " + metabolites_1h_averaged.columns[i] + " != " + metabolites_names.columns[i])
            columns_match = False
            break
    if columns_match:
        metabolites_1h_averaged.columns = metabolites_names.values[0]
        metabolites_6h_averaged.columns = metabolites_names.values[0]
        metabolites_7h_averaged.columns = metabolites_names.values[0]
        metabolites_10h_averaged.columns = metabolites_names.values[0]
else:
    print("The columns don't match on the size") 
### reading the name_map
namemap = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/name_map_full_mouse_fixed.csv", sep=',')
### change the column names in the metabolites dataframes to the real original names (the ids are the column called Query, and the original names are in Match)
original_name_list = []
for i in range(len(metabolites_1h_averaged.columns)):
    metaboliteToChange = metabolites_1h_averaged.columns[i]
    real_name_pd = namemap[namemap["Query"] == metaboliteToChange]["Match"]
    real_name = metaboliteToChange
    if len(real_name_pd) > 0:
        real_name = real_name_pd.values[0]
    else:
        print("The metabolite " + metaboliteToChange + " doesn't have a real name")
    original_name_list.append(real_name)

metabolites_1h_averaged.columns = original_name_list
metabolites_6h_averaged.columns = original_name_list
metabolites_7h_averaged.columns = original_name_list
metabolites_10h_averaged.columns = original_name_list

# example plotting for the AT1-metabolites iteration matrix
## get the iteration matrix for the AT1-metabolites
AT1_metabolites_iterationMatrix = iterationMatrices['AT1-metabolites'].copy()
## adding the 1h real values to the simulated data, since it's the first timepoint
### shifting all the timepoints by the minimum interval between them (that is the intra timestep)
timestep = AT1_metabolites_iterationMatrix.index.astype(float)[1] - AT1_metabolites_iterationMatrix.index.astype(float)[0]
AT1_metabolites_iterationMatrix.index = AT1_metabolites_iterationMatrix.index.astype(float) + timestep
### adding the 1h real values to the simulated data
pd_row_to_add = pd.DataFrame(index=[0], columns=AT1_metabolites_iterationMatrix.columns)
metabolites_for_AT1_1h = metabolites_1h_averaged[metabolites_1h_averaged.index == 'AT1']
if len(metabolites_for_AT1_1h) > 0:
    for i in range(len(AT1_metabolites_iterationMatrix.columns)):
        metaboliteHMDBName = AT1_metabolites_iterationMatrix.columns[i]
        metaboliteOriginalName = namemap[namemap['HMDB'] == metaboliteHMDBName]['Match']
        if len(metaboliteOriginalName) > 0:
            metaboliteName = metaboliteOriginalName.values[0]
        else:
            metaboliteName = metaboliteHMDBName
        ## check if the metabolite is in the metabolites_1h_averaged dataframe
        ## if it is, add the value to the pd_row_to_add dataframe
        ## if it is not, add 0
        if metaboliteName in metabolites_for_AT1_1h.columns:
            pd_row_to_add[metaboliteHMDBName] = metabolites_for_AT1_1h[metaboliteName].values[0]
        else:
            pd_row_to_add[metaboliteHMDBName] = 0
    pd_row_to_add.index = [0]
    pd_row_to_add.index.name = 'time'
    AT1_metabolites_iterationMatrix = pd.concat([pd_row_to_add, AT1_metabolites_iterationMatrix], axis=0)
# plot the iteration matrix
## plot it in different subplots of 3 rows and 3 columns to show all the nodes
# current_completed_plots = 0
# while current_completed_plots < len(AT1_metabolites_iterationMatrix.columns):
#     fig, axs = plt.subplots(3, 3)
#     fig.suptitle('AT1-metabolites')
#     for i in range(3):
#         for j in range(3):
#             if current_completed_plots < len(AT1_metabolites_iterationMatrix.columns):
#                 axs[i, j].plot(AT1_metabolites_iterationMatrix.index, AT1_metabolites_iterationMatrix[AT1_metabolites_iterationMatrix.columns[current_completed_plots]])
#                 HMDB_id = AT1_metabolites_iterationMatrix.columns[current_completed_plots]
#                 ## get the original name of the metabolite
#                 original_name_Series = namemap[namemap['HMDB'] == HMDB_id]['Match']
#                 if len(original_name_Series) > 0:
#                     original_name = original_name_Series.values[0]
#                 else:
#                     original_name = HMDB_id
#                 axs[i, j].set_title(original_name)
#                 # change the orientation of the x labels to be slighly rotated
#                 axs[i, j].tick_params(axis='x', rotation=90)
#                 current_completed_plots += 1
#             else:
#                 break
#     plt.show()

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
                ## plot the real metabolite value as a single point in the plot, specifically for AT1-metabolites
                valueForMetabolite_1h = 0
                valueForMetabolite_6h = 0
                valueForMetabolite_7h = 0
                valueForMetabolite_10h = 0
                metabolites_for_AT1_1h = metabolites_1h_averaged[metabolites_1h_averaged.index == 'AT1']
                if original_name in metabolites_1h_averaged.columns:
                    valueForMetabolite_1h = metabolites_1h_averaged[original_name].values[0]
                metabolites_for_AT1_6h = metabolites_6h_averaged[metabolites_6h_averaged.index == 'AT1']
                if original_name in metabolites_6h_averaged.columns:
                    valueForMetabolite_6h = metabolites_6h_averaged[original_name].values[0]
                metabolites_for_AT1_7h = metabolites_7h_averaged[metabolites_7h_averaged.index == 'AT1']
                if original_name in metabolites_7h_averaged.columns:
                    valueForMetabolite_7h = metabolites_7h_averaged[original_name].values[0]
                metabolites_for_AT1_10h = metabolites_10h_averaged[metabolites_10h_averaged.index == 'AT1']
                if original_name in metabolites_10h_averaged.columns:
                    valueForMetabolite_10h = metabolites_10h_averaged[original_name].values[0]
                if not( "v-in" in original_name or "v-out" in original_name ):
                    fig.add_scatter(x=[0, 5, 6, 9], y=[valueForMetabolite_1h, valueForMetabolite_6h, valueForMetabolite_7h, valueForMetabolite_10h], mode='markers', name='Real for ' + original_name)
                
                current_completed_plots += 1
            else:
                break
    fig.show()