import pandas as pd
import numpy as np
import os
from dotenv import load_dotenv
# read the environment variables
load_dotenv()
# get the path from the environment variables
type = os.environ.get("TYPE")
outputPath_matrices_all_experiments = os.environ.get("ITERATION-MATRICES-PATH") 

# read the different experiments names (the names of the folders inside the outputPath_matrices_all_experiments environment variable)
experiments = os.listdir(outputPath_matrices_all_experiments)

# reading the shared data between the different experiments ( the real data and some other things to make it consistent)
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
source_namemap = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/python/temporalSingleCell/KEGG-to-original-names-all.tsv", sep='\t')


### change the column names in the metabolites dataframes to the HMDB names, there needs to be a translation from (source names) -> KEGG ids -> HMDB
node_name_list = []
for i in range(len(metabolites_1h_averaged.columns)):
    metaboliteToChange = metabolites_1h_averaged.columns[i]
    # kegg_name_pd = source_namemap[source_namemap["name"] == metaboliteToChange]["KEGG"]
    #real_name_pd = namemap[namemap["Query"] == metaboliteToChange]["Match"]
    final_name = metaboliteToChange
    # if len(kegg_name_pd) > 0:
    HMDB_name_pd = namemap[namemap["KEGG"] == metaboliteToChange]["HMDB"].dropna()
    if(len(HMDB_name_pd) > 0):
        final_name = HMDB_name_pd.values[0]
    else:
        print("The metabolite " + metaboliteToChange + " doesn't have a HMDB identifier")
    # else:
    #     print("The metabolite " + metaboliteToChange + " doesn't have a real name")
    node_name_list.append(final_name)

metabolites_1h_averaged.columns = node_name_list
metabolites_6h_averaged.columns = node_name_list
metabolites_7h_averaged.columns = node_name_list
metabolites_10h_averaged.columns = node_name_list

# creating the dataframe that store the MSEs for the experiments connected to a specific timepoint
MSE_6h_df = pd.DataFrame(columns = ["MSE","experiment"])
# start the loop going through the different experiments, and for each experiment, select the type from the environment variable, and the different timepoints, then compute the MSE for each timepoint (sum of the squared differences between the predicted and the real values for each node in the selected graph)
for experiment in experiments:
    print("Computing MSE for experiment: " + experiment)
    ## control if the file for the iteration matrix exists
    iterationMatrixSelected_path = os.path.join(outputPath_matrices_all_experiments, experiment,"iterationMatrices", type + "-metabolites.tsv")
    iterationMatrixSelected = pd.DataFrame()
    if( not os.path.exists(iterationMatrixSelected_path)):
        print("The file " + iterationMatrixSelected_path + " doesn't exist")
    else:
        temp_iterationMatrix = pd.read_csv(iterationMatrixSelected_path, sep='\t')
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
        ## get the iteration matrix 
        iterationMatrixSelected = temp_iterationMatrix
    print("The iteration for "+ experiment+ " matrix has been read")
    ## get the corresponding timepoints for the iteration matrix
    timepoints = iterationMatrixSelected.index
    ## get the timepoints for the metabolites, the 1h timepoint is not set since it's the input of the simulation, the timepoints are shifted by one hour and 0.2h (to account for the missing first time point)
    timepoints_metabolites = [ "4.8", "5.8", "8.8"]
    ## create the variable to store the MSE for each timepoint
    mse_6h = 0
    mse_7h = 0
    mse_10h = 0
    ### computing for the 6h timepoint
    metabolites_6h_forSelectedType = metabolites_6h_averaged[metabolites_6h_averaged.index == type]
    iterationMatrix_6h = iterationMatrixSelected[iterationMatrixSelected.index == timepoints_metabolites[0]]
    metabolites_count = 0
    for i in range(len(metabolites_6h_averaged.columns)):
        if metabolites_6h_averaged.columns[i] in iterationMatrix_6h.columns:
            # the column is present in the simulated data, so we can compute the error for it
            metabolites_count = metabolites_count + 1
            error = metabolites_6h_averaged[metabolites_6h_averaged.columns[i]].values[0] - iterationMatrix_6h[metabolites_6h_averaged.columns[i]].values[0]
            mse_6h = mse_6h + np.square(error)
        else:
            # the column is not present in the simulated data, so no error can be computed
            print("There is no metabolite in the simulation for the (" + metabolites_6h_averaged.columns[i] + ") metabolite")
    mse_6h = mse_6h / metabolites_count
    MSE_6h_df = pd.concat([MSE_6h_df,pd.DataFrame({"MSE": [mse_6h],"experiment":experiment})])
    #print("MSE_6h for the experiment " + experiment + " is: " + str(mse_6h))


print(MSE_6h_df)
        
