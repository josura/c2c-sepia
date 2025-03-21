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

# read the fluxes data with the different time-points

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