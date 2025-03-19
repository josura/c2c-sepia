import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

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

