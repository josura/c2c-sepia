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
        

