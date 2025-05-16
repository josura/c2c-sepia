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