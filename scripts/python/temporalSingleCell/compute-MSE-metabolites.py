import pandas as pd
import numpy as np
import os
import pyenv

# read the environment variables
env = pyenv.readEnvFile("/home/josura/Projects/ccc/c2c-sepia/scripts/bash/cluster/performanceAnalysis/localPerformanceTimes.tsv")
# get the path from the environment variables
type = env["TYPE"]
outputPath = 