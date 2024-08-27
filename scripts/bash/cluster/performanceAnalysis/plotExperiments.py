import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
#import seaborn as sns
import numpy as np
import os

# read the timeseries data, the data is in the format of:
# graphNodesNumber	numberProcesses	numberTypes	numberIterations	time

# read the data
data = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/bash/cluster/performanceAnalysis/performanceTimes.tsv", sep="\t")

# add a column with the average time, group by the number of nodes and the number of processes
data["averageTime"] = data.groupby(["graphNodesNumber", "numberProcesses"])["time"].transform("mean")

# plot in 3d the number of processes, the number of nodes and the time
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data["graphNodesNumber"], data["numberProcesses"], data["averageTime"])
ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('Time')
plt.show()