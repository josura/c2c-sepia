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
data = pd.read_csv("/home/josura/Projects/ccc/c2c-sepia/scripts/bash/cluster/performanceAnalysis/localPerformanceTimes.tsv", sep="\t")

# create a new dataframe with the average time, group by the number of nodes and the number of processes
data = data.groupby(["graphNodesNumber", "numberProcesses"]).agg({"time": "mean"}).reset_index()

# plot in 3d the number of processes, the number of nodes and the time
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data["graphNodesNumber"], data["numberProcesses"], data["averageTime"])
ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('Time')
plt.show()

# add column with the speedup per number of nodes and number of processes
# speedup = avg_time(1 processor) / avg_time(n processors)
data["speedup"] = data.groupby(["graphNodesNumber", "numberProcesses"])["averageTime"].transform("min") / data["averageTime"]

# plot the speedup
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data["graphNodesNumber"], data["numberProcesses"], data["speedup"])
ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('Speedup')
plt.show()