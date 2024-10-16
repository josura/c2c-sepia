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
ax.scatter(data["graphNodesNumber"], data["numberProcesses"], data["time"])
ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('average Time')
plt.show()

# also plot as wireframe
ax = plt.figure().add_subplot(projection='3d')
# Create data
nodes = data["graphNodesNumber"].unique()
processes = data["numberProcesses"].unique()
X, Y = np.meshgrid(nodes, processes)
Z = np.zeros_like(X, dtype=np.float64)

for i, node in enumerate(nodes):
    for j, process in enumerate(processes):
        Z[j, i] = data[(data['graphNodesNumber'] == node) & (data['numberProcesses'] == process)]['time'].values[0]

# Plot the wireframe
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, Z)

ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('Average Time')
plt.show()

# add column with the speedup per number of nodes and number of processes
# speedup(nProcessors,Nnodes) = avg_time(1processor,Nnodes) / avg_time(nProcessors, Nnodes)
speedup = []
for i in range(len(data)):
    timeForOneProcessor = data[(data["numberProcesses"] == 1) & (data["graphNodesNumber"] == data["graphNodesNumber"][i])]["time"].values[0]
    speedup.append(timeForOneProcessor / data["time"][i])

data["speedup"] = speedup


# plot the speedup
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data["graphNodesNumber"], data["numberProcesses"], data["speedup"])
ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('Speedup')
plt.show()

# also plot as wireframe
ax = plt.figure().add_subplot(projection='3d')
# Create data
nodes = data["graphNodesNumber"].unique()
processes = data["numberProcesses"].unique()
X, Y = np.meshgrid(nodes, processes)
Z = np.zeros_like(X, dtype=np.float64)

for i, node in enumerate(nodes):
    for j, process in enumerate(processes):
        Z[j, i] = data[(data['graphNodesNumber'] == node) & (data['numberProcesses'] == process)]['speedup'].values[0]

# Plot the wireframe
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(X, Y, Z)

ax.set_xlabel('Number of Nodes')
ax.set_ylabel('Number of Processes')
ax.set_zlabel('Speedup')

plt.show()

# save aggregated data into a file
