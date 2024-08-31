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
x = data["graphNodesNumber"].unique()
y = data["numberProcesses"].unique()
X, Y = np.meshgrid(x, y)
Z = np.array(data["time"]).reshape(len(y), len(x))
# Plot the 3D surface
ax.plot_surface(X, Y, Z, edgecolor='royalblue', lw=0.5, rstride=8, cstride=8,
                alpha=0.3)

# Plot projections of the contours for each dimension.  By choosing offsets
# that match the appropriate axes limits, the projected contours will sit on
# the 'walls' of the graph
ax.contourf(X, Y, Z, zdir='z', offset=-100, cmap='coolwarm')
ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap='coolwarm')
ax.contourf(X, Y, Z, zdir='y', offset=40, cmap='coolwarm')

ax.set(xlim=(-40, 40), ylim=(-40, 40), zlim=(-100, 100),
       xlabel='X', ylabel='Y', zlabel='Z')

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