import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np
import os

# get different time series from different folder in the specified folder

# get the folders in the specified folder
specified_folder = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes'
folders = [f for f in os.listdir(specified_folder) if os.path.isdir(os.path.join(specified_folder, f))]

edges_data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/syntheticGraphs/100Nodes/edge_data.tsv', sep='\t')

nodes_data = []
# iterate over the folders and get the time series
for folder in folders:
    nodes_data.append(pd.read_csv(specified_folder + '/' + folder + '/allFiles/fullGraph_output.tsv', sep='\t'))

#cast Start and End columns to Strings
edges_data['Start'] = edges_data['Start'].astype(str)
edges_data['End'] = edges_data['End'].astype(str)

# Initialize the graph
G = nx.Graph()

# Add edges to the graph
for _, row in edges_data.iterrows():
    G.add_edge(row['Start'], row['End'], weight=row['Weight'])


# create another list of dataframes with thresholded values, all the columns after the first one should be thresholded since they are the time series, 0 values means that the node is not infected and 1 means that it is infected
thresholded_data_sorted = nodes_data.copy()
for i in range(len(thresholded_data_sorted)):
    thresholded_data_sorted[i].iloc[:, 1:] = np.where(thresholded_data_sorted[i].iloc[:, 1:] > 0.5, 1, 0)
    # order the rows by the iteration column
    thresholded_data_sorted[i] = thresholded_data_sorted[i].sort_values(by=['time'])

    


# create a list of lists with the percentage of infected at each time step for different folders, the index of the list is the time step (iteration column in the dataframe)
infected_nodes = []
for i in range(len(thresholded_data_sorted)):
    infected_nodes.append([])
    for _, row in thresholded_data_sorted[i].iterrows():
        infected_nodes[i].append(len(row[row == 1])/len(row))



# plot the percentages through time, on the x axis we have the time steps and on the y axis we have the percentage of infected nodes, each line represents a different folder
fig, ax = plt.subplots()
ax.set_xlim(0, len(infected_nodes[0]))
ax.set_ylim(0, 1)
ax.set_xlabel('Time')
ax.set_ylabel('Percentage of Infected Nodes')
ax.set_title('Percentage of Infected Nodes, Different Parameters (Dissipation scale factor, Propagation scale factor)')
for i in range(len(infected_nodes)):
    folder = folders[i]
    folderNoPath = folder.split('/')[-1]
    # the file has the following format: dissipationScaleFactor<float>_propagationScaleFactor<float>
    dissipationScaleFactor = folderNoPath.split('_')[0]
    dissipationScaleFactor = dissipationScaleFactor.split('dissipationScaleFactor')[1]
    propagationScaleFactor = folderNoPath.split('_')[1]
    propagationScaleFactor = propagationScaleFactor.split('propagationScaleFactor')[1]
    label = '( ' + dissipationScaleFactor + ' , ' + propagationScaleFactor + ' )'
    ax.plot(infected_nodes[i], label=label)
# augment dimensions of the legend box
ax.legend(title="(\u03BB (n), \u03C9 (n))", fontsize=12)
plt.show()

