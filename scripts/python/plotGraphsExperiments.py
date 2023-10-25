import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np

# get different time series from different folder in the specified folder

# get the folders in the specified folder
folders = [f for f in os.listdir('/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes') if os.path.isdir(os.path.join('/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes', f))]

edges_data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/syntheticGraphs/100Nodes/edge_data.tsv', sep='\t')

nodes_data = []
# iterate over the folders and get the time series
for folder in folders:
    nodes_data.append(pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes/' + folder + '/allfiles/fullGraph_output.tsv', sep='\t'))

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
    thresholded_data_sorted[i] = thresholded_data_sorted[i].sort_values(by=['iteration'])

    


# create a list of lists with the percentage of infected at each time step for different folders, the index of the list is the time step (iteration column in the dataframe)
infected_nodes = []
for i in range(len(thresholded_data_sorted)):
    infected_nodes.append([])
    for _, row in thresholded_data_sorted[i].iterrows():
        infected_nodes[i].append(len(row[row == 1])/len(row))


for _, row in thresholded_data_sorted.iterrows():
    infected_nodes.append(row[row == 1].index.tolist())



# plot the percentages through time, on the x axis we have the time steps and on the y axis we have the percentage of infected nodes, each line represents a different folder
fig, ax = plt.subplots()
ax.set_xlim(0, len(infected_nodes[0]))
ax.set_ylim(0, 1)
ax.set_xlabel('Time')
ax.set_ylabel('Percentage of Infected Nodes')
ax.set_title('Percentage of Infected Nodes')
for i in range(len(infected_nodes)):
    folder = folders[i]
    folderNoPath = folder.split('/')[-1]
    ax.plot(infected_nodes[i], label=folderNoPath)
ax.legend()
plt.show()

# plot the number of susceptible nodes
ax.plot(susceptible_nodes_stats, label='Susceptible Nodes')
ax.legend()
plt.show()

