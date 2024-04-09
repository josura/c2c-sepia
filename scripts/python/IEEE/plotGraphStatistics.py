import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np

# Sample data for nodes and edges
nodes_data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes/fullGraph_output.tsv', sep='\t')
edges_data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/syntheticGraphs/100Nodes/edge_data.tsv', sep='\t')

#cast Start and End columns to Strings
edges_data['Start'] = edges_data['Start'].astype(str)
edges_data['End'] = edges_data['End'].astype(str)

# Initialize the graph
G = nx.Graph()

# Add edges to the graph
for _, row in edges_data.iterrows():
    G.add_edge(row['Start'], row['End'], weight=row['Weight'])

# create another dataframe with thresholded values, all the columns after the first one should be thresholded since they are the time series, 0 values means that the node is not infected and 1 means that it is infected
thresholded_data_sorted = nodes_data.copy()
thresholded_data_sorted.iloc[:, 1:] = np.where(thresholded_data_sorted.iloc[:, 1:] > 0.5, 1, 0)

# order the rows by the iteration column
thresholded_data_sorted = thresholded_data_sorted.sort_values(by=['iteration'])

# create a list of lists with the nodes that are infected at each time step, the index of the list is the time step (iteration column in the dataframe)
infected_nodes = []
for _, row in thresholded_data_sorted.iterrows():
    infected_nodes.append(row[row == 1].index.tolist())

# create a list of lists with the nodes that are not infected at each time step
susceptible_nodes = []
for _, row in thresholded_data_sorted.iterrows():
    susceptible_nodes.append(row[row == 0].index.tolist())


# get the statistics(number of infected vs number of susceptible) of the states of the nodes at each time step
infected_nodes_stats = []
susceptible_nodes_stats = []
for i in range(len(infected_nodes)):
    infected_nodes_stats.append(len(infected_nodes[i]))
    susceptible_nodes_stats.append(len(susceptible_nodes[i]))



# plot the statistics, on the x axis we have the time steps and on the y axis we have the number of infected and susceptible nodes
fig, ax = plt.subplots()
ax.set_xlim(0, len(infected_nodes))
ax.set_ylim(0, len(G.nodes))
ax.set_xlabel('Time')
ax.set_ylabel('Number of Nodes')
ax.set_title('Number of Infected vs Number of Susceptible Nodes')
ax.plot(infected_nodes_stats, label='Infected Nodes')

# plot the number of susceptible nodes
ax.plot(susceptible_nodes_stats, label='Susceptible Nodes')
ax.legend()
plt.show()

