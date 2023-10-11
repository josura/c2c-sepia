import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np

# Sample data for nodes and edges
nodes_data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes/fullGraph_output.tsv', sep='\t')
edges_data = pd.read_csv('/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/syntheticGraphs/100Nodes/edge_data.tsv', sep='\t')

# Initialize the graph
G = nx.Graph()

# Add edges to the graph
for _, row in edges_data.iterrows():
    G.add_edge(row['Start'], row['End'], weight=row['Weight'])

# Function to update node colors based on their values
def update(iteration):
    plt.cla()  # Clear previous plot
    # search for the row where iteration is equal to the current iteration
    # rowIteration = nodes_data[nodes_data['iteration'] == iteration]
    node_colors = [nodes_data[node][nodes_data['iteration'] == iteration] for node in nodes_data.columns[1:]]
    nx.draw(G, with_labels=True, node_color=node_colors, cmap=plt.cm.RdYlBu, node_size=500, font_size=12)
    plt.title(f'Iteration {iteration}')
    plt.show()

# Function to update node colors based on their values, maintains node positions
def updateTrue(iteration):
    plt.cla()  # Clear previous plot
    # search for the row where iteration is equal to the current iteration
    # rowIteration = nodes_data[nodes_data['iteration'] == iteration]
    node_colors = [nodes_data[node][nodes_data['iteration'] == iteration] for node in nodes_data.columns[1:]]
    # interpolate values from 0 to 1 to node colors, from white to blue
    nx.draw(G, with_labels=True, node_color=node_colors, cmap=plt.cm.cool , node_size=500, font_size=12, pos=nx.spring_layout(G, seed=42))
    plt.title(f'Iteration {iteration}')
    plt.show()


# Generate animation
num_iterations = len(nodes_data)
ani = FuncAnimation(plt.gcf(), updateTrue, frames=num_iterations, repeat=False)
plt.show()
