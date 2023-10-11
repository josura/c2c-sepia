import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd

# Sample data for nodes and edges
nodes_data = pd.read_csv('nodes_data.tsv', sep='\t')
edges_data = pd.read_csv('edges_data.tsv', sep='\t')

# Initialize the graph
G = nx.Graph()

# Add edges to the graph
for _, row in edges_data.iterrows():
    G.add_edge(row['Start'], row['End'], weight=row['Weight'])

# Function to update node colors based on their values
def update(iteration):
    plt.cla()  # Clear previous plot
    node_colors = [nodes_data[node][iteration] for node in nodes_data.columns[1:]]
    nx.draw(G, with_labels=True, node_color=node_colors, cmap=plt.cm.RdYlBu, node_size=5000, font_size=12)
    plt.title(f'Iteration {iteration}')
    plt.show()

# Generate animation
num_iterations = len(nodes_data)
ani = FuncAnimation(plt.gcf(), update, frames=num_iterations, repeat=False)
plt.show()
