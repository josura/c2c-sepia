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

# create another dataframe with thresholded values, all the columns after the first one should be thresholded since they are the time series
thresholded_data = nodes_data.copy()
thresholded_data.iloc[:, 1:] = np.where(thresholded_data.iloc[:, 1:] > 0.5, 1, 0)