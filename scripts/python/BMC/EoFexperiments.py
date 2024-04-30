import networkx as nx
import matplotlib.pyplot as plt
import EoN
import pandas as pd
import numpy as np

# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFile = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/preferentialAttachment/10000nodes/1/edge_data.tsv'
G = nx.Graph()
with open(graphEdgesFile, 'r') as f:
    # skip the header
    f.readline()
    for line in f:
        start, end, weight = line.strip().split('\t')
        G.add_edge(start, end, weight=float(weight))

# set the parameters for the simulation
gamma = 0.2
tau = 0.5
# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFile = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/preferentialAttachment/10000nodes/1/node_conditions.tsv'
initialConditions = {}
with open(initialConditionFile, 'r') as f:
    for line in f:
        node, condition = line.strip().split('\t')
        initialConditions[node] = condition

# get the initial infecteds
initial_infecteds = []
for node in G.nodes():
    if initialConditions[node] == 'Infectious':
        initial_infecteds.append(node)


# run the SIS simulation in the graphs
t_sis, S_sis, I_sis = EoN.fast_SIS(G, tau, gamma, initial_infecteds=initial_infecteds, tmax=100)

# run the SIR simulation in the graph
t_sir, S_sir, I_sir, R_sir = EoN.fast_SIR(G, tau, gamma, initial_infecteds=initial_infecteds, tmax=100)

# read timeseries from the file
timeseriesFile = 'BMC/timeseries_prefAtt_10000Nodes.tsv'

timeSeries_dataframe = pd.read_csv(timeseriesFile, sep='\t')
thresholded_timeseries = timeSeries_dataframe.copy()

thresholded_timeseries.iloc[:, 1:] = np.where(thresholded_timeseries.iloc[:, 1:] > 0.5, 1, 0)
# sort the rows by the time column
thresholded_timeseries = thresholded_timeseries.sort_values(by='time')

# get the number of infecteds at each time
infecteds = thresholded_timeseries.iloc[:, 1:].sum(axis=1)

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t_sis, I_sis, label='EoN SIS simulation')
plt.plot(t_sir, I_sir, label='EoN SIR simulation')
plt.plot(thresholded_timeseries['time'], infecteds, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()

plt.show()