import networkx as nx
import matplotlib.pyplot as plt
import EoN
import pandas as pd
import numpy as np
from collections import defaultdict

# set the parameters for the simulation
gamma = 0.2
tau = 0.5

# PREFERENTIAL ATTACHMENT

# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFile = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/preferentialAttachment/10000nodes/1/edge_data.tsv'
G_prefAtt = nx.Graph()
with open(graphEdgesFile, 'r') as f:
    # skip the header
    f.readline()
    for line in f:
        start, end, weight = line.strip().split('\t')
        G_prefAtt.add_edge(start, end, weight=float(weight))

# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFilePrefAtt = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/preferentialAttachment/10000nodes/1/node_conditions.tsv'
initialConditionsPrefAtt = {}
with open(initialConditionFilePrefAtt, 'r') as f:
    for line in f:
        node, condition = line.strip().split('\t')
        initialConditionsPrefAtt[node] = condition

# get the initial infecteds
initial_infecteds_prefAtt = []
for node in G_prefAtt.nodes():
    if initialConditionsPrefAtt[node] == 'Infectious':
        initial_infecteds_prefAtt.append(node)


# run the SIS simulation in the graphs
t_sis, S_sis, I_sis = EoN.fast_SIS(G_prefAtt, tau, gamma, initial_infecteds=initial_infecteds_prefAtt, tmax=100)

# run the SIR simulation in the graph
t_sir, S_sir, I_sir, R_sir = EoN.fast_SIR(G_prefAtt, tau, gamma, initial_infecteds=initial_infecteds_prefAtt, tmax=100)

# SIRS simulation setup

H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
H.add_edge('I', 'R', rate = 0.2)   # recovery rate
H.add_edge('R', 'S', rate = 0.1)   # immunity loss rate

J = nx.DiGraph()    #DiGraph showing transition that does require an interaction.
J.add_edge(('I', 'S'), ('I', 'I'), rate = 0.5)  #IS->II

IC_prefAtt = defaultdict(lambda: 'S')
for node in initial_infecteds_prefAtt:
    IC_prefAtt[node] = 'I'

return_statuses = ('S', 'I', 'R')

t_sirs, S_sirs, I_sirs, R_sirs = EoN.Gillespie_simple_contagion(G_prefAtt, H, J, IC_prefAtt, return_statuses, tmax = 100)


# read timeseries from the file
timeseriesFilePrefAtt = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC/preferentialAttachment/10000Nodes/1/allFiles/fullGraph_output.tsv'

timeSeries_dataframe_prefAtt = pd.read_csv(timeseriesFilePrefAtt, sep='\t')
thresholded_timeseries_prefAtt = timeSeries_dataframe_prefAtt.copy()

thresholded_timeseries_prefAtt.iloc[:, 1:] = np.where(thresholded_timeseries_prefAtt.iloc[:, 1:] > 0.5, 1, 0)
# sort the rows by the time column
thresholded_timeseries_prefAtt = thresholded_timeseries_prefAtt.sort_values(by='time')

# get the number of infecteds at each time
infecteds = thresholded_timeseries_prefAtt.iloc[:, 1:].sum(axis=1)

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t_sis, I_sis, label='EoN SIS simulation')
plt.plot(t_sir, I_sir, label='EoN SIR simulation')
plt.plot(t_sirs, I_sirs, label='EoN SIRS simulation')
plt.plot(thresholded_timeseries_prefAtt['time'], infecteds, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the preferential attachment network of 10000 nodes')

plt.show()
