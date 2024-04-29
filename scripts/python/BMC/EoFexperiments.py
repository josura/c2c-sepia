import networkx as nx
import matplotlib.pyplot as plt
import EoN
import pandas as pd

# get the network from one of the generated files that contains the edge data (source, target, weight)
graphEdgesFile = 'scripts/R/epidemics/significanceAnalysys/preferentialAttachment/10000nodes/1/edge_data.tsv'
G = nx.read_edgelist(graphEdgesFile, delimiter='\t', data=(('weight', float),))

# set the parameters for the simulation
gamma = 0.2
tau = 0.5
# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFile = 'scripts/R/epidemics/significanceAnalysys/preferentialAttachment/10000nodes/1/node_condititions.tsv'
initialConditions = {}
with open(initialConditionFile, 'r') as f:
    for line in f:
        node, condition = line.strip().split('\t')
        initialConditions[node] = condition

# get the initial infecteds
initial_infecteds = []
for node in G.nodes():
    if initialConditions[node] == 'Infected':
        initial_infecteds.append(node)


# run the SIS simulation in the graphs
t, S, I = EoN.fast_SIS(G, tau, gamma, initial_infecteds=initial_infecteds, tmax=100)

# read timeseries from the file
timeseriesFile = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC/preferentialAttachment/10000nodes/1/allFiles/fullGraph_output.tsv'

timeSeries_dataframe = pd.read_csv(timeseriesFile, sep='\t')
thresholded_timeseries = timeSeries_dataframe.copy()

thresholded_timeseries.iloc[:, 1:] = np.where(thresholded_timeseries.iloc[:, 1:] > 0.5, 1, 0)
# sort the rows by the time column
thresholded_timeseries = thresholded_timeseries.sort_values(by='time')

# get the number of infecteds at each time
infecteds = thresholded_timeseries.iloc[:, 1:].sum(axis=1)

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t, I, label='EoN simulation')
plt.plot(thresholded_timeseries['time'], infecteds, label='Thresholded timeseries')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()

plt.show()