import networkx as nx
import matplotlib.pyplot as plt
import EoN
import pandas as pd
import numpy as np

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

# read timeseries from the file
timeseriesFilePrefAtt = 'BMC/timeseries_prefAtt_10000Nodes.tsv'

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
plt.plot(thresholded_timeseries_prefAtt['time'], infecteds, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the preferential attachment network of 10000 nodes')

plt.show()

# PREFERENTIAL ATTACHMENT WITH AGING

# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFilePrefAttAging = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/preferentialAttachmentAging/10000nodes/1/edge_data.tsv'
G_prefAttAging = nx.Graph()
with open(graphEdgesFile, 'r') as f:
    # skip the header
    f.readline()
    for line in f:
        start, end, weight = line.strip().split('\t')
        G_prefAttAging.add_edge(start, end, weight=float(weight))



# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFilePrefAttAging = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/preferentialAttachmentAging/10000nodes/1/node_conditions.tsv'
initialConditionsPrefAttAging = {}
with open(initialConditionFilePrefAttAging, 'r') as f:
    for line in f:
        node, condition = line.strip().split('\t')
        initialConditionsPrefAttAging[node] = condition

# get the initial infecteds
initial_infecteds_prefAttAging = []
for node in G_prefAttAging.nodes():
    if initialConditionsPrefAttAging[node] == 'Infectious':
        initial_infecteds_prefAttAging.append(node)

# run the SIS simulation in the graphs
t_sis, S_sis, I_sis = EoN.fast_SIS(G_prefAttAging, tau, gamma, initial_infecteds=initial_infecteds_prefAttAging, tmax=100)

# run the SIR simulation in the graph
t_sir, S_sir, I_sir, R_sir = EoN.fast_SIR(G_prefAttAging, tau, gamma, initial_infecteds=initial_infecteds_prefAttAging, tmax=100)


# read timeseries from the file
timeseriesFilePrefAttAging = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC/preferentialAttachmentAging/10000Nodes/1/allFiles/fullGraph_output.tsv'

timeSeries_dataframe_prefAttAging = pd.read_csv(timeseriesFilePrefAttAging, sep='\t')
thresholded_timeseries_prefAttAging = timeSeries_dataframe_prefAttAging.copy()

thresholded_timeseries_prefAttAging.iloc[:, 1:] = np.where(thresholded_timeseries_prefAttAging.iloc[:, 1:] > 0.5, 1, 0)
# sort the rows by the time column
thresholded_timeseries_prefAttAging = thresholded_timeseries_prefAttAging.sort_values(by='time')

# get the number of infecteds at each time
infecteds = thresholded_timeseries_prefAttAging.iloc[:, 1:].sum(axis=1)

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t_sis, I_sis, label='EoN SIS simulation')
plt.plot(t_sir, I_sir, label='EoN SIR simulation')
plt.plot(thresholded_timeseries_prefAttAging['time'], infecteds, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the preferential attachment network with aging of 10000 nodes')

plt.show()

# ERDOS-RENYI

# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFileErdosRenyi = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/erdosRenyi/10000nodes/1/edge_data.tsv'
G_erdosRenyi = nx.Graph()
with open(graphEdgesFileErdosRenyi, 'r') as f:
    # skip the header
    f.readline()
    for line in f:
        start, end, weight = line.strip().split('\t')
        G_erdosRenyi.add_edge(start, end, weight=float(weight))

# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFileErdosRenyi = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/erdosRenyi/10000nodes/1/node_conditions.tsv'
initialConditionsErdosRenyi = {}
with open(initialConditionFileErdosRenyi, 'r') as f:
    for line in f:
        node, condition = line.strip().split('\t')
        initialConditionsErdosRenyi[node] = condition

# get the initial infecteds
initial_infecteds_erdosRenyi = []
for node in G_erdosRenyi.nodes():
    if initialConditionsErdosRenyi[node] == 'Infectious':
        initial_infecteds_erdosRenyi.append(node)

# run the SIS simulation in the graphs
t_sis, S_sis, I_sis = EoN.fast_SIS(G_erdosRenyi, tau, gamma, initial_infecteds=initial_infecteds_erdosRenyi, tmax=100)

# run the SIR simulation in the graph
t_sir, S_sir, I_sir, R_sir = EoN.fast_SIR(G_erdosRenyi, tau, gamma, initial_infecteds=initial_infecteds_erdosRenyi, tmax=100)

# read timeseries from the file
timeseriesFileErdosRenyi = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC/erdosRenyi/10000Nodes/1/allFiles/fullGraph_output.tsv'

timeSeries_dataframe_erdosRenyi = pd.read_csv(timeseriesFileErdosRenyi, sep='\t')
thresholded_timeseries_erdosRenyi = timeSeries_dataframe_erdosRenyi.copy()

thresholded_timeseries_erdosRenyi.iloc[:, 1:] = np.where(thresholded_timeseries_erdosRenyi.iloc[:, 1:] > 0.5, 1, 0)
# sort the rows by the time column
thresholded_timeseries_erdosRenyi = thresholded_timeseries_erdosRenyi.sort_values(by='time')

# get the number of infecteds at each time
infecteds = thresholded_timeseries_erdosRenyi.iloc[:, 1:].sum(axis=1)

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t_sis, I_sis, label='EoN SIS simulation')
plt.plot(t_sir, I_sir, label='EoN SIR simulation')
plt.plot(thresholded_timeseries_erdosRenyi['time'], infecteds, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the Erdos-Renyi network of 10000 nodes')

plt.show()

# REGULAR

# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFileRegular = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/regular/10000nodes/edge_data.tsv'

G_regular = nx.Graph()
with open(graphEdgesFileRegular, 'r') as f:
    # skip the header
    f.readline()
    for line in f:
        start, end, weight = line.strip().split('\t')
        G_regular.add_edge(start, end, weight=float(weight))

# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFileRegular = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/regular/10000nodes/node_conditions.tsv'
initialConditionsRegular = {}
with open(initialConditionFileRegular, 'r') as f:
    for line in f:
        node, condition = line.strip().split('\t')
        initialConditionsRegular[node] = condition

# get the initial infecteds
initial_infecteds_regular = []
for node in G_regular.nodes():
    if initialConditionsRegular[node] == 'Infectious':
        initial_infecteds_regular.append(node)

# run the SIS simulation in the graphs
t_sis, S_sis, I_sis = EoN.fast_SIS(G_regular, tau, gamma, initial_infecteds=initial_infecteds_regular, tmax=100)

# run the SIR simulation in the graph
t_sir, S_sir, I_sir, R_sir = EoN.fast_SIR(G_regular, tau, gamma, initial_infecteds=initial_infecteds_regular, tmax=100)

# read timeseries from the file
timeseriesFileRegular = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC/regular/10000Nodes/allFiles/fullGraph_output.tsv'

timeSeries_dataframe_regular = pd.read_csv(timeseriesFileRegular, sep='\t')
thresholded_timeseries_regular = timeSeries_dataframe_regular.copy()

thresholded_timeseries_regular.iloc[:, 1:] = np.where(thresholded_timeseries_regular.iloc[:, 1:] > 0.5, 1, 0)

# sort the rows by the time column
thresholded_timeseries_regular = thresholded_timeseries_regular.sort_values(by='time')

# get the number of infecteds at each time
infecteds = thresholded_timeseries_regular.iloc[:, 1:].sum(axis=1)

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t_sis, I_sis, label='EoN SIS simulation')
plt.plot(t_sir, I_sir, label='EoN SIR simulation')
plt.plot(thresholded_timeseries_regular['time'], infecteds, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the regular network of 10000 nodes')

plt.show()