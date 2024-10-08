import networkx as nx
import matplotlib.pyplot as plt
import EoN
import pandas as pd
import numpy as np
from collections import defaultdict
import sys

# set the parameters for the simulation
gamma = 0.2
tau = 0.5
scale_factor = 1

# control if the user has provided the arguments of the command line, that is the graph folder and the time series file generated by the MASFENON simulation
if len(sys.argv) < 3:
    print('Usage: python generalEoFparameterFitting.py <graph_folder> <time_series_file> scaling_factor(optional)') # TODO expand with another time series file for the SIRS simulation
    sys.exit(1)

# read the graph folder from the command line
graphFolder = sys.argv[1]


# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFile = graphFolder + '/edge_data.tsv'
# control if the file for the graph edges exists
try:
    with open(graphEdgesFile, 'r') as f:
        pass
except FileNotFoundError:
    print('The file ' + graphEdgesFile + ' does not exist')
    sys.exit(1)

G = nx.Graph()
with open(graphEdgesFile, 'r') as f:
    # skip the header
    f.readline()
    for line in f:
        start, end, weight = line.strip().split('\t')
        G.add_edge(start, end, weight=float(weight))

# read the initial condition of every node from a file
# the file should have the following format:
# node_name  condition(Susceptible, Infected)
initialConditionFile = graphFolder + '/node_conditions.tsv'
# control if the file for the node conditions exists
try:
    with open(initialConditionFile, 'r') as f:
        pass
except FileNotFoundError:
    print('The file ' + initialConditionFile + ' does not exist')
    sys.exit(1)

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

# SIRS simulation setup

H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
H.add_edge('I', 'R', rate = 0.2)   # recovery rate
H.add_edge('R', 'S', rate = 0.1)   # immunity loss rate

J = nx.DiGraph()    #DiGraph showing transition that does require an interaction.
J.add_edge(('I', 'S'), ('I', 'I'), rate = 0.5)  #IS->II

IC = defaultdict(lambda: 'S')
for node in initial_infecteds:
    IC[node] = 'I'

return_statuses = ('S', 'I', 'R')

t_sirs, S_sirs, I_sirs, R_sirs = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = 100)


# read timeseries from the file
timeSeriesParameterSIS = sys.argv[2]
# control if the file for the time series exists
try:
    with open(timeSeriesParameterSIS, 'r') as f:
        pass
except FileNotFoundError:
    print('The file ' + timeSeriesParameterSIS + ' does not exist')
    sys.exit(1)

timeSeries_dataframeSIS = pd.read_csv(timeSeriesParameterSIS, sep='\t')
thresholded_timeseriesSIS = timeSeries_dataframeSIS.copy()

thresholded_timeseriesSIS.iloc[:, 1:] = np.where(thresholded_timeseriesSIS.iloc[:, 1:] > 0.5, 1, 0)
# sort the rows by the time column
thresholded_timeseriesSIS = thresholded_timeseriesSIS.sort_values(by='time')

# get the number of infecteds at each time
infectedsSIS = thresholded_timeseriesSIS.iloc[:, 1:].sum(axis=1)

# add initial infected quantity as the first time point of the MASFENON simulation
infectedsSIS = pd.concat([pd.Series([len(initial_infecteds)]), infectedsSIS], ignore_index=True)
# remove the last time point from the MASFENON simulation
infectedsSIS = infectedsSIS[:-1]

# scale the time points of MASFENON to simulate lesser timesteps
if len(sys.argv) == 4:
    scale_factor = float(sys.argv[3])    
scaled_times = thresholded_timeseriesSIS['time'].values
# convert the time points to double values
scaled_times = scaled_times * scale_factor

# print the least squares error between the EoN simulations and the MASFENON simulation for every time point
# SIS simulation
I_sis_filtered = []
for time in scaled_times:
    #I_sis_filtered.append(I_sis[np.where(t_sis == scaled_time)[0]]) # doesn't work since the time points in the SIS simulation are double values and are increased, the comparison should be done with the closest time point
    closestTimePoint = min(t_sis, key=lambda x:abs(x-time))
    I_sis_filtered.append(I_sis[np.where(t_sis == closestTimePoint)[0]])

# flatten and sort the list based on the time points
I_sis_filtered = np.array(I_sis_filtered).flatten()
I_sis_filtered = I_sis_filtered[np.argsort(scaled_times)]

leastSquaresError = np.sum(np.square(I_sis_filtered - infectedsSIS))
print('Least squares error between EoN SIS simulation and MASFENON simulation: ' + str(leastSquaresError))
# SIR simulation
I_sir_filtered = []
for time in scaled_times:
    closestTimePoint = min(t_sir, key=lambda x:abs(x-time))
    I_sir_filtered.append(I_sir[np.where(t_sir == closestTimePoint)[0]])

I_sir_filtered = np.array(I_sir_filtered).flatten()
I_sir_filtered = I_sir_filtered[np.argsort(scaled_times)]

leastSquaresError = np.sum(np.square(I_sir_filtered - infectedsSIS))
print('Least squares error between EoN SIR simulation and MASFENON simulation: ' + str(leastSquaresError))
# SIRS simulation
I_sirs_filtered = []
for time in scaled_times:
    closestTimePoint = min(t_sirs, key=lambda x:abs(x-time))
    I_sirs_filtered.append(I_sirs[np.where(t_sirs == closestTimePoint)[0]])

I_sirs_filtered = np.array(I_sirs_filtered).flatten()
I_sirs_filtered = I_sirs_filtered[np.argsort(scaled_times)]

leastSquaresError = np.sum(np.square(I_sirs_filtered - infectedsSIS))
print('Least squares error between EoN SIRS simulation and MASFENON simulation: ' + str(leastSquaresError))

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
plt.plot(t_sis, I_sis, label='EoN SIS simulation')
plt.plot(t_sir, I_sir, label='EoN SIR simulation')
plt.plot(t_sirs, I_sirs, label='EoN SIRS simulation')
plt.plot(scaled_times, infectedsSIS, label='MASFENON simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the network of 10000 nodes')

plt.show()
