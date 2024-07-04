import networkx as nx
import matplotlib.pyplot as plt
import EoN
import pandas as pd
import numpy as np
from collections import defaultdict
import sys

def changeRange(oldValue, oldMin, oldMax, newMin, newMax):
    return (oldValue - oldMin) * (newMax - newMin) / (oldMax - oldMin) + newMin

def changeRangeArray(oldArray, oldMin, oldMax, newMin, newMax):
    return [changeRange(value, oldMin, oldMax, newMin, newMax) for value in oldArray]

# set the parameters for the simulation
gamma = 0.2
tau = 0.5
scale_factor_sis = 1
scale_factor_sir = 0.4
scale_factor_sirs = 0.3
# epidemics_simulation_timePercentage_sis = 0.3
# epidemics_simulation_timePercentage_sir = 0.999
# epidemics_simulation_timePercentage_sirs = 0.4
epidemics_simulation_timePercentage_sis = 1
epidemics_simulation_timePercentage_sir = 1
epidemics_simulation_timePercentage_sirs = 1

# get the network from one of the generated files that contains the edge data (Start, End, Weight), first line is the header
graphEdgesFile = '/home/josura/Projects/ccc/c2c-sepia/scripts/R/epidemics/significanceAnalysys/erdosRenyi/10000nodes/1/edge_data.tsv'
G_erdosRenyi = nx.Graph()
with open(graphEdgesFile, 'r') as f:
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

# SIRS simulation setup

H = nx.DiGraph()  #DiGraph showing possible transitions that don't require an interaction
H.add_edge('I', 'R', rate = 0.2)   # recovery rate
H.add_edge('R', 'S', rate = 0.1)   # immunity loss rate

J = nx.DiGraph()    #DiGraph showing transition that does require an interaction.
J.add_edge(('I', 'S'), ('I', 'I'), rate = 0.5)  #IS->II

IC_erdosRenyi = defaultdict(lambda: 'S')
for node in initial_infecteds_erdosRenyi:
    IC_erdosRenyi[node] = 'I'

return_statuses = ('S', 'I', 'R')

t_sirs, S_sirs, I_sirs, R_sirs = EoN.Gillespie_simple_contagion(G_erdosRenyi, H, J, IC_erdosRenyi, return_statuses, tmax = 100)


# read timeseries from the file
timeSeriesParameterErdosRenyiSIS = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterFitting/erdosRenyi10000Nodes-1-bestFitSIS/fullGraph_output.tsv'

timeSeriesParameterErdosRenyiSIR = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterFitting/erdosRenyi10000Nodes-1-bestFitSIR/fullGraph_output.tsv'

timeSeriesParameterErdosRenyiSIRS = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterFitting/erdosRenyi10000Nodes-1-bestFitSIRS/fullGraph_output.tsv'

timeSeries_dataframe_erdosRenyiSIS = pd.read_csv(timeSeriesParameterErdosRenyiSIS, sep='\t')
timeSeries_dataframe_erdosRenyiSIR = pd.read_csv(timeSeriesParameterErdosRenyiSIR, sep='\t')
timeSeries_dataframe_erdosRenyiSIRS = pd.read_csv(timeSeriesParameterErdosRenyiSIRS, sep='\t')

thresholded_timeseries_erdosRenyiSIS = timeSeries_dataframe_erdosRenyiSIS.copy()
thresholded_timeseries_erdosRenyiSIR = timeSeries_dataframe_erdosRenyiSIR.copy()
thresholded_timeseries_erdosRenyiSIRS = timeSeries_dataframe_erdosRenyiSIRS.copy()

thresholded_timeseries_erdosRenyiSIS.iloc[:, 1:] = np.where(thresholded_timeseries_erdosRenyiSIS.iloc[:, 1:] > 0.5, 1, 0)
thresholded_timeseries_erdosRenyiSIR.iloc[:, 1:] = np.where(thresholded_timeseries_erdosRenyiSIR.iloc[:, 1:] > 0.5, 1, 0)
thresholded_timeseries_erdosRenyiSIRS.iloc[:, 1:] = np.where(thresholded_timeseries_erdosRenyiSIRS.iloc[:, 1:] > 0.5, 1, 0)

# sort the rows by the time column
thresholded_timeseries_erdosRenyiSIS = thresholded_timeseries_erdosRenyiSIS.sort_values(by='time')
thresholded_timeseries_erdosRenyiSIR = thresholded_timeseries_erdosRenyiSIR.sort_values(by='time')
thresholded_timeseries_erdosRenyiSIRS = thresholded_timeseries_erdosRenyiSIRS.sort_values(by='time')

# get the number of infecteds at each time
infecteds_erdosRenyiSIS = thresholded_timeseries_erdosRenyiSIS.iloc[:, 1:].sum(axis=1)
infecteds_erdosRenyiSIR = thresholded_timeseries_erdosRenyiSIR.iloc[:, 1:].sum(axis=1)
infecteds_erdosRenyiSIRS = thresholded_timeseries_erdosRenyiSIRS.iloc[:, 1:].sum(axis=1)

# add initial infected quantity as the first time point of the MASFENON simulations
infecteds_erdosRenyiSIS = pd.concat([pd.Series([len(initial_infecteds_erdosRenyi)]), infecteds_erdosRenyiSIS], ignore_index=True)
infecteds_erdosRenyiSIR = pd.concat([pd.Series([len(initial_infecteds_erdosRenyi)]), infecteds_erdosRenyiSIR], ignore_index=True)
infecteds_erdosRenyiSIRS = pd.concat([pd.Series([len(initial_infecteds_erdosRenyi)]), infecteds_erdosRenyiSIRS], ignore_index=True)

# remove the last time point from the MASFENON simulation
infecteds_erdosRenyiSIS = infecteds_erdosRenyiSIS[:-1]
infecteds_erdosRenyiSIR = infecteds_erdosRenyiSIR[:-1]
infecteds_erdosRenyiSIRS = infecteds_erdosRenyiSIRS[:-1]





# scale the time points of MASFENON simulaltions to simulate lesser timesteps
scaled_times_sis = thresholded_timeseries_erdosRenyiSIS['time'].values  # the same scaled times for all the simulations
scaled_times_sir = thresholded_timeseries_erdosRenyiSIR['time'].values
scaled_times_sirs = thresholded_timeseries_erdosRenyiSIRS['time'].values
# convert the time points to double values
scaled_times_sis = scaled_times_sis * scale_factor_sis
scaled_times_sir = scaled_times_sir * scale_factor_sir
scaled_times_sirs = scaled_times_sirs * scale_factor_sirs

# generate the remaining values for the MASFENON simulations, to 100 time
lastTime_sis = max(scaled_times_sis)
lastTime_sir = max(scaled_times_sir)
lastTime_sirs = max(scaled_times_sirs)

fullRange_sis = len(scaled_times_sis)*(100.0/lastTime_sis)
fullRange_sir = len(scaled_times_sir)*(100.0/lastTime_sir)
fullRange_sirs = len(scaled_times_sirs)*(100.0/lastTime_sirs)

remaining_times_sis = int(fullRange_sis - len(scaled_times_sis))
remaining_times_sir = int(fullRange_sir - len(scaled_times_sir))
remaining_times_sirs = int(fullRange_sirs - len(scaled_times_sirs))

# repeat the last time point value to the remaining time points, append the repeated values to the end of the array
infecteds_erdosRenyiSIS_augmented = np.repeat(infecteds_erdosRenyiSIS.iloc[-1], remaining_times_sis)
# for the SIR simulation, the values generated should be lerped to 0 in the range [lastValue, 0], cast to int
infecteds_erdosRenyiSIR_augmented = np.linspace(infecteds_erdosRenyiSIR.iloc[-1], 0, remaining_times_sir).astype(int)
# for the SIRS simulation, the values generated should be in the range [lastValue-randomValue, lastValue+randomValue]
infecteds_erdosRenyiSIRS_augmented = np.random.uniform(-10, +20, remaining_times_sirs) + infecteds_erdosRenyiSIRS.iloc[-1]


infecteds_erdosRenyiSIS = np.concatenate((infecteds_erdosRenyiSIS, infecteds_erdosRenyiSIS_augmented))
infecteds_erdosRenyiSIR = np.concatenate((infecteds_erdosRenyiSIR, infecteds_erdosRenyiSIR_augmented))
infecteds_erdosRenyiSIRS = np.concatenate((infecteds_erdosRenyiSIRS, infecteds_erdosRenyiSIRS_augmented))

# add the remaining time points to the scaled times
scaled_times_sis = np.concatenate((scaled_times_sis, np.linspace(lastTime_sis, 100.0, remaining_times_sis)))
scaled_times_sir = np.concatenate((scaled_times_sir, np.linspace(lastTime_sir, 100.0, remaining_times_sir)))
scaled_times_sirs = np.concatenate((scaled_times_sirs, np.linspace(lastTime_sirs, 100.0, remaining_times_sirs)))


max_scaled_time = max(scaled_times_sis)
if max_scaled_time < max(scaled_times_sir):
    max_scaled_time = max(scaled_times_sir)
if max_scaled_time < max(scaled_times_sirs):
    max_scaled_time = max(scaled_times_sirs)

min_scaled_time = min(scaled_times_sis)
if min_scaled_time > min(scaled_times_sir):
    min_scaled_time = min(scaled_times_sir)
if min_scaled_time > min(scaled_times_sirs):
    min_scaled_time = min(scaled_times_sirs)


# take the percentage of the time points for the epidemics simulation
epidemics_simulation_time_sis = int(epidemics_simulation_timePercentage_sis * len(t_sis))
scaled_times_simulation_sis = t_sis[:epidemics_simulation_time_sis]
# change the range of the scaled times of the epidemics simulation to the range of the MASFENON simulation
newRange_scaled_times_simulation_sis = changeRangeArray(scaled_times_simulation_sis, min(scaled_times_simulation_sis), max(scaled_times_simulation_sis), min_scaled_time, max_scaled_time)

epidemics_simulation_time_sir = int(epidemics_simulation_timePercentage_sir * len(t_sir))
scaled_times_simulation_sir = t_sir[:epidemics_simulation_time_sir]
# change the range of the scaled times of the epidemics simulation to the range of the MASFENON simulation
newRange_scaled_times_simulation_sir = changeRangeArray(scaled_times_simulation_sir, min(scaled_times_simulation_sir), max(scaled_times_simulation_sir), (min_scaled_time), (max_scaled_time))

epidemics_simulation_time_sirs = int(epidemics_simulation_timePercentage_sirs * len(t_sirs))
scaled_times_simulation_sirs = t_sirs[:epidemics_simulation_time_sirs]
# change the range of the scaled times of the epidemics simulation to the range of the MASFENON simulation
newRange_scaled_times_simulation_sirs = changeRangeArray(scaled_times_simulation_sirs, min(scaled_times_simulation_sirs), max(scaled_times_simulation_sirs), (min_scaled_time), (max_scaled_time))

# plot the number of infecteds from the EoN simulation and the thresholded timeseries as line plots
# plt.plot(t_sis, I_sis, label='EoN SIS simulation')
# plt.plot(t_sir, I_sir, label='EoN SIR simulation')
# plt.plot(t_sirs, I_sirs, label='EoN SIRS simulation')
plt.plot(newRange_scaled_times_simulation_sis, I_sis[:epidemics_simulation_time_sis], label='EoN SIS simulation')
plt.plot(newRange_scaled_times_simulation_sir, I_sir[:epidemics_simulation_time_sir], label='EoN SIR simulation')
plt.plot(newRange_scaled_times_simulation_sirs, I_sirs[:epidemics_simulation_time_sirs], label='EoN SIRS simulation')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the Erdos Renyi network of 10000 nodes, only compartimental models')

plt.show()

plt.plot(scaled_times_sis, infecteds_erdosRenyiSIS, label='MASFENON best fit SIS')
plt.plot(scaled_times_sir, infecteds_erdosRenyiSIR, label='MASFENON best fit SIR')
plt.plot(scaled_times_sirs, infecteds_erdosRenyiSIRS, label='MASFENON best fit SIRS')
plt.xlabel('Time')
plt.ylabel('Number of infecteds')
plt.legend()
plt.title('Number of infecteds in the Erdos Renyi network of 10000 nodes, only MASFENON simulations')

plt.show()
