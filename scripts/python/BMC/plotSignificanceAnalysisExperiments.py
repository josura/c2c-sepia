import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import seaborn as sns
import numpy as np
import os

# get different time series from different folder in the specified folder

# get the folders in the specified folder
specified_folder = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC'
networks_types = [f for f in os.listdir(specified_folder) if os.path.isdir(os.path.join(specified_folder, f))]
nodesNumber_list = [100, 1000, 10000]

# every folder specifiedfolder/networks_types[i]/nodesNumber_list[j]Nodes contains thirty folders that represent the different graphs, every graph folder contains the time series in a file allfiles/fullGraph_output.tsv

# list of lists of lists, the first list is the network type, the second list is the time series of the different samples, the third list is the percentage of infected nodes at each time step
all_infected_nodes_100 = []
all_infected_nodes_1000 = []
all_infected_nodes_10000 = []

# the graph data that is needed is only the nodes data, and is already available in the time series folder, so no graph needs to be loaded
for network_type in networks_types:
    for nodesNumber in nodesNumber_list:
        nodes_data = []
        network_type_nodes_folder = specified_folder + '/' + network_type + '/' + str(nodesNumber) + 'Nodes/'
        folders = [f for f in os.listdir(network_type_nodes_folder) if os.path.isdir(os.path.join(network_type_nodes_folder, f))]
        # iterate over the folders and get the time series
        if(len(folders) > 2):
            for folder in folders:
                nodes_data.append(pd.read_csv(network_type_nodes_folder + folder + '/allFiles/fullGraph_output.tsv', sep='\t'))
                # print(network_type_nodes_folder + folder + '/allFiles/fullGraph_output.tsv')
        else:  #regular graph with only the allFiles folder
            nodes_data.append(pd.read_csv(network_type_nodes_folder + 'allFiles/fullGraph_output.tsv', sep='\t'))
            # print(network_type_nodes_folder + 'allFiles/fullGraph_output.tsv')
        # create another list of dataframes with thresholded values, all the columns after the first one should be thresholded since they are the time series, 0 values means that the node is not infected and 1 means that it is infected
        thresholded_data_sorted = nodes_data.copy()
        for i in range(len(thresholded_data_sorted)):
            thresholded_data_sorted[i].iloc[:, 1:] = np.where(thresholded_data_sorted[i].iloc[:, 1:] > 0.5, 1, 0)
            # order the rows by the iteration column
            thresholded_data_sorted[i] = thresholded_data_sorted[i].sort_values(by=['time'])

        # create a list of lists with the percentage of infected at each time step for different folders, the index of the list is the time step (iteration column in the dataframe)
        infected_nodes = []
        for i in range(len(thresholded_data_sorted)):
            infected_nodes.append([])
            for _, row in thresholded_data_sorted[i].iterrows():
                infected_nodes[i].append(len(row[row == 1])/len(row))
        

        # create another list of lists with the length of one of the infected node lists(the timesteps that are all the same) , effectively transposing the infected nodes list
        infected_nodes_transposed = []
        for i in range(len(infected_nodes[0])):
            infected_nodes_transposed.append([])
            for j in range(len(infected_nodes)):
                infected_nodes_transposed[i].append(infected_nodes[j][i])

        # print the infected nodes transposed list sizes
        for i in range(len(infected_nodes_transposed)):
            print('timestep ' + str(i) + ' has ' + str(len(infected_nodes_transposed[i])) + ' samples')
        
        if nodesNumber == 100:
            all_infected_nodes_100.append (infected_nodes_transposed)
        elif nodesNumber == 1000:
            all_infected_nodes_1000.append (infected_nodes_transposed)
        elif nodesNumber == 10000:
            all_infected_nodes_10000.append (infected_nodes_transposed)


        # # add the violin plot of the percentage of infected nodes at a specific time, the different folders are the samples, only one violin plot per network type, nodes number and time step
        # fig, ax = plt.subplots()
        # ax.violinplot(infected_nodes_transposed, showmeans=False, showmedians=True)
        # ax.set_xlabel('Time')
        # ax.set_ylabel('Density of Infected Nodes')
        # ax.set_title('Density of Infected Nodes, ' + network_type + ' Network, ' + str(nodesNumber) + ' Nodes')
        # # pause the execution to see the plot, close the plot to continue with the next on
        # plt.show()
        
# save the data into a pandas dataframe for further analysis, the rows are the samples and the columns are the timesteps, the values are the percentage of infected nodes, an additional column is added to specify the network type

pd_all_infected_nodes_100 = pd.DataFrame(columns=['networkType', 'nodesNumber', 'timeStep', 'sample', 'percentageInfected'])
pd_all_infected_nodes_1000 = pd.DataFrame(columns=['networkType', 'nodesNumber', 'timeStep', 'sample', 'percentageInfected'])
pd_all_infected_nodes_10000 = pd.DataFrame(columns=['networkType', 'nodesNumber', 'timeStep', 'sample', 'percentageInfected'])

for networkTypeIndex in range(len(all_infected_nodes_100)):
    for timeStepIndex in range(len(all_infected_nodes_100[networkTypeIndex])):
        for sampleIndex in range(len(all_infected_nodes_100[networkTypeIndex][timeStepIndex])):
            #ignore regular graphs
            if (networks_types[networkTypeIndex] == 'regular'):
                continue
            new_row_df = pd.DataFrame({'networkType': networks_types[networkTypeIndex], 'nodesNumber': 100, 'timeStep': timeStepIndex, 'sample': sampleIndex, 'percentageInfected': all_infected_nodes_100[networkTypeIndex][timeStepIndex][sampleIndex]}, index=[0])
            pd_all_infected_nodes_100 = pd.concat([pd_all_infected_nodes_100, new_row_df], ignore_index=True)

for networkTypeIndex in range(len(all_infected_nodes_1000)):
    for timeStepIndex in range(len(all_infected_nodes_1000[networkTypeIndex])):
        for sampleIndex in range(len(all_infected_nodes_1000[networkTypeIndex][timeStepIndex])):
            if (networks_types[networkTypeIndex] == 'regular'):
                continue
            new_row_df = pd.DataFrame({'networkType': networks_types[networkTypeIndex], 'nodesNumber': 1000, 'timeStep': timeStepIndex, 'sample': sampleIndex, 'percentageInfected': all_infected_nodes_1000[networkTypeIndex][timeStepIndex][sampleIndex]}, index=[0])
            pd_all_infected_nodes_1000 = pd.concat([pd_all_infected_nodes_1000, new_row_df], ignore_index=True)
            
for networkTypeIndex in range(len(all_infected_nodes_10000)):
    for timeStepIndex in range(len(all_infected_nodes_10000[networkTypeIndex])):
        for sampleIndex in range(len(all_infected_nodes_10000[networkTypeIndex][timeStepIndex])):
            # if (networks_types[networkTypeIndex] == 'regular'):
            #     continue
            new_row_df = pd.DataFrame({'networkType': networks_types[networkTypeIndex], 'nodesNumber': 10000, 'timeStep': timeStepIndex, 'sample': sampleIndex, 'percentageInfected': all_infected_nodes_10000[networkTypeIndex][timeStepIndex][sampleIndex]}, index=[0])
            pd_all_infected_nodes_10000 = pd.concat([pd_all_infected_nodes_10000, new_row_df], ignore_index=True)
            
fig, ax = plt.subplots()

#sns.set(rc={"figure.figsize":(8, 4)}) #width=8, height=4

sns.violinplot(ax = ax,
               data = pd_all_infected_nodes_100,
               x = 'timeStep',
               y = 'percentageInfected',
               hue = 'networkType',
               split = False)
#set title
plt.title('Percentage of Infected Nodes, 100 Nodes')
plt.show()

fig, ax = plt.subplots()

sns.violinplot(ax = ax,
                data = pd_all_infected_nodes_1000,
                x = 'timeStep',
                y = 'percentageInfected',
                hue = 'networkType',
                split = False)
#set title
plt.title('Percentage of Infected Nodes, 1000 Nodes')
plt.show()

fig, ax = plt.subplots()

sns.violinplot(ax = ax,
                data = pd_all_infected_nodes_10000,
                x = 'timeStep',
                y = 'percentageInfected',
                hue = 'networkType',
                split = False)
#set title
plt.title('Percentage of Infected Nodes, 10000 Nodes')
plt.show()

# boxplots
fig, ax = plt.subplots()


sns.boxplot(ax = ax,
               data = pd_all_infected_nodes_100,
               x = 'timeStep',
               y = 'percentageInfected',
               hue = 'networkType')
#set title
plt.title('Percentage of Infected Nodes, 100 Nodes')
plt.show()

fig, ax = plt.subplots()

sns.boxplot(ax = ax,
                data = pd_all_infected_nodes_1000,
                x = 'timeStep',
                y = 'percentageInfected',
                hue = 'networkType')
#set title
plt.title('Percentage of Infected Nodes, 1000 Nodes')
plt.show()

fig, ax = plt.subplots()

sns.boxplot(ax = ax,
                data = pd_all_infected_nodes_10000,
                x = 'timeStep',
                y = 'percentageInfected',
                hue = 'networkType')
#set title
plt.title('Percentage of Infected Nodes, 10000 Nodes')
plt.show()
          
# Shapiro-Wilk test for normality for every combination of timestep and network type
from scipy.stats import shapiro

pd_shapiro_100 = pd.DataFrame(columns=['networkType', 'nodesNumber', 'timeStep', 'p-value'])
for networkTypeIndex in range(len(all_infected_nodes_100)):
    for timeStepIndex in range(len(all_infected_nodes_100[networkTypeIndex])):
        if (networks_types[networkTypeIndex] == 'regular'):
            continue
        stat, p = shapiro(all_infected_nodes_100[networkTypeIndex][timeStepIndex])
        print('network type: ' + networks_types[networkTypeIndex] + ', timestep: ' + str(timeStepIndex) + ', p-value: ' + str(p))
        new_row_df = pd.DataFrame({'networkType': networks_types[networkTypeIndex], 'nodesNumber': 100, 'timeStep': timeStepIndex, 'p-value': p}, index=[0])
        pd_shapiro_100 = pd.concat([pd_shapiro_100, new_row_df], ignore_index=True)

pd_shapiro_1000 = pd.DataFrame(columns=['networkType', 'nodesNumber', 'timeStep', 'p-value'])
for networkTypeIndex in range(len(all_infected_nodes_1000)):
    for timeStepIndex in range(len(all_infected_nodes_1000[networkTypeIndex])):
        if (networks_types[networkTypeIndex] == 'regular'):
            continue
        stat, p = shapiro(all_infected_nodes_1000[networkTypeIndex][timeStepIndex])
        print('network type: ' + networks_types[networkTypeIndex] + ', timestep: ' + str(timeStepIndex) + ', p-value: ' + str(p))
        new_row_df = pd.DataFrame({'networkType': networks_types[networkTypeIndex], 'nodesNumber': 1000, 'timeStep': timeStepIndex, 'p-value': p}, index=[0])
        pd_shapiro_1000 = pd.concat([pd_shapiro_1000, new_row_df], ignore_index=True)

pd_shapiro_10000 = pd.DataFrame(columns=['networkType', 'nodesNumber', 'timeStep', 'p-value'])
for networkTypeIndex in range(len(all_infected_nodes_10000)):
    for timeStepIndex in range(len(all_infected_nodes_10000[networkTypeIndex])):
        if (networks_types[networkTypeIndex] == 'regular'):
            continue
        stat, p = shapiro(all_infected_nodes_10000[networkTypeIndex][timeStepIndex])
        print('network type: ' + networks_types[networkTypeIndex] + ', timestep: ' + str(timeStepIndex) + ', p-value: ' + str(p))
        new_row_df = pd.DataFrame({'networkType': networks_types[networkTypeIndex], 'nodesNumber': 10000, 'timeStep': timeStepIndex, 'p-value': p}, index=[0])
        pd_shapiro_10000 = pd.concat([pd_shapiro_10000, new_row_df], ignore_index=True)

# plot the p-values of the Shapiro-Wilk test
fig, ax = plt.subplots()
sns.barplot(ax = ax,
            data = pd_shapiro_100,
            x = 'timeStep',
            y = 'p-value',
            hue = 'networkType')
plt.title('Shapiro-Wilk Test p-values, 100 Nodes')
# add significance line
plt.axhline(0.05, color='r', linestyle='--')
plt.show()

fig, ax = plt.subplots()
sns.barplot(ax = ax,
            data = pd_shapiro_1000,
            x = 'timeStep',
            y = 'p-value',
            hue = 'networkType')
plt.title('Shapiro-Wilk Test p-values, 1000 Nodes')
plt.axhline(0.05, color='r', linestyle='--')
plt.show()

fig, ax = plt.subplots()
sns.barplot(ax = ax,
            data = pd_shapiro_10000,
            x = 'timeStep',
            y = 'p-value',
            hue = 'networkType')
plt.title('Shapiro-Wilk Test p-values, 10000 Nodes')
plt.axhline(0.05, color='r', linestyle='--')
plt.show()