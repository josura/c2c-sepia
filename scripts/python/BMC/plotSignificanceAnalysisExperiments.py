import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np
import os

# get different time series from different folder in the specified folder

# get the folders in the specified folder
specified_folder = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics-BMC'
networks_types = [f for f in os.listdir(specified_folder) if os.path.isdir(os.path.join(specified_folder, f))]
nodesNumber_list = [100, 1000, 10000]

# every folder specifiedfolder/networks_types[i]/nodesNumber_list[j]Nodes contains thirty folders that represent the different graphs, every graph folder contains the time series in a file allfiles/fullGraph_output.tsv

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


        # add the violin plot of the percentage of infected nodes at a specific time, the different folders are the samples, only one violin plot per network type, nodes number and time step
        fig, ax = plt.subplots()
        ax.violinplot(infected_nodes_transposed, showmeans=False, showmedians=True)
        ax.set_xlabel('Time')
        ax.set_ylabel('Density of Infected Nodes')
        ax.set_title('Density of Infected Nodes, ' + network_type + ' Network, ' + str(nodesNumber) + ' Nodes')
        # pause the execution to see the plot, close the plot to continue with the next on
        plt.show()
        


    





# plot the percentages through time, on the x axis we have the time steps and on the y axis we have the percentage of infected nodes, each line represents a different folder
# fig, ax = plt.subplots()
# ax.set_xlim(0, len(infected_nodes[0]))
# ax.set_ylim(0, 1)
# ax.set_xlabel('Time')
# ax.set_ylabel('Percentage of Infected Nodes')
# ax.set_title('Percentage of Infected Nodes, Different Parameters (Dissipation scale factor, Propagation scale factor)')
# for i in range(len(infected_nodes)):
#     folder = folders[i]
#     folderNoPath = folder.split('/')[-1]
#     # the file has the following format: dissipationScaleFactor<float>_propagationScaleFactor<float>
#     dissipationScaleFactor = folderNoPath.split('_')[0]
#     dissipationScaleFactor = dissipationScaleFactor.split('dissipationScaleFactor')[1]
#     propagationScaleFactor = folderNoPath.split('_')[1]
#     propagationScaleFactor = propagationScaleFactor.split('propagationScaleFactor')[1]
#     label = '( ' + dissipationScaleFactor + ' , ' + propagationScaleFactor + ' )'
#     ax.plot(infected_nodes[i], label=label)
# # augment dimensions of the legend box
# ax.legend(title="(\u03BB (n), \u03C9 (n))", fontsize=12)
# plt.show()

