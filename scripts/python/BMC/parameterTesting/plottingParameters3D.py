import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.animation import FuncAnimation
import pandas as pd
import numpy as np
import os

# get different time series from different folder in the specified folder

network_types = ['erdosRenyi', 'preferentialAttachment', 'preferentialAttachmentAging', 'regular']

for network_type in network_types:
    # get the folders in the specified folder
    # specified_folder_changingDissipations = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterTesting/erdosRenyi/dissipations'
    # specified_folder_changingPropagations = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterTesting/erdosRenyi/propagations'
    specified_folder_changingDissipations = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterTesting/' + network_type + '/dissipations' 
    specified_folder_changingPropagations = '/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/parameterTesting/' + network_type + '/propagations'
    # DISSIPATIONS PLOTTING
    folders = [f for f in os.listdir(specified_folder_changingDissipations) if os.path.isdir(os.path.join(specified_folder_changingDissipations, f))]

    nodes_data = []
    # iterate over the folders and get the time series
    for folder in folders:
        nodes_data.append(pd.read_csv(specified_folder_changingDissipations + '/' + folder + '/allFiles/fullGraph_output.tsv', sep='\t'))

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

    dissipationScaleFactors = []
    for i in range(len(infected_nodes)):
        folder = folders[i]
        folderNoPath = folder.split('/')[-1]
        # the file has the following format: dissipationScaleFactor<float>_propagationScaleFactor<float>
        dissipationScaleFactor = folderNoPath.split('_')[0]
        dissipationScaleFactor = dissipationScaleFactor.split('dissipationScaleFactor')[1]
        propagationScaleFactor = folderNoPath.split('_')[1]
        propagationScaleFactor = propagationScaleFactor.split('propagationScaleFactor')[1]
        dissipationScaleFactors.append(float(dissipationScaleFactor))


    # plot the percentages on the z axis through time as a wireframe, on the x axis we have the time steps and on the y axis we have the changing Dissipation scale factor
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, len(infected_nodes[0]))
    ax.set_ylim(min(dissipationScaleFactors), max(dissipationScaleFactors))
    ax.set_zlim(0, 1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Dissipation Scale Factor')
    ax.set_zlabel('Percentage of Infected Nodes')
    #ax.set_title('Percentage of Infected Nodes for '+ network_type + ' network, Different Dissipation Scale Factors')
    for i in range(len(infected_nodes)):
        ax.plot(range(len(infected_nodes[i])), [dissipationScaleFactors[i]]*len(infected_nodes[i]), infected_nodes[i])

    # save the plot
    plt.savefig('/home/josura/Projects/ccc/c2c-sepia/docs/plots/epidemics/BMC/3dplot' + network_type + 'Dissipations.png')

    # PROPAGATIONS PLOTTING
    folders = [f for f in os.listdir(specified_folder_changingPropagations) if os.path.isdir(os.path.join(specified_folder_changingPropagations, f))]

    nodes_data = []
    # iterate over the folders and get the time series
    for folder in folders:
        nodes_data.append(pd.read_csv(specified_folder_changingPropagations + '/' + folder + '/allFiles/fullGraph_output.tsv', sep='\t'))

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

    propagationScaleFactors = []
    for i in range(len(infected_nodes)):
        folder = folders[i]
        folderNoPath = folder.split('/')[-1]
        # the file has the following format: dissipationScaleFactor<float>_propagationScaleFactor<float>
        dissipationScaleFactor = folderNoPath.split('_')[0]
        dissipationScaleFactor = dissipationScaleFactor.split('dissipationScaleFactor')[1]
        propagationScaleFactor = folderNoPath.split('_')[1]
        propagationScaleFactor = propagationScaleFactor.split('propagationScaleFactor')[1]
        propagationScaleFactors.append(float(propagationScaleFactor))

    # plot the percentages on the z axis through time as a wireframe, on the x axis we have the time steps and on the y axis we have the changing Propagation scale factor
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, len(infected_nodes[0]))
    ax.set_ylim(min(propagationScaleFactors), max(propagationScaleFactors))
    ax.set_zlim(0, 1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Propagation Scale Factor')
    ax.set_zlabel('Percentage of Infected Nodes')
    #ax.set_title('Percentage of Infected Nodes for '+ network_type + ' network, Different Propagation Scale Factors')
    for i in range(len(infected_nodes)):
        ax.plot(range(len(infected_nodes[i])), [propagationScaleFactors[i]]*len(infected_nodes[i]), infected_nodes[i])

    # save the plot
    plt.savefig('/home/josura/Projects/ccc/c2c-sepia/docs/plots/epidemics/BMC/3dplot' + network_type + 'Propagations.png')
