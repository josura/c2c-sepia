import networkx as nx
import matplotlib.pyplot as plt
import EoN

# get the network from one of the generated files that contains the edge data (source, target, weight)
graphEdgesFile = 'scripts/R/epidemics/significanceAnalysys/preferentialAttachmentAging/10000nodes/1/edge_data.tsv'
G = nx.read_edgelist(graphEdgesFile, delimiter='\t', data=(('weight', float),))

# set the parameters for the simulation
gamma = 1
tau = 0.2
rho = 0.005
# read the initial condition of every node from a file
initialConditionFile = 'scripts/R/epidemics/significanceAnalysys/preferentialAttachmentAging/10000nodes/1/node_condititions.tsv'