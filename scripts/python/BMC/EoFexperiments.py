import networkx as nx
import matplotlib.pyplot as plt
import EoN

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


# run the simulation in the graphs
t, S, I, R = EoN.fast_SIR(G, tau, gamma, initial_infecteds, tmax=float('inf'), return_full_data=True)