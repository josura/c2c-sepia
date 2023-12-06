# c2c-sepia
Cell-to-cell communication with enriched pathways project, using scRNA-seq data.

A framework to work with scRNA-seq profiles is presented here, the algorithm will take as an input:
- The scRNA-seq profiles and log2fold-change
- The type for every cell.

## DEPENDENCIES
- armadillo
- MPI implementation (tested on OpenMPI)

## BUILD
```bash
mkdir build
cmake .
make
```

## PARAMETER CUSTOMIZATION
The custom functions used for the scaling of the perturbation values (dissipation and conservation) are defined in <src/CustomScalingFunctions> in the source folder


## USAGE
```bash
./build/c2c-sepia --fUniqueGraph [graph].tsv --fInitialPerturbationPerType [matrix].tsv --typeInteractionFolder [typesInteractionFolder]
```

The options are the following:
- **--help**  => print help section
- **--fUniqueGraph** (obligatory if there is only a single topology structure for the agents, otherwise see) => graph filename, for an example graph see in resources or in test data. NOTE: if this option is chosen, graphsFilesFolder cannot be used")
- **--graphsFilesFolder** (obligatory if --fUniqueGraph is not specified) => graphs file folder
- **--fInitialPerturbationPerType** (obligatory if the folder with the initial inputs is not specified, this option is only available if a unique graph was specified for every agent, since all the agents will have the same topology and the inputs can be specified as a matrix in a single file) => initialPerturbationPerType matrix filename, for an example see in data
- **--initialPerturbationPerTypeFolder** (obligatory if fInitialPerturbationPerType was not specified) => initial inputs folder, every file name follows the names of the types defined 
- **--subtypes** => subtypes filename, every row one string, representing the types that should be used, the intersection of the subtypes defined and the types obtained from the graph folder or from the initialValues matrix should be non-empty
- **--typeInteractionFolder** => directory for the type interactions
- **--ensembleGeneNames**" (no parameter) => use ensemble gene names, since the metapathway used in resources have entrez_ids, a map will be done from ensemble to entrez, the map is available in resources, if the metapathway is consistent with the data used for the log-fold changes and the interactions (the have the same gene names), this parameter should not be used. Only use in case of external data sources that have ensemble ids and the metapathway used is the one in the resources.
- **--sameTypeCommunication** (no parameter) => "use same type communication, since it is not permitted as the standard definition of the model, this adds a virtual node for the same type")
- **--outputFolder** (obligatory) => output folder for output of the algorithm at each iteration
- **--intertypeIterations** => number of iterations for intertype communication
- **--intratypeIterations** => number of iterations for intratype communication
- **--timestep** => timestep to use for the iteration, the final time is iterationIntracell\*iterationIntercell\*timestep
- **--dissipationModel** => (string) the dissipation model for the computation, available models are: 'none (default)','power','random','periodic','scaled' and 'custom'
- **--dissipationModelParameters** => parameters for the dissipation model, for the power dissipation indicate the base, for the random dissipation indicate the min and max value, for the periodic dissipation indicate the period
- **--graphsFilesFolder** => graphs (pathways or other types of graphs) file folder IMPORTANT not yet implemented loading for these graphs into the different computations
- **--conservationModel** => (string) the conservation model used for the computation, available models are: 'none (default)','scaled','random' and 'custom'
- **--conservationModelParameters** => parameters for the dissipation model, for the scaled parameter the constant used to scale the conservation final results, in the case of random the upper and lower limit (between 0 and 1)
- **--saturation** => use saturation of values, default to 1, if another value is needed, use the saturationTerm
- **--saturationTerm** => defines the limits of the saturation [-saturationTerm,saturationTerm]
- **--conservateInitialNorm** => conservate the initial euclidean norm of the perturbation values, that is ||Pn|| <= ||Initial||, default to false
- **undirectedEdges** => edges in the graphs are undirected
- **undirectedTypeEdges** => edges between types are undirected
- **loggingOptions** => (string) logging options, available options are: 'all','none', default to all
    

## INPUTS 
For the structure of the input see the following reference example data in the repository:
- Graph: 
    - Generic graph: data/testdata/testGraph/edges-Graph1-general.tsv
    - entrez-id graph(the name of the nodes are entrez-ids): data/testdata/edges.tsv
- Graphs folder (the name of the files should be the types, with a tsv format and extension): data/testdata/testHeterogeneousGraph
- Initial values of the nodes for every type:
    - Matrix (only with same graph for every type): data/testdata/testGraph/initialValues-general.tsv
    - Single vector (usable in both cases of homogeneous and heterogeneous graphs): data/testdata/testHeterogeneousGraph/initialValues/t0.tsv
- Initial values files folder: data/testdata/testHeterogeneousGraph/initialValues
- Type interactions(see the file in the folder for the structure, you can create more files that will be seen during the execution in the directory specified as a parameter): data/testdata/testHeterogeneousGraph/interactions

## EXAMPLES

### different graphs, saturation at 1, dissipation scaled at 0.2, conservation scaled at 0.5
```bash
./c2c-sepia --graphsFilesFolder data/testdata/testHeterogeneousGraph/graphs \
            --initialPerturbationPerTypeFolder data/testdata/testHeterogeneousGraph/initialValuesPartialTypes \
            --typeInteractionFolder data/testdata/testHeterogeneousGraph/interactions \
            --subtypes data/testdata/testGraph/subcelltypes.txt \
            --dissipationModel scaled \
            --dissipationModelParameters 0.2 \
            --saturation \
            --conservationModel scaled \
            --conservationModelParameters 0.5 \
            --outputFolder outputs
```

## USE CASES
The framework can be used in every situation with the following structure:
- network of networks, where every network is a weighted graph and the networks are interacting with each other via some connection of two nodes inside the two networks
- every network can be abstracted to a single type (like in the biological case, every meta-pathway or pathway is related to a single cell or a single place)
- the interactions of the types are done via a node in one graph(related to a type) that has a weighted edge to a node on the other graph (related to the same type or another type), the difference between the interaction inside the graphs and outside the graph is that perturbation inside the graph is done more quickly (more iteration to propagate the values intra-network) while the perturbation outside of the type (inter-type propagation) is done more slowly and after a number of propagation intra-type.

## things to consider
The map file that is used for the mapping of ensemble gene names to entrez ids (required if the metapathway in the __resourse__ folder are used) can b e found in the metapathways folders themselves. The one used in the tool is found at __resources/graphs/metapathwayReactome/nodes.txt__

The mapping file is only available to use entrez ids for the metapathway used and loaded into this repository, if the graph used for the computation is different there is no need to the mapping file and it will not be even considered during the computation

If the custom  functions and scaling functions are used, the project need to be rebuilt with make


Under development

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
