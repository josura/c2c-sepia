# steps taken to generate the data

## changing some HMDB IDs to match the ones in the database
Alpha-D-glucose HMDB0003345 -> D-Glucose HMDB0000122   https://hmdb.ca/metabolites/HMDB0000122

## All metabolites are not in the namemap, so we should still generate the rest of the metabolites
- changed pyrimidine since it was not in the namemap
- CC00000 is a placeholder for the modules that create the same metabolite as the source module. If these modules are considered, loops are created in the graph, but this type of behaviour is too unstable, at least until we have a coarser granularity in the metabolite layer.
- some metabolites are not really defined and also use placeholders like CC00001, for example Steroid_hormone or Fatty Acid. I don't know how to treat these metabolites.

## coarser VS finer granularity in the creation of the metabolite layer
coarser granularity makes the metabolites in every module collapse into one single metabolite, while finer granularity makes the metabolites in every module remain separate so they will have their names considered to add them in the database used to find the interactions (metalinks should be expanded with these metabolites)

## consistency in the generated inputs
consistency in the generated inputs should always be paramount since the algorithm requires it and it is very important to define the data used:
- generated graphs for the metabolite layer should be the same for every celltype, the only thing that changes are the weights of the edges(given byu the fluxes)
- generated graphs for the gene layer could be different for every celltype, the structure and weights differ
- generated graphs for one single celltype should have a consistent strategy to integrate the different time points:
  - only the first time point is considered for the generation of the graphs for every celltype
  - every time point is considered for the generation of the graphs for every celltype. This means that all the genes should be unioned or intersected to generate the gene layer.

## flux rate are considered to generate the weights for the metabolite layer
the flux rate is used to generate the weights for the metabolite layer, the higher the flux rate, the higher the weight. Since the metabolite layer could be dependent on the type of granularity used, for same interactions modules, the weights could be collapsed into one single weight (for coarser granularity) or remain separate (for finer granularity). Also, flux rates are given for every cell, so some kind of aggregation should be done to generate the weights for the metabolite layer for the single celltype.

## selection of highly variable cells and treatment of 7h data or every time point separately
there are some problems with the 7h data that seems to be related to the low quality of it (there are some Nan or infinite values during the rank aggregate function to establish the interaction between metabolites and genes), the possible approaches are:
- remove the 7h data, only consider it for training
- consider the 7h data, treat it separately from the other time points to make the code work
- consider the 7h data, treat it as the other time points, standardize a selection of the genes that are highly variable.
- consider the 7h data, impute the data (high level of rework for the code), this should be considered as the last option because almost all of the code used to generate the data should be changed, and an additional analysis on the significance of the imputation should be done, since the imputation could introduce high levels of noise in the data.