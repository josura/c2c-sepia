
# MPI examples

This page contains a collection of MPI examples of running the framework on test data.

## Homogeneous graph
    
### Example 1: custom propagation model (neighbors if the code in src/propagationModelCustom.cxx is unmodified), scaled dissipation and conservation models set at constants, saturation at default value (1.0)
```bash
mpirun -np 2 ./build/c2c-sepia-MPI --fUniqueGraph data/testdata/testGraph/edges-Graph1-general.tsv \
            --fInitialPerturbation data/testdata/testGraph/initialValues-general.tsv \
            --fInteractions data/testdata/testGraph/interactions.tsv \
            --dissipationModel scaled \
            --dissipationModelParameters 0.2 \
            --conservationModel scaled \
            --conservationModelParameters 0.2 \
            --propagationModel customPropagation \
            --propagationModelParameters 0.2 \
            --saturation \
            --outputFolder outputs/testingMPIsingle
```

## Heterogeneous temporal graph

### Example 1: custom propagation model (neighbors if the code in src/propagationModelCustom.cxx is unmodified), scaled dissipation and conservation models set at constants, saturation at default value (1.0)
```bash
mpirun -np 2 ./build/c2c-sepia-MPI --graphsFilesFolder data/testdata/testHeterogeneousGraph/graphs \
            --initialPerturbationPerTypeFolder data/testdata/testHeterogeneousTemporalGraph/initialValuesPartialTypes \
            --typeInteractionFolder data/testdata/testHeterogeneousTemporalGraph/interactions \
            --dissipationModel scaled \
            --dissipationModelParameters 0.2 \
            --conservationModel scaled \
            --conservationModelParameters 0.2 \
            --propagationModel customPropagation \
            --propagationModelParameters 0.2 \
            --virtualNodesGranularity typeAndNode \
            --saturation \
            --outputFolder outputs/testingMPIgranularDifferentProcessors
```