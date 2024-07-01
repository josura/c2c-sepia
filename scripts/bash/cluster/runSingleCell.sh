# submit mouse experiment with MASFENON data
fUniqueGraph="../datiIdo/metapathwayEdges.tsv"
initialPerturbationFolder="../datiIdo/inputValues"
typeInteractionsFolder="../datiIdo/interactions"
nodesDescriptionFile="../datiIdo/nodesInfo.tsv"
outputFolder="output/cluster/sepia/ido"
# run the simulation
srun ./build/c2c-sepia-MPI --fUniqueGraph $fUniqueGraph \
        --initialPerturbationPerTypeFolder $initialPerturbationFolder \
        --typeInteractionFolder $typeInteractionsFolder \
        --nodeDescriptionFile $nodesDescriptionFile \
        --dissipationModel scaled \
        --dissipationModelParameters 0.2 \
        --propagationModel neighbors \
        --propagationModelParameters 0.5 \
        --intertypeIterations 20 \
        --intratypeIterations 5 \
        --virtualNodesGranularity typeAndNode \
        --saturation \
        --undirectedEdges \
        --undirectedTypeEdges \
        --outputFolder $outputFolder