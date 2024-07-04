export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=32

# generate 10 scale factors for the dissipation, from 0 to 1
listDissipationScaleFactors=$(seq 0.0 0.10 1.0)
# generate 10 scale factors for the propagation, from 0 to 10
listPropagationScaleFactors=$(seq 0.0 1.0 10.0)


# submit mouse experiment with MASFENON data
fUniqueGraph="../datiIdo/metapathwayEdges.tsv"
initialPerturbationFolder="../datiIdo/inputValues"
typeInteractionsFolder="../datiIdo/interactions"
nodesDescriptionFile="../datiIdo/nodesInfo.tsv"
outputFolder="/leonardo_scratch/fast/CNHPC_1452526/datiIdo_results_changing"
# run the simulations

# run the simulation for different numbers of dissipation scale factors
propagationScaleFactor=0.5
for dissipationScaleFactor in ${listDissipationScaleFactors[@]}; do
        outputFolderName="$outputFolder/changingDissipations/dissipationScaleFactor${dissipationScaleFactor}_propagationScaleFactor${propagationScaleFactor}"
        if [ ! -d $outputFolderName ]; then
                mkdir -p $outputFolderName
        fi
        srun ./build/c2c-sepia-MPI --fUniqueGraph $fUniqueGraph \
                --initialPerturbationPerTypeFolder $initialPerturbationFolder \
                --typeInteractionFolder $typeInteractionsFolder \
                --nodeDescriptionFile $nodesDescriptionFile \
                --dissipationModel scaled \
                --dissipationModelParameters $dissipationScaleFactor \
                --propagationModel neighbors \
                --propagationModelParameters $propagationScaleFactor \
                --intertypeIterations 20 \
                --intratypeIterations 5 \
                --virtualNodesGranularity typeAndNode \
                --saturation \
                --outputFolder $outputFolderName
done


# run the simulation for different numbers of propagation scale factors
dissipationScaleFactor=0.5
for propagationScaleFactor in ${listPropagationScaleFactors[@]}; do
        outputFolderName="$outputFolder/changingPropagations/dissipationScaleFactor${dissipationScaleFactor}_propagationScaleFactor${propagationScaleFactor}"
        if [ ! -d $outputFolderName ]; then
                mkdir -p $outputFolderName
        fi
        srun ./build/c2c-sepia-MPI --fUniqueGraph $fUniqueGraph \
                --initialPerturbationPerTypeFolder $initialPerturbationFolder \
                --typeInteractionFolder $typeInteractionsFolder \
                --nodeDescriptionFile $nodesDescriptionFile \
                --dissipationModel scaled \
                --dissipationModelParameters $dissipationScaleFactor \
                --propagationModel neighbors \
                --propagationModelParameters $propagationScaleFactor \
                --intertypeIterations 20 \
                --intratypeIterations 5 \
                --virtualNodesGranularity typeAndNode \
                --saturation \
                --outputFolder $outputFolderName
done