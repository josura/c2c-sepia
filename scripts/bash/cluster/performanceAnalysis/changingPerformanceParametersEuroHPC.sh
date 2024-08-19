# fill this part with the configuration for slurm

# control if parameters are passed, if only three arguments are passed, the maximum number of processors is also passed
if [ $# -ne 3 ]; then
    echo $# arguments passed
    echo "Please pass the folder with the inputs and the folder where to save the outputs and the maximum number of processors"
    exit 1
fi

inputsFolder=$1
outputFolder=$2
processors=$3

srun -n $processors ./build/c2c-sepia-MPI --fUniqueGraph $inputsFolder/metapathwayEdges.tsv \
    --initialPerturbationPerTypeFolder $inputsFolder/inputValues \
    --typeInteractionFolder $inputsFolder/interactions \
    --nodeDescriptionFile $inputsFolder/nodesInfo.tsv \
    --dissipationModel scaled \
    --dissipationModelParameters 0.5 \
    --propagationModel neighbors \
    --propagationModelParameters 0.5 \
    --intertypeIterations 20 \
    --intratypeIterations 5 \
    --virtualNodesGranularity typeAndNode \
    --saturation \
    --outputFolder $outputFolder
