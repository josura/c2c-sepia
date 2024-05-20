
# control if parameters are passed, if only three arguments are passed, the number of processors is set in the main script
if [ $# -ne 2 ]; then
    echo $# arguments passed
    echo "Please pass the folder with the inputs and the folder where to save the outputs"
    exit 1
fi



#example usage: sh scripts/bash/experiments/epidemicsExperimentsBMC.sh wlan0 scripts/R/epidemics/syntheticGraphs/graphtype/numberOfNodes outputs/epidemics/graphtype/numberOfNodes 
# get the directory of the inputs from the first argument
inputsFolder=$1
outputsFolder=$2


#get the graph file, the initial perturbation and the type interactions from the folder
graphsFolderName=$(ls $inputsFolder | grep "graph")
initialPerturbationFolderName=$(ls $inputsFolder | grep "node_conditions_discr")
typeInteractionsFolderName=$(ls $inputsFolder | grep "interactions")

# get the full path for the graph file, the initial perturbation and the type interactions 
graphsFolder=$inputsFolder/$graphsFolderName
initialPerturbationFolder=$inputsFolder/$initialPerturbationFolderName
typeInteractionsFolder=$inputsFolder/$typeInteractionsFolderName
nodesFolder=$inputsFolder/communities

echo "Graph file: $graphsFolderName"
echo "Initial perturbation file: $initialPerturbationFolderName"
echo "Type interactions file: $typeInteractionsFolderName"
echo "Nodes file: $nodesFolderName"

# generate 50 scale factors for the dissipation, from 0 to 1
listDissipationScaleFactors=$(seq 0.0 0.02 1.0)
# generate 50 scale factors for the propagation, from 0 to 10
listPropagationScaleFactors=$(seq 0.0 0.2 10.0)

# loop through the propagation scale factors
dissipationScaleFactor=0.5
for propagationScaleFactor in ${listPropagationScaleFactors[@]}; do
    # echo "Dissipation scale factor: $dissipationScaleFactor"
    # echo "Propagation scale factor: $propagationScaleFactor"
    # get the output folder
    outputFolder="$outputsFolder/changingDissipations/dissipationScaleFactor${dissipationScaleFactor}_propagationScaleFactor${propagationScaleFactor}"
    # run the simulation
    echo "mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include wlan0 -np 12 ./build/c2c-sepia-MPI --graphsFilesFolder $graphsFolder \
        --initialPerturbationPerTypeFolder $initialPerturbationFolder \
        --typeInteractionFolder $typeInteractionsFolder \
        --nodeDescriptionFolder $nodesFolder \
        --dissipationModel scaled \
        --dissipationModelParameters $dissipationScaleFactor \
        --propagationModel neighbors \
        --propagationModelParameters $propagationScaleFactor \
        --intertypeIterations 20 \
        --intratypeIterations 5 \
        --virtualNodesGranularity typeAndNode \
        --saturation \
        --undirectedEdges \
        --undirectedTypeEdges \
        --outputFolder $outputFolder"
done

# loop through the dissipation scale factors
propagationScaleFactor=1.0
for dissipationScaleFactor in ${listDissipationScaleFactors[@]}; do
    # echo "Dissipation scale factor: $dissipationScaleFactor"
    # echo "Propagation scale factor: $propagationScaleFactor"
    # get the output folder
    outputFolder="$outputsFolder/changingPropagations/dissipationScaleFactor${dissipationScaleFactor}_propagationScaleFactor${propagationScaleFactor}"
    # run the simulation
    echo "mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include wlan0 -np 12 ./build/c2c-sepia-MPI --graphsFilesFolder $graphsFolder \
        --initialPerturbationPerTypeFolder $initialPerturbationFolder \
        --typeInteractionFolder $typeInteractionsFolder \
        --nodeDescriptionFolder $nodesFolder \
        --dissipationModel scaled \
        --dissipationModelParameters $dissipationScaleFactor \
        --propagationModel neighbors \
        --propagationModelParameters $propagationScaleFactor \
        --intertypeIterations 20 \
        --intratypeIterations 5 \
        --virtualNodesGranularity typeAndNode \
        --saturation \
        --undirectedEdges \
        --undirectedTypeEdges \
        --outputFolder $outputFolder"
done
