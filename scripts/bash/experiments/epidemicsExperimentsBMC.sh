
# control if parameters are passed, if only three arguments are passed, the number of processors is set to 4
if [ $# -ne 3 ] && [ $# -ne 4 ]; then
    echo $# arguments passed
    echo "Please pass the interface, the folder with the inputs and the folder where to save the outputs, and optionally the number of processors to use"
    exit 1
fi



# interface is passed as an argument
interface=$1

#example usage: sh scripts/bash/experiments/epidemicsExperimentsBMC.sh wlan0 scripts/R/epidemics/syntheticGraphs/graphtype/numberOfNodes outputs/epidemics/graphtype/numberOfNodes 
# get the directory of the inputs from the first argument
inputsFolder=$2
outputsFolder=$3

numberOfProcessors=4

if [ $# -eq 4 ]; then
    numberOfProcessors=$4
fi

dissipationScaleFactor=0.2
propagationScaleFactor=0.5

# from the input folder where the graphs are stored, select the inputs for every graph and echo the command to run the simulation, the name of the graphs folders go from 1 to 30
for i in {1..30} 
do
    #get the graph file, the initial perturbation and the type interactions from the folder
    graphsFolderName=$(ls $inputsFolder/$i | grep "graph") 
    initialPerturbationFolderName=$(ls $inputsFolder/$i | grep "node_conditions_discr")
    typeInteractionsFolderName=$(ls $inputsFolder/$i | grep "interactions")
    

    # get the full path for the graph file, the initial perturbation and the type interactions 
    graphsFolder=$inputsFolder/$i/$graphsFolderName
    initialPerturbationFolder=$inputsFolder/$i/$initialPerturbationFolderName
    typeInteractionsFolder=$inputsFolder/$i/$typeInteractionsFolderName
    nodesFolder=$inputsFolder/$i/communities

    echo "Graph file: $graphsFolderName"
    echo "Initial perturbation file: $initialPerturbationFolderName"
    echo "Type interactions file: $typeInteractionsFolderName"
    echo "Nodes file: $nodesFolderName"

    # get the output folder
    #outputFolder="$outputsFolder/dissipationScaleFactor${dissipationScaleFactor}_propagationScaleFactor${propagationScaleFactor}"
    outputFolder="$outputsFolder/$i"
    # run the simulation
    mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np $numberOfProcessors ./build/c2c-sepia-MPI --graphsFilesFolder $graphsFolder \
        --initialPerturbationPerTypeFolder $initialPerturbationFolder \
        --typeInteractionFolder $typeInteractionsFolder \
        --nodeDescriptionFolder $nodesFolder \
        --dissipationModel scaled \
        --dissipationModelParameters $dissipationScaleFactor \
        --propagationModel neighbors \
        --propagationModelParameters $propagationScaleFactor \
        --intertypeIterations 10 \
        --intratypeIterations 5 \
        --virtualNodesGranularity typeAndNode \
        --saturation \
        --undirectedEdges \
        --undirectedTypeEdges \
        --outputFolder $outputFolder
done