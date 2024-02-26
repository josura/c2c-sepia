
# control if parameters are passed
if [ $# -ne 3 ]; then
    echo $# arguments passed
    echo "Please pass the interface, the folder with the inputs and the folder where to save the outputs"
    exit 1
fi

# interface is passed as an argument
interface=$1

#example usage: sh scripts/bash/experiments/epidemicsExperimentsBMC.sh wlan0 scripts/R/epidemics/syntheticGraphs/ outputs/epidemics/ 
# get the directory of the inputs from the first argument
inputsFolder=$1
outputsFolder=$2

#get the graph file, the initial perturbation and the type interactions from the folder
graphsFolderName=$(ls $inputsFolder | grep "graph")
initialPerturbationFolderName=$(ls $inputsFolder | grep "node_conditions_discr")
typeInteractionsFolderName=$(ls $inputsFolder | grep "interactions")

# get the full path for the graph file, the initial perturbation and the type interactions 
graphsFolder=$inputsFolder$graphsFolderName
initialPerturbationFolder=$inputsFolder$initialPerturbationFolderName
typeInteractionsFolder=$inputsFolder$typeInteractionsFolderName

echo "Graph file: $graphsFolderName"
echo "Initial perturbation file: $initialPerturbationFolderName"
echo "Type interactions file: $typeInteractionsFolderName"

listDissipationScaleFactors=(0 0.2 0.4 0.6)
listPropagationScaleFactors=(0.5 1.0 2.0 4.0)

# loop through the permutation of scale factors (propagation vs dissipation)
for dissipationScaleFactor in ${listDissipationScaleFactors[@]}; do
    for propagationScaleFactor in ${listPropagationScaleFactors[@]}; do
        echo "Dissipation scale factor: $dissipationScaleFactor"
        echo "Propagation scale factor: $propagationScaleFactor"
        # get the output folder
        outputFolder="$outputsFolder/dissipationScaleFactor${dissipationScaleFactor}_propagationScaleFactor${propagationScaleFactor}"
        # run the simulation
        mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np 4 ./build/c2c-sepia-MPI --graphsFilesFolder $graphsFolder \
            --initialPerturbationPerTypeFolder $initialPerturbationFolder \
            --typeInteractionFolder $typeInteractionsFolder \
            --dissipationModel scaled \
            --dissipationModelParameters $dissipationScaleFactor \
            --propagationModel neighbors \
            --propagationModelParameters $propagationScaleFactor \
            --saturation \
            --undirectedEdges \
            --undirectedTypeEdges \
            --outputFolder $outputFolder
    done
done
