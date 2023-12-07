# interface is passed as an argument
interface=$1


#example usage: sh scripts/bash/experiments/100NodesEpidemicsExperiments.sh scripts/R/epidemics/syntheticGraphs/100Nodes/ outputs/100NodesEpidemics/ 
# get the directory of the inputs from the first argument
inputsFolder=$2
outputsFolder=$3

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

        mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np 4 ./build/c2c-sepia-MPI --graphsFilesFolder data/testdata/testHeterogeneousGraph/graphs \
                    --initialPerturbationPerTypeFolder data/testdata/testHeterogeneousGraph/initialValuesPartialTypes \
                    --typeInteractionFolder data/testdata/testHeterogeneousGraph/interactions \
                    --dissipationModel scaled \
                    --dissipationModelParameters 0.2 \
                    --conservationModel scaled \
                    --conservationModelParameters 0.2 \
                    --propagationModel customPropagation \
                    --propagationModelParameters 0.2 \
                    --saturation \
                    --undirectedEdges \
                    --undirectedTypeEdges \
                    --outputFolder outputs/testingMPI
    done
done