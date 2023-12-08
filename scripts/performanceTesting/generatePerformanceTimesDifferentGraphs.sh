# example usage: sh generatePerformanceTimesDifferentGraphs.sh <interface> graphsFolder/ outputs/times/ 

# control if parameters are passed
if [ $# -ne 3 ]; then
    echo "Please pass the interface, the folder with the inputs and the folder where to save the outputs"
    exit 1
fi

# interface is passed as an argument
interface=$1

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



# loop through the files (the name of the file should be <numnodes>Nodes.tsv), the number of nodes start from 1000 and increase by 1000 until 30000
for graphFile in $(ls $graphsFolder | grep "Nodes.tsv"); do
    echo "Graph file: $graphFile"
    # get the number of nodes from the file name
    numNodes=$(echo $graphFile | cut -d'N' -f 1)
    echo "$numNodes\t" >> $outputsFolder/times.txt
    # run the simulation
    echo "mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np 4 ./build/c2c-sepia-MPI \ 
                --graphsFilesFolder $graphsFolder \
                --initialPerturbationPerTypeFolder $initialPerturbationFolder \
                --typeInteractionFolder $typeInteractionsFolder \
                --dissipationModel scaled \
                --dissipationModelParameters 0.2 \
                --conservationModel scaled \
                --conservationModelParameters 0.2 \
                --propagationModel scaled \
                --propagationModelParameters 0.2 \
                --saturation \
                --undirectedEdges \
                --undirectedTypeEdges \
                --outputFolder $outputFolder \
                --savePerformance scripts/performanceTesting/timesDifferentGraphs.tsv"

done