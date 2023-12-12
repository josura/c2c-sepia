# example usage: sh generatePerformanceTimesDifferentGraphs.sh <interface> graphsFolder/ outputs/times/ 

# control if parameters are passed
if [ $# -ne 3 ]; then
    echo $# arguments passed
    echo "Please pass the interface, the folder with the singleGraph and the folder where to save the outputs"
    exit 1
fi

# interface is passed as an argument
interface=$1

# get the directory of the inputs from the first argument
inputFolder=$2
outputsFolder=$3

echo "Graph folder: $inputFolder"

#get the graph file, the initial perturbation and the type interactions from the folder
graphsFolderName=$(ls $inputFolder | grep "graph")
initialPerturbationFolderName=$(ls $inputFolder | grep "node_conditions_discr")
typeInteractionsFolderName=$(ls $inputFolder | grep "interactions")

# get the full path for the graph file, the initial perturbation and the type interactions 
graphsFolder=$inputFolder$graphsFolderName
initialPerturbationFolder=$inputFolder$initialPerturbationFolderName
typeInteractionsFolder=$inputFolder$typeInteractionsFolderName

echo "Graph file: $graphsFolderName"
echo "Initial perturbation file: $initialPerturbationFolderName"
echo "Type interactions file: $typeInteractionsFolderName"


#change the number of processors from 1 to 16 (the processor cores available on my laptop)
for numProcessors in $(seq 1 16); do
    echo "Number of processors: $numProcessors"
    # create the folder where to save the outputs
    mkdir -p $outputsFolder/$numProcessors

    # run the simulation
    mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np $numProcessors ./build/c2c-sepia-MPI --graphsFilesFolder $graphsFolder \
                --initialPerturbationPerTypeFolder $initialPerturbationFolder \
                --typeInteractionFolder $typeInteractionsFolder \
                --dissipationModel scaled \
                --dissipationModelParameters 0.2 \
                --conservationModel scaled \
                --conservationModelParameters 0.2 \
                --propagationModel scaled \
                --propagationModelParameters 0.2 \
                --intertypeIterations 3 \
                --intratypeIterations 3 \
                --saturation \
                --undirectedEdges \
                --undirectedTypeEdges \
                --outputFolder $outputsFolder/$numProcessors \
                --savePerformance scripts/performanceTesting/times.tsv
done
