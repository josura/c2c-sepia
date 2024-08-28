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
interface="wlan0"

mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np $processors ./build/c2c-sepia-MPI --graphsFilesFolder $inputsFolder/graphs \
    --initialPerturbationPerTypeFolder $inputsFolder/node_conditions_discr \
    --typeInteractionFolder $inputsFolder/interactions \
    --nodeDescriptionFolder $inputsFolder/communities \
    --dissipationModel scaled \
    --dissipationModelParameters 0.5 \
    --propagationModel neighbors \
    --propagationModelParameters 0.5 \
    --intertypeIterations 20 \
    --intratypeIterations 5 \
    --timestep 5 \
    --virtualNodesGranularity typeAndNode \
    --saturation \
    --savePerformance /home/josura/Projects/ccc/c2c-sepia/scripts/bash/cluster/performanceAnalysis/localPerformanceTimes.tsv \
    --outputFolder $outputFolder
