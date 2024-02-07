# interface is passed as an argument
interface=$1

# control if the interface is passed as an argument
if [ -z "$interface" ]
then
    echo "Interface not passed as argument"
    exit 1
fi

mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np 2 ./build/c2c-sepia-MPI --graphsFilesFolder data/testdata/testHeterogeneousGraph/graphs \
            --initialPerturbationPerTypeFolder data/testdata/testHeterogeneousTemporalGraph/initialValuesPartialTypes \
            --typeInteractionFolder data/testdata/testHeterogeneousTemporalGraph/interactions \
            --dissipationModel scaled \
            --dissipationModelParameters 0.2 \
            --conservationModel scaled \
            --conservationModelParameters 0.2 \
            --propagationModel customPropagation \
            --propagationModelParameters 0.2 \
            --virtualNodesGranularity typeAndNode \
            --saturation \
            --outputFolder outputs/testingMPIgranularDifferentProcessors
