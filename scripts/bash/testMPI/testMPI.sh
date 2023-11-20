# interface is passed as an argument
interface=$1

mpirun --mca pml ob1 --mca btl tcp,self --mca btl_tcp_if_include $interface -np 2 ./build/c2c-sepia-MPI --graphsFilesFolder data/testdata/testHeterogeneousGraph/graphs \
            --initialPerturbationPerTypeFolder data/testdata/testHeterogeneousGraph/initialValuesPartialTypes \
            --typeInteractionFolder data/testdata/testHeterogeneousGraph/interactions \
            --dissipationModel scaled \
            --dissipationModelParameters 0.2 \
            --propagationModel neighbors \
            --propagationModelParameters 0.2 \
            --saturation \
            --undirectedEdges \
            --undirectedTypeEdges \
            --outputFolder outputs/testingMPI