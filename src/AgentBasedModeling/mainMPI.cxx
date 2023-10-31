#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char** argv) {
    // starting data for each process
    
    // these types will be used for indexing the single processes, rank 0 is the master process, rank 1 will get the first type, rank 2 the second and so on
    std::vector<std::string> types = {"type1", "type2", "type3", "type4"}; //TODO get from files
    std::map<std::string, int> typeToIndex;
    for (int i = 0; i < types.size(); ++i) {
        typeToIndex[types[i]] = i;
    }

    // initialize MPI
    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int workloadPerProcess = types.size() / numProcesses;
    int startIdx = rank * workloadPerProcess;
    int endIdx = (rank == numProcesses - 1) ? types.size() : (rank + 1) * workloadPerProcess;

    // Each process works on its assigned workload
    for (int i = startIdx; i < endIdx; ++i) {
        // Compute for types[i]
        std::cout << "Process " << rank << " is computing for type: " << types[i] << "\n";
        // ... Perform computation for types[i]

        // Send results to master process
        //MPI_Send(&types[i], 1, MPI_STRING, 0, 0, MPI_COMM_WORLD);

        // Send virtual outputs in the current process to the target process taken from the typeToIndex map
        // MPI_Send(&virtualOutputs, 1, MPI_STRING, typeToIndex[types[i]], 0, MPI_COMM_WORLD);

        // Receive results from master process

        // MPI_Recv(&types[i], 1, MPI_STRING, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Finalize();
    return 0;
}
