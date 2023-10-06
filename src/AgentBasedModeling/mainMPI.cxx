#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Assuming types is the workload to be distributed
    std::vector<std::string> types = {"type1", "type2", "type3", "type4"};

    int workloadPerProcess = types.size() / numProcesses;
    int startIdx = rank * workloadPerProcess;
    int endIdx = (rank == numProcesses - 1) ? types.size() : (rank + 1) * workloadPerProcess;

    // Each process works on its assigned workload
    for (int i = startIdx; i < endIdx; ++i) {
        // Compute for types[i]
        std::cout << "Process " << rank << " is computing for type: " << types[i] << "\n";
        // ... Perform computation for types[i]

        // Send results to master process
        MPI_Send(&types[i], 1, MPI_STRING, 0, 0, MPI_COMM_WORLD);

        // Receive results from master process

        // MPI_Recv(&types[i], 1, MPI_STRING, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Finalize();
    return 0;
}
