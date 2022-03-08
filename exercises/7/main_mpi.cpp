#include <iostream>
#include <mpi.h>

int main(int argc, char* argv[]) {
    int s, r;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);

    // Check we have exactly 2 processes for this demo
    if (s != 2) {
        std::cout << "Only 2 processes allowed with this program!" << std::endl;
        MPI_Finalize();
        return 1;
    }

    // Generate a random number
    srand(time(0) + r*100);
    double x = (double)(rand())/RAND_MAX;
    double y = 0.0f;

    std::cout << "Performing exchange... " << std::endl;
    std::cout << "Rank " << r << " sends " << x << std::endl;

    if (r == 0) {
        // Only execute on process 0
        // MPI_Send (&var, count, type, dest, tag, comm)
        MPI_Send(&x, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);

        // MPI_Recv (&var, count, type, src, tag, comm, status)
        MPI_Recv(&y, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        // Execute on all other processes (1)
        // Recv from 0
        MPI_Recv(&y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Send to 0s
        MPI_Send(&x, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }

    std::cout << "Complete!" << std::endl;
    std::cout << "Rank: " << r << " receives: " << y << std::endl;

    MPI_Finalize();
    return 0;
}