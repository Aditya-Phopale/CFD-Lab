#include "Communication.hpp"

void Communication::init_parallel(int argn, char **args, int &rank, int &size) {
  MPI_Init(&argn, &args);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
}

void Communication::finalize() { MPI_Finalize(); }

void Communication::communicate() {}
