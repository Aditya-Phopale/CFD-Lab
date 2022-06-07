/*
In this file, we define our main function and ensure that a valid input
date file is provided.
*/
#include <iostream>
#include <string>
#include<mpi.h>

#include "Case.hpp"

int main(int argn, char **args) {
  int rank, size;
  MPI_Init(&argn, &args);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (argn > 1) {
    std::string file_name{args[1]};
    Case problem(file_name, argn, args, rank);
    problem.simulate();
  } else {
    std::cout << "Error: No input file is provided to fluidchen." << std::endl;
    std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat"
              << std::endl;
  }
  MPI_Finalize();
}
