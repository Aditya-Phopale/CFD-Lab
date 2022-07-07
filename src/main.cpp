/*
In this file, we define our main function and ensure that a valid input
date file is provided.
*/
#include <mpi.h>

#include <chrono>
#include <iostream>
#include <string>

#include "Case.hpp"
#include "Communication.hpp"

int main(int argn, char **args) {
  auto start = std::chrono::steady_clock::now();
  Communication::init_parallel(argn, args);
  argn = 2;
  if (argn > 1) {
    // std::string file_name{args[1]};
    std::string file_name{
        "/home/gaurav/Desktop/CFDLAB/project/ws2-configs-and-geometries/"
        "configs-and-geometries/ShearFlow/ShearFlowSurface.dat"};
    Case problem(file_name, argn, args);
    problem.simulate();
  } else {
    std::cout << "Error: No input file is provided to fluidchen." << std::endl;
    std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat"
              << std::endl;
  }
  Communication::finalize();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  std::cout << "\nElapsed time: " << elapsed_time.count() << "\n";
}
