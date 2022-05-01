#include <iostream>
#include <string>

#include "Case.hpp"

int main(int argn, char **args) {
  argn = 2;
  if (argn > 1) {
    // std::string file_name{args[1]};
    std::string file_name{
        "/home/gaurav/Desktop/CFDLAB/worksheet1/group-gap-cfd-lab/"
        "example_cases/LidDrivenCavity/LidDrivenCavity.dat"};
    Case problem(file_name, argn, args);
    problem.simulate();
  } else {
    std::cout << "Error: No input file is provided to fluidchen." << std::endl;
    std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat"
              << std::endl;
  }
}
