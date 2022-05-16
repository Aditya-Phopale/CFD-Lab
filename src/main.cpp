/*
In this file, we define our main function and ensure that a valid input
date file is provided.
*/
#include <iostream>
#include <string>

#include "Case.hpp"

int main(int argn, char **args) {
  argn = 2;
  if (argn > 1) {
    // std::string file_name{args[1]};
    std::string file_name{
        "/home/gaurav/Desktop/CFDLAB/worksheet2_new/group-gap-cfd-lab/"
        "ws2-configs-and-geometries/configs-and-geometries/ChannelWithBFS/"
        "ChannelWithBFS.dat"};
    Case problem(file_name, argn, args);
    problem.simulate();
  } else {
    std::cout << "Error: No input file is provided to fluidchen." << std::endl;
    std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat"
              << std::endl;
  }
}
