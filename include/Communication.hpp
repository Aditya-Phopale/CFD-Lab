#pragma once
#include <mpi.h>
#include <vector>

#include "Domain.hpp"
#include "Fields.hpp"
class Communication {
 public:
  static int rank;
  static int root_rank;
  static int size;
  static void init_parallel(int argn, char **args);
  static void finalize();
  static void communicate(Matrix<double> &A, Domain domain);
  static double reduce_min(double &dt);
  static double reduce_sum();
};