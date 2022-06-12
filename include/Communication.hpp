#pragma once
#include <mpi.h>

#include "Domain.hpp"
#include "Fields.hpp"
class Communication {
 public:
  static int rank;
  static int size;
  static void init_parallel(int argn, char **args);
  static void finalize();
  static void communicate(Fields &A, Domain domain);
  static double reduce_min();
  static double reduce_sum();
};