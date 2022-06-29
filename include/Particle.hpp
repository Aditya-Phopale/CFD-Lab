#pragma once
#include <cmath>
#include <vector>

#include "Grid.hpp"

struct particle {
  double x;
  double y;
};

std::vector<particle> initialize_particles(int ppc);
