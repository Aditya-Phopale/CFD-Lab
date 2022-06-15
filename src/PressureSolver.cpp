/*
In this file, we solve for pressure and calculate the residual to check whether
the solution is converging or not.
*/
#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid,
                  const std::vector<std::unique_ptr<Boundary>> &boundaries) {
  double dx = grid.dx();
  double dy = grid.dy();

  double coeff =
      _omega /
      (2.0 * (1.0 / (dx * dx) +
              1.0 / (dy * dy)));  // = _omega * h^2 / 4.0, if dx == dy == h

  for (auto currentCell : grid.fluid_cells()) {
    int i = currentCell->i();
    int j = currentCell->j();

    field.p(i, j) =
        (1.0 - _omega) * field.p(i, j) +
        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) -
                 field.rs(i, j));
  }

  double res = 0.0;
  double rloc = 0.0;

  // Using squared value of difference to calculate residual

  for (auto currentCell : grid.fluid_cells()) {
    int i = currentCell->i();
    int j = currentCell->j();

    double val =
        Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
    rloc += (val * val);
  }
  //std::cout<<rloc<<"\n";
  rloc = Communication::reduce_sum(rloc);
  int fluid_cells_size = grid.fluid_cells().size();
  fluid_cells_size = Communication::reduce_sum_integer(fluid_cells_size);
  {
    res = rloc / fluid_cells_size;
    res = std::sqrt(res);
  }

  return res;
}
