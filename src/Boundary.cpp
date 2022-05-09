/*
In this file, we will apply boundary conditions for the walls of our domain. Our
domain has 4 walls - 1 infinitely long lid at the top moving with a fixed
velocity and 3 other stationary walls. Based on their position they will each be
described a specific boundary condition.
*/
#include "Boundary.hpp"

#include <cmath>
#include <iostream>

/*
In the following code section, you will see 2 constructors for each boundary
type. For this worksheet, we will use the first one as we do not need wall
temperature.
*/

//   For the 3 fixed walls
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells)
    : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells,
                                     std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

/*
In the following code section, we will apply boundary conditions to the fixed
walls of our domain based on where the fluid lies with respect to that wall. For
instance, for the bottom wall; the fluid lies to the top of it. Similarly, we
implement the same for the two other fixed walls
*/
void FixedWallBoundary::apply(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::TOP)) {
      field.u(i, j) = -field.u(i, j + 1);
      field.v(i, j) = 0.0;
      field.p(i, j) = field.p(i, j + 1);
      field.g(i, j) = field.v(i, j);
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) = -field.v(i + 1, j);
      field.p(i, j) = field.p(i + 1, j);
      field.f(i, j) = field.u(i, j);
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) = -field.v(i - 1, j);
      field.p(i, j) = field.p(i - 1, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) = -field.u(i, j - 1);
      field.v(i, j) = 0.0;
      field.p(i, j) = field.p(i, j - 1);
      field.g(i, j) = field.v(i, j);
      continue;
    }
  }
}
//  For the moving wall
MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells,
                                       double wall_velocity)
    : _cells(cells) {
  _wall_velocity.insert(
      std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

/*
For the moving wall, i.e., the lid the fluid will always be below it
*/

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells,
                                       std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells),
      _wall_velocity(wall_velocity),
      _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.u(i, j - 1);
      field.v(i, j - 1) = 0.0;
      field.p(i, j) = field.p(i, j - 1);
      field.g(i, j - 1) = field.v(i, j - 1);
      continue;
    }
    if (_cells[i]->is_border(border_position::TOP)) {
      field.u(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.u(i, j + 1);
      field.v(i, j - 1) = 0.0;
      field.p(i, j) = field.p(i, j + 1);
      field.g(i, j) = field.v(i, j);
      continue;
    }
    if (_cells[i]->is_border(border_position::RIGHT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.v(i + 1, j);
      field.p(i, j) = field.p(i, j - 1);
      field.f(i, j) = field.u(i, j);
      continue;
    }
    if (_cells[i]->is_border(border_position::LEFT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.v(i - 1, j);
      field.p(i, j) = field.p(i, j - 1);
      field.f(i, j) = field.u(i, j);
      continue;
    }
  }
}
