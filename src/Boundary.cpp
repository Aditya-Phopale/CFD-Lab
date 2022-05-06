#include "Boundary.hpp"

#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells)
    : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells,
                                     std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
  for (auto cells : _cells) {
    int i = cells->i();
    int j = cells->j();
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

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells,
                                       double wall_velocity)
    : _cells(cells) {
  _wall_velocity.insert(
      std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells,
                                       std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells),
      _wall_velocity(wall_velocity),
      _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {
  for (auto cells : _cells) {
    int i = cells->i();
    int j = cells->j();
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.u(i, j - 1);
      field.v(i, j - 1) = 0.0;
      field.p(i, j) = field.p(i, j - 1);
      field.g(i, j - 1) = field.v(i, j - 1);
      continue;
    }
    // TO check
    //   if (_cells[i]->is_border(border_position::TOP)) {
    //     field.u(posx, posy) =
    //         (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
    //         field.u(posx, posy + 1);
    //     field.v(posx, posy - 1) = 0.0;
    //     field.p(posx, posy) = field.p(posx, posy + 1);
    //     field.g(posx, posy) = field.v(posx, posy);
    //     continue;
    //   }
    //   if (_cells[i]->is_border(border_position::RIGHT)) {
    //     field.u(posx, posy) = 0.0;
    //     field.v(posx, posy) =
    //         (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
    //         field.v(posx + 1, posy);
    //     field.p(posx, posy) = field.p(posx, posy - 1);
    //     field.f(posx, posy) = field.u(posx, posy);
    //     continue;
    //   }
    //   if (_cells[i]->is_border(border_position::LEFT)) {
    //     field.u(posx, posy) = 0.0;
    //     field.v(posx, posy) =
    //         (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
    //         field.v(posx - 1, posy);
    //     field.p(posx, posy) = field.p(posx, posy - 1);
    //     field.f(posx, posy) = field.u(posx, posy);
    //     continue;
    //   }
    // }
  }
}
