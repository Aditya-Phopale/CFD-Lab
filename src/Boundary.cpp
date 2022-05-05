#include "Boundary.hpp"

#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells)
    : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells,
                                     std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {
  int posx, posy;
  for (int i = 0; i < _cells.size(); i++) {
    posx = _cells[i]->i();
    posy = _cells[i]->j();
    if (_cells[i]->is_border(border_position::TOP)) {
      field.u(posx, posy) = -field.u(posx, posy + 1);
      field.v(posx, posy) = 0.0;
      field.p(posx, posy) = field.p(posx, posy + 1);
      field.g(posx, posy) = field.v(posx, posy);
      continue;
    }
    if (_cells[i]->is_border(border_position::RIGHT)) {
      field.u(posx, posy) = 0.0;
      field.v(posx, posy) = -field.v(posx + 1, posy);
      field.p(posx, posy) = field.p(posx + 1, posy);
      field.f(posx, posy) = field.u(posx, posy);
      continue;
    }
    if (_cells[i]->is_border(border_position::LEFT)) {
      field.u(posx - 1, posy) = 0.0;
      field.v(posx, posy) = -field.v(posx - 1, posy);
      field.p(posx, posy) = field.p(posx - 1, posy);
      field.f(posx - 1, posy) = field.u(posx - 1, posy);
      continue;
    }
    if (_cells[i]->is_border(border_position::BOTTOM)) {
      field.u(posx, posy) = -field.u(posx, posy - 1);
      field.v(posx, posy) = 0.0;
      field.p(posx, posy) = field.p(posx, posy - 1);
      field.g(posx, posy) = field.v(posx, posy);
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
  int posx, posy;
  for (int i = 0; i < _cells.size(); i++) {
    posx = _cells[i]->i();
    posy = _cells[i]->j();
    if (_cells[i]->is_border(border_position::BOTTOM)) {
      field.u(posx, posy) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.u(posx, posy - 1);
      field.v(posx, posy - 1) = 0.0;
      field.p(posx, posy) = field.p(posx, posy - 1);
      field.g(posx, posy - 1) = field.v(posx, posy - 1);
      continue;
    }
    // TO check
    if (_cells[i]->is_border(border_position::TOP)) {
      field.u(posx, posy) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.u(posx, posy + 1);
      field.v(posx, posy - 1) = 0.0;
      field.p(posx, posy) = field.p(posx, posy + 1);
      field.g(posx, posy) = field.v(posx, posy);
      continue;
    }
    if (_cells[i]->is_border(border_position::RIGHT)) {
      field.u(posx, posy) = 0.0;
      field.v(posx, posy) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.v(posx + 1, posy);
      field.p(posx, posy) = field.p(posx, posy - 1);
      field.f(posx, posy) = field.u(posx, posy);
      continue;
    }
    if (_cells[i]->is_border(border_position::LEFT)) {
      field.u(posx, posy) = 0.0;
      field.v(posx, posy) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.v(posx - 1, posy);
      field.p(posx, posy) = field.p(posx, posy - 1);
      field.f(posx, posy) = field.u(posx, posy);
      continue;
    }
  }
}
