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
    double temp = 0.0;
    bool energy_flag = field.energy_eq();
    if (energy_flag) {
      if (cells->type() == cell_type::FIXED_WALL3) {
        temp = _wall_temperature[cellID::fixed_wall_3];
      } else {
        temp = _wall_temperature[cellID::fixed_wall_4];
      }
    }
    if (cells->is_border(border_position::NORTHEAST)) {
      field.u(i, j) = 0.0;
      field.v(i, j) = 0.0;
      // field.p(i, j) = (field.p(i, j + 1) + field.p(i + 1, j)) / 2;
      field.g(i, j) = field.v(i, j);
      field.f(i, j) = field.u(i, j);
      if (energy_flag) {
        field.t(i, j) =
            (2.0) * (temp) - (field.t(i, j + 1) + field.t(i + 1, j)) / 2;
      }
      continue;
    }
    if (cells->is_border(border_position::NORTHWEST)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) = 0.0;
      // field.p(i, j) = (field.p(i, j + 1) + field.p(i - 1, j)) / 2;
      field.g(i, j) = field.v(i, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      if (energy_flag) {
        field.t(i, j) =
            (2.0) * (temp) - (field.t(i, j + 1) + field.t(i - 1, j)) / 2;
      }
      continue;
    }
    if (cells->is_border(border_position::SOUTHWEST)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j - 1) = 0.0;
      // field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2;
      field.g(i, j - 1) = field.v(i, j - 1);
      field.f(i - 1, j) = field.u(i - 1, j);
      if (energy_flag) {
        field.t(i, j) =
            (2.0) * (temp) - (field.t(i, j - 1) + field.t(i - 1, j)) / 2;
      }
      continue;
    }
    if (cells->is_border(border_position::SOUTHEAST)) {
      field.u(i, j) = 0.0;
      field.v(i, j - 1) = 0.0;
      // field.p(i, j) = (field.p(i, j - 1) + field.p(i + 1, j)) / 2;
      field.g(i, j - 1) = field.v(i, j - 1);
      field.f(i, j) = field.u(i, j);
      if (energy_flag) {
        field.t(i, j) =
            (2.0) * (temp) - (field.t(i, j - 1) + field.t(i + 1, j)) / 2;
      }
      continue;
    }

    if (cells->is_border(border_position::TOP)) {
      field.u(i, j) = -field.u(i, j + 1);
      field.v(i, j) = 0.0;
      // field.p(i, j) = field.p(i, j + 1);
      field.g(i, j) = field.v(i, j);
      if (energy_flag) {
        field.t(i, j) = (2.0) * (temp)-field.t(i, j + 1);
      }
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) = -field.u(i, j - 1);
      field.v(i, j - 1) = 0.0;
      // field.p(i, j) = field.p(i, j - 1);
      field.g(i, j - 1) = field.v(i, j - 1);
      if (energy_flag) {
        field.t(i, j) = (2.0) * (temp)-field.t(i, j - 1);
      }
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) = -field.v(i + 1, j);
      // field.p(i, j) = field.p(i + 1, j);
      field.f(i, j) = field.u(i, j);
      if (energy_flag) {
        field.t(i, j) = (2.0) * (temp)-field.t(i + 1, j);
      }
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) = -field.v(i - 1, j);
      // field.p(i, j) = field.p(i - 1, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      if (energy_flag) {
        field.t(i, j) = (2.0) * (temp)-field.t(i - 1, j);
      }
      continue;
    }
  }
}

void FixedWallBoundary::apply_pressure(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();

    if (cells->is_border(border_position::NORTHEAST)) {
      field.p(i, j) = (field.p(i, j + 1) + field.p(i + 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::NORTHWEST)) {
      field.p(i, j) = (field.p(i, j + 1) + field.p(i - 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::SOUTHWEST)) {
      field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2;
      continue;
    }
    if (cells->is_border(border_position::SOUTHEAST)) {
      field.p(i, j) = (field.p(i, j - 1) + field.p(i + 1, j)) / 2;
      continue;
    }

    if (cells->is_border(border_position::TOP)) {
      field.p(i, j) = field.p(i, j + 1);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.p(i, j) = field.p(i, j - 1);
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.p(i, j) = field.p(i + 1, j);
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.p(i, j) = field.p(i - 1, j);
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
      // field.p(i, j) = field.p(i, j - 1);
      field.g(i, j - 1) = field.v(i, j - 1);
      continue;
    }
    if (_cells[i]->is_border(border_position::TOP)) {
      field.u(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.u(i, j + 1);
      field.v(i, j) = 0.0;
      // field.p(i, j) = field.p(i, j + 1);
      field.g(i, j) = field.v(i, j);
      continue;
    }
    if (_cells[i]->is_border(border_position::RIGHT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.v(i + 1, j);
      // field.p(i, j) = field.p(i + 1, j);
      field.f(i, j) = field.u(i, j);
      continue;
    }
    if (_cells[i]->is_border(border_position::LEFT)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) =
          (2.0) * (_wall_velocity[LidDrivenCavity::moving_wall_id]) -
          field.v(i - 1, j);
      // field.p(i, j) = field.p(i - 1, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      continue;
    }
  }
}

void MovingWallBoundary::apply_pressure(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::BOTTOM)) {
      field.p(i, j) = field.p(i, j - 1);
      continue;
    }
    if (_cells[i]->is_border(border_position::TOP)) {
      field.p(i, j) = field.p(i, j + 1);
      continue;
    }
    if (_cells[i]->is_border(border_position::RIGHT)) {
      field.p(i, j) = field.p(i + 1, j);
      continue;
    }
    if (_cells[i]->is_border(border_position::LEFT)) {
      field.p(i, j) = field.p(i - 1, j);
      continue;
    }
  }
}

OutletBoundary::OutletBoundary(std::vector<Cell *> cells) : _cells(cells) {}

OutletBoundary::OutletBoundary(std::vector<Cell *> cells,
                               std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void OutletBoundary::apply(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::LEFT)) {
      field.u(i - 1, j) = field.u(i - 2, j);
      field.v(i, j) = field.v(i - 1, j);
      // field.p(i, j) = -field.p(i - 1, j);
      field.f(i - 1, j) = field.u(i - 1, j);  //- 2*field.dt()*field.p(i,j);
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) = -field.v(i + 1, j);
      // field.p(i, j) = field.p(i + 1, j);
      field.f(i, j) = field.u(i, j);
      continue;
    }
    if (cells->is_border(border_position::TOP)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) = -field.v(i - 1, j);
      // field.p(i, j) = field.p(i - 1, j);
      field.g(i, j) = field.v(i, j);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) = -field.u(i, j - 1);
      field.v(i, j) = 0.0;
      // field.p(i, j) = field.p(i, j - 1);
      field.g(i, j) = field.v(i, j);
      continue;
    }
  }
}

void OutletBoundary::apply_pressure(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::LEFT)) {
      field.p(i, j) = -field.p(i - 1, j);
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.p(i, j) = field.p(i + 1, j);
      continue;
    }
    if (cells->is_border(border_position::TOP)) {
      field.p(i, j) = field.p(i - 1, j);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.p(i, j) = field.p(i, j - 1);
      continue;
    }
  }
}

InletBoundary::InletBoundary(std::vector<Cell *> cells, double uin, double vin)
    : _cells(cells), _uin(uin), _vin(vin) {}

InletBoundary::InletBoundary(std::vector<Cell *> cells, double uin, double vin,
                             std::map<int, double> wall_temperature)
    : _cells(cells),
      _uin(uin),
      _vin(vin),
      _wall_temperature(wall_temperature) {}

void InletBoundary::apply(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::RIGHT)) {
      field.u(i, j) = _uin;
      field.v(i, j) = (2.0 * _vin) - field.v(i + 1, j);
      field.p(i, j) = field.p(i + 1, j);
      field.f(i, j) = field.u(i, j);
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.u(i - 1, j) = _uin;
      field.v(i, j) = (2.0 * _vin) - field.v(i - 1, j);
      field.p(i, j) = field.p(i - 1, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      continue;
    }
    if (cells->is_border(border_position::TOP)) {
      field.u(i, j) = (2.0 * _uin) - field.u(i, j + 1);
      field.v(i, j) = _vin;
      field.p(i, j) = field.p(i, j + 1);
      field.g(i, j) = field.v(i, j);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) = (2.0 * _uin) - field.u(i, j - 1);
      field.v(i, j - 1) = _vin;
      field.p(i, j) = field.p(i, j - 1);
      field.g(i, j - 1) = field.v(i, j - 1);
      continue;
    }
  }
}

void InletBoundary::apply_pressure(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::RIGHT)) {
      field.p(i, j) = field.p(i + 1, j);
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.p(i, j) = field.p(i - 1, j);
      continue;
    }
    if (cells->is_border(border_position::TOP)) {
      field.p(i, j) = field.p(i, j + 1);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.p(i, j) = field.p(i, j - 1);
      continue;
    }
  }
}

AdiabaticBoundary::AdiabaticBoundary(std::vector<Cell *> cells)
    : _cells(cells) {}

void AdiabaticBoundary::apply(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::NORTHEAST)) {
      field.u(i, j) = 0.0;
      field.v(i, j) = 0.0;
      // field.p(i, j) = (field.p(i, j + 1) + field.p(i + 1, j)) / 2;
      field.g(i, j) = field.v(i, j);
      field.f(i, j) = field.u(i, j);
      field.t(i, j) = (field.t(i, j + 1) + field.t(i + 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::NORTHWEST)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) = 0.0;
      // field.p(i, j) = (field.p(i, j + 1) + field.p(i - 1, j)) / 2;
      field.g(i, j) = field.v(i, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      field.t(i, j) = (field.t(i, j + 1) + field.t(i - 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::SOUTHWEST)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j - 1) = 0.0;
      // field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2;
      field.g(i, j - 1) = field.v(i, j - 1);
      field.f(i - 1, j) = field.u(i - 1, j);
      field.t(i, j) = (field.t(i, j - 1) + field.t(i - 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::SOUTHEAST)) {
      field.u(i, j) = 0.0;
      field.v(i, j - 1) = 0.0;
      // field.p(i, j) = (field.p(i, j - 1) + field.p(i + 1, j)) / 2;
      field.g(i, j - 1) = field.v(i, j - 1);
      field.f(i, j) = field.u(i, j);
      field.t(i, j) = (field.t(i, j - 1) + field.t(i + 1, j)) / 2;
      continue;
    }

    if (cells->is_border(border_position::TOP)) {
      field.u(i, j) = -field.u(i, j + 1);
      field.v(i, j) = 0.0;
      // field.p(i, j) = field.p(i, j + 1);
      field.g(i, j) = field.v(i, j);
      field.t(i, j) = field.t(i, j + 1);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.u(i, j) = -field.u(i, j - 1);
      field.v(i, j - 1) = 0.0;
      // field.p(i, j) = field.p(i, j - 1);
      field.g(i, j - 1) = field.v(i, j - 1);
      field.t(i, j) = field.t(i, j - 1);
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.u(i, j) = 0.0;
      field.v(i, j) = -field.v(i + 1, j);
      // field.p(i, j) = field.p(i + 1, j);
      field.f(i, j) = field.u(i, j);
      field.t(i, j) = field.t(i + 1, j);
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.u(i - 1, j) = 0.0;
      field.v(i, j) = -field.v(i - 1, j);
      // field.p(i, j) = field.p(i - 1, j);
      field.f(i - 1, j) = field.u(i - 1, j);
      field.t(i, j) = field.t(i - 1, j);
      continue;
    }
  }
}

void AdiabaticBoundary::apply_pressure(Fields &field) {
  int i, j;
  for (auto cells : _cells) {
    i = cells->i();
    j = cells->j();
    if (cells->is_border(border_position::NORTHEAST)) {
      field.p(i, j) = (field.p(i, j + 1) + field.p(i + 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::NORTHWEST)) {
      field.p(i, j) = (field.p(i, j + 1) + field.p(i - 1, j)) / 2;
      continue;
    }
    if (cells->is_border(border_position::SOUTHWEST)) {
      field.p(i, j) = (field.p(i - 1, j) + field.p(i, j - 1)) / 2;
      continue;
    }
    if (cells->is_border(border_position::SOUTHEAST)) {
      field.p(i, j) = (field.p(i, j - 1) + field.p(i + 1, j)) / 2;
      continue;
    }

    if (cells->is_border(border_position::TOP)) {
      field.p(i, j) = field.p(i, j + 1);
      continue;
    }
    if (cells->is_border(border_position::BOTTOM)) {
      field.p(i, j) = field.p(i, j - 1);
      continue;
    }
    if (cells->is_border(border_position::RIGHT)) {
      field.p(i, j) = field.p(i + 1, j);
      continue;
    }
    if (cells->is_border(border_position::LEFT)) {
      field.p(i, j) = field.p(i - 1, j);
      continue;
    }
  }
}