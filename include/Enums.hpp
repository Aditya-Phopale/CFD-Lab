#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 4: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 4;
const double wall_velocity = 1.0;
}  // namespace LidDrivenCavity

namespace cellID {
const int fluid = 0;
const int inflow = 1;
const int outflow = 2;
const int fixed_wall_3 = 3;
const int fixed_wall_4 = 4;
const int fixed_wall_5 = 5;
const int empty = 7;
}  // namespace cellID

enum class border_position {
  TOP,
  BOTTOM,
  LEFT,
  RIGHT,
  NORTHWEST,
  SOUTHEAST,
  NORTHEAST,
  SOUTHWEST,
  NORTHSOUTH,
  EASTWEST,
  NORTHEASTWEST,
  NORTHEASTSOUTH,
  NORTHWESTSOUTH,
  EASTWESTSOUTH,
  NORTHEASTWESTSOUTH
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
}  // namespace border

enum class cell_type {
  FLUID,
  EMPTY,
  SURFACE,
  FIXED_WALL3,
  FIXED_WALL4,
  MOVING_WALL,
  INLET,
  OUTLET,
  ADIABATIC_WALL,
  DEFAULT
};
