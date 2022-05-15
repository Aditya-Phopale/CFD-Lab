#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
 public:
  /**
   * @brief Main method to patch the boundary conditons to given field and
   * grid
   *
   * @param[in] Field to be applied
   */
  virtual void apply(Fields &field) = 0;
  virtual ~Boundary() = default;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
 public:
  FixedWallBoundary(std::vector<Cell *> cells);
  FixedWallBoundary(std::vector<Cell *> cells,
                    std::map<int, double> wall_temperature);
  virtual ~FixedWallBoundary() = default;
  virtual void apply(Fields &field);

 private:
  std::vector<Cell *> _cells;
  std::map<int, double> _wall_temperature;
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
 public:
  MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
  MovingWallBoundary(std::vector<Cell *> cells,
                     std::map<int, double> wall_velocity,
                     std::map<int, double> wall_temperature);
  virtual ~MovingWallBoundary() = default;
  virtual void apply(Fields &field);

 private:
  std::vector<Cell *> _cells;
  std::map<int, double> _wall_velocity;
  std::map<int, double> _wall_temperature;
};

class InletBoundary : public Boundary {
 public:
  InletBoundary(std::vector<Cell *> cells, double uin, double vin);
  InletBoundary(std::vector<Cell *> cells,
                     double uin, double vin,
                     std::map<int, double> wall_temperature);
  virtual ~InletBoundary() = default;
  virtual void apply(Fields &field);

 private:
  std::vector<Cell *> _cells;
  double _uin;
  double _vin;
  std::map<int, double> _wall_temperature;
};

class OutletBoundary : public Boundary {
 public:
  OutletBoundary(std::vector<Cell *> cells);
  OutletBoundary(std::vector<Cell *> cells,
                     std::map<int, double> wall_temperature);
  virtual ~OutletBoundary() = default;
  virtual void apply(Fields &field);

 private:
  std::vector<Cell *> _cells;
  std::map<int, double> _wall_temperature;
};
