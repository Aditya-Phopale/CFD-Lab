/*
In this file, we calculate the velocity, and the time by which
the current timestep will be advanced to the next one
*/
#include "Fields.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

Fields::Fields(double nu, double alpha, double beta, double dt, double tau,
               int imax, int jmax, double UI, double VI, double PI, double TI,
               double GX, double GY, bool energy_eq)
    : _nu(nu),
      _dt(dt),
      _tau(tau),
      _alpha(alpha),
      _beta(beta),
      _gx(GX),
      _gy(GY),
      _energy_eq(energy_eq) {
  _U = Matrix<double>(imax + 2, jmax + 2,
                      UI);  // Matrix for velocity along the X-direction
  _V = Matrix<double>(imax + 2, jmax + 2,
                      VI);  // Matrix for velocity along the Y-direction
  _P =
      Matrix<double>(imax + 2, jmax + 2,
                     PI);  // Matrix for the pressure values in the cell centers
  if (energy_eq) {
    _T = Matrix<double>(
        imax + 2, jmax + 2,
        TI);  // Matrix for the temperature values in the cell centers
  }
  _F = Matrix<double>(imax + 2, jmax + 2,
                      0.0);  // Matrix containing discretized differential data
  _G =
      Matrix<double>(imax + 2, jmax + 2,
                     0.0);  // of the momentum equation for U and V respectively
  _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);

  _VOF = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

// Calculating differential data for Explicit Euler Scheme

void Fields::calculate_fluxes(Grid &grid) {
  int i, j;

  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();

    if (cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
      _F(i, j) = _U(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                   Discretization::convection_u(_U, _V, i, j));
      if (_energy_eq) {
        _F(i, j) =
            _F(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i + 1, j)) * _gx;
      }
    }
    if (cell->neighbour(border_position::TOP)->type() == cell_type::FLUID) {
      _G(i, j) = _V(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                   Discretization::convection_v(_U, _V, i, j));

      if (_energy_eq) {
        _G(i, j) =
            _G(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i, j + 1)) * _gy;
      }
    }
  }
}

void Fields::calculate_temperature(Grid &grid) {
  int i, j;
  Matrix<double> T_old = _T;
  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();

    _T(i, j) =
        T_old(i, j) + _dt * (_alpha * Discretization::diffusion(T_old, i, j) -
                             Discretization::convection_T(_U, _V, T_old, i, j));
  }
}

void Fields::calculate_rs(Grid &grid) {
  int i, j;
  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();
    double term1 = (_F(i, j) - _F(i - 1, j)) / grid.dx();
    double term2 = (_G(i, j) - _G(i, j - 1)) / grid.dy();
    _RS(i, j) = (term1 + term2) / _dt;
  }
}

// Applying explicit Euler method

void Fields::calculate_velocities(Grid &grid) {
  int i, j;
  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();

    if (cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
      _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();
    }
    if (cell->neighbour(border_position::TOP)->type() == cell_type::FLUID) {
      _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();
    }
  }
}

void Fields::calculate_vof(Grid &grid) {
  int i, j;

  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();
    // std::vector<double> fluxes(4);  // top,right,bottom,left

    // if (_U(i, j) > 0)
    //   fluxes.at(1) = _VOF(i, j) * _U(i, j);
    // else
    //   fluxes.at(1) = _VOF(i + 1, j) * _U(i, j);

    // if (_U(i - 1, j) > 0)
    //   fluxes.at(3) = _VOF(i - 1, j) * _U(i - 1, j);
    // else
    //   fluxes.at(3) = _VOF(i, j) * _U(i - 1, j);

    // if (_V(i, j) > 0)
    //   fluxes.at(0) = _VOF(i, j) * _V(i, j);
    // else
    //   fluxes.at(0) = _VOF(i, j + 1) * _V(i, j);

    // if (_V(i, j - 1) > 0)
    //   fluxes.at(2) = _VOF(i, j - 1) * _V(i, j - 1);
    // else
    //   fluxes.at(2) = _VOF(i, j) * _V(i, j - 1);

    // if (_VOF(i, j) > 0) {
    //   double n_x = (1 / grid.dx()) * (_VOF(i + 1, j + 1) + 2 * _VOF(i + 1, j)
    //   +
    //                                   _VOF(i + 1, j - 1) - _VOF(i - 1, j + 1)
    //                                   - 2 * _VOF(i - 1, j) - _VOF(i - 1, j -
    //                                   1));

    //   double n_y = (1 / grid.dy()) * (_VOF(i + 1, j + 1) + 2 * _VOF(i, j + 1)
    //   +
    //                                   _VOF(i - 1, j + 1) - _VOF(i + 1, j - 1)
    //                                   - 2 * _VOF(i, j - 1) - _VOF(i - 1, j -
    //                                   1));

    //   double angle_beta = atan(-n_x / n_y);
    //   double angle_alpha = atan(grid.dx() * tan(angle_beta) / grid.dy());

    //   int cases = set_case(angle_alpha, _VOF(i, j));

    //   std::vector<double> side_fractions =
    //       set_side_vf(cases, angle_alpha, _VOF(i, j));

    //   double u_left = _U(i - 1, j);
    //   double u_right = _U(i, j);
    //   double u_top = _V(i, j);
    //   double u_bottom = _V(i, j - 1);

    //   calculate_vf_fluxes(fluxes, u_left, u_right, u_top, u_bottom, _VOF(i,
    //   j),
    //                       side_fractions, angle_beta, cases, grid);

    _VOF(i, j) =
        _VOF(i, j) - _dt * (Discretization::convection_T(_U, _V, _VOF, i, j));
  }

  // _VOF(i, j) = _VOF(i, j) - _dt * ((fluxes.at(1) - fluxes.at(3)) /
  // grid.dx() +
  //                                  (fluxes.at(0) - fluxes.at(2)) /
  //                                  grid.dy());

  // std::cout << _VOF(50, 50) << "\n";
}

void Fields::initialise_vof(Grid &grid) {
  double radius = 5;
  int i, j;
  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();

    if ((i - grid.domain().imax / 2) * (i - grid.domain().imax / 2) +
            (j - grid.domain().jmax / 2) * (j - grid.domain().jmax / 2) <
        radius * radius) {
      _VOF(i, j) = 1.0;
    }
  }
}

// Calculating dt based on CFL conditions

double Fields::calculate_dt(Grid &grid) {
  double CFLu = 0.0;
  double CFLv = 0.0;
  double CFLnu = 0.0;
  double CFLt = 0.0;

  double dx2 = grid.dx() * grid.dx();
  double dy2 = grid.dy() * grid.dy();

  double u_max = 0;
  double v_max = 0;

  for (auto cell : grid.fluid_cells()) {
    int i = cell->i();
    int j = cell->j();
    u_max = std::max(u_max, fabs(_U(i, j)));
    v_max = std::max(v_max, fabs(_V(i, j)));
  }

  u_max = Communication::reduce_max(u_max);
  v_max = Communication::reduce_max(v_max);

  CFLu = grid.dx() / u_max;
  CFLv = grid.dy() / v_max;
  CFLnu = (0.5 / _nu) * (1.0 / (1.0 / dx2 + 1.0 / dy2));
  _dt = std::min({
      CFLnu,
      CFLu,
      CFLv,
  });

  if (_energy_eq) {
    CFLt = (0.5 / _alpha) * (1.0 / (1.0 / dx2 + 1.0 / dy2));
    _dt = std::min({_dt, CFLt});
  }
  _dt = _tau * _dt;

  return _dt;
}

// Functions to return the corresponding values

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }
double &Fields::vof(int i, int j) { return _VOF(i, j); }
bool Fields::energy_eq() { return _energy_eq; }

Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::t_matrix() { return _T; }
Matrix<double> &Fields::g_matrix() { return _G; }
Matrix<double> &Fields::rs_matrix() { return _RS; }
Matrix<double> &Fields::vof_matrix() { return _VOF; }

double Fields::dt() const { return _dt; }

int Fields::set_case(double angle_alpha, double vof) {
  if (angle_alpha < M_PI / 4) {
    if (vof <= 0.5 * tan(angle_alpha))
      return 1;
    else if (vof <= 1 - 0.5 * tan(angle_alpha))
      return 2;
    else
      return 4;
  } else {
    if (vof <= 0.5 / tan(angle_alpha))
      return 1;
    else if (vof <= 1 - 0.5 / tan(angle_alpha))
      return 3;
    else
      return 4;
  }
}

std::vector<double> Fields::set_side_vf(int cases, double angle_alpha,
                                        double vf) {
  std::vector<double> side_fractions(4);  // top,right,bottom,left
  if (cases == 1) {
    side_fractions.at(0) = 0;
    side_fractions.at(1) = std::sqrt(2 * vf * tan(angle_alpha));
    side_fractions.at(2) = std::sqrt(2 * vf / tan(angle_alpha));
    side_fractions.at(3) = 0;
  } else if (cases == 2) {
    side_fractions.at(0) = 0;
    side_fractions.at(1) = vf + 0.5 * tan(angle_alpha);
    side_fractions.at(2) = 1;
    side_fractions.at(3) = vf - 0.5 * tan(angle_alpha);
  } else if (cases == 3) {
    side_fractions.at(0) = vf - 0.5 / tan(angle_alpha);
    side_fractions.at(1) = 1;
    side_fractions.at(2) = vf + 0.5 / tan(angle_alpha);
    side_fractions.at(3) = 0;
  } else {
    side_fractions.at(0) = 1 - std::sqrt(2 * vf / tan(angle_alpha));
    side_fractions.at(1) = 1;
    side_fractions.at(2) = 1;
    side_fractions.at(3) = 1 - std::sqrt(2 * vf * tan(angle_alpha));
  }

  return side_fractions;
}

void Fields::calculate_vf_fluxes(std::vector<double> &fluxes, double u_left,
                                 double u_right, double u_top, double u_bottom,
                                 double vf, std::vector<double> side_fractions,
                                 double angle_beta, int cases, Grid &grid) {
  if (cases == 1) {
    if (u_top > 0) {
      if (u_top * _dt >= (1 - side_fractions.at(1)) * grid.dy())
        fluxes.at(0) = 0;
      else
        fluxes.at(0) =
            0.5 *
            pow((u_top * _dt - (1 - side_fractions.at(1)) * grid.dy()), 2) /
            tan(angle_beta);
    }
    if (u_right > 0) {
      if (u_right * _dt >= side_fractions.at(2) * grid.dy())
        fluxes.at(1) = vf * grid.dx() * grid.dy();
      else
        fluxes.at(1) = 0.5 * u_right * _dt *
                       (2 - u_right * _dt * grid.dx() / side_fractions.at(2)) *
                       side_fractions.at(1) * grid.dy();
    }
    if (u_bottom < 0) {
      if (u_bottom * _dt >= side_fractions.at(1) * grid.dy())
        fluxes.at(2) = -vf * grid.dx() * grid.dy();
      else
        fluxes.at(2) = -0.5 * u_bottom * _dt *
                       (2 - u_bottom * _dt * grid.dy() / side_fractions.at(1)) *
                       side_fractions.at(2) * grid.dx();
    }
    if (u_left < 0) {
      if (u_left * _dt <= (1 - side_fractions.at(2)) * grid.dx())
        fluxes.at(3) = 0;
      else
        fluxes.at(3) =
            -0.5 *
            pow((u_left * _dt - (1 - side_fractions.at(3)) * grid.dx()), 2) *
            tan(angle_beta);
    }
  }

  if (cases == 2) {
    if (u_top > 0) {
      if (u_top * _dt >= (1 - side_fractions.at(1)) * grid.dy())
        fluxes.at(0) = 0;
      else if (u_top * _dt <= (1 - side_fractions.at(3)) * grid.dy())
        fluxes.at(0) =
            0.5 *
            pow((u_top * _dt - (1 - side_fractions.at(1)) * grid.dy()), 2) /
            tan(angle_beta);
      else
        fluxes.at(0) = u_right * _dt *
                       (side_fractions.at(1) * grid.dy() -
                        0.5 * u_right * _dt * tan(angle_beta));
    }
    if (u_right > 0) {
      fluxes.at(1) = u_right * _dt *
                     (side_fractions.at(1) * grid.dy() -
                      0.5 * u_right * _dt * tan(angle_beta));
    }
    if (u_bottom < 0) {
      if (u_bottom * _dt >= side_fractions.at(3) * grid.dy())
        fluxes.at(2) = -u_bottom * _dt * grid.dx();
      else if (u_bottom * _dt <= side_fractions.at(1) * grid.dy())
        fluxes.at(2) =
            -u_bottom * _dt * grid.dx() -
            0.5 * pow((u_bottom * _dt - side_fractions.at(3) * grid.dy()), 2) /
                tan(angle_beta);
      else
        fluxes.at(2) = -vf * grid.dx() * grid.dy();
    }
    if (u_left < 0) {
      fluxes.at(3) = -u_left * _dt *
                     (side_fractions.at(3) * grid.dy() +
                      0.5 * u_left * _dt * tan(angle_beta));
    }
  }
  if (cases == 3) {
    if (u_top > 0) {
      fluxes.at(0) = u_top * _dt *
                     (side_fractions.at(0) * grid.dx() +
                      0.5 * u_top * _dt / tan(angle_beta));
    }
    if (u_right > 0) {
      if (u_right * _dt >= side_fractions.at(0) * grid.dx())
        fluxes.at(1) = u_right * _dt * grid.dy();
      else if (u_right * _dt <= side_fractions.at(2) * grid.dx())
        fluxes.at(1) =
            u_right * _dt * grid.dy() -
            0.5 * pow((u_right * _dt - side_fractions.at(0) * grid.dx()), 2) *
                tan(angle_beta);
      else
        fluxes.at(1) = vf * grid.dx() * grid.dy();
    }
    if (u_bottom < 0) {
      fluxes.at(2) = -u_bottom * _dt *
                     (side_fractions.at(2) * grid.dx() -
                      0.5 * u_bottom * _dt / tan(angle_beta));
    }
    if (u_left < 0) {
      if (u_left * _dt <= side_fractions.at(2) * grid.dx())
        fluxes.at(3) = 0;
      else if (u_left * _dt <= side_fractions.at(1) * grid.dx())
        fluxes.at(3) =
            -0.5 *
            pow((u_left * _dt - (1 - side_fractions.at(3)) * grid.dx()), 2) *
            tan(angle_beta);
      else
        fluxes.at(3) =
            -u_left * _dt * grid.dy() - (1 - vf) * grid.dx() * grid.dy();
    }
  }
  if (cases == 4) {
    if (u_top > 0) {
      if (u_top * _dt >= (1 - side_fractions.at(3) * grid.dy()))
        fluxes.at(0) =
            u_top * _dt * grid.dx() - (1 - vf) * grid.dx() * grid.dy();
      else
        fluxes.at(0) =
            u_top * _dt *
            (side_fractions.at(0) * _dt + 0.5 * u_top * _dt / tan(angle_beta));
    }
    if (u_right > 0) {
      if (u_right * _dt >= side_fractions.at(0) * grid.dx())
        fluxes.at(1) = u_right * _dt * grid.dy();
      else
        fluxes.at(1) =
            u_right * _dt * grid.dy() -
            0.5 * tan(angle_beta) *
                pow((u_right * _dt - side_fractions.at(0) * grid.dx()), 2);
    }
    if (u_bottom < 0) {
      if (u_bottom * _dt <= side_fractions.at(3) * grid.dy())
        fluxes.at(2) = -u_bottom * _dt * grid.dx();
      else
        fluxes.at(2) =
            -u_bottom * _dt * grid.dx() -
            0.5 * pow((u_bottom * _dt - side_fractions.at(3) * grid.dy()), 2) /
                tan(angle_beta);
    }
    if (u_left < 0) {
      if (u_left * _dt <= (1 - side_fractions.at(0)) * grid.dx())
        fluxes.at(3) =
            -u_left * _dt * grid.dy() - (1 - vf) * grid.dx() * grid.dy();
      else
        fluxes.at(3) = -u_left * _dt *
                       (side_fractions.at(3) * grid.dy() +
                        0.5 * u_left * _dt * tan(angle_beta));
    }
  }
}
