/*
In this file, we calculate the velocity, and the time by which
the current timestep will be advanced to the next one
*/
#include "Fields.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

Fields::Fields(double nu, double Re, double alpha, double beta, double dt,
               double tau, int imax, int jmax, double UI, double VI, double PI,
               double TI, double GX, double GY, bool energy_eq)
    : _nu(nu),
      _Re(Re),
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
}

// Calculating differential data for Explicit Euler Scheme

void Fields::calculate_fluxes(Grid &grid) {
  int i, j;

  for (auto cell : grid.fluid_cells()) {
    i = cell->i();
    j = cell->j();

    if (cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID) {
      _F(i, j) = _U(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                   Discretization::convection_u(_U, _V, i, j) + _gx);
      if (_energy_eq) {
        _F(i, j) =
            _F(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i + 1, j)) * _gx;
      }
    }
    if (cell->neighbour(border_position::TOP)->type() == cell_type::FLUID) {
      _G(i, j) = _V(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                   Discretization::convection_v(_U, _V, i, j) + _gy);

      if (_energy_eq) {
        _G(i, j) =
            _G(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i, j + 1)) * _gy;
      }
    }
  }
  for (auto cell : grid.surface_cells()) {
    i = cell->i();
    j = cell->j();
    if (cell->neighbour(border_position::BOTTOM)->type() == cell_type::EMPTY ||
        cell->neighbour(border_position::SOUTHWEST)->type() ==
            cell_type::EMPTY ||
        cell->neighbour(border_position::LEFT)->type() == cell_type::EMPTY) {
      _F(i, j) = _U(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                   Discretization::convection_u(_U, _V, i, j) + _gx);
      if (_energy_eq) {
        _F(i, j) =
            _F(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i + 1, j)) * _gx;
      }
      _G(i, j) = _V(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                   Discretization::convection_v(_U, _V, i, j) + _gy);

      if (_energy_eq) {
        _G(i, j) =
            _G(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i, j + 1)) * _gy;
      }
    } else if (cell->neighbour(border_position::NORTHWEST)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::TOP)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::NORTHWESTSOUTH)->type() ==
                   cell_type::EMPTY) {
      _F(i, j) = _U(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                   Discretization::convection_u(_U, _V, i, j) + _gx);
      if (_energy_eq) {
        _F(i, j) =
            _F(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i + 1, j)) * _gx;
      }
    } else if (cell->neighbour(border_position::SOUTHEAST)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::RIGHT)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::EASTWESTSOUTH)->type() ==
                   cell_type::EMPTY) {
      _G(i, j) = _V(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                   Discretization::convection_v(_U, _V, i, j) + _gy);

      if (_energy_eq) {
        _G(i, j) =
            _G(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i, j + 1)) * _gy;
      }
      // } else if (cell->neighbour(border_position::NORTHWESTSOUTH)->type() ==
      // cell_type::EMPTY) {
      //     _F(i, j) = _U(i, j) +
      //                _dt * (_nu * Discretization::diffusion(_U, i, j) -
      //                Discretization::convection_u(_U, _V, i, j));
      //     if (_energy_eq) {
      //         _F(i, j) = _F(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i + 1,
      //         j)) * _gx;
      //     }
      // } else if (cell->neighbour(border_position::EASTWESTSOUTH)->type() ==
      // cell_type::EMPTY) {
      //     _G(i, j) = _V(i, j) +
      //                _dt * (_nu * Discretization::diffusion(_V, i, j) -
      //                Discretization::convection_v(_U, _V, i, j));

      //     if (_energy_eq) {
      //         _G(i, j) = _G(i, j) - (_beta * _dt / 2) * (_T(i, j) + _T(i, j +
      //         1)) * _gy;
      //     }
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

  for (auto cell : grid.surface_cells()) {
    i = cell->i();
    j = cell->j();
    if (cell->neighbour(border_position::BOTTOM)->type() == cell_type::EMPTY ||
        cell->neighbour(border_position::SOUTHWEST)->type() ==
            cell_type::EMPTY ||
        cell->neighbour(border_position::LEFT)->type() == cell_type::EMPTY) {
      _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();
      _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();

    } else if (cell->neighbour(border_position::NORTHWEST)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::TOP)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::NORTHWESTSOUTH)->type() ==
                   cell_type::EMPTY) {
      _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();

    } else if (cell->neighbour(border_position::SOUTHEAST)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::RIGHT)->type() ==
                   cell_type::EMPTY ||
               cell->neighbour(border_position::EASTWESTSOUTH)->type() ==
                   cell_type::EMPTY) {
      _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();

      // } else if (cell->neighbour(border_position::NORTHWESTSOUTH)->type() ==
      // cell_type::EMPTY) {
      //     _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();

    }  // } else if (cell->neighbour(border_position::EASTWESTSOUTH)->type() ==
       // cell_type::EMPTY) {
       //     _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();
       // }
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
bool Fields::energy_eq() { return _energy_eq; }

Matrix<double> &Fields::p_matrix() { return _P; }
Matrix<double> &Fields::u_matrix() { return _U; }
Matrix<double> &Fields::v_matrix() { return _V; }
Matrix<double> &Fields::f_matrix() { return _F; }
Matrix<double> &Fields::t_matrix() { return _T; }
Matrix<double> &Fields::g_matrix() { return _G; }
Matrix<double> &Fields::rs_matrix() { return _RS; }

double Fields::dt() const { return _dt; }

double Fields::Re() const { return _Re; }

double Fields::gx() const { return _gx; };

double Fields::gy() const { return _gy; };
