#include "Fields.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI,
               double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
  _U = Matrix<double>(imax + 2, jmax + 2, UI);
  _V = Matrix<double>(imax + 2, jmax + 2, VI);
  _P = Matrix<double>(imax + 2, jmax + 2, PI);

  _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
  _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
  _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}
//_F(i,j) = Discretization::convection_u(U,i,j)
void Fields::calculate_fluxes(Grid &grid) {
  // int i;
  // int j;
  // for(auto cell : grid.fluid_cells()){
  //   i = cell->i();
  //   j = cell->j();
  //   _F(i,j) = _U(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
  //                                  Discretization::convection_u(_U, _V, i, j));
  //   _G(i, j) = _V(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
  //                                  Discretization::convection_u(_U, _V, i, j));
  // }




  for (int i{1}; i < grid.imax(); i++) {
    for (int j{1}; j < grid.jmax() + 1; j++) {
      _F(i, j) = _U(i, j) + _dt * (_nu * Discretization::diffusion(_U, i, j) -
                                   Discretization::convection_u(_U, _V, i, j));
    }
  }

  for (int i{1}; i < grid.imax() + 1; i++) {
    for (int j{1}; j < grid.jmax(); j++) {
      _G(i, j) = _V(i, j) + _dt * (_nu * Discretization::diffusion(_V, i, j) -
                                   Discretization::convection_u(_U, _V, i, j));
    }
  }
}

void Fields::calculate_rs(Grid &grid) {
  int i;
  int j;
  for(auto cell : grid.fluid_cells()){
    i = cell->i();
    j = cell->j();
    double term1 = (_F(i, j) - _F(i - 1, j)) / grid.dx();
    double term2 = (_G(i, j) - _G(i, j - 1)) / grid.dy();
    _RS(i, j) = (term1 + term2) / _dt;
  }



//   for (int i{1}; i < grid.imax() + 1; i++) {
//     for (int j{1}; j < grid.jmax() + 1; j++) {
//       double term1 = (_F(i, j) - _F(i - 1, j)) / grid.dx();
//       double term2 = (_G(i, j) - _G(i, j - 1)) / grid.dy();
//       _RS(i, j) = (term1 + term2) / _dt;
//     }
//   }
}

void Fields::calculate_velocities(Grid &grid) {
  // int i;
  // int j;
  // for(auto cell : grid.fluid_cells()){
  //   i = cell->i();
  //   j = cell->j();
  //   _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();
  //   _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();

  // }



  for (int i{1}; i < grid.imax(); i++) {
    for (int j{1}; j < grid.jmax() + 1; j++) {
      _U(i, j) = _F(i, j) - _dt * (_P(i + 1, j) - _P(i, j)) / grid.dx();
    }
  }

  for (int i{1}; i < grid.imax() + 1; i++) {
    for (int j{1}; j < grid.jmax(); j++) {
      _V(i, j) = _G(i, j) - _dt * (_P(i, j + 1) - _P(i, j)) / grid.dy();
    }
  }
}

double Fields::calculate_dt(Grid &grid) {
  double CFLu = 0.0;
  double CFLv = 0.0;
  double CFLnu = 0.0;

  double dx2 = grid.dx() * grid.dx();
  double dy2 = grid.dy() * grid.dy();

  double u_max = 0.0;
  double v_max = 0.0;

  int i,j;
  for(auto cell : grid.fluid_cells()){
    i = cell->i();
    j = cell->j();
    u_max = std::max(u_max, fabs(_U(i, j)));
    v_max = std::max(v_max, fabs(_V(i, j)));
  }

  // for (int i{1}; i < grid.imax(); i++) {
  //   for (int j{1}; j < grid.jmax() + 1; j++) {
  //     u_max = std::max(u_max, fabs(_U(i, j)));
  //   }
  // }
  // for (int i{1}; i < grid.imax() + 1; i++) {
  //   for (int j{1}; j < grid.jmax(); j++) {
  //     v_max = std::max(v_max, fabs(_V(i, j)));
  //   }
  // }

  CFLu = grid.dx() / u_max;
  CFLv = grid.dy() / v_max;
  CFLnu = (0.5 / _nu) * (1.0 / (1 / dx2 + 1 / dy2));

  _dt = _tau * (std::min({CFLnu, CFLu, CFLv}));

  return _dt;
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
