/*
In this file, we discretize our Navier-Stokes equation in space with finite
differences We discretize the convection terms; horizontal velocity U, vertical
velocity V, the diffusion term, laplacian, and the SOR helper function.
*/
#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
  _dx = dx;
  _dy = dy;
  _gamma = gamma;
}

// Calculating the value of convective part of U
double Discretization::convection_u(const Matrix<double> &U,
                                    const Matrix<double> &V, int i, int j) {
  double term1 =
      (1 / _dx) * (interpolate(U, i, j, 1, 0) * interpolate(U, i, j, 1, 0) -
                   interpolate(U, i, j, -1, 0) * interpolate(U, i, j, -1, 0)) +
      _gamma / _dx *
          (fabs(interpolate(U, i, j, 1, 0)) * (U(i, j) - U(i + 1, j)) / 2 -
           fabs(interpolate(U, i, j, -1, 0)) * (U(i - 1, j) - U(i, j)) / 2);

  double term2 =
      (1 / _dy) *
          (interpolate(V, i, j, 1, 0) * interpolate(U, i, j, 0, 1) -
           interpolate(V, i, j - 1, 1, 0) * interpolate(U, i, j, 0, -1)) +
      (_gamma / _dy) *
          (fabs(interpolate(V, i, j, 1, 0)) * (U(i, j) - U(i, j + 1)) / 2 -
           fabs(interpolate(V, i, j - 1, 1, 0)) * (U(i, j - 1) - U(i, j)) / 2);

  return term1 + term2;
}

// Calculating the value of convective part of V

double Discretization::convection_v(const Matrix<double> &U,
                                    const Matrix<double> &V, int i, int j) {
  double term1 =
      (1 / _dy) * (interpolate(V, i, j, 0, 1) * interpolate(V, i, j, 0, 1) -
                   interpolate(V, i, j, 0, -1) * interpolate(V, i, j, 0, -1)) +
      (_gamma / _dy) *
          (fabs(interpolate(V, i, j, 0, 1)) * (V(i, j) - V(i, j + 1)) / 2 -
           fabs(interpolate(V, i, j, 0, -1)) * (V(i, j - 1) - V(i, j)) / 2);

  double term2 =
      (1 / _dx) *
          (interpolate(U, i, j, 0, 1) * interpolate(V, i, j, 1, 0) -
           interpolate(U, i - 1, j, 0, 1) * interpolate(V, i, j, -1, 0)) +
      (_gamma / _dx) *
          (fabs(interpolate(U, i, j, 0, 1)) * (V(i, j) - V(i + 1, j)) / 2 -
           fabs(interpolate(U, i - 1, j, 0, 1)) * (V(i - 1, j) - V(i, j)) / 2);

  return term1 + term2;
}

// Calculating the value of convective part of T
double Discretization::convection_T(const Matrix<double> &U,
                                    const Matrix<double> &V,
                                    const Matrix<double> &T, int i, int j) {
  double term1 =
      (1 / (2 * _dx)) * (U(i, j) * (T(i, j) + T(i + 1, j)) -
                         U(i - 1, j) * (T(i - 1, j) + T(i, j))) +
      (1 / (2 * _dx)) * (fabs(U(i, j)) * (T(i, j) - T(i + 1, j)) -
                         fabs(U(i - 1, j)) * (T(i - 1, j) - T(i, j)));

  double term2 =
      (1 / (2 * _dy)) * (V(i, j) * (T(i, j) + T(i, j + 1)) -
                         V(i, j - 1) * (T(i, j - 1) + T(i, j))) +
      (1 / (2 * _dy)) * (fabs(V(i, j)) * (T(i, j) - T(i, j + 1)) -
                         fabs(V(i, j - 1)) * (T(i, j - 1) - T(i, j)));

  return term1 + term2;
};

// Using the same for calculating diffusive part of U and V

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
  double term1 = (A(i + 1, j) - 2 * A(i, j) + A(i - 1, j)) / (_dx * _dx);
  double term2 = (A(i, j + 1) - 2 * A(i, j) + A(i, j - 1)) / (_dy * _dy);

  return term1 + term2;
}

// Calculating the laplacian part of the equation

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
  double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                  (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
  return result;
}

// Calculating the SOR Helper

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
  double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) +
                  (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
  return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j,
                                   int i_offset, int j_offset) {
  return (A(i, j) + A(i + i_offset, j + j_offset)) / 2;
}