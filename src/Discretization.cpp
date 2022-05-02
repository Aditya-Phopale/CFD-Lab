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

// page number 6
double Discretization::convection_u(const Matrix<double> &U,
                                    const Matrix<double> &V, int i, int j) {
  double term1 =
      (pow(0.5 * (U(i, j) + U(i + 1, j)), 2) -
       pow(0.5 * (U(i - 1, j) + U(i, j)), 2)) /
          _dx +
      _gamma *
          (0.5 * abs(U(i, j) + U(i + 1, j)) * 0.5 * (U(i, j) - U(i + 1, j)) -
           0.5 * abs(U(i - 1, j) + U(i, j)) * 0.5 * (U(i - 1, j) - U(i, j))) /
          _dx;

  double term2 =
      (0.5 * (U(i, j) + U(i, j + 1) * 0.5 * (V(i, j) + V(i + 1, j))) -
       0.5 * (U(i, j - 1) + U(i, j)) * 0.5 * (V(i, j - 1) + V(i + 1, j - 1))) /
          _dy +
      _gamma *
          (0.5 * abs(V(i, j) + V(i + 1, j)) * 0.5 * (U(i, j) - U(i, j + 1)) -
           0.5 * abs(V(i, j - 1) + V(i + 1, j - 1)) * 0.5 *
               (U(i, j - 1) - U(i, j))) /
          _dy;

  return term1 + term2;
}

double Discretization::convection_v(const Matrix<double> &U,
                                    const Matrix<double> &V, int i, int j) {
  double term1 =
      (pow(0.5 * (V(i, j) + V(i, j + 1)), 2) -
       pow(0.5 * (V(i, j - 1) + V(i, j)), 2)) /
          _dy +
      _gamma *
          (0.5 * abs(V(i, j) + V(i, j + 1)) * 0.5 * (V(i, j) - V(i, j + 1)) -
           0.5 * abs(V(i, j - 1) + V(i, j)) * 0.5 * (V(i, j - 1) - V(i, j))) /
          _dy;

  double term2 =
      (0.5 * (U(i, j) + U(i, j + 1) * 0.5 * (V(i, j) + V(i + 1, j))) -
       0.5 * (U(i - 1, j) + U(i - 1, j + 1)) * 0.5 * (V(i - 1, j) + V(i, j))) /
          _dy +
      _gamma *
          (0.5 * abs(U(i, j) + U(i, j + 1)) * 0.5 * (V(i, j) - V(i + 1, j)) -
           0.5 * abs(U(i - 1, j - 1) + U(i - 1, j + 1)) * 0.5 *
               (U(i - 1, j) - U(i, j))) /
          _dy;

  return term1 + term2;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {
  // std::cout << A(i - 1, j) << "\n";
  double term1 = (A(i + 1, j) - 2 * A(i, j) + A(i - 1, j)) / (_dx * _dx);
  double term2 = (A(i, j + 1) - 2 * A(i, j) + A(i, j - 1)) / (_dy * _dy);

  return term1 + term2;
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
  double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                  (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
  return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
  double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) +
                  (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
  return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j,
                                   int i_offset, int j_offset) {}