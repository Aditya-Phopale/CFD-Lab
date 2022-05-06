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
  
  double term1 = (1/_dx) * (((U(i,j) + U(i+1,j))*(U(i,j) + U(i+1,j))/4) - ((U(i-1,j) + U(i,j))*(U(i-1,j) + U(i,j))/4)) +
	    _gamma/(4*_dx)*(fabs(U(i,j) + U(i+1,j))*(U(i,j) - U(i+1,j)) -
				fabs(U(i-1,j) + U(i,j))*(U(i-1,j) - U(i,j)));
  
  double term2 = (1/(4*_dy)) * ((V(i,j) + V(i+1,j))*(U(i,j) + U(i,j+1)) - (V(i,j-1) + V(i+1,j-1))*(U(i,j-1) + U(i,j))) +
	    _gamma/(4*_dy)*(fabs(V(i,j) + V(i+1,j))*(U(i,j) - U(i,j+1)) -
				fabs(V(i,j-1) + V(i+1,j-1))*(U(i,j-1) - U(i,j)));
  
  return term1 + term2;
}

double Discretization::convection_v(const Matrix<double> &U,
                                    const Matrix<double> &V, int i, int j) {
  double term1 = (1/_dy) * (((V(i,j) + V(i,j+1))*(V(i,j) + V(i,j+1))/4) - ((V(i,j-1) + V(i,j))*(V(i,j-1) + V(i,j))/4)) +
	    _gamma/(4*_dy)*(fabs(V(i,j) + V(i,j+1))*(V(i,j) - V(i,j+1)) -
				fabs(V(i,j-1) + V(i,j))*(V(i,j-1) - V(i,j)));

  double term2 = (1/(4*_dx)) * ((U(i,j) + U(i,j+1))*(V(i,j) + V(i+1,j)) - (U(i-1,j) + U(i-1,j+1))*(V(i-1,j) + V(i,j))) +
	    _gamma/(4 * _dx)*(fabs(U(i,j) + U(i,j+1))*(V(i,j) - V(i+1,j)) -
				fabs (U(i-1,j) + U(i-1,j+1))*(V(i-1,j) - V(i,j)));

  return term1 + term2;
}

double Discretization::diffusion(const Matrix<double> &A, int i, int j) {

  double term1 = ( A(i+1, j) - 2 * A(i, j) + A(i-1, j) ) / (_dx * _dx);
  double term2 = ( A(i, j+1) - 2 * A(i, j) + A(i, j-1) ) / (_dy * _dy);

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