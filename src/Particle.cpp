#include "Particle.hpp"

#include <vector>

void Particle::advance_particle(double &dt) {
    this->x += this->vel_u * dt;
    this->y += this->vel_v * dt;
}

void Particle::calculate_velocities(double &dx, double &dy, Matrix<double> &u , Matrix<double> &v) {
    int i = x / dx + 1;
    int j = (y + dy / 2) / dy + 1;

    double x1, x2, y1, y2;
    x1 = (i - 1) * dx;
    x2 = i * dx;
    y1 = ((j - 1) - 0.5) * dy;
    y2 = (j - 0.5) * dy;

    vel_u = ((x2 - x) * (y2 - y) * u(i - 1, j - 1) + (x - x1) * (y2 - y) * u(i, j - 1) +
             (x2 - x) * (y - y1) * u(i - 1, j) + (x - x1) * (y - y1) * u(i, j)) /
            (dx * dy);

    i = (x + dx/2)/dx + 1;
    j = y/dy + 1;

    x1 = ((i-1)-0.5)*dx;
    x2 = (i-0.5)*dx;
    y1 = (j-1)*dy;
    y2 = j*dy; 

    vel_v = ((x2 - x) * (y2 - y) * v(i - 1, j - 1) + (x - x1) * (y2 - y) * v(i, j - 1) +
             (x2 - x) * (y - y1) * v(i - 1, j) + (x - x1) * (y - y1) * v(i, j)) /
            (dx * dy);
}

double& Particle::x_pos(){
    return x;
}

double& Particle::y_pos(){
    return y;
}

double& Particle::u(){
    return vel_u;
}

double &Particle::v(){
    return vel_v;
}

