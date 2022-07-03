#pragma once
#include "Datastructures.hpp"

class Particle {
  public:

    void advance_particle(double &dt);

    void calculate_velocities(double &dx, double &dy, Matrix<double> &u , Matrix<double> &v);

    //void set_particles(Grid &grid, Fields &field);

    double& x_pos();

    double& y_pos();

    double& u();

    double& v();

  private:
    double x;
    double y;
    double vel_u;
    double vel_v;
};

// std::vector<particle> initialize_particles(int ppc);
