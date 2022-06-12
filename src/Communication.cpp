#include "Communication.hpp"

int Communication::rank = 0;
int Communication::size = 0;

void Communication::init_parallel(int argn, char **args) {
  MPI_Init(&argn, &args);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
}

void Communication::finalize() { MPI_Finalize(); }

void Communication::communicate(Fields &A, Domain domain) {
  if (domain.domain_neighbors.at(0) != -1) {
    int count = domain.jmax - domain.jmin - 2;
    std::vector<double> buffer_e(count);
    std::vector<double> buffer_recv(count);
    for (int k = domain.jmin + 1; k < domain.jmax - 1; k++) {
      buffer_e.at(k) = A.at(domain.imax - 1, k);
    }
    MPI_Sendrecv(buffer_e.data(), count, MPI_DOUBLE,
                 domain.domain_neighbors.at(0), 0, buffer_recv.data(), count,
                 MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
    for (int k = domain.jmin; k < domain.jmax; k++) {
      A.at(domain.imax - 2, k) = buffer_recv.at(k);
    }
  }

  if (domain.domain_neighbors.at(1) != -1) {
    std::vector < double buffer_n(domain.imax - domain.imin - 1);
    for (int k = domain.imin; k < domain.imax; k++) {
      buffer_n.at(k) = A.at(k, jmax + 1);
    }
  }

  if (domain.domain_neighbors.at(2) != -1) {
    std::vector < double buffer_w(domain.jmax - domain.jmin - 1);
    for (int k = domain.jmin; k < domain.jmax; k++) {
      buffer_w.at(k) = A.at(imin, k);
    }
  }

  if (domain.domain_neighbors.at(3) != -1) {
    std::vector < double buffer_s(domain.imax - domain.imin - 1);
    for (int k = domain.imin; k < domain.imax; k++) {
      buffer_s.at(k) = A.at(k, jmin);
    }
  }
}
