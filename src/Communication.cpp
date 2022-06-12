#include "Communication.hpp"

int Communication::rank = 0;
int Communication::size = 0;
int Communication::root_rank = 0;

void Communication::init_parallel(int argn, char **args) {
  MPI_Init(&argn, &args);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
}

void Communication::finalize() { MPI_Finalize(); }

double Communication::reduce_min(double &dt) {
  double result = 0;
  MPI_Reduce(&dt, &result, 1, MPI_DOUBLE, MPI_MIN, Communication::root_rank,
             MPI_COMM_WORLD);
  return result;
}

double Communication::reduce_max(double &vel) {
  double result = 0;
  MPI_Reduce(&vel, &result, 1, MPI_DOUBLE, MPI_MAX, Communication::root_rank,
             MPI_COMM_WORLD);
  return result;
}

double Communication::reduce_sum(double &res) {
  double result = 0;
  MPI_Reduce(&res, &result, 1, MPI_DOUBLE, MPI_SUM, Communication::root_rank,
             MPI_COMM_WORLD);
  return result;
}

void Communication::communicate(Matrix<double> &A, Domain domain) {
  int count_v = domain.jmax - domain.jmin - 2;
  int count_h = domain.imax - domain.imin - 2;
  std::vector<double> buffer_send_v(count_v);
  std::vector<double> buffer_send_h(count_h);
  std::vector<double> buffer_recv_v(count_v);
  std::vector<double> buffer_recv_h(count_h);

  if (domain.domain_neighbors.at(2) != -1) {
    for (int k = 1; k <= count_v; k++) {
      buffer_send_v.at(k - 1) = A(1, k);
    }
    if (domain.domain_neighbors.at(0) != -1) {
      MPI_Sendrecv(buffer_send_v.data(), count_v, MPI_DOUBLE,
                   domain.domain_neighbors.at(2), 0, buffer_recv_v.data(),
                   count_v, MPI_DOUBLE, domain.domain_neighbors.at(0), 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (int k = 1; k <= count_v; k++) {
        A(domain.size_x + 1, k) = buffer_recv_v.at(k - 1);
      }
    } else {
      MPI_Send(buffer_send_v.data(), count_v, MPI_DOUBLE,
               domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD);
    }
  }

  if (domain.domain_neighbors.at(0) != -1) {
    for (int k = 1; k <= count_v; k++) {
      buffer_send_v.at(k - 1) = A(domain.size_x, k);
    }
    if (domain.domain_neighbors.at(2) != -1) {
      MPI_Sendrecv(buffer_send_v.data(), count_v, MPI_DOUBLE,
                   domain.domain_neighbors.at(0), 0, buffer_recv_v.data(),
                   count_v, MPI_DOUBLE, domain.domain_neighbors.at(2), 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (int k = 1; k <= count_v; k++) {
        A(0, k) = buffer_recv_v.at(k - 1);
      }
    } else {
      MPI_Send(buffer_send_v.data(), count_v, MPI_DOUBLE,
               domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD);
    }
  }

  if (domain.domain_neighbors.at(1) != -1) {
    for (int k = 1; k <= count_h; k++) {
      buffer_send_v.at(k - 1) = A(k, domain.size_y);
    }
    if (domain.domain_neighbors.at(3) != -1) {
      MPI_Sendrecv(buffer_send_v.data(), count_h, MPI_DOUBLE,
                   domain.domain_neighbors.at(1), 0, buffer_recv_v.data(),
                   count_h, MPI_DOUBLE, domain.domain_neighbors.at(3), 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (int k = 1; k <= count_h; k++) {
        A(k, 0) = buffer_recv_v.at(k - 1);
      }
    } else {
      MPI_Send(buffer_send_v.data(), count_h, MPI_DOUBLE,
               domain.domain_neighbors.at(1), 0, MPI_COMM_WORLD);
    }
  }

  if (domain.domain_neighbors.at(3) != -1) {
    for (int k = 1; k <= count_h; k++) {
      buffer_send_v.at(k - 1) = A(k, 1);
    }
    if (domain.domain_neighbors.at(1) != -1) {
      MPI_Sendrecv(buffer_send_v.data(), count_h, MPI_DOUBLE,
                   domain.domain_neighbors.at(3), 0, buffer_recv_v.data(),
                   count_h, MPI_DOUBLE, domain.domain_neighbors.at(1), 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for (int k = 1; k <= count_h; k++) {
        A(k, domain.size_y + 1) = buffer_recv_v.at(k - 1);
      }
    } else {
      MPI_Send(buffer_send_v.data(), count_h, MPI_DOUBLE,
               domain.domain_neighbors.at(3), 0, MPI_COMM_WORLD);
    }
  }
}
