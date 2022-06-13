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
  MPI_Allreduce(&dt, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  return result;
}

double Communication::reduce_max(double &vel) {
  double result = 0;
  MPI_Allreduce(&vel, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return result;
}

double Communication::reduce_sum(double &res) {
  double result = 0;

  MPI_Allreduce(&res, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return result;
}

void Communication::communicate(Matrix<double> &A, Domain domain) {
  // int count_v = domain.jmax - domain.jmin - 2;
  // int count_h = domain.imax - domain.imin - 2;

  if (domain.domain_neighbors.at(2) != -1) {
    // for (int k = 1; k <= count_v; k++) {
    //   buffer_send_v.at(k - 1) = A(1, k);
    // }
    // std::vector<double> buffer_send_v;
    std::vector<double> buffer_recv_v(A.get_col(domain.size_x).size());
    // buffer_send_v = A.get_col(1);
    if (domain.domain_neighbors.at(0) != -1) {
      MPI_Status status;
      MPI_Sendrecv(A.get_col(1).data(), A.get_col(1).size(), MPI_DOUBLE,
                   domain.domain_neighbors.at(2), 0, buffer_recv_v.data(),
                   buffer_recv_v.size(), MPI_DOUBLE,
                   domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD, &status);
      // for (int k = 1; k <= count_v; k++) {
      //   A(domain.size_x + 1, k) = buffer_recv_v.at(k - 1);
      // }
      A.set_col(buffer_recv_v, domain.size_x + 1);
    } else {
      MPI_Send(A.get_col(1).data(), A.get_col(1).size(), MPI_DOUBLE,
               domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD);
      MPI_Recv(buffer_recv_v.data(), buffer_recv_v.size(), MPI_DOUBLE,
               domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      // for (int k = 1; k <= count_v; k++) {
      //   A(0, k) = buffer_recv_v.at(k - 1);
      // }
      A.set_col(buffer_recv_v, 0);
    }
  }

  if (domain.domain_neighbors.at(0) != -1) {
    // for (int k = 1; k <= count_v; k++) {
    //   buffer_send_v.at(k - 1) = A(domain.size_x, k);
    // }
    // std::vector<double> buffer_send_v;
    std::vector<double> buffer_recv_v(A.get_col(1));
    // buffer_send_v = A.get_col(domain.size_x);
    if (domain.domain_neighbors.at(2) != -1) {
      MPI_Status status;
      MPI_Sendrecv(A.get_col(domain.size_x).data(),
                   A.get_col(domain.size_x).size(), MPI_DOUBLE,
                   domain.domain_neighbors.at(0), 0, buffer_recv_v.data(),
                   buffer_recv_v.size(), MPI_DOUBLE,
                   domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD, &status);
      // for (int k = 1; k <= count_v; k++) {
      //   A(0, k) = buffer_recv_v.at(k - 1);
      // }
      A.set_col(buffer_recv_v, 0);
    } else {
      MPI_Recv(buffer_recv_v.data(), buffer_recv_v.size(), MPI_DOUBLE,
               domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Send(A.get_col(domain.size_x).data(), A.get_col(domain.size_x).size(),
               MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD);
      // for (int k = 1; k <= count_v; k++) {
      //   A(domain.size_x + 1, k) = buffer_recv_v.at(k - 1);
      // }
      A.set_col(buffer_recv_v, domain.size_x + 1);
    }
  }

  if (domain.domain_neighbors.at(1) != -1) {
    // for (int k = 1; k <= count_h; k++) {
    //   buffer_send_v.at(k - 1) = A(k, domain.size_y);
    // }
    // std::vector<double> buffer_send_h;
    std::vector<double> buffer_recv_h(A.get_row(1));
    // buffer_send_h = A.get_row(domain.size_y);
    if (domain.domain_neighbors.at(3) != -1) {
      MPI_Status status;
      MPI_Sendrecv(
          A.get_row(domain.size_y).data(), A.get_row(domain.size_y).size(),
          MPI_DOUBLE, domain.domain_neighbors.at(1), 0, buffer_recv_h.data(),
          buffer_recv_h.size(), MPI_DOUBLE, domain.domain_neighbors.at(3), 0,
          MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // for (int k = 1; k <= count_h; k++) {
      //   A(k, 0) = buffer_recv_v.at(k - 1);
      // }
      A.set_row(buffer_recv_h, 0);
    } else {
      MPI_Send(A.get_row(domain.size_y).data(), A.get_row(domain.size_y).size(),
               MPI_DOUBLE, domain.domain_neighbors.at(1), 0, MPI_COMM_WORLD);
      MPI_Recv(buffer_recv_h.data(), buffer_recv_h.size(), MPI_DOUBLE,
               domain.domain_neighbors.at(1), 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      // for (int k = 1; k <= count_h; k++) {
      //   A(k, domain.size_y + 1) = buffer_recv_v.at(k - 1);
      // }
      A.set_row(buffer_recv_h, domain.size_y + 1);
    }
  }

  if (domain.domain_neighbors.at(3) != -1) {
    // for (int k = 1; k <= count_h; k++) {
    //   buffer_send_v.at(k - 1) = A(k, 1);
    // }
    // std::vector<double> buffer_send_h;
    std::vector<double> buffer_recv_h(A.get_row(domain.size_y));
    // buffer_send_h = A.get_row(1);
    if (domain.domain_neighbors.at(1) != -1) {
      MPI_Status status;
      MPI_Sendrecv(A.get_row(1).data(), A.get_row(1).size(), MPI_DOUBLE,
                   domain.domain_neighbors.at(3), 0, buffer_recv_h.data(),
                   buffer_recv_h.size(), MPI_DOUBLE,
                   domain.domain_neighbors.at(1), 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);

      // for (int k = 1; k <= count_h; k++) {
      //   A(k, domain.size_y + 1) = buffer_recv_v.at(k - 1);
      // }
      A.set_row(buffer_recv_h, domain.size_y + 1);
    } else {
      MPI_Recv(buffer_recv_h.data(), buffer_recv_h.size(), MPI_DOUBLE,
               domain.domain_neighbors.at(3), 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      MPI_Send(A.get_row(1).data(), A.get_row(1).size(), MPI_DOUBLE,
               domain.domain_neighbors.at(3), 0, MPI_COMM_WORLD);
      // for (int k = 1; k <= count_h; k++) {
      //   A(k, 0) = buffer_recv_v.at(k - 1);
      // }
      A.set_row(buffer_recv_h, 0);
    }
  }
}
