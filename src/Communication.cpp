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

double Communication::reduce_min(double dt) {
   return MPI_Reduce(&Communication::rank, &dt, 1, MPI_DOUBLE, MPI_MIN, Communication::root_rank, MPI_COMM_WORLD);
}

void Communication::communicate(Matrix<double> &A, Domain domain) {
    int count_v = domain.jmax - domain.jmin - 2;
    int count_h = domain.imax - domain.imin - 2;
    std::vector<double> buffer_send_v(count_v);
    std::vector<double> buffer_send_h(count_h);
    std::vector<double> buffer_recv_v(count_v);
    std::vector<double> buffer_recv_h(count_h);
    
    if(domain.domain_neighbors.at(2) != -1){ //&& domain.domain_neighbors.at(0) != -1){
      for (int k = 1; k <= count_v; k++) {
        buffer_send_v.at(k-1) = A(domain.imin + 1, k + domain.jmin);
      }
      if(domain.domain_neighbors.at(0) != -1){
        MPI_Sendrecv(buffer_send_v.data(), count_v, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, buffer_recv_v.data(), count_v,
                 MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (int k = 1; k <= count_v; k++) {
          A(domain.imax - 1, k + domain.jmin) = buffer_recv_v.at(k-1);
        }
      }
      else{
        MPI_Send(buffer_send_v.data(), count_v, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD);
      }
    }
    
    if(domain.domain_neighbors.at(0) != -1){ //&& domain.domain_neighbors.at(0) != -1){
      for (int k = 1; k <= count_v; k++) {
        buffer_send_v.at(k-1) = A(domain.imax - 2, k + domain.jmin);
      }
      if(domain.domain_neighbors.at(2) != -1){
        MPI_Sendrecv(buffer_send_v.data(), count_v, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, buffer_recv_v.data(), count_v,
                 MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (int k = 1; k <= count_v; k++) {
          A(domain.imin, k + domain.jmin) = buffer_recv_v.at(k-1);
        }
      }
      else{
        MPI_Send(buffer_send_v.data(), count_v, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD);
      }
    }
    // for (int k = 1; k <= count_h; k++) {
    //   buffer_send_v.at(k-1) = A(k + domain.imin, domain.imax);
    // }
    if(domain.domain_neighbors.at(1) != -1){ //&& domain.domain_neighbors.at(0) != -1){
      for (int k = 1; k <= count_h; k++) {
        buffer_send_v.at(k-1) = A(k + domain.imin, domain.jmax - 2);
      }
      if(domain.domain_neighbors.at(3) != -1){
        MPI_Sendrecv(buffer_send_v.data(), count_h, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, buffer_recv_v.data(), count_h,
                 MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (int k = 1; k <= count_h; k++) {
          A(k + domain.imin, domain.jmin) = buffer_recv_v.at(k-1);
        }
      }
      else{
        MPI_Send(buffer_send_v.data(), count_h, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD);
      }
    }

    if(domain.domain_neighbors.at(3) != -1){ //&& domain.domain_neighbors.at(0) != -1){
      for (int k = 1; k <= count_h; k++) {
        buffer_send_v.at(k-1) = A(k + domain.imin, domain.jmin + 1);
      }
      if(domain.domain_neighbors.at(1) != -1){
        MPI_Sendrecv(buffer_send_v.data(), count_h, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, buffer_recv_v.data(), count_h,
                 MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        for (int k = 1; k <= count_h; k++) {
          A(k + domain.imin, domain.jmax - 1) = buffer_recv_v.at(k-1);
        }
      }
      else{
        MPI_Send(buffer_send_v.data(), count_h, MPI_DOUBLE,
                 domain.domain_neighbors.at(2), 0, MPI_COMM_WORLD);
      }
    }
    
  //   MPI_Sendrecv(buffer_w.data(), count, MPI_DOUBLE,
  //                domain.domain_neighbors.at(0), 0, buffer_recv.data(), count,
  //                MPI_DOUBLE, domain.domain_neighbors.at(0), 0, MPI_COMM_WORLD,
  //                MPI_STATUS_IGNORE);
  //   for (int k = domain.jmin; k < domain.jmax; k++) {
  //     A(domain.imax - 1, k) = buffer_recv.at(k);
  //   }

  // if (domain.domain_neighbors.at(1) != -1) {
  //   std::vector <double> buffer_n(domain.imax - domain.imin - 1);
  //   for (int k = domain.imin; k < domain.imax; k++) {
  //     buffer_n.at(k) = A(k, domain.jmax + 1);
  //   }
  // }

  // if (domain.domain_neighbors.at(2) != -1) {
  //   //std::vector <double> buffer_w(domain.jmax - domain.jmin - 1);
  //   for (int k = domain.jmin; k < domain.jmax; k++) {
  //     buffer_w.at(k) = A(domain.imin, k);
  //   }
  // }

  // if (domain.domain_neighbors.at(3) != -1) {
  //   std::vector <double> buffer_s(domain.imax - domain.imin - 1);
  //   for (int k = domain.imin; k < domain.imax; k++) {
  //     buffer_s.at(k) = A(k, domain.jmin);
  //   }
  // }
}
