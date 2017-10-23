#include "mpi.h"
#include <iostream>
#include <cmath>

/*
  Solution of u - u_xx = f with homogeneous Dirichlet boundary conditions
 */

const int N = 10;   // Number of points per CPU

using namespace std;

double f(double x) {
  if (x<0.5) return 1;
  else return -1;
};

struct Requests {
  MPI_Request send_right, send_left;
  MPI_Request recv_left, recv_right;
};

void update_ghost_zones_begin(double* u, int N, int rank, int ncpu, Requests* req) {

  int right = rank+1;
  int left  = rank-1;

  // 0 -> 1 -> 2 -> 3 -> ...
  if (right < ncpu) 
    MPI_Isend( &(u[N]), 1, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &(req->send_right));
  if (left>=0)
    MPI_Irecv( &(u[0]), 1, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &(req->recv_left));

  // 0 <- 1 <- 2 <- 3 <- ...
  if (left>=0) 
    MPI_Isend( &(u[1]), 1, MPI_DOUBLE, left, 1, MPI_COMM_WORLD, &(req->send_left));
  if (right < ncpu) 
    MPI_Irecv( &(u[N+1]), 1, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, &(req->recv_right));
}

void update_ghost_zones_end(int rank, int ncpu, Requests* req) {
  MPI_Status s;
  int right = rank+1;
  int left  = rank-1;

  // 0 -> 1 -> 2 -> 3 -> ...
  if (right < ncpu) 
    MPI_Wait(&(req->send_right), &s);
  if (left>=0)
    MPI_Wait(&(req->recv_left), &s);

  // 0 <- 1 <- 2 <- 3 <- ...
  if (left>=0) 
    MPI_Wait(&(req->send_left), &s);
  if (right < ncpu) 
    MPI_Wait(&(req->recv_right), &s);
}

int main(int argc, char **argv) {

    double u[N+2], un[N+2];

    int ncpu, rank;

    Requests req;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double h = 1.0 / ncpu / N;

    // Initial condition
    for (int i=0; i<N+2; i++) 
      u[i] = 0;


    for (int iter=0; iter<100; iter++) {

      update_ghost_zones_begin(u, N, rank, ncpu, &req);

      for (int i=2; i<=N-1; i++) {
	double x = (i-1 + rank*N) * h;
	un[i] = (h*h*f(x) + u[i-1] + u[i+1]) / (h*h + 2); 
      }

      update_ghost_zones_end(rank, ncpu, &req);
      
      {
	int i = 1;
	double x = (i-1 + rank*N) * h;
	un[i] = (h*h*f(x) + u[i-1] + u[i+1]) / (h*h + 2); 
      }
      {
	int i = N;
	double x = (i-1 + rank*N) * h;
	un[i] = (h*h*f(x) + u[i-1] + u[i+1]) / (h*h + 2); 
      }

      double rez = 0;
      for (int i=1; i<=N; i++) {
	rez += fabs(u[i]-un[i]);
	u[i] = un[i];
      }

      if (rank==0) std::cout << iter << " " << rez << std::endl;

    }


    MPI_Finalize();
    return 0;
}
