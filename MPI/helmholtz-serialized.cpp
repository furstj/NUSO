#include "mpi.h"
#include <iostream>
#include <cmath>

/*
  Solution of u - u_xx = f with homogeneous Dirichlet boundary conditions

  NOTE: Bad communication pattern - implicit serialization
 */

const int N = 10;   // Number of points per CPU

using namespace std;

double f(double x) {
  if (x<0.5) return 1;
  else return -1;
};

void update_ghost_zones(double* u, int N, int rank, int ncpu) {
  MPI_Status s;

  int right = rank+1;
  int left  = rank-1;

  // 0 -> 1 -> 2 -> 3 -> ...
  if (right < ncpu) 
    MPI_Send( &(u[N]), 1, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);
  if (left>=0)
    MPI_Recv( &(u[0]), 1, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &s);

  // 0 <- 1 <- 2 <- 3 <- ...
  if (left>=0) 
    MPI_Send( &(u[1]), 1, MPI_DOUBLE, left, 1, MPI_COMM_WORLD);
  if (right < ncpu) 
    MPI_Recv( &(u[N+1]), 1, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, &s);
}

int main(int argc, char **argv) {

    double u[N+2], un[N+2];

    int ncpu, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double h = 1.0 / ncpu / N;

    // Initial condition
    for (int i=0; i<N+2; i++) 
      u[i] = 0;


    for (int iter=0; iter<100; iter++) {

      update_ghost_zones(u, N, rank, ncpu);

      for (int i=1; i<=N; i++) {
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
