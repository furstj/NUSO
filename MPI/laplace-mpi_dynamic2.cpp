#include "mpi.h"
#include <iostream>
#include <stdlib.h>

#include "timer.hpp"

const int N = 500;

// Alokace 2D pole o indexech N1 <= i < N2 a M1 <= j < M2 ... black magic
double** matrix(int N1, int N2, int M1, int M2) {
    double **a;
    a = (double**) malloc((N2-N1)*sizeof(double*))-N1;
    a[N1] = (double*) malloc((N2-N1)*(M2-M1)*sizeof(double))-M1;
    for (int i=N1+1; i<N2; i++)
	a[i] = a[i-1] + (M2-M1);
    return a;
}

int min(int a, int b) {
    if (a<b) return a;
    else return b;
}

int max(int a, int b) {
    if (a>b) return a;
    else return b;
}

void update_ghost_zones(double** u,int beg,int end, int rank, int ncpu) {
    MPI_Status s;

    if (rank < ncpu-1)
	MPI_Send(&(u[end-1][0]), N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
    if (rank > 0)
	MPI_Recv(&u[beg-1][0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &s);

    if (rank > 0)
	MPI_Send(&u[beg][0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
    if (rank < ncpu-1)
	MPI_Recv(&u[end][0], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &s);
}

int main(int argc, char **argv) {

    int ncpu, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int ncols = N / ncpu;
    int beg = rank*ncols + min(rank, N % ncpu);
    int end = (rank+1)*ncols + min(rank+1, N % ncpu);
    
    double** u  = matrix(beg-1, end+1, 0, N);
    double** un = matrix(beg-1, end+1, 0, N);

    std::cout << rank << "\t" << beg << "\t" << end << "\n";
    

    // Pocatecni podminka
    for (int i=beg; i<end; i++)
	for (int j=0; j<N; j++)
	    u[i][j] = 0;

    // Okrajova podminka na dolni stene
    for (int i=beg; i<end; i++) {
	double x = (double)i / (N-1);
	u[i][0] = 4*x*(1-x);
    }

    Timer stopky("vypocet");

    // Iterace
    stopky.start();
    for (int n=0; n<10000; n++) {

	// Vymena dat na okrajich pasu
	update_ghost_zones(u,beg,end,rank,ncpu);

	// Vypocet pro vnitrni body
	for (int i=max(1,beg); i<min(N-1,end); i++)
	    for (int j=1; j<N-1; j++)
		un[i][j] = (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]) / 4;
		

	// Kopirovani un do u
	for (int i=max(1,beg); i<min(N-1,end); i++)
	    for (int j=1; j<N-1; j++)
		u[i][j] = un[i][j];

	// Tisk na obrazovku
	if (beg <= N/2 && N/2 < end)
	    if (n%100 == 0) std::cout << "n=" << n 
				      << "\t u(0.5,0.5)="<<u[N/2][N/2]<<"\n";

    }
    stopky.stop();
	
    
    MPI_Finalize();
    return 0;
}
