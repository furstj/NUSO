#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (size!=2) {
    cout << "Code is for 2 cpus only!" << endl;
    MPI_Finalize();
    return 0;
  }


  if (rank==0) 
    {
      double x = 42;
      MPI_Send( &x, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD );
    }
  else
    {
      double z;
      MPI_Status stat;

      MPI_Recv( &z, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
      cout << "Recieved " << z << endl;
    }

  MPI_Finalize();

  return 0;
}

