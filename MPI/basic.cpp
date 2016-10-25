#include <iostream>
#include "mpi.h"

using namespace std;

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  cout << "I'm " << rank << " from " << size << " cpus." << endl;

  MPI_Finalize();

  return 0;
}

