#include <petscdmplex.h>
#include <iostream>

using namespace std;

void PrintMyCoordinates(const DM& dm, int printRank) {
  int myRank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &myRank);

  int dim;
  DMGetCoordinateDim(dm, &dim);

  Vec coord;
  DMGetCoordinates(dm, &coord);

  if (myRank==printRank) {
    int istart, iend;
    VecGetOwnershipRange(coord, &istart, &iend);
    int n = (iend - istart) / dim;

    cout << "Nodal coordinates at cpu " << myRank << ", dim=" << dim << ", n=" << n << endl;

    double* x;
    VecGetArray(coord, &x);
    for (int i=0; i<n; i++) {
      for (int d=0; d<dim; d++) cout << x[dim*i+d] << "\t";
      cout << endl;
    }
    VecRestoreArray(coord, &x);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
  int myRank;
  PetscInitialize(&argc, &argv, NULL, 0);
  MPI_Comm_rank(PETSC_COMM_WORLD, &myRank);

  DM dm;
  DMPlexCreateGmshFromFile(PETSC_COMM_WORLD, "rectangle.msh", PETSC_TRUE, &dm);
  DMSetFromOptions(dm);


  //DMView(dm, PETSC_VIEWER_STDOUT_WORLD);
  { DM dmDist;
    DMPlexDistribute(dm, 1, NULL, &dmDist);
    if (dmDist) { DMDestroy(&dm); dm = dmDist; }
  }
  DMView(dm, PETSC_VIEWER_STDOUT_WORLD);

  int ncpus;
  MPI_Comm_size(PETSC_COMM_WORLD, &ncpus);

  // Nodal coordinates
  for (int cpu=0; cpu<ncpus; cpu++)
    PrintMyCoordinates(dm, cpu);

  DMDestroy(&dm);
  PetscFinalize();
  return 0;
}
