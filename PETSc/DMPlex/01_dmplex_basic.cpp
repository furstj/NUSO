#include "petscdmplex.h"

int main(int argc, char** argv) {
  int myRank;
  PetscInitialize(&argc, &argv, NULL, 0);
  MPI_Comm_rank(PETSC_COMM_WORLD, &myRank);

  // Vytvoreni strukturovane site o velikosti 4x2 
  int nx = 4;
  int ny = 2;

  DM dm;
  DMPlexCreate(PETSC_COMM_WORLD, &dm);

  int ncells = nx*ny;
  int nverts = (nx+1)*(ny+1);
  int nfaces = nx*(ny+1) + (nx+1)*ny;

  if (myRank==0) {

    DMPlexSetChart(dm, 0, ncells + nverts + nfaces);
    for (int c=0; c<ncells; c++) 
      DMPlexSetConeSize(dm, c, 4);
    for (int e=ncells+nverts; e<ncells+nverts+nfaces; e++)
      DMPlexSetConeSize(dm, e, 2);
  }

  DMSetUp(dm);

  if (myRank==0) {

    { // Steny definuici bunky
      int fid = ncells + nverts;

      int fi[nx+1][ny];
      for (int j=0; j<ny; j++)
	for (int i=0; i<nx+1; i++) 
	  fi[i][j] = fid++;

      int fj[nx+1][ny];
      for (int j=0; j<ny+1; j++)
	for (int i=0; i<nx; i++) 
	  fj[i][j] = fid++;
      
      for (int j=0; j<ny; j++)
	for (int i=0; i<nx; i++) {
	  int v[] = { fi[i][j], fj[i+1][j], fi[i][j+1], fj[i][j] };
	  DMPlexSetCone(dm, i+j*nx, v);
	}
    }

    { // Vrcholy definujici steny
      int fid = ncells + nverts;

      for (int j=0; j<ny; j++)
	for (int i=0; i<nx+1; i++) {
	  int v0 = ncells + i + j*(nx+1);
	  int v[] = {v0, v0+nx+1};
	  DMPlexSetCone(dm, fid++, v);
	}

      for (int j=0; j<ny+1; j++)
	for (int i=0; i<nx; i++) {
	  int v0 = ncells + i + j*(nx+1);
	  int v[] = {v0, v0+1};
	  DMPlexSetCone(dm, fid++, v);
	}
    }
  }

  DMPlexSymmetrize(dm);
  DMPlexStratify(dm);
  DMSetDimension(dm, 2);

  DMView(dm, PETSC_VIEWER_STDOUT_WORLD);

  { DM dmDist;
    DMPlexDistribute(dm, 0, NULL, &dmDist);
    if (dmDist) { DMDestroy(&dm); dm = dmDist; }
  }

  DMView(dm, PETSC_VIEWER_STDOUT_WORLD);

  
  DMDestroy(&dm);
  PetscFinalize();
  return 0;
}
