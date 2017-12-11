#include "petscdmplex.h"
#include <iostream>

/*
  Strukturovana sit o rozmerech 4x2

  10----11----12----13----14
  |      |     |     |     |
  |  4   |  5  |  6  |  7  |
  |      |     |     |     |
  5------6-----7-----8-----9
  |      |     |     |     |
  |  0   |  1  |  2  |  3  |
  |      |     |     |     |
  0------1-----2-----3-----4


 */

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

    { // Steny definujici bunky
      int fid = ncells + nverts;

      int fi[nx+1][ny];
      for (int j=0; j<ny; j++)
	for (int i=0; i<nx+1; i++) 
	  fi[i][j] = fid++;

      int fj[nx][ny+1];
      for (int j=0; j<ny+1; j++)
	for (int i=0; i<nx; i++) 
	  fj[i][j] = fid++;
      
      for (int j=0; j<ny; j++)
	for (int i=0; i<nx; i++) {
	  int v[] = { fi[i][j], fj[i][j], fi[i+1][j], fj[i][j+1] };
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
