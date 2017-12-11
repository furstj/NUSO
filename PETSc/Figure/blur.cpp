#include <iostream>
#include <fstream>

#include <petscdmda.h>

const int width  = 256;
const int height = 256;

int max(int a, int b) {
  if (a>b) return a;
  else return b;
}

int min(int a, int b) {
  if (a<b) return a;
  else return b;
}

using namespace std;

int main(int argc, char** argv) {

  PetscInitialize( &argc, &argv, (char*)0, 0);

  // Vytvoreni distribuovaneho obrazku
  DM da;
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
             DMDA_STENCIL_STAR,
	     width, height, PETSC_DECIDE, PETSC_DECIDE, 
	     1, 1, PETSC_NULL, PETSC_NULL, &da);
  DMSetUp(da);
  
  // Lokalni obraz pole s pridanymi bunkami
  Vec l;
  DMCreateLocalVector(da, &l);    
  
  // Globalni vektor na ulozeni obrazku
  Vec g;
  DMCreateGlobalVector(da, &g);

  // Nastaveni pocatecnich hodnot obrazku
  double **u, **un;
  DMDAVecGetArray(da, g, &u);

  int i1, j1, n, m;
  DMDAGetCorners(da, &i1, &j1, PETSC_NULL, &n, &m, PETSC_NULL);
  cout << i1 << "\t" << j1 << "\t" << n << "\t" << m << endl;
  for (int i=i1; i<i1+n; i++)
    for (int j=j1; j<j1+m; j++)
      {
	double x = (i-0.5*width) * 2.0 / width;
	double y = (j-0.5*height) * 2.0 / height;
	if (x*x + y*y < 1.0/16.) u[j][i] = 0;
	else u[j][i] = 1;
      }

  DMDAVecRestoreArray(da, g, &u);

  // Vyhlazeni obrazku (explicitni metoda)
  for (int iter=0; iter<1000; iter++) 
    {
      // Preposlani hranic
      DMGlobalToLocalBegin(da, g, INSERT_VALUES, l);
      DMGlobalToLocalEnd(da, g, INSERT_VALUES, l);

      DMDAVecGetArray(da, l, &u);
      DMDAVecGetArray(da, g, &un);
      
      for (int i=max(i1,1); i<min(i1+n,width-1); i++)
	for (int j=max(j1,1); j<min(j1+m,height-1); j++)
	  un[j][i] = u[j][i] + 0.1*(u[j+1][i]+u[j-1][i]+u[j][i+1]+u[j][i-1]
				    -4*u[j][i]);
      
      DMDAVecRestoreArray(da, g, &un);
      DMDAVecRestoreArray(da, l, &u);

    }

  
  // Ulozeni obrazku do pgm formatu
  // Preposlani obrazku na procesor 0
  {
    Vec nat;
    DMDACreateNaturalVector(da, &nat);
    DMDAGlobalToNaturalBegin(da, g, INSERT_VALUES, nat);
    DMDAGlobalToNaturalEnd(da, g, INSERT_VALUES, nat);
    
    Vec u_loc;
    VecCreateSeq(PETSC_COMM_SELF, width*height, &u_loc);
    
    VecScatter ctx;
    VecScatterCreateToZero(nat,&ctx,&u_loc);
    VecScatterBegin(ctx, nat, u_loc, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, nat, u_loc, INSERT_VALUES, SCATTER_FORWARD);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank==0) {
      double* image;
      VecGetArray(u_loc, &image);
      ofstream pgm("obrazek.pgm");
      pgm << "P2" << endl;
      pgm << "# output of blur.cpp" << endl;
      pgm << width << " " << height << endl;
      pgm << 256 << endl;
      for (int i=0; i<width*height; i++)
	pgm << (int)(256*image[i]) << endl;
      VecRestoreArray(u_loc, &image);
      
    }
    
  }

  PetscFinalize();  
  return 0;
}
