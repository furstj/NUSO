static char help[] = "Solves Poisson equation - Delta u = 1 on unit square using DMDA.\n\n";

#include <iostream>
#include <fstream>

#include "mpi.h"
#include <petscdmda.h>
#include <petscksp.h>

Mat CreateLaplacianMatrix(DM& da) {
  Mat A;
  DMCreateMatrix(da,&A);
  
  DMDALocalInfo  info;
  DMDAGetLocalInfo(da,&info);

  PetscReal hx  = 1./info.mx, hy = 1./info.my;

  for (int j=info.ys; j<info.ys+info.ym; j++) {
    for (int i=info.xs; i<info.xs+info.xm; i++) {

      MatStencil  row = {0},col[5] = {{0}};
      PetscScalar v[5];
      PetscInt    ncols = 0;
      row.j = j; row.i = i;

      col[ncols].j = j; col[ncols].i = i; v[ncols++] = 2*(hx/hy + hy/hx);

      /* boundaries */
      if (i>0)         {col[ncols].j = j;   col[ncols].i = i-1; v[ncols++] = -hy/hx;}
      if (i<info.mx-1) {col[ncols].j = j;   col[ncols].i = i+1; v[ncols++] = -hy/hx;}
      if (j>0)         {col[ncols].j = j-1; col[ncols].i = i;   v[ncols++] = -hx/hy;}
      if (j<info.my-1) {col[ncols].j = j+1; col[ncols].i = i;   v[ncols++] = -hx/hy;}
      MatSetValuesStencil(A,1,&row,ncols,col,v,INSERT_VALUES);
    }
  }
  
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  return A;
}


void SaveSolution(DM& da, Vec& x) {
  DMDALocalInfo  info;
  DMDAGetLocalInfo(da,&info);

  Vec x_loc;
  VecCreateSeq(PETSC_COMM_SELF, info.mx*info.my, &x_loc);

  VecScatter ctx;
  VecScatterCreateToZero(x,&ctx,&x_loc);
  VecScatterBegin(ctx,x,x_loc,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x,x_loc,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  AO ao;
  DMDAGetAO(da, &ao);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  if (rank==0) {
    double* data;
    VecGetArray(x_loc, &data);

    std::ofstream output("laplace2d.dat");

    for (int i=0; i<info.mx; i++) {
      for (int j=0; j<info.my; j++) {
	int ia = j + i*info.my;
	AOApplicationToPetsc(ao, 1, &ia);
	output << i << "\t" << j << "\t" << data[ia] << std::endl;
      }
      output << std::endl;
    }
    VecRestoreArray(x_loc, &data);
  }

}


int main(int argc, char** argv) {

  PetscInt n = 100;

  PetscInitialize(&argc,&argv,(char*)0,help);

  // Vytvoreni distribuce dat
  DM da;
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
	       DMDA_STENCIL_STAR,n,n,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&da);
  DMSetUp(da);
  
  // Sestaveni matice
  Mat A = CreateLaplacianMatrix(da);
  
  // Prava strana a vektor pro reseni
  Vec u, b;
  DMCreateGlobalVector(da,&u);
  VecDuplicate(u,&b);
  VecSet(b,1.);

  // Reseni soustavy
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,A,A);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,b,u);
  
  SaveSolution(da, u);

  VecDestroy(&b);
  VecDestroy(&u);
  MatDestroy(&A);
  DMDestroy(&da);
  PetscFinalize();
  return 0;
}
