/*    
    Solution of 
    
    - (u_xx + u_yy) = f(x,y)
    u|Gamma = 0
    
    via finite differences.
    
    Note: the code is sequential in order to keep myMatMul as simple as possible!
    */

#include <iostream>

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <fstream>

struct MyContext
{
    std::size_t nx;
    std::size_t ny;
};

PetscErrorCode myMatMul(Mat mfMatrix, Vec x, Vec y)
{
    MyContext* ctx;
    MatShellGetContext(mfMatrix, (void**) &ctx);
    size_t nx = ctx->nx;
    size_t ny = ctx->ny;

    const PetscScalar *xp;
    PetscScalar *yp;
    
    VecGetArrayRead(x, &xp);
    VecGetArray(y, &yp);

    // Laplacian in internal points:
    for (size_t i=0; i<nx; i++)
        for (size_t j=0; j<ny; j++)
        {
            yp[i+j*nx] = 4*xp[i+j*nx];
            if (i>0)    yp[i+j*nx] -= xp[i-1+j*nx];
            if (i<nx-1) yp[i+j*nx] -= xp[i+1+j*nx];
            if (j>0)    yp[i+j*nx] -= xp[i+(j-1)*nx];
            if (j<ny-1) yp[i+j*nx] -= xp[i+(j+1)*nx];
        };

    VecRestoreArrayRead(x, &xp);
    VecRestoreArray(y, &yp);
    
    return 0;
}



int main(int argc,char **args)
{
  int         proc_size, rank;
  
  /* Inicializace */
  CHKERRQ( PetscInitialize( &argc , &args , (char *)0 , 0 ) );

  /* Informace o poctu procesoru a mem cisle */
  MPI_Comm_size( PETSC_COMM_WORLD , &proc_size );
  if (proc_size > 1) {
      std::cerr << "Example is only for single processor!" << std::endl;
      return 1;
  }


  /* Vytvoreni ulohy o velikosti 100x100 bodu */
  MyContext ctx {100, 100};
    
  Vec x,f;
  VecCreateMPI(PETSC_COMM_WORLD, ctx.nx*ctx.ny, PETSC_DETERMINE, &x);
  VecSet(x, 0.0);
  VecDuplicate(x, &f);

  // Vypocet prave strany
  VecSet(f, 1.0/(ctx.nx*ctx.nx));

  // Vytvoreni matrix-free matice
  Mat A;
  MatCreateShell(
      PETSC_COMM_WORLD, ctx.nx*ctx.ny, ctx.nx*ctx.ny,
      PETSC_DETERMINE, PETSC_DETERMINE, &ctx, &A);

  MatSetUp(A);
  MatShellSetOperation(A, MATOP_MULT, (void(*)(void))myMatMul);

  // Vytvoreni predpodminovace
  Mat P;
  MatCreate(PETSC_COMM_WORLD, &P);
  MatSetSizes(P, ctx.nx*ctx.ny, ctx.nx*ctx.ny, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSetFromOptions(P);
  MatSetUp(P);
  for (size_t i=0; i<ctx.nx*ctx.ny; i++)
      MatSetValue(P,i,i,4.0,INSERT_VALUES);
  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

  // Reseni soustavy s matrix-free matici
  KSP solver;
  KSPCreate(PETSC_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, P);

  KSPSetType(solver, KSPCG);
  PC prec;
  KSPGetPC(solver,&prec);
  PCSetType(prec,PCJACOBI);
  
  KSPSetFromOptions(solver);
  KSPSetUp(solver);

  KSPSolve(solver, f, x);

  std::ofstream of("13_matrix_free.dat");
  const PetscScalar *xp;
  VecGetArrayRead(x, &xp);
  for (int i=0; i<ctx.nx; i++) {
      for (int j=0; j<ctx.ny; j++)
          of << "\t" << xp[i+j*ctx.nx] << std::endl;
      of << std::endl;
  }
  VecRestoreArrayRead(x, &xp);
      
  
  KSPDestroy(&solver);
  MatDestroy(&P);
  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&f);

  CHKERRQ( PetscFinalize() );
  
  return 0;
}
