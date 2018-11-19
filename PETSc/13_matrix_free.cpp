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

struct MyContext
{
    std::size_t nx;
    std::size_t ny;
};

PetscErrorCode myMatMul(Mat mfMatrix, Vec x, Vec y)
{
    MyContext* ctx;
    MatShellGetContext(mvMatrix, (void**) &ctx);
    size_t nx = ctx->nx;
    size_t ny = ctx->ny;

    PetsScalar *xp, *yp;
    VecGetArray(x, &xp);
    VecGetArray(y, &yp);

    // Internal points:
    for (size_t i=1; i<nx-1; i++)
        for (size_t j=1; i<ny-1; j++)
            yp[i+j*nx] = 4*xp[i+j*nx]
                - xp[i+1+j*nx] - xp[i-1+j*nx] - xp[i+(j+1)*nx] - xp[i+(j-1)*nx];

    // Boundary points
    for (size_t i=0; i<nx; i++) {
        yp[i+0*nx] = 0.0;
        yp[i+(ny-1)*nx] = 0.0;
    }
    for (size_t j=0; j<ny; j++) {
        yp[0+j*nx] = 0.0;
        yp[nx-1+j*nx] = 0.0;
    }

    VecRestoreArray(x, &xp);
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
    
  Vec x;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ctx->nx*ctx->ny, &x);
  
  /* Nastaveni vsech prvku x na 1.0 */
  VecSet(x, 1.0);

  // Vytvoreni matice o velkosti 10x10
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, 10, 10);
  MatSetFromOptions(A);
  MatSetUp(A);

  // Nastaveni prvku na 1, -2, 1
  PetscInt istart, iend;
  MatGetOwnershipRange(A, &istart, &iend);
  for (int i=istart; i<iend; i++) {
    if (i>0) MatSetValue(A, i, i-1, 1, INSERT_VALUES);
    MatSetValue(A, i, i, -2.0, INSERT_VALUES);
    if (i<9) MatSetValue(A, i, i+1, 1, INSERT_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  // Vypocet y=A*x
  Vec y;
  VecDuplicate(x, &y);
  MatMult(A, x, y);
  VecView(y, PETSC_VIEWER_STDOUT_WORLD);
 
  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&y);

  CHKERRQ( PetscFinalize() );
  
  return 0;
}
