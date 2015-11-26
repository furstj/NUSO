#include <iostream>

#include <petscvec.h>
#include <petscmat.h>

int main(int argc,char **args)
{
  int         proc_size, rank;
  
  /* Inicializace */
  CHKERRQ( PetscInitialize( &argc , &args , (char *)0 , 0 ) );

  /* Informace o poctu procesoru a mem cisle */
  MPI_Comm_size( PETSC_COMM_WORLD , &proc_size );
  MPI_Comm_rank( PETSC_COMM_WORLD , &rank );
  
  /* Vytvoreni vektoru a nastaveni jeho celkove velikosti na 10 */
  Vec x;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 10, &x);

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
