#include <stdio.h>
#include <petscvec.h>

int main(int argc,char **args)
{
  
  PetscErrorCode   ierr ;
  int         proc_size, rank, i, n;
  Vec x, y, z ;
  PetscScalar a;
  
  /* Inicializace */
  ierr = PetscInitialize( &argc , &args , (char *)0 , 0 ) ;
  CHKERRQ( ierr ) ;       

  /* Informace o poctu procesoru a mem cisle */
  MPI_Comm_size( PETSC_COMM_WORLD , &proc_size );
  MPI_Comm_rank( PETSC_COMM_WORLD , &rank );
  
  //if (rank==0) printf("Using %d processors.\n", proc_size);
  PetscPrintf(PETSC_COMM_WORLD, "Using %d processors.\n", proc_size);


  /* Vytvoreni vektoru a nastaveni jeho celkove velikosti na 10 */
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 10, &x);

  /* Vytvoreni stejne velkych a stejne distribovanych vektoru */
  VecDuplicate(x, &y);
  VecDuplicate(x, &z);

  /* Nastaveni vsech prvku x na 1.0 */
  VecSet(x, 1.0);

  PetscPrintf(PETSC_COMM_WORLD, "x = ");
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  

  /* Nastaveni prvku vektoru y na y[i]=i, praci provadi pouze "master" */
  VecGetSize(y, &n);
  if (rank==0) 
    for (i=0; i<n; i++) VecSetValue(y, i, (double)i, INSERT_VALUES);
  
  /* Provedeni komunikaci */
  VecAssemblyBegin(y); 
  VecAssemblyEnd(y);

  PetscPrintf(PETSC_COMM_WORLD, "y = ");
  VecView(y, PETSC_VIEWER_STDOUT_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  
  /* Vypocet z = a*x + y */
  a = 2;
  VecWAXPY(z, a, x, y);

  PetscPrintf(PETSC_COMM_WORLD, "z = ");
  VecView(z, PETSC_VIEWER_STDOUT_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  

  VecDestroy(&x);
  VecDestroy(&y);
  VecDestroy(&z);

  ierr = PetscFinalize() ;
  CHKERRQ( ierr ) ;       
  
  return 0;
}
