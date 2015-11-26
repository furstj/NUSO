#include <stdio.h>
#include <petscvec.h>

int main(int argc,char **args)
{
  
  PetscErrorCode   ierr ;
  int         proc_size, rank, i, istart, iend, n;
  Vec x, y, z ;
  PetscScalar a;
  PetscScalar *px, *py, *pz;

  /* Inicializace */
  ierr = PetscInitialize( &argc , &args , (char *)0 , 0 ) ;
  CHKERRQ( ierr ) ;       

  /* Informace o poctu procesoru a mem cisle */
  MPI_Comm_size( PETSC_COMM_WORLD , &proc_size );
  MPI_Comm_rank( PETSC_COMM_WORLD , &rank );
  PetscPrintf(PETSC_COMM_WORLD, "Using %d processors.\n", proc_size);

  /* Vytvoreni vektoru a nastaveni jeho celkove velikosti na 10 */
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, 10, &x);

  /* Vytvoreni stejne velkych a stejne distribovanych vektoru */
  VecDuplicate(x, &y);
  VecDuplicate(x, &z);

  /* Nastaveni vsech prvku x na 1.0 */
  VecSet(x, 1.0);

  /* Nastaveni prvku vektoru y na y[i]=i, praci provadi ten, komu dany prvek patri */
  VecGetOwnershipRange(y, &istart, &iend);
  for (i=istart; i<iend; i++) VecSetValue(y, i, (double)i, INSERT_VALUES);
  
  /* Provedeni komunikaci */
  VecAssemblyBegin(y); 
  VecAssemblyEnd(y);

  PetscPrintf(PETSC_COMM_WORLD, "y = ");
  VecView(y, PETSC_VIEWER_STDOUT_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  
  /* Vypocet z = a*x + y */
  VecGetArray(x, &px);
  VecGetArray(y, &py);
  VecGetArray(z, &pz);

  VecGetLocalSize(x, &n); 
  a = 2;
  for (i=0; i<n; i++) 
    pz[i] = a * px[i] + py[i]; 
  
  VecRestoreArray(x, &px); 
  VecRestoreArray(y, &py);
  VecRestoreArray(z, &pz);

  PetscPrintf(PETSC_COMM_WORLD, "z = ");
  VecView(z, PETSC_VIEWER_STDOUT_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  

  ierr = PetscFinalize() ;
  CHKERRQ( ierr ) ;       
  
  return 0;
}
