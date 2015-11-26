#include <stdio.h>
#include <petscvec.h>

int main(int argc,char **args)
{
  
  PetscErrorCode   ierr ;
  int         proc_size, rank, i, istart, iend, n;
  Vec x, y, z ;
  PetscScalar a, xx, yy;
  Vec xloc, yloc, zloc;

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
  VecGetOwnershipRange(x, &istart, &iend);
  for (i=istart; i<iend; i++) VecSetValue(y, i, (double)i, INSERT_VALUES);
  
  /* Provedeni komunikaci */
  VecAssemblyBegin(y); 
  VecAssemblyEnd(y);

  /* Vypocet z = a*x + y */
  VecDuplicate(x, &xloc);
  VecDuplicate(x, &yloc);
  VecDuplicate(x, &zloc);

  VecGetLocalVectorRead(x, xloc);
  VecGetLocalVectorRead(y, yloc);
  VecGetLocalVector(z, zloc);
  
  a = 2;
  VecWAXPY(zloc, a, xloc, yloc);  

  VecRestoreLocalVectorRead(x, xloc);
  VecRestoreLocalVectorRead(y, yloc);
  VecRestoreLocalVector(z, zloc);
  

  PetscPrintf(PETSC_COMM_WORLD, "z = ");
  VecView(z, PETSC_VIEWER_STDOUT_WORLD);
  PetscPrintf(PETSC_COMM_WORLD, "\n");
  

  ierr = PetscFinalize() ;
  CHKERRQ( ierr ) ;       
  
  return 0;
}
