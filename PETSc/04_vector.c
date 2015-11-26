#include <stdio.h>
#include <petscvec.h>

int main(int argc,char **args)
{
  
  PetscErrorCode   ierr ;
  int         proc_size, rank, i, istart, iend, n;
  Vec x, y, z ;
  PetscScalar a;
  const PetscScalar *pz;
  Vec zloc;
  VecScatter ctx;

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
  a = 2;
  VecWAXPY(z, a, x, y);  

  /* Tisk vektoru procesorem c. 0 */
  VecGetSize(z, &n);
  VecCreateSeq(PETSC_COMM_SELF, n, &zloc);
  
  VecScatterCreateToZero(z, &ctx, &zloc);
  VecScatterBegin(ctx, z, zloc, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, z, zloc, INSERT_VALUES, SCATTER_FORWARD);
  
  VecGetArrayRead(zloc, &pz);
  if (rank==0) 
    for (i=0; i<n; i++) printf("z[%i]=%lf\n", i, pz[i]);
  VecRestoreArrayRead(zloc, &pz);

	   
  VecScatterDestroy(&ctx);
  VecDestroy(&zloc);
  VecDestroy(&x);  VecDestroy(&y);  VecDestroy(&z);

  ierr = PetscFinalize() ;
  CHKERRQ( ierr ) ;       
  
  return 0;
}
