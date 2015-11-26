#include <iostream>
#include <petscvec.h>
#include <petscao.h>

const int n = 10;

int main(int argc,char **args)
{
  CHKERRQ( PetscInitialize( &argc , &args , (char *)0 , 0 ) );
  
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Vytvoreni vektoru a nastaveni jeho celkove velikosti na 10
  Vec x;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &x);

  // Precislovani vektoru tak, ze v petsc vektoru jsou data
  // permutovana takto: [0, 2, 4, 6, 8, 1, 3, 5, 7, 9]
  AO ao;
  {
    int i1, i2;
    VecGetOwnershipRange(x, &i1, &i2);
    int l = i2 - i1;
    int a[l], p[l];
    for (int i=i1; i<i2; i++) {
      p[i-i1] = i;
      a[i-i1] = i / 2 + n/2 * (i%2);
    }
    AOCreateBasic(PETSC_COMM_WORLD, l, a, p, &ao);  
  }
  
  AOView(ao, PETSC_VIEWER_STDOUT_WORLD);

  // Nastaveni vektoru pomoci VecSetValue
  {
    int i1, i2;
    VecGetOwnershipRange(x, &i1, &i2);
    for (int i=i1; i<i2; i++) {
      int ip = i;
      AOApplicationToPetsc(ao, 1, &ip);
      VecSetValue(x, ip, i, INSERT_VALUES);
    }
  }
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  VecDestroy(&x);

  CHKERRQ( PetscFinalize() );
  
  return 0;
}
