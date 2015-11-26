/*
  Reseni nelinearni soustavy rovnic

  x^2 + y^2 - 1 = 0
  sin(x) - y    = 0
  
  hledame koren s x>0

  Tipy:
  -snes_fd
  -snes_mf

  Pozn.:
  priklad je pouze pro beh na jednom procesoru

 */

#include <petscsnes.h>

PetscErrorCode CalculateFunc(SNES snes, Vec x, Vec f, void* ctx) {
  PetscPrintf(PETSC_COMM_WORLD, "CalculateFunc\n");

  const PetscScalar *xl;
  PetscScalar *fl;
  VecGetArrayRead(x, &xl);
  VecGetArray(f, &fl);

  fl[0] = xl[0]*xl[0] + xl[1]*xl[1] - 1;
  fl[1] = sin(xl[0])  - xl[1];

  VecRestoreArrayRead(x, &xl);
  VecRestoreArray(f, &fl);

  return 0;
}


PetscErrorCode CalculateJacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  PetscPrintf(PETSC_COMM_WORLD, "CalculateJacobian\n");
  int ncpus;
  MPI_Comm_size(PETSC_COMM_WORLD,&ncpus);
  if (ncpus > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
			"Example is only for sequential runs");
  const PetscScalar *xl;
  VecGetArrayRead(x, &xl);

  MatSetValue(P, 0, 0, 2*xl[0], INSERT_VALUES);
  MatSetValue(P, 0, 1, 2*xl[1], INSERT_VALUES);
  MatSetValue(P, 1, 0, cos(xl[0]), INSERT_VALUES);
  MatSetValue(P, 1, 1, -1.0, INSERT_VALUES);

  MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

  VecRestoreArrayRead(x, &xl);

  return 0;
}


int main(int argc, char** argv) {
  PetscInitialize( &argc, &argv, (char *)0, 0 ) ;

  int ncpus;
  MPI_Comm_size(PETSC_COMM_WORLD,&ncpus);
  if (ncpus > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
			"Example is only for sequential runs");

  Vec x, r;
  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, 2);
  VecSetFromOptions(x);
  VecDuplicate(x, &r);

  Mat J;
  MatCreate(PETSC_COMM_WORLD, &J);
  MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, 2, 2);
  MatSetFromOptions(J);
  MatSetUp(J);


  SNES snes;
  SNESCreate(PETSC_COMM_WORLD, &snes);

  SNESSetFunction(snes, r, CalculateFunc, NULL);
  SNESSetJacobian(snes, J, J, CalculateJacobian, NULL);
  SNESSetFromOptions(snes);

  VecSet(x, 1.0);
  SNESSolve(snes, NULL, x);
  VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  SNESGetFunction(snes,&r,0,0);
  VecView(r,PETSC_VIEWER_STDOUT_WORLD);

  SNESDestroy(&snes);
  MatDestroy(&J);
  VecDestroy(&r);
  VecDestroy(&x);

  PetscFinalize();
  return 0;
}


