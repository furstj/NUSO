/*
  Reseni soustavy ODR

  dk/dt = - beta_star * k * omega
  domega/dt = -beta * omega^2

  s pocatecni podminkou k(0) = 1, omega(0) = 1000

  Pozn.:
  priklad je pouze pro beh na jednom procesoru

  Tipy:
  -ts_monitor_lg_solution -draw_pause 0.1
  -ts_monitor -ts_type rk

 */

#include <petscts.h>

struct MyContext {
  double beta_star;
  double beta;
};

PetscErrorCode CalculateFunc(TS ts, PetscReal t, Vec x, Vec f, void* ctx) {
  const PetscScalar *xl;
  PetscScalar *fl;
  VecGetArrayRead(x, &xl);
  VecGetArray(f, &fl);

  double beta_star = static_cast<MyContext*>(ctx) -> beta_star;
  double beta      = static_cast<MyContext*>(ctx) -> beta;

  fl[0] = - beta_star * xl[0] * xl[1];
  fl[1] = -beta       * xl[1] * xl[1];

  VecRestoreArrayRead(x, &xl);
  VecRestoreArray(f, &fl);

  return 0;
}



int main(int argc, char** argv) {
  PetscInitialize( &argc, &argv, (char *)0, 0 ) ;

  int ncpus;
  MPI_Comm_size(PETSC_COMM_WORLD,&ncpus);
  if (ncpus > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,
			"Example is only for sequential runs");

  MyContext ctx{0.09, 0.075}; 

  Vec x;
  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, 2);
  VecSetFromOptions(x);

  VecSetValue(x, 0, 1.0, INSERT_VALUES);
  VecSetValue(x, 1, 100.0, INSERT_VALUES);
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  TS ts;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetProblemType(ts, TS_NONLINEAR);
  
  TSSetRHSFunction(ts, NULL, CalculateFunc, &ctx);

  TSSetSolution(ts, x);  // Pocatecni podminka
  TSSetTimeStep(ts, 0.01);
  TSSetDuration(ts, 1000000, 1.0);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
  TSSetFromOptions(ts);

  TSSolve(ts, x);

  VecView(x,PETSC_VIEWER_STDOUT_WORLD);

  TSView(ts,PETSC_VIEWER_STDOUT_WORLD);

  TSDestroy(&ts);
  VecDestroy(&x);

  PetscFinalize();
  return 0;
}


