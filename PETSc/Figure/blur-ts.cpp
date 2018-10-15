#include <iostream>
#include <fstream>

#include <petscdmda.h>
#include <petscts.h>


/***********************************************************************
 *  Rozostreni obrazku ziskane resenim rovnice
 *
 *  u_t = u_xx + u_yy
 *
 *
 *  blur_ts -ts_type euler -ts_dt 0.1 -ts_final_time 10
 *
 *  beuler stabilni pro dt 0.5
 *  euler stabilni pro dt 0.25
 ***********************************************************************/

const int width  = 256;
const int height = 256;

int max(int a, int b) {
  if (a>b) return a;
  else return b;
}

int min(int a, int b) {
  if (a<b) return a;
  else return b;
}

// Context 
struct AppCtx {
  DM da;
  Vec local;

  AppCtx() {
    // Vytvoreni distribuovaneho obrazku
    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
             DMDA_STENCIL_STAR,
	     width, height, PETSC_DECIDE, PETSC_DECIDE, 
	     1, 1, PETSC_NULL, PETSC_NULL, &da);
    DMSetUp(da);
    
    // Vytvoreni lokalniho pracovniho pole 
    DMCreateLocalVector(da, &local);    
  }

};


// Aproximace laplacianu na siti
PetscErrorCode FormFunction(TS ts, PetscReal t, Vec in, Vec out, 
			    void *ctx) {
  AppCtx *appCtx = (AppCtx*) ctx;
  DM da = appCtx->da;
  Vec l = appCtx->local;

  VecZeroEntries(out);

  // Presunuti do lokalniho vektoru s hranicnimi body
  DMGlobalToLocalBegin(da, in, INSERT_VALUES, l);
  DMGlobalToLocalEnd(da, in, INSERT_VALUES, l);

  // Zjisteni indexu me casti poli
  int i1, j1, n, m;
  DMDAGetCorners(da, &i1, &j1, PETSC_NULL, &n, &m, PETSC_NULL);

  // Vypocet prave strany
  double **u, **un;

  DMDAVecGetArray(da, l, &u);
  DMDAVecGetArray(da, out, &un);

  for (int i=max(i1,1); i<min(i1+n,width-1); i++)
    for (int j=max(j1,1); j<min(j1+m,height-1); j++)
      un[j][i] = u[j+1][i]+u[j-1][i]+u[j][i+1]+u[j][i-1]-4*u[j][i]; 

  DMDAVecRestoreArray(da,out,&un);
  DMDAVecRestoreArray(da,l,&u);

  return 0;
}

// Vypocet priblizneho Jakobianu prave strany
PetscErrorCode FormJacobian(TS ts,PetscReal t,Vec in,Mat J,Mat B, void *ctx) {
  int i1, i2;
  MatGetOwnershipRange(B, &i1, &i2);
  for (int i=i1; i<i2; i++) {
    MatSetValue(B, i, i, -4, INSERT_VALUES);
    if (i+1<i2) MatSetValue(B, i, i+1, 1, INSERT_VALUES);
    if (i-1>=i1) MatSetValue(B, i, i-1, 1, INSERT_VALUES);
  }
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

  if (J != B) {
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
  }

  return 0;
}



void InitialCondition(DM da, Vec g) {
  double **u;
  int i1, j1, n, m;

  DMDAVecGetArray(da, g, &u);
  DMDAGetCorners(da, &i1, &j1, PETSC_NULL, &n, &m, PETSC_NULL);
  for (int i=i1; i<i1+n; i++)
    for (int j=j1; j<j1+m; j++)
      {
	double x = (i-0.5*width) * 2.0 / width;
	double y = (j-0.5*height) * 2.0 / height;
	if (x*x + y*y < 1.0/16.) u[j][i] = 0;
	else u[j][i] = 1;
      }

  DMDAVecRestoreArray(da, g, &u);
}



using namespace std;

int main(int argc, char** argv) {

  PetscInitialize( &argc, &argv, (char*)0, 0);

  // Kontext aplikace (DA, pracovni pole, ...)
  AppCtx appCtx;

  // Globalni vektor na ulozeni obrazku
  Vec g;
  DMCreateGlobalVector(appCtx.da, &g);

  // Nastaveni pocatecnich hodnot obrazku
  InitialCondition(appCtx.da, g);

  // Vytvoreni resice pro "nelinearni" ulohu
  TS ts;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetProblemType(ts, TS_NONLINEAR);

  // Nastaveni funkce pro vypocet prave strany
  TSSetRHSFunction(ts, NULL, FormFunction, (void*)&appCtx);

  // Nastaveni funkce pro vypocet Jakobiho matice prave strany
  Mat J; 
  MatCreate(PETSC_COMM_WORLD, &J);
  MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, width*height, width*height);
  MatSetUp(J);
  MatSetFromOptions(J);

  TSSetRHSJacobian(ts, J, J, FormJacobian, &appCtx);

  // Nastaveni parametru resice
  double dt = 0.1, t_end=2;
  
  TSSetType(ts, TSBEULER);

  TSSetTime(ts, 0.0);
  TSSetTimeStep(ts, dt);
  TSSetMaxSteps(ts, 1000);
  TSSetMaxTime(ts, t_end);
  TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);

  // Nacteni parametru z prikazove radky
  TSSetFromOptions(ts);

  // Nastaveni pocatecni podminky a vyreseni
  TSSolve(ts, g);

  // Ulozeni obrazku do pgm formatu
  // Preposlani obrazku na procesor 0
  {
    Vec nat;
    DMDACreateNaturalVector(appCtx.da, &nat);
    DMDAGlobalToNaturalBegin(appCtx.da, g, INSERT_VALUES, nat);
    DMDAGlobalToNaturalEnd(appCtx.da, g, INSERT_VALUES, nat);
    
    Vec u_loc;
    VecCreateSeq(PETSC_COMM_SELF, width*height, &u_loc);
    
    VecScatter ctx;
    VecScatterCreateToZero(nat,&ctx,&u_loc);
    VecScatterBegin(ctx, nat, u_loc, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, nat, u_loc, INSERT_VALUES, SCATTER_FORWARD);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (rank==0) {
      double* image;
      VecGetArray(u_loc, &image);
      ofstream pgm("obrazek.pgm");
      pgm << "P2" << endl;
      pgm << "# output of blur.cpp" << endl;
      pgm << width << " " << height << endl;
      pgm << 256 << endl;
      for (int i=0; i<width*height; i++)
	pgm << (int)(256*image[i]) << endl;
      VecRestoreArray(u_loc, &image);
      
    }
    
  }

  PetscFinalize();  
  return 0;
}
