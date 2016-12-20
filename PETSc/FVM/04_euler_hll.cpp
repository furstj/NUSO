/*

  Program pro reseni Eulerovych rovnic v oblasti ve tvaru kanalu

  Tip:
  time mpirun -np 4 ./04_euler_hll -ts_monitor -ts_final_time 50 -ts_type pseudo -snes_mf -ts_pseudo_max_dt 1.
*/

#include "ChannelGrid.hpp"
#include "DecomposeGrid.hpp"

#include <set>
#include <vector>
#include <cmath>
#include <functional>
#include <iostream>
#include <fstream>
#include <petscvec.h>
#include <petscts.h>

using namespace std;

const int Ny = 50;   // Pocety bunek v jednotlivych smerech
const int Nx = 3 * Ny;

struct Model {
  double rho;
  double rhoU;
  double rhoV;
  double rhoE;
  static double kappa;
};

double Model::kappa = 1.4;

double Pressure(const Model& W) {
  return (W.kappa-1) * (W.rhoE - 0.5*(pow(W.rhoU,2)+pow(W.rhoV,2))/W.rho); 
};

double SoundSpeed(const Model& W) {
  return sqrt(W.kappa * Pressure(W) / W.rho);
};

Model& operator+=(Model& me, const Model& other) {
  me.rho += other.rho;
  me.rhoU += other.rhoU;
  me.rhoV += other.rhoV;
  me.rhoE += other.rhoE;
  return me;
};

Model& operator-=(Model& me, const Model& other) {
  me.rho -= other.rho;
  me.rhoU -= other.rhoU;
  me.rhoV -= other.rhoV;
  me.rhoE -= other.rhoE;
  return me;
};

Model& operator/=(Model& me, double other) {
  me.rho /= other;
  me.rhoU /= other;
  me.rhoV /= other;
  me.rhoE /= other;
  return me;
};

Model& operator*=(Model& me, double other) {
  me.rho *= other;
  me.rhoU *= other;
  me.rhoV *= other;
  me.rhoE *= other;
  return me;
};

Model operator+(const Model& first, const Model& second) {
  Model m(first);
  return m+=second;
};

Model operator-(const Model& first, const Model& second) {
  Model m(first);
  return m-=second;
};

Model operator*(double a, const Model& first) {
  Model m(first);
  return m*=a;
};

Model operator/(const Model& first, double a) {
  Model m(first);
  return m/=a;
};

typedef std::function<Model(const Model&, const Vector&)> BocoFunction;
typedef std::map<int,BocoFunction> BocoMap;

Model WallBoco(const Model& W, const Vector& s) {
  double q = (W.rhoU*s[0] + W.rhoV*s[1]) / (s[0]*s[0] + s[1]*s[1]);
  return Model{ W.rho, W.rhoU - 2*q*s[0], W.rhoV - 2*q*s[1], W.rhoE };
}

Model OutletBoco(const Model& W, const Vector& s) {
  const double pOut = 0.737;
  double p2 = Pressure(W);
  double rhoE = W.rhoE + (pOut - p2) / W.kappa;
  return Model{ W.rho, W.rhoU, W.rhoV, rhoE };
}

Model InletBoco(const Model& W, const Vector& s) {
  const double p0 = 1;
  const double rho0 = 1;

  double p = min(Pressure(W), p0);
  double rho = rho0 * pow(p/p0, 1/W.kappa);
  double a2 = W.kappa * p / rho;
  double M2 = 2/(W.kappa-1) * ( (W.kappa * p0 / rho0) / a2 - 1 );
  double U = sqrt(M2*a2);
  double rhoE = p/(W.kappa-1) + 0.5*rho*U*U;
  return Model{ rho, rho*U, 0.0, rhoE };
}


double VelocityMag(const Model& W) {
  return sqrt( pow(W.rhoU,2) + pow(W.rhoV,2) ) / W.rho;
};







Vec CreateGhostedBlockVector(Grid& g, int blockSize) {
  auto dist = boost::any_cast<const Decomposition&>(g.userData());

  int n;
  ISLocalToGlobalMappingGetSize(dist.locToGlobMap, &n);
  const PetscInt *idx;
  ISLocalToGlobalMappingGetIndices(dist.locToGlobMap, &idx);

  Vec v;
  PetscInt nLoc = dist.cellEnd - dist.cellStart;
  std::cout << n-nLoc << std::endl;

  VecCreateGhostBlock(PETSC_COMM_WORLD, blockSize, g.cells().size()*blockSize, PETSC_DECIDE, 
  		 n - nLoc, idx + nLoc, &v);
  ISLocalToGlobalMappingRestoreIndices(dist.locToGlobMap, &idx);

  return v;
}


// TODO: zapisuji jen hodnotu v bunkach a ne jejich souradnice
void Save(Vec u, const char* filename) {
  int n, rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  VecGetSize(u, &n);

  Vec uloc;
  VecCreateSeq(PETSC_COMM_SELF, n, &uloc);

  VecScatter ctx;
  VecScatterCreateToZero(u, &ctx, &uloc);
  VecScatterBegin(ctx, u, uloc, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, u, uloc, INSERT_VALUES, SCATTER_FORWARD);

  int bs = sizeof(Model)/sizeof(double);
  cout << "bs = " << bs << endl;
  const PetscScalar *pu;
  VecGetArrayRead(uloc, &pu);
  if (rank==0) {
    std::ofstream out(filename);
    for (size_t i=0; i<n; i+=bs) {
      for (size_t ii=i; ii<i+bs; ii++) out << pu[ii] << "\t";
      out << std::endl;
      if ( (i+bs) % (Nx*bs) == 0) out << std::endl;
    }
  }
  VecRestoreArrayRead(uloc, &pu);
  
  VecScatterDestroy(&ctx);
  VecDestroy(&uloc);
}


Model Flux(const Model& Wl, const Model& Wr, const Vector& s) {
  double sMag = sqrt(s[0]*s[0] + s[1]*s[1]);
  double pl = Pressure(Wl);
  double pr = Pressure(Wr);
  double ql = (Wl.rhoU*s[0] + Wl.rhoV*s[1]) / Wl.rho;
  double qr = (Wr.rhoU*s[0] + Wr.rhoV*s[1]) / Wr.rho;

  Model Fl{Wl.rho * ql, Wl.rhoU*ql + pl*s[0], Wl.rhoV*ql + pl*s[1], (Wl.rhoE + pl) * ql};
  Model Fr{Wr.rho * qr, Wr.rhoU*qr + pl*s[0], Wr.rhoV*qr + pr*s[1], (Wr.rhoE + pr) * qr};

  double sr = max( ql + sMag*SoundSpeed(Wl), qr + sMag*SoundSpeed(Wr) ); 
  double sl = min( ql - sMag*SoundSpeed(Wl), qr - sMag*SoundSpeed(Wr) ); 

  if (sl > 0) return Fl;
  else if (sr < 0) return Fr;
  
  return (sr*Fl - sl*Fr + sl*sr*(Wr-Wl)) / (sr-sl);
};

/**/
void CalculateRezidual(Grid& g, Vec u, 
		       Vec r, const BocoMap& bocos) {
  Vec uLoc, rLoc;
  const double *ul;
  double *rl;
  VecGhostGetLocalForm(u, &uLoc);  
  VecGetArrayRead(uLoc, &ul);

  VecGhostGetLocalForm(r, &rLoc);
  VecSet(rLoc, 0.0);
  VecGetArray(rLoc, &rl);

  const Model *W = (Model*)(ul);
  Model *R = (Model*)(rl);
  int blockSize = sizeof(Model)/sizeof(double);

  for (auto f: g.internalFaces()) {
    auto flux = Flux( W[f.owner], W[f.neighbour], f.s);
    R[f.owner] -= flux;
    R[f.neighbour] += flux;
  }

  for (auto bnd: g.boundaryPatches()) {
    auto bcf = bocos.at(bnd.first);
    for (auto f: bnd.second) {
      R[f.owner] -= Flux( W[f.owner], bcf(W[f.owner],f.s), f.s);
    }
  }

  VecRestoreArray(rLoc, &rl);
  VecGhostRestoreLocalForm(r, &rLoc);
  VecRestoreArrayRead(uLoc, &ul);
  VecGhostRestoreLocalForm(u, &uLoc);

  VecGhostUpdateBegin(r,ADD_VALUES,SCATTER_REVERSE);
  VecGhostUpdateEnd(r,ADD_VALUES,SCATTER_REVERSE);
  
  VecGetArray(r, &rl);
  R = (Model*)(rl);
  for (size_t i=0; i<g.cells().size(); i++) R[i] /= g.cells()[i].vol;
  VecRestoreArray(r, &rl);
}

struct MyContext {
  Grid* gptr;
  std::map<int,BocoFunction> bocos;
};


PetscErrorCode CalculateRHS(TS ts, PetscReal t, Vec u, Vec r, void* ctx) {
  Grid& g = *( static_cast<MyContext*>(ctx) -> gptr );
  BocoMap& bocos = static_cast<MyContext*>(ctx) -> bocos;
  VecGhostUpdateBegin(u, INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(u, INSERT_VALUES,SCATTER_FORWARD);
  CalculateRezidual(g, u, r, bocos);
}


int main(int argc, char **argv) {
  
  PetscInitialize( &argc, &argv, (char*)0, 0);

  Grid g = DecomposeGrid( ChannelGrid(Nx, Ny) );

  std::cout << g.cells().size() << std::endl;

  Vec u,r;
  u = CreateGhostedBlockVector(g, sizeof(Model)/sizeof(double));

  VecDuplicate(u, &r);
  double dt = 0.1 / (2.5*Nx/3);
  
  // Pocatecni podminka (TODO: nastavit jen mou lokalni cast);
  for (int i=0; i<Nx*Ny; i++) {
    double W0[] = {1.0, 0.1, 0.0, 3.0 };
    VecSetValuesBlocked(u, 1, &i, W0, INSERT_VALUES);
  };
  VecAssemblyBegin(u);
  VecAssemblyEnd(u);
  
  MyContext ctx;
  ctx.gptr = &g;
  // Okrajove podminky
  ctx.bocos[BND_LEFT]   = InletBoco;
  ctx.bocos[BND_RIGHT]  = OutletBoco;
  ctx.bocos[BND_TOP]    = WallBoco;
  ctx.bocos[BND_BOTTOM] = WallBoco;

  double t = 0;
  double tEnd = 30.;

  TS ts;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetProblemType(ts, TS_NONLINEAR);
  TSSetSolution(ts, u);
  TSSetRHSFunction(ts, NULL, CalculateRHS, &ctx);
  TSSetType(ts, TSEULER);
  TSSetInitialTimeStep(ts, 0.0, dt);
  TSSetDuration(ts, 10000000, tEnd);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
  TSSetFromOptions(ts);
  
  TSSolve(ts, u);

  Save(u, "W.dat");
  
  VecDestroy(&r);
  VecDestroy(&u);
  PetscFinalize();

  return 0;
}
