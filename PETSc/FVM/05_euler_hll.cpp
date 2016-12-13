/*

  Program pro reseni Eulerovych rovnic v oblasti ve tvaru kanalu

  Pridan vypocet Jakobianu

  Tip:
  time mpirun -np 4 ./05_euler_hll -ts_dt 1.0 -ts_monitor 
  time mpirun -np 4 ./05_euler_hll -ts_monitor -ts_type pseudo -ts_pseudo_max_dt 5 
  time mpirun -np 4 ./05_euler_hll -ts_monitor -ts_type pseudo -ts_pseudo_max_dt 1. -snes_mf
*/

#include "ChannelGrid.hpp"
#include "DecomposeGrid.hpp"

#include <set>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
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
  static constexpr size_t size() { return 4; }
  double operator[](int i) const { return ((double*)this)[i]; }
  double& operator[](int i) { return ((double*)this)[i]; }
  const double* begin() const { return &this->rho; }
  const double* end() const { return begin()+size(); }
  double* begin() { return &this->rho; }
  double* end()   { return begin()+size(); }
};

double Model::kappa = 1.4;

double Pressure(const Model& W) {
  return (W.kappa-1) * (W.rhoE - 0.5*(pow(W.rhoU,2)+pow(W.rhoV,2))/W.rho); 
};

double SoundSpeed(const Model& W) {
  return sqrt(W.kappa * Pressure(W) / W.rho);
};

Model& operator+=(Model& me, const Model& other) {
  for (size_t i=0; i<me.size(); i++) 
    me[i] += other[i];
  return me;
};

Model& operator-=(Model& me, const Model& other) {
  for (size_t i=0; i<me.size(); i++) 
    me[i] -= other[i];
  return me;
};

Model& operator/=(Model& me, double other) {
  for (size_t i=0; i<me.size(); i++) 
    me[i] /= other;
  return me;
};

Model& operator*=(Model& me, double other) {
  for (size_t i=0; i<me.size(); i++) 
    me[i] *= other;
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

  int bs = Model::size();
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
  Vec vol;
  std::map<int,BocoFunction> bocos;
};


PetscErrorCode CalculateRHS(TS ts, PetscReal t, Vec u, Vec r, void* ctx) {
  Grid& g = *( static_cast<MyContext*>(ctx) -> gptr );
  BocoMap& bocos = static_cast<MyContext*>(ctx) -> bocos;
  VecGhostUpdateBegin(u, INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(u, INSERT_VALUES,SCATTER_FORWARD);
  CalculateRezidual(g, u, r, bocos);
}


void FluxJacobian(std::function<Model(const Model&)> f, const Model& W, Model Df[] ) {
  Model f0 = f(W);
  const double eps = 1.e-8;
  for (size_t j=0; j<Model::size(); j++) {
    Model Weps(W);
    Weps[j] += eps;
    for (size_t i=0; i<Model::size(); i++) {
      Model fEps = f(Weps);
      Df[i][j] = (fEps[i] - f0[i]) / eps;
    }
  }
}

void UpdateJacobianBlock(Mat J, PetscInt bi, PetscInt bj, Model A[], double factor) {
  const int n = Model::size();
  double aa[n*n];
  int ai[n], aj[n];
  
  int idx = 0;
  for (int i=0; i<Model::size(); i++) 
    for (int j=0; j<Model::size(); j++) 
      aa[idx++] = factor * A[i][j];
  
  for (int i=0; i<Model::size(); i++) {
    ai[i] = bi*n + i;
    aj[i] = bj*n + i;
  }
  MatSetValues(J, n, ai, n, aj, aa, ADD_VALUES);
}


PetscErrorCode  CalculateJacobian(TS ts,PetscReal t,Vec u, Mat J, Mat B, void* ctx) {
  //cout << "CALCULATE_JACOBIAN ...";
  Grid& g = *( static_cast<MyContext*>(ctx) -> gptr );
  auto dist = boost::any_cast<const Decomposition&>(g.userData());
  BocoMap& bocos = static_cast<MyContext*>(ctx) -> bocos;
  Vec vol = static_cast<MyContext*>(ctx) -> vol;

  Mat jac = B;

  MatZeroEntries(jac);

  VecGhostUpdateBegin(u, INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(u, INSERT_VALUES,SCATTER_FORWARD);

  Vec uLoc, volLoc;
  const double* ul;
  VecGhostGetLocalForm(u, &uLoc);  
  VecGetArrayRead(uLoc, &ul);

  const double *volume;
  VecGhostGetLocalForm(vol, &volLoc);  
  VecGetArrayRead(volLoc, &volume);

  const Model *W = (Model*)(ul);
  int blockSize = Model::size();

  for (auto f: g.internalFaces()) {
    Model Al[Model::size()], Ar[Model::size()];
   
    int go = dist.GlobalIndex(f.owner);
    int gn = dist.GlobalIndex(f.neighbour);

    auto FWl = [&](const Model& Wl) { return Flux(Wl, W[f.neighbour], f.s); };
    FluxJacobian(FWl, W[f.owner], Al);

    auto FWr = [&](const Model& Wr) { return Flux(W[f.owner], Wr, f.s); };
    FluxJacobian(FWr, W[f.neighbour], Ar);

    double vol = 1.0/(Ny*Ny);

    UpdateJacobianBlock(jac, go, go, Al, -1/volume[f.owner]);
    UpdateJacobianBlock(jac, go, gn, Ar, -1/volume[f.owner]);

    UpdateJacobianBlock(jac, gn, go, Al, 1/volume[f.neighbour]);
    UpdateJacobianBlock(jac, gn, gn, Ar, 1/volume[f.neighbour]);
    
  }

  for (auto bnd: g.boundaryPatches()) {
    auto bcf = bocos.at(bnd.first);
    for (auto f: bnd.second) {
      Model Al[Model::size()];
      auto FWl = [&](const Model& Wl) { return Flux(Wl, bcf(Wl,f.s), f.s); };
      FluxJacobian(FWl, W[f.owner], Al);
      int go = dist.GlobalIndex(f.owner);
      UpdateJacobianBlock(jac, go, go, Al, -1/volume[f.owner]);
    }
  }
  VecRestoreArrayRead(uLoc, &ul);
  VecGhostRestoreLocalForm(u, &uLoc);

  VecRestoreArrayRead(volLoc, &volume);
  VecGhostRestoreLocalForm(vol, &volLoc);

  MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);

  //MatView(J, PETSC_VIEWER_STDOUT_WORLD);
  //cout << "END" << endl;
  return 0;
}

int main(int argc, char **argv) {
  
  PetscInitialize( &argc, &argv, (char*)0, 0);

  const int modelSize = sizeof(Model) / sizeof(double);

  Grid g = DecomposeGrid( ChannelGrid(Nx, Ny) );

  std::cout << g.cells().size() << std::endl;
  
  Mat J;
  MatCreate(PETSC_COMM_WORLD, &J);
  MatSetSizes(J, g.cells().size()*modelSize, g.cells().size()*modelSize, 
	      PETSC_DECIDE, PETSC_DECIDE);
  MatSetFromOptions(J);
  MatSetUp(J);

  Vec u,r;
  u = CreateGhostedBlockVector(g, sizeof(Model)/sizeof(double));

  VecDuplicate(u, &r);
  double dt = 0.1 / (2.5*Nx/3);
  
  MyContext ctx;
  ctx.gptr = &g;
  ctx.vol = CreateGhostedBlockVector(g, 1);
  // Okrajove podminky
  ctx.bocos[BND_LEFT]   = InletBoco;
  ctx.bocos[BND_RIGHT]  = OutletBoco;
  ctx.bocos[BND_TOP]    = WallBoco;
  ctx.bocos[BND_BOTTOM] = WallBoco;

  // Pocatecni podminka 
  int istart, iend;
  VecGetOwnershipRange(ctx.vol, &istart, &iend);
  for (int i=istart; i<iend; i++) {
    double W0[] = {1.0, 0.1, 0.0, 3.0 };
    VecSetValuesBlocked(u, 1, &i, W0, INSERT_VALUES);
    VecSetValue(ctx.vol, i, g.cells()[i-istart].vol, INSERT_VALUES);
  };
  VecAssemblyBegin(u);
  VecAssemblyEnd(u);
  VecAssemblyBegin(ctx.vol);
  VecAssemblyEnd(ctx.vol);
  VecGhostUpdateBegin(ctx.vol, INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(ctx.vol, INSERT_VALUES,SCATTER_FORWARD);
  
  double t = 0;
  double tEnd = 50.;

  TS ts;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetProblemType(ts, TS_NONLINEAR);
  TSSetSolution(ts, u);
  TSSetRHSFunction(ts, NULL, CalculateRHS, &ctx);
  TSSetRHSJacobian(ts, J, J, CalculateJacobian, &ctx);
  TSSetType(ts, TSEULER);
  TSSetInitialTimeStep(ts, 0.0, dt);
  TSSetDuration(ts, 10000000, tEnd);
  TSSetType(ts, TSBEULER);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
  TSSetFromOptions(ts);
  
  TSSolve(ts, u);

  Save(u, "W.dat");
  
  VecDestroy(&r);
  VecDestroy(&u);
  PetscFinalize();

  return 0;
}
