/*

  Program pro reseni ulohy

  u_t + a u_x + b u_y = 0

  na ctverci [0,1]x[0,1] s nulovou pocatecni podminkou a okrajovymi podminkami 

  u(x,y) = 1 pro y=0
  u(x,y) = 0 pro x=0

  Tipy:
  

  03_scalar_upwind_ts -ts_monitor
  03_scalar_upwind_ts -ts_monitor -ts_type rk
  03_scalar_upwind_ts -ts_monitor -ts_type beuler -snes_mf -ts_dt 0.2 -ts_final_time 2

*/

#include "CartesianGrid.hpp"

#include <set>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <petscvec.h>
#include <petscts.h>

const int N = 50;   // Pocet bunek v jednom smeru
double a = 1;      // Rychlosti v jednotlivych smerech
double b = 2.0/3.0;
const double tEnd = 2.0/3.0;


struct Decomposition {
  PetscInt cellStart;
  PetscInt cellEnd;
  ISLocalToGlobalMapping locToGlobMap;

  Decomposition(PetscInt istart, PetscInt iend, const Grid::FaceList& faces):
    cellStart(istart),
    cellEnd(iend)
  {
    std::set<PetscInt> ghosts;
    for (auto f: faces)
      if ((f.neighbour < istart || iend <= f.neighbour) && f.neighbour >= 0) 
	ghosts.insert(f.neighbour);

    std::vector<PetscInt> glob(cellEnd - cellStart + ghosts.size());

    int id=0;
    for (int i=cellStart; i<cellEnd; i++) glob[id++] = i;
    for (auto g: ghosts) glob[id++] = g;

    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, 1, glob.size(), glob.data(),
				 PETSC_COPY_VALUES, &locToGlobMap);
  }

  ~Decomposition() {}
};

/* 
   Rozdeleni site mezi procesory takto:
   1) vrcholy jsou nakopirovany na vsechny procesory
   2) bunky jsou rozdeleny rovnmerne dle petsc
   3) steny jsou rozdeleny tak, ze cpu(f) = cpu(cell(f.owner))

   Indexy bunek u sten jsou zmeneny na lokalni cislovani
 */
Grid DecomposeGrid(const Grid& g) {
  PetscInt istart, iend;
  Vec tmp;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, g.cells().size(), &tmp);
  VecGetOwnershipRange(tmp, &istart, &iend);
  VecDestroy(&tmp);

  Grid::CellList cells(0);
  for (auto c: g.cells())
    if (istart <= c.id && c.id < iend) cells.push_back(c);
    
  Grid::FaceList faces;
  for (auto f: g.internalFaces())
    if (istart <= f.owner && f.owner < iend) faces.push_back(f);

  for (auto bnd: g.boundaryPatches())
    for (auto f: bnd.second)
      if (istart <= f.owner && f.owner < iend) faces.push_back(f);


  Decomposition dist(istart, iend, faces);
  for (auto& f: faces) {
    int glob[1], loc[1], one;
    glob[0] = f.owner;
    ISGlobalToLocalMappingApply(dist.locToGlobMap, IS_GTOLM_MASK, 1, glob, &one, loc);
    f.owner = loc[0];

    if (f.neighbour >=0) {
      glob[0] = f.neighbour;
      ISGlobalToLocalMappingApply(dist.locToGlobMap, IS_GTOLM_MASK, 1, glob, &one, loc);
      f.neighbour = loc[0];
    }
  }
  
  Grid dg(g.points(), faces, cells);
  
  dg.userData() = dist;

  return dg;
}


Vec CreateGhostedVector(Grid& g) {
  auto dist = boost::any_cast<const Decomposition&>(g.userData());

  int n;
  ISLocalToGlobalMappingGetSize(dist.locToGlobMap, &n);
  const PetscInt *idx;
  ISLocalToGlobalMappingGetIndices(dist.locToGlobMap, &idx);

  Vec v;
  PetscInt nLoc = dist.cellEnd - dist.cellStart;
  std::cout << n-nLoc << std::endl;

  VecCreateGhost(PETSC_COMM_WORLD, g.cells().size(), PETSC_DECIDE, 
  		 n - nLoc, idx + nLoc, &v);
  ISLocalToGlobalMappingRestoreIndices(dist.locToGlobMap, &idx);

  return v;
}

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

  const PetscScalar *pu;
  VecGetArrayRead(uloc, &pu);
  if (rank==0) {
    int sq =  sqrt(n);
    std::ofstream out(filename);
    for (size_t i=0; i<n; i++) {
      out << pu[i] << std::endl;
      if ((i+1)%sq == 0) out << std::endl;
    }
  }
  VecRestoreArrayRead(uloc, &pu);
  
  VecScatterDestroy(&ctx);
  VecDestroy(&uloc);
}


double Flux(double ul, double ur, const Vector& s) {
  double phi = a*s[0] + b*s[1];
  if (phi>0) return phi*ul;
  else return phi*ur;
}


/**/
void CalculateRezidual(Grid& g, Vec u, 
		       Vec r) {
  Vec uLoc, rLoc;
  const double *ul;
  double *rl;
  VecGhostGetLocalForm(u, &uLoc);  
  VecGetArrayRead(uLoc, &ul);

  VecGhostGetLocalForm(r, &rLoc);
  VecSet(rLoc, 0.0);
  VecGetArray(rLoc, &rl);


  for (auto f: g.internalFaces()) {
    double flux = Flux( ul[f.owner], ul[f.neighbour], f.s);
    rl[f.owner] -= flux;
    rl[f.neighbour] += flux;
  }

  for (auto bnd: g.boundaryPatches()) {
    double ub = 0;
    if (bnd.first == BND_BOTTOM) ub = 1;
    
    for (auto f: bnd.second) {
      rl[f.owner] -= Flux( ul[f.owner], ub, f.s);
    }
  }


  VecRestoreArray(rLoc, &rl);
  VecGhostRestoreLocalForm(r, &rLoc);
  VecRestoreArrayRead(uLoc, &ul);
  VecGhostRestoreLocalForm(u, &uLoc);

  VecGhostUpdateBegin(r,ADD_VALUES,SCATTER_REVERSE);
  VecGhostUpdateEnd(r,ADD_VALUES,SCATTER_REVERSE);
  
  
  VecGetArray(r, &rl);
  for (size_t i=0; i<g.cells().size(); i++) rl[i] /= g.cells()[i].vol;
  VecRestoreArray(r, &rl);
}

struct MyContext {
  Grid* gptr;
};

PetscErrorCode CalculateRHS(TS ts, PetscReal t, Vec u, Vec r, void* ctx) {
  Grid& g = *( static_cast<MyContext*>(ctx) -> gptr );
  VecGhostUpdateBegin(u, INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(u, INSERT_VALUES,SCATTER_FORWARD);
  CalculateRezidual(g, u, r);
}


int main(int argc, char **argv) {
  
  PetscInitialize( &argc, &argv, (char*)0, 0);

  Grid g = DecomposeGrid( CartesianGrid(N, N) );

  std::cout << g.cells().size() << std::endl;

  Vec u,r;
  u = CreateGhostedVector(g);

  VecDuplicate(u, &r);
  double dt = 1.0 / (N*(fabs(a)+fabs(b)));
  
  // Pocatecni podminka
  VecSet(u, 0.0);
  
  MyContext ctx;
  ctx.gptr = &g;
 
  double t = 0;
  
  TS ts;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetProblemType(ts, TS_NONLINEAR);
  TSSetSolution(ts, u);
  TSSetRHSFunction(ts, NULL, CalculateRHS, &ctx);
  TSSetType(ts, TSEULER);
  TSSetInitialTimeStep(ts, 0.0, dt);
  TSSetDuration(ts, 10000000, tEnd);
  TSSetFromOptions(ts);
  
  TSSolve(ts, u);

  Save(u, "u.dat");
  
  VecDestroy(&r);
  VecDestroy(&u);
  PetscFinalize();

  return 0;
}
