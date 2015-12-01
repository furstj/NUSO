#include "DecomposeGrid.hpp"
#include <set>

Decomposition::Decomposition(PetscInt istart, PetscInt iend, const Grid::FaceList& faces):
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
