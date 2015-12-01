#include "Grid.hpp"

#include <petscvec.h>

struct Decomposition {
  PetscInt cellStart;
  PetscInt cellEnd;
  ISLocalToGlobalMapping locToGlobMap;
  Decomposition(PetscInt istart, PetscInt iend, const Grid::FaceList& faces);
  ~Decomposition() {};
};

Grid DecomposeGrid(const Grid& g);

