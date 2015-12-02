#include "Grid.hpp"

#include <petscvec.h>

struct Decomposition {
  PetscInt cellStart;
  PetscInt cellEnd;
  ISLocalToGlobalMapping locToGlobMap;

  Decomposition(PetscInt istart, PetscInt iend, const Grid::FaceList& faces);
  ~Decomposition() {};
  
  PetscInt GlobalIndex(PetscInt local) const;
};

Grid DecomposeGrid(const Grid& g);

