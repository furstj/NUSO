#include "Grid.hpp"
#include <map>

Grid::Grid(const Grid::PointList& points, 
	   const Grid::FaceList& faces,
	   const Grid::CellList& cells):
  points_(points),
  internalFaces_(0),
  cells_(cells),
  boundaryPatches_(),
  userData_(0)
{

  auto reversed = [](const Face& f)
    { return Face{f.id, f.neighbour, f.owner, {-f.s[0], -f.s[1]}}; };


  for (auto f : faces) {
    if (f.neighbour < f.owner) f = reversed(f); 

    if (f.owner>=0) 
      internalFaces_.push_back(f);
    else
      boundaryPatches_[f.owner].push_back(reversed(f));
  }
}
