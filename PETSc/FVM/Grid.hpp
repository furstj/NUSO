#pragma once

#include <array>
#include <vector>
#include <map>
#include <boost/any.hpp>

typedef int Label;
typedef std::array<double,2> Point;
typedef std::array<double,2> Vector;

struct Face {
  Label  id;
  Label  owner;
  Label  neighbour;
  Vector s; 
};

struct Cell {
  Label  id;
  double vol;
};


struct Grid {      

public:

  typedef std::vector<Point>       PointList;
  typedef std::vector<Face>        FaceList;
  typedef std::vector<Cell>        CellList;
  typedef std::map<Label,FaceList> PatchList;
  

  Grid(const PointList& points, const FaceList& faces, const CellList& cells);

  virtual ~Grid() {};
  
  const PointList&  points() const { return points_; }
  const CellList&   cells()  const { return cells_; }
  const FaceList&   internalFaces()  const { return internalFaces_; }
  const PatchList&  boundaryPatches() const { return boundaryPatches_; }

  boost::any& userData() { return userData_; }

private:

  PointList points_;
  FaceList  internalFaces_;
  CellList  cells_;
  PatchList boundaryPatches_;
  boost::any userData_;
};
