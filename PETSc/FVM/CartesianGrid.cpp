#include "CartesianGrid.hpp"
#include <iostream>

Grid CartesianGrid(std::size_t Nx, std::size_t Ny) {
  
  double dx = 1.0 / Nx;
  double dy = 1.0 / Ny;

  Label id = 0;
  Grid::PointList pts( (Nx+1)*(Ny+1) );
  for (int j=0; j<Ny+1; j++)
    for (int i=0; i<Nx+1; i++)
      pts[id++] = Point{ i*dx, j*dy };

  id = 0;
  Grid::CellList cells( Nx*Ny );
  for (int j=0; j<Ny; j++)
    for (int i=0; i<Nx; i++) {
      cells[id] = Cell{id, dx*dy};
      id++;
    }
  
  Grid::FaceList faces( (Nx+1)*Ny + Nx*(Ny+1) );
  id = 0;
  for (int j=0; j<Ny; j++) 
    for (int i=0; i<Nx+1; i++) {
      faces[id].id = id;
      faces[id].s = {dy, 0.0};
      faces[id].owner     = (i>0 ? i-1 + j*Nx : BND_LEFT);
      faces[id].neighbour = (i<Nx ? i + j*Nx  : BND_RIGHT);
      id++;
    }
  for (int j=0; j<Ny+1; j++) 
    for (int i=0; i<Nx; i++) {
      faces[id].id = id;
      faces[id].s = {0.0, dx};
      faces[id].owner = (j>0 ? i + (j-1)*Nx : BND_BOTTOM);
      faces[id].neighbour = (j<Ny ? i + j*Nx  : BND_TOP);
      id++;
    }

  return Grid(pts, faces, cells);
}

