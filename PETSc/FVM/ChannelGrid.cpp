#include "ChannelGrid.hpp"
#include <iostream>

Grid ChannelGrid(std::size_t Nx, std::size_t Ny) {
  
  double dx = 3.0 / Nx;
  double dy = 1.0 / Ny;

  Label id = 0;
  Grid::PointList pts( (Nx+1)*(Ny+1) );
  for (int j=0; j<Ny+1; j++)
    for (int i=0; i<Nx+1; i++) {
      double x = i*dx;
      double y = j*dy;
      double R = 1.3;
      if ( 1.0 < x &&  x < 2.0) 
	y = y * 1 + (1-y) * (sqrt(pow(R,2) - pow(x-1.5,2)) - 1.2);
      pts[id++] = Point{ x, y };
    }

  id = 0;
  Grid::CellList cells( Nx*Ny );
  for (int j=0; j<Ny; j++)
    for (int i=0; i<Nx; i++) {
      auto A = pts[i + j*(Nx+1)];
      auto B = pts[i+1 + j*(Nx+1)];
      auto C = pts[i+1 + (j+1)*(Nx+1)];
      auto D = pts[i + (j+1)*(Nx+1)];
      Vector e { C[0]-A[0], C[1]-A[1] };
      Vector f { D[0]-B[0], D[1]-B[1] };
      double v = (e[0]*f[1] - e[1]*f[0])/2;
      cells[id] = Cell{id, v};
      id++;
    }
  
  Grid::FaceList faces( (Nx+1)*Ny + Nx*(Ny+1) );
  id = 0;
  for (int j=0; j<Ny; j++) 
    for (int i=0; i<Nx+1; i++) {
      auto A = pts[i + j*(Nx+1)];
      auto B = pts[i + (j+1)*(Nx+1)];
      faces[id].id = id;
      faces[id].s = {B[1]-A[1], A[0]-B[0]};
      faces[id].owner     = (i>0 ? i-1 + j*Nx : BND_LEFT);
      faces[id].neighbour = (i<Nx ? i + j*Nx  : BND_RIGHT);
      id++;
    }
  for (int j=0; j<Ny+1; j++) 
    for (int i=0; i<Nx; i++) {
      auto A = pts[i + j*(Nx+1)];
      auto B = pts[i+1 + j*(Nx+1)];
      faces[id].id = id;
      faces[id].s = {A[1]-B[1], B[0]-A[0]};
      faces[id].owner = (j>0 ? i + (j-1)*Nx : BND_BOTTOM);
      faces[id].neighbour = (j<Ny ? i + j*Nx  : BND_TOP);
      id++;
    }

  return Grid(pts, faces, cells);
}

