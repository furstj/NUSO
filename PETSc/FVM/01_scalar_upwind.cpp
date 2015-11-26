/*

  Program pro reseni ulohy

  u_t + a u_x + b u_y = 0

  na ctverci [0,1]x[0,1] s nulovou pocatecni podminkou a okrajovymi podminkami 

  u(x,y) = 1 pro y=0
  u(x,y) = 0 pro x=0
*/

#include "CartesianGrid.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

const int N = 50;   // Pocet bunek v jednom smeru
double a = 1;      // Rychlosti v jednotlivych smerech
double b = 2.0/3.0;
const double tEnd = 2;


void Save(const std::vector<double>& u, const char* filename) {
  auto nSqr = u.size();
  int n =  sqrt(nSqr);
  std::ofstream out(filename);
  for (size_t i=0; i<nSqr; i++) {
    out << u[i] << std::endl;
    if ((i+1)%n == 0) out << std::endl;
  }
}


double Flux(double ul, double ur, const Vector& s) {
  double phi = a*s[0] + b*s[1];
  if (phi>0) return phi*ul;
  else return phi*ur;
}


void CalculateRezidual(const Grid& g, const std::vector<double>& u, 
 		       std::vector<double>& r) {
  
  for (auto& ri : r) ri = 0;

  for (auto f: g.internalFaces()) {
    double flux = Flux( u[f.owner], u[f.neighbour], f.s);
    r[f.owner] -= flux;
    r[f.neighbour] += flux;
  }

  for (auto bnd: g.boundaryPatches()) {
    double ub = 0;
    if (bnd.first == BND_BOTTOM) ub = 1;
    
    for (auto f: bnd.second) 
      r[f.owner] -= Flux( u[f.owner], ub, f.s);
  }
  
  for (auto c : g.cells()) r[c.id] /= c.vol;
}


int main(int argc, char **argv) {
  
  Grid g = CartesianGrid(N, N);

  std::vector<double> u(g.cells().size());
  std::vector<double> r(g.cells().size());
  
  double dt = 1.0 / (N*(fabs(a)+fabs(b)));
  
  // Pocatecni podminka
  for (auto& ui : u) ui = 0;
  
  double t = 0;
  
  while (t<tEnd)  {
    
    CalculateRezidual(g, u, r);

    for (size_t i = 0; i<u.size(); i++)
      u[i] += dt * r[i];
    
    t += dt;
  };

  Save(u, "u.dat");

   return 0;
}

  
