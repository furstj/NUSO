// FEM 
// sestaveni matice na prvnim procesoru, reseni soustavy paralelne

#include <iostream>
#include <fstream>
#include <string>
#include <valarray>

#include <petscksp.h>

using namespace std;

struct Point {
  double x;          // Souradnice vrcholu
  double y;
};

struct Element {
  int id[3];         // Indexy 3 vrcholu tvorici trojuhelnik
};

// ======================================================================

valarray<Point>    node;  // Uzly site
valarray<Element>  elem;  // Trojuhelniky
valarray<int>      flag;  // Znacky vrcholu

// ======================================================================

void nacti_sit(const string& jmeno) {
  int itmp;
    
  // Nacteni souradnic a znacek vrcholu
  {
    string jmeno_souboru = jmeno+".node";
    ifstream input(jmeno_souboru.c_str());
    int nnodes;
    input >> nnodes >> itmp >> itmp >> itmp;
    cout << "Nacitam " << nnodes << " uzlu" << endl;

    node.resize(nnodes);
    flag.resize(nnodes);
    for (int i=0; i<nnodes; i++)
      input >> itmp >> node[i].x >> node[i].y >> flag[i];
  }

  // Nacteni elementu
  {
    string jmeno_souboru = jmeno+".ele";
    ifstream input(jmeno_souboru.c_str());
    int nelems;
    input >> nelems >> itmp >> itmp;
    cout << "Nacitam " << nelems << " trojuhelniku" << endl;
    elem.resize(nelems);
    for (int i=0; i<nelems; i++) {
      int i0, i1, i2;
      input >> itmp >> i0 >> i1 >> i2; 
      elem[i].id[0] = i0-1;
      elem[i].id[1] = i1-1;
      elem[i].id[2] = i2-1;
    }
	
  }
}

// Vypocet plochy trojuhelnika 
double volume(const Point& a, const Point& b, const Point& c) {
  return ( (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x) ) / 2.0 ;
}

// Vypocet gradientu linearni funkce dane hodnotami ve trech bodech
Point grad(const Point& a, double ua, 
	   const Point& b, double ub, 
	   const Point& c, double uc ) {
  Point grad;
  double detA = 2*volume(a,b,c);
    
  grad.x = ( (c.y-a.y)*(ub-ua) - (b.y-a.y)*(uc-ua) ) / detA;
  grad.y = (-(c.x-a.x)*(ub-ua) + (b.x-a.x)*(uc-ua) ) / detA;

  return grad;
}

// ======================================================================


int main(int argc, char **argv) {

  int rank;
  
  PetscInitialize( &argc, &argv, (char *)0, 0 );
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  if (rank==0) 
    nacti_sit("ctverec.2");

  int n = node.size();
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  valarray<double> u(node.size());
  valarray<double> f(node.size());

  // Prava strana rovnice f = -1
  for (size_t i=0; i<f.size(); ++i) f[i] = -1.0;

  // Alokace matice a prave strany
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetUp(A);
  MatSetFromOptions(A);

  Vec b;
  VecCreate(PETSC_COMM_WORLD, &b);
  VecSetFromOptions(b);
  VecSetSizes(b, PETSC_DECIDE, n);

  if (rank==0) {

    for (size_t e=0; e<elem.size(); ++e) {   // Sestaveni matice
      Point p1 = node[elem[e].id[0]];
      Point p2 = node[elem[e].id[1]];
      Point p3 = node[elem[e].id[2]];
      double vol = volume(p1, p2, p3);
      
      for (int ii=0; ii<3; ++ii) {       // Doplnim radek i
	int i = elem[e].id[ii];
	
	if (flag[i]==0) {              // Vrchol i je vnitrni bod
	  VecSetValue(b, i, f[i]*vol/3, ADD_VALUES);
	  
	  double phi[3] = {0, 0, 0};
	  phi[ii] = 1;
	  Point gradI = grad(p1,phi[0], p2,phi[1], p3,phi[2]);
	  
	  for (int jj=0; jj<3; ++jj) {
	    int j = elem[e].id[jj];
	    if (flag[j]==0) {
	      double phi[3] = {0, 0, 0};
	      phi[jj] = 1;
	      
	      Point gradJ = grad(p1,phi[0], p2,phi[1], p3,phi[2]);
	      double prod = vol*(gradI.x*gradJ.x+gradI.y*gradJ.y);
	      MatSetValue(A, i, j, -prod, ADD_VALUES);
	    }
	  }
	}
	else {                      // Vrchol i je hranicni bod
	  VecSetValue(b, i, 0.0, ADD_VALUES);
	  MatSetValue(A,i,i,1.0,ADD_VALUES);
	}
      }
    }
  }
   
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  
  // Mam soustavu, ted uz ji jen vyresim
  KSP ksp;
  Vec x;
  VecDuplicate(b, &x);
    
  KSPCreate( PETSC_COMM_WORLD , &ksp ) ;
  KSPSetOperators( ksp , A , A );
  KSPSolve( ksp, b , x ) ;

  // Ted musim z PETSc vymamit reseni (scatter to zero)
  
  PetscFinalize();
  return 0;
}


    
