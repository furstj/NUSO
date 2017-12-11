#include <iostream>
#include <fstream>
#include <string>
#include <valarray>
#include <iomanip>

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

  nacti_sit("ctverec.2");

  int n = node.size();

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
  VecSetSizes(b,PETSC_DECIDE, n);
  VecSetFromOptions(b);

  // Sestaveni matice a prave strany 
  cout << "Sestavuji matici" << endl;

  int i1, i2; // Tyto radky patri mne
  MatGetOwnershipRange(A, &i1, &i2);

  for (size_t e=0; e<elem.size(); ++e) {   // Sestaveni bez ohledu na okraje
    Point p1 = node[elem[e].id[0]];
    Point p2 = node[elem[e].id[1]];
    Point p3 = node[elem[e].id[2]];
    double vol = volume(p1, p2, p3);
    
    for (int ii=0; ii<3; ++ii) {       // Doplnim radek i
      int i = elem[e].id[ii];
      if (i1<=i && i<i2) {
	
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

  // Rozeslani matice 
  cout << "Komunikace" << endl;
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // Mam soustavu, ted uz ji jen vyresim
  KSP ksp;
  Vec x;
  VecDuplicate(b, &x);
    
  cout << "Reseni" << endl;
  KSPCreate( PETSC_COMM_WORLD , &ksp ) ;
  KSPSetOperators( ksp , A , A ) ;
  KSPSetFromOptions( ksp ) ;
  KSPSolve( ksp, b , x ) ;

  // Ted musim z PETSc vymamit reseni, udelam to tak, ze vektor necham poslat na rank==0
  cout << "Ulozeni vysledku" << endl;
  Vec x_loc;
  VecCreateSeq(PETSC_COMM_SELF, n, &x_loc);

  VecScatter ctx;
  VecScatterCreateToZero(x,&ctx,&x_loc);
  VecScatterBegin(ctx,x,x_loc,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x,x_loc,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  if (rank==0) {
    double *solution;
    VecGetArray(x_loc, &solution);

    // Vypis reseni ve formatu vtk 
    ofstream vtk("solution.vtk");
    vtk << "# vtk DataFile Version 3.1" << endl;
    vtk << "Solution to Poisson equation" << endl;
    vtk << "ASCII" << endl;
    vtk << "DATASET UNSTRUCTURED_GRID" << endl;
    vtk << endl;
    vtk << scientific << setprecision(6);

    vtk << "POINTS " << node.size() << " FLOAT" << endl;
    for (int i=0; i<n; i++) 
      vtk << node[i].x << "\t" << node[i].y << "\t" << 0.0 << endl;
    vtk << "CELLS " << elem.size() << " " << 4*elem.size() << endl;
    for (size_t j=0; j<elem.size(); j++) 
      vtk << "3" << "\t" << elem[j].id[0] << "\t" << elem[j].id[1] 
	  << "\t" << elem[j].id[2] << endl;
    vtk << "CELL_TYPES " << elem.size() << endl;
    for (size_t j=0; j<elem.size(); j++) 
      vtk << "5" << endl;

    vtk << "POINT_DATA " << n<< endl << endl;
    vtk << "SCALARS U DOUBLE" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for (int i=0; i<n; ++i)
      vtk << solution[i] << endl;
    
    VecRestoreArray(x_loc, &solution);
  }

  PetscFinalize();
  return 0;
}


    
