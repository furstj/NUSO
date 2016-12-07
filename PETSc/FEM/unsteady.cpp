#include <iostream>
#include <fstream>
#include <string>
#include <valarray>
#include <iomanip>
#include <cassert>
#include <petscksp.h>
#include <petscts.h>

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

// Vypocet \int_T grad phi_i phi_j 
double intGradGrad(const double phi_i[], const double phi_j[], 
		   const Point& a, const Point& b, const Point& c) 
{
  double vol = volume(a, b, c);
  assert(vol>0);
  Point grad_i = grad(a, phi_i[0], b, phi_i[1], c, phi_i[2]);
  Point grad_j = grad(a, phi_j[0], b, phi_j[1], c, phi_j[2]);
  return (grad_i.x*grad_j.x + grad_i.y*grad_j.y)*vol;
}

double intPhiPhi(const double phi_i[], const double phi_j[], 
		   const Point& a, const Point& b, const Point& c) 
{
  // 3-point integration
  double beta1[3] = {2./3., 1./3.,  1./3.};
  double beta2[3] = {1./3., 2./3.,  1./3.};
  double beta3[3] = {1./3., 1./3.,  2./3.};
  double w[3] = {1./3., 1./3., 1./3. };
  double f;
  f = 0;
  for (int i=0; i<3; i++) {
    double p_i = phi_i[0]*beta1[i]+phi_i[1]*beta2[i]+phi_i[2]*beta3[i];
    double p_j = phi_j[0]*beta1[i]+phi_j[1]*beta2[i]+phi_j[2]*beta3[i];
    f += w[i]*p_i*p_j;
  }
  return f*volume(a,b,c);
}

// ======================================================================


int main(int argc, char **argv) {

  int rank;
  
  PetscInitialize( &argc, &argv, (char *)0, 0 );
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //nacti_sit("ctverec.2");
  nacti_sit("ctverec.1");

  int n = node.size();
  // Alokace matice tuhosti
  Mat A;
  //MatCreateAIJ(PETSC_COMM_WORLD, 
  //	       PETSC_DECIDE, PETSC_DECIDE, n, n, 
  //	       0, PETSC_NULL, 0, PETSC_NULL, &A);
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetUp(A);
  MatSetFromOptions(A);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // Alokace matice hmotnosti
  Mat B;
  MatCreate(PETSC_COMM_WORLD, &B);
  MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &B);

  // Alokace reseni
  Vec u;
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u,PETSC_DECIDE, n);
  VecSetFromOptions(u);

  // Sestaveni matic a pocatecni podminky
  cout << "Sestavuji matici" << endl;
  
  int i1, i2; // Tyto radky patri mne
  MatGetOwnershipRange(A, &i1, &i2);

  for (size_t e=0; e<elem.size(); ++e) {   // Sestaveni matic po elementech
    Point p1 = node[elem[e].id[0]];
    Point p2 = node[elem[e].id[1]];
    Point p3 = node[elem[e].id[2]];
    
    for (int ii=0; ii<3; ++ii) {       // Doplnim radek i
      int i = elem[e].id[ii];
      if (i1<=i && i<i2) {
	
	if (flag[i]==0) {              // Vrchol i je vnitrni bod
	  double phi_i[3] = {0, 0, 0};
	  phi_i[ii] = 1;
	  
	  for (int jj=0; jj<3; ++jj) {
	    int j = elem[e].id[jj];
	    if (flag[j]==0) {
	      double phi_j[3] = {0, 0, 0};
	      phi_j[jj] = 1;

	      double aij = - intGradGrad(phi_i, phi_j, p1, p2, p3);
	      double bij = intPhiPhi(phi_i, phi_j, p1, p2, p3);

	      MatSetValue(A, i, j, aij, ADD_VALUES);
	      MatSetValue(B, i, j, bij, ADD_VALUES);
	    }
	  }
	}
	else {                      // Vrchol i je hranicni bod
	  MatSetValue(B,i,i, 1.0, ADD_VALUES);
	  MatSetValue(A,i,i, 0.0, ADD_VALUES);
	}
      }
    }
  }


  // Rozeslani matice 
  cout << "Zacatek rozesilani matic" << endl;
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);

  cout << "Vypocet pocatecni podminky" << endl;
  for (int i=i1; i<i2; i++) {
    if (pow(node[i].x-0.5,2)+pow(node[i].y-0.5,2) < 0.125) 
      VecSetValue(u, i, 1.0, INSERT_VALUES);
    else
      VecSetValue(u, i, 0.0, INSERT_VALUES);
  }

  cout << "Zacatek rozesilani pocatecni podminky" << endl;
  VecAssemblyBegin(u);

  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  cout << "Dokonceno rozesilani matic" << endl;

  VecAssemblyEnd(u);
  cout << "Dokonceno rozesilani pocatecni podminky" << endl;

  /* ==================================================================
     Ted mam sestavenou ulohu
     
     B u_t = A u

     s pocatecni podminkou ulozenou v u

     Tuto ulohu budu resit pomoci modulu TS
     ==================================================================*/
  TS ts;
  TSCreate(PETSC_COMM_WORLD, &ts);
  TSSetProblemType(ts, TS_LINEAR);

  //TSSetMatrices(ts,A, PETSC_NULL, B, PETSC_NULL, 
  //		SAME_NONZERO_PATTERN, PETSC_NULL);
  //TSSetMatrices(ts,A, PETSC_NULL, PETSC_NULL, PETSC_NULL, 
  //		SAME_NONZERO_PATTERN, PETSC_NULL);
  TSSetRHSFunction(ts, NULL, TSComputeRHSFunctionLinear, NULL);
  TSSetRHSJacobian(ts, A, A, TSComputeRHSJacobianConstant, NULL);

  TSSetIFunction(ts,NULL,TSComputeIFunctionLinear,NULL);
  TSSetIJacobian(ts,B,B,TSComputeIJacobianConstant, NULL);


  double dt = 1.e-2;
  double t_end = 1.0;
  TSSetInitialTimeStep(ts, 0.0, dt);
  TSSetDuration(ts, 1000000, t_end);
  TSSetType(ts, TSBEULER);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
  TSSetFromOptions(ts);

  cout << "Reseni ODE" << endl; 
  TSSolve(ts, u);


  // Ted musim z PETSc vymamit reseni, udelam to tak, ze vektor necham poslat na rank==0
  Vec x_loc;
  VecCreateSeq(PETSC_COMM_SELF, n, &x_loc);

  VecScatter ctx;
  VecScatterCreateToZero(u,&ctx,&x_loc);
  VecScatterBegin(ctx,u,x_loc,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,u,x_loc,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  if (rank==0) {
    double *solution;
    VecGetArray(x_loc, &solution);
    
    // Vypis reseni ve formatu vtk 
    ofstream vtk("unsteady.vtk");
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


    
