#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>

#include <mkl_dss.h>
#include <mkl_spblas.h>

using namespace std;

struct SparseMatrix {
  int n, m;               // Rozmery matice
  int nnz;                // Pocet nenulovych clenu
  int* cptr;
  int* cidx;
  double* val; 

  ~SparseMatrix() {
    delete[] cptr;
    delete[] cidx;
    delete[] val;
  }
};

// Nacteni ridke matice ze souboru ve formatu MatrixMarket
void loadMatrixMM(const char* filename, SparseMatrix& A) {
  string line;
  ifstream file(filename);
  if (!file.good()) {cerr << "Chyba\n"; exit(1);};

  do {                    // Preskocime hlavicku
    getline(file, line);
  } while (line[0]=='%');

  // Nacti rozmery matice
  stringstream(line) >> A.n >> A.m >> A.nnz;
  
  // Nacti matici ve formatu i, j, A(i,j) s upravou indexovani od 0
  int *ii = new int[A.nnz];
  int *jj = new int[A.nnz];
  double *AA = new double[A.nnz];
  for (int l = 0; l<A.nnz; l++) {
    file >> ii[l] >> jj[l] >> AA[l];
    ii[l]--;
    jj[l]--;
  };

  // Preved matici na CSR format pomoci funkci MKL
  A.cptr = new int[A.n+1];
  A.cidx = new int[A.nnz];
  A.val  = new double[A.nnz];

  MKL_INT job[6] = {2, 0, 0, 0, A.nnz, 0};
  MKL_INT info;
  mkl_dcsrcoo(job, &A.n, A.val, A.cidx, A.cptr, &A.nnz, AA, ii, jj, &info);
  cout << info << endl;;
  delete[] ii;
  delete[] jj;
  delete[] AA;
}

//======================================================================

int main() {
  SparseMatrix A;
  loadMatrixMM("poisson3Da.mtx",A); // Nacteni matice

  double *b = new double[A.n];    // Prava strana
  double *x = new double[A.n];    // Reseni

  // Vypocet prave strany (reseni je vektor x = [1, 1, 1, ... 1]
  for (int i=0; i<A.n; i++) b[i]=0;
  for (int j=0; j<A.m; j++)
    for (int k=A.cptr[j]; k<A.cptr[j+1]; k++)
      b[j] += A.val[A.cidx[k]];

  _INTEGER_t opt = MKL_DSS_DEFAULTS + MKL_DSS_ZERO_BASED_INDEXING;
  _MKL_DSS_HANDLE_t handle;

  dss_create_(&handle, &opt);

  opt = MKL_DSS_NON_SYMMETRIC;
  dss_define_structure_(&handle, &opt, A.cptr, &A.n, &A.m, A.cidx, &A.nnz);

  _INTEGER_t perm[A.n];
  opt = MKL_DSS_AUTO_ORDER;
  dss_reorder_(&handle, &opt, perm); 

  opt = MKL_DSS_INDEFINITE;
  dss_factor_real_(&handle, &opt, A.val);
  
  opt = 0;
  int one = 1;
  dss_solve_real_(&handle, &opt, b, &one, x);

  double stat[8];
  dss_statistics_(&handle, &opt, "FactorTime", stat);
  cout << "Time " << stat[0] << endl;

  dss_delete_(&handle, &opt);

  delete[] b;
  delete[] x;

  return 0;
}
