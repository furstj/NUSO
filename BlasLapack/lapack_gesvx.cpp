// g++ -g lapack_gesv.cpp -lopenblas
#include <iostream>
#include <lapacke.h>
#include <cmath>

using namespace std;

class Matrix {
public:

  Matrix(size_t rows, size_t cols):
    rows_(rows),
    cols_(cols)
  { data_ = new double[rows*cols]; }

  ~Matrix() { delete[] data_; }

  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  size_t lda() const { return cols_; }

  double* data() { return data_; }

  double* operator[](int row) { return data_ + row*lda(); }

private:
  size_t rows_;
  size_t cols_;
  double* data_;
};


int main() {
  int n = 1000;

  double* x = new double[n];
  double* y = new double[n];
  double alpha = 1.0;
  double beta = 1.0;

  Matrix A(n, n);
  
  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<n; j++) 
      A[i][j] = 1.0;
    A[i][i] = n+1;
    y[i] = 2.0;
  }
    
  // x = A \ y
  int* ipiv = new int[n];
  char equed[2] = "N"; 
  double r[n], c[n], rcond, ferr, berr, rpivot;

  int info = LAPACKE_dgesvx(LAPACK_ROW_MAJOR, 'E', 'N', 
			    A.rows(), 1, A.data(), A.lda(),
			    A.data(), A.lda(),
			    ipiv, equed, r, c,
			    y, 1, 
			    x, 1,
			    &rcond, &ferr, &berr, &rpivot);


  cout << "info = " << info << endl;
  cout << "ferr = " << ferr << endl;
  cout << "berr = " << berr << endl;
  cout << "rcond= " << rcond << endl;
  cout << endl;

  int i = 0;
  //for (size_t i=0; i<A.rows(); i++) 
    cout << "x[i] = " << x[i] << endl;

  double err = 0.0;
  for (size_t i=0; i<A.rows(); i++) 
    err += fabs( x[i] - 1.0/n );
  cout << "err  = " << err << endl;
  
    
  delete[] ipiv;
  delete[] y;
  delete[] x;

  return 0;
}
