// g++ -g lapack_gesv.cpp -lopenblas
#include <iostream>
#include <lapacke.h>

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
  int n = 10;

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

  int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, A.rows(), 1, A.data(), A.lda(),
			   ipiv, y, 1);

  cout << "info = " << info << endl;

  for (size_t i=0; i<A.rows(); i++) 
    cout << "y[i] = " << y[i] << endl;


  delete[] ipiv;
  delete[] y;

  return 0;
}
