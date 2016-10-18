// g++ cblas1.cpp -lcblas -lblas

#include <iostream>
#include <cblas.h>

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

  double alpha = 1.0;
  double beta = 1.0;

  Matrix A(n, n);
  Matrix B(n, n);
  Matrix C(n, n);

  for (size_t i=0; i<n; i++) {
    for (size_t j=0; j<n; j++) 
      A[i][j] = 1.0;
    A[i][i] = n+0.5;

    for (size_t j=0; j<n; j++) {
      B[i][j] = A[i][j];
      C[i][j] = 0;
    }
  }
    
  // C = alpha * A * B + beta * C
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
	      A.rows(), A.cols(), B.rows(), alpha, A.data(), A.lda(),
	      B.data(), B.lda(), 
	      beta, C.data(), C.lda());

  for (size_t i=0; i<C.rows(); i++) {
    for (size_t j=0; j<C.cols(); j++) 
      cout << C[i][j] << " ";
    cout << endl;
  }

  return 0;
}
