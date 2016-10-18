// g++ cblas1.cpp -lcblas -lblas

#include <iostream>
#include <cblas.h>

using namespace std;

class BandedMatrix {
public:

  BandedMatrix(size_t rows, size_t cols, size_t kl, size_t ku):
    rows_(rows),
    cols_(cols),
    kl_(kl),
    ku_(ku)
  { data_ = new double[rows*(kl+ku+1)]; }

  ~BandedMatrix() { delete[] data_; }

  size_t rows() const { return rows_; }
  size_t cols() const { return cols_; }
  size_t kl() const { return kl_; }
  size_t ku() const { return ku_; }
  size_t lda() const { return 1 + kl() + ku(); }
  double* data() { return data_; }

  double* operator[](int row) { return data_ + kl_ + row*lda() - row; }

private:
  size_t rows_;
  size_t cols_;
  size_t kl_;
  size_t ku_;
  double* data_;
};

int main() {
  size_t n = 10;

  double* x = new double[n];
  double* y = new double[n];
  double alpha = 1.0;
  double beta = 0.0;

  BandedMatrix A(n, n, 2, 2);
  
  for (size_t i=0; i<n; i++) {
    for (size_t j=max(0, int(i-2)); j<min(n, i+3); j++) 
      A[i][j] = 1.0;
    A[i][i] = 4;
    x[i] = 1.0;
    y[i] = 3.0;
  }

  // y = alpha * A + beta * y
  cblas_dgbmv(CblasRowMajor, CblasNoTrans,
	      A.rows(), A.cols(), A.kl(), A.ku(), alpha, A.data(), A.lda(),
	      x, 1, 
	      beta, y, 1);

  for (size_t i=0; i<A.rows(); i++) 
    cout << "y[i] = " << y[i] << endl;

  delete[] y;
  delete[] x;

  return 0;
}
