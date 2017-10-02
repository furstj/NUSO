// g++ cblas1.cpp -lcblas -lblas

#include <iostream>
#include <cblas.h>

using namespace std;

int main() {
  int n = 1000;

  double* x = new double[n];
  double* y = new double[n];
  double a = 3;

  for (int i=0; i<n; i++) {
    y[i] = 1;
    x[i] = 2;
  }

  // y = y + a*x
  cblas_daxpy(n, a, x, 1, y, 1);
  cout << "y[0] = " << y[0] << endl;

  // sum(y)
  double s = cblas_dasum(n, y, 1);
  cout << "sum(y) = " << s << endl;

  // ||y||
  double nrm2 = cblas_dnrm2(n, y, 1);
  cout << "||y|| = " << nrm2 << endl;

  // (x.y)
  double dot = cblas_ddot(n, x, 1, y, 1);
  cout << "(x.y) = " << dot << endl;

  // y *= 1/2
  cblas_dscal(n, 0.5, y, 1);
  cout << "y[0] = " << y[0] << endl;
 
  delete[] y;
  delete[] x;

  return 0;
}
