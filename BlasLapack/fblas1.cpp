// g++ cblas1.cpp -lblas

#include <iostream>
#include <openblas/f77blas.h>

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
  int one = 1;
  daxpy_(&n, &a, x, &one, y, &one);

  cout << "y[0] = " << y[0] << endl;

  // sum(y)
  double s = dasum_(&n, y, &one);
  cout << "sum(y) = " << s << endl;

  // ||y||
  double nrm2 = dnrm2_(&n, y, &one);
  cout << "||y|| = " << nrm2 << endl;

  // (x.y)
  double dot = ddot_(&n, x, &one, y, &one);
  cout << "(x.y) = " << dot << endl;

  // y *= 1/2
  double coeff = 0.5;
  dscal_(&n, &coeff, y, &one);
  cout << "y[0] = " << y[0] << endl;
 
  delete[] y;
  delete[] x;

  return 0;
}
