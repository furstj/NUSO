#include <iostream>
#include <cmath>
#include <chrono>
#include <cblas.h>

using namespace std;


int main() {
  const int lda = 1000;
  const int n = 1000;
  
  double* a = new double[lda*lda];
  double* b = new double[lda*lda];
  double* c = new double[lda*lda];

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      a[i+lda*j] = i+j;
      b[i+lda*j] = i-j;
    };

  auto t1 = chrono::high_resolution_clock::now();
  
  //one = 1.0;
  //zero = 0.0;
  //dgemm_("N", "N", &n, &n, &n, &one, a, &lda, b, &lda, &zero, c, &lda);

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
  	      n, n, n, (double)1, a, lda, b, lda, (double)0, c, lda);

  auto t2 = chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_time = t2 - t1;

  double mflops = 2*pow(double(n),3)/elapsed_time.count()/1000000.0;

  cout << "time   = " << elapsed_time.count() << endl;
  cout << "mflops = " << mflops << endl;

  return 0;
}
