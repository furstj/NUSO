#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

int main() {
  //const int lda = 1024;
  const int lda = 1001;
  const int n = 1000;
  const int bsize = 51;
  
  double* a = new double[lda*lda];
  double* b = new double[lda*lda];
  double* c = new double[lda*lda];

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      a[i+lda*j] = i+j;
      b[i+lda*j] = i-j;
      c[i+lda*j] = 0.0;
    };

  auto t1 = chrono::system_clock::now();

  for (int k=0; k<n; k+=bsize) {
    int kend = min(k+bsize, n);

    for (int j=0; j<n; j+=bsize) {
      int jend = min(j+bsize, n);
      
      for (int i=0; i<n; i+=bsize) {
	int iend = min(i+bsize, n);

	
	for (int kk=k; kk<kend; kk++)
	  for (int jj=j; jj<jend; jj++)
	    for (int ii=i; ii<iend; ii++)	 
	      c[ii+lda*jj] += a[ii+lda*kk] * b[kk+lda*jj];
      }
    }
  }

  auto t2 = chrono::system_clock::now();

  std::chrono::duration<double> elapsed_time = t2 - t1;

  double mflops = 2*pow(double(n),3)/elapsed_time.count()/1000000.0;

  cout << "time   = " << elapsed_time.count() << endl;
  cout << "mflops = " << mflops << endl;

  return 0;
}
