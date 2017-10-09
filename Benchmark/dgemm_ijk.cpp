#include <iostream>
#include <chrono>
#include<cmath>

using namespace std;

int main() {
  //const int lda = 1024;
  const int lda = 1001;
  const int n = 1000;
  
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
  
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) 
      for (int k=0; k<n; k++)
	c[i+lda*j] += a[i+lda*k] * b[k+lda*j];
    

  auto t2 = chrono::system_clock::now();

  std::chrono::duration<double> elapsed_time = t2 - t1;

  double mflops = 2*pow(double(n),3)/elapsed_time.count()/1000000.0;

  cout << "time   = " << elapsed_time.count() << endl;
  cout << "mflops = " << mflops << endl;

  return 0;
}
