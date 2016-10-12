#include <iostream>
#include <chrono>
#include<cmath>

using namespace std;

double** allocate_matrix(int n, int m) {
  double** A = new double*[n];
  for (int i=0; i<n; i++)
    A[i] = new double[m];
  return A;
}

int main() {
  const int lda = 1001;
  const int n = 1000;
  
  auto a = allocate_matrix(n, n);
  auto b = allocate_matrix(n, n);
  auto c = allocate_matrix(n, n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      a[i][j] = i+j;
      b[i][j] = i-j;
    };

  auto t1 = chrono::high_resolution_clock::now();
  
    for (int j=0; j<n; j++) 
  for (int i=0; i<n; i++)
      {
	c[i][j] = 0.0;
	for (int k=0; k<n; k++)
	  c[i][j] += a[i][k] * b[k][j];
      }
  
  auto t2 = chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_time = t2 - t1;

  double mflops = 2*pow(double(n),3)/elapsed_time.count()/1000000.0;

  cout << "time   = " << elapsed_time.count() << endl;
  cout << "mflops = " << mflops << endl;

  return 0;
}
