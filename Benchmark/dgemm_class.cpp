#include <iostream>
#include <chrono>
#include <cmath>

using namespace std;

class Matrix {
public:
  Matrix() = delete;
  Matrix(int n, int m): 
    n_(n), m_(m)
  {
    values_ = new double[n*m];
  }
  ~Matrix() 
  {
    delete[] values_;
  }

  inline double* operator[](int i) { return &(values_[i*n_]); }
  //inline double& operator()(int i, int j) { return values_[i*n_+j]; }

private:
  int n_;
  int m_;
  double* values_;
};


int main() {
  const int n = 1000;
  //const int n = 1024;

  auto a = Matrix(n, n);
  auto b = Matrix(n, n);
  auto c = Matrix(n, n);

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
