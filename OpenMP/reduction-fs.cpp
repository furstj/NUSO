#include <iostream>
#include <omp.h>

using namespace std;

int main() {

  double a[100];
  double* mm;

  int nthreads = omp_get_max_threads();

  mm = new double[nthreads*32];

  for (int i=0; i<100; i++)
    a[i] = i;

  for (int i=0; i<nthreads; i++) 
    mm[i] = a[0];


#pragma omp parallel 
  {
    int myid = 32*omp_get_thread_num();

    #pragma omp for
    for (int i=0; i<100; i++)
      mm[myid] = max(mm[myid], a[i]);
  }

  double m = mm[0];
  for (int i=0; i<nthreads; i++) 
    m = max(mm[32*i],m);

  cout << "max=" << m << endl;

  return 0;
}
