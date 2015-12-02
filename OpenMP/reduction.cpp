#include <iostream>

using namespace std;

int main() {

  double a[100];
  double m;

  for (int i=0; i<100; i++)
    a[i] = i;

  m = a[0];
#pragma omp parallel for reduction(max:m)
  for (int i=0; i<100; i++)
    m = max(m, a[i]);

  cout << "max=" << m << endl;

  return 0;
}
