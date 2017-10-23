#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  double a[4];
  int idx[4] = {0, 2, 1, 3};
  int i;

#pragma omp parallel for shared(a,idx)
  for (i=0; i<4; i++)
    a[idx[i]] = i;

  return 0;
}
