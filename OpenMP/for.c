#include <stdio.h>

#define N 1000

int main() {
    double x[N], y[N], s;
    int i;
    
    for (i=0; i<N; i++) {
	x[i] = i;
	y[i] = i;
    };

    s = 0;
    #pragma omp parallel for reduction(+:s)
    for (i=0; i<N; i++)
	s += x[i]*y[i];

    printf("Dot product = %lf\n", s);

    return 0;
}
