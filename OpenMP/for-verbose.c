#include <stdio.h>

#include <omp.h>

#define N 20

int main() {
    double x[N], y[N], s;
    int i,myproc;
    
    for (i=0; i<N; i++) {
	x[i] = i;
	y[i] = i;
    };

    s = 0;
    #pragma omp parallel for reduction(+:s) private(myproc) schedule(dynamic,2)
    for (i=0; i<N; i++) {
	s += x[i]*y[i];
	myproc=omp_get_thread_num();
	printf("myporc=%i,  i=%i, s=%lf\n", myproc, i, s);
    }

    printf("Dot product = %lf\n", s);

    return 0;
}
