#include <stdio.h>

#define N 1000

int main() {
    double *a;
    double *b;
    double *c;
    double s, t;
    
    int i, j, k;

    a = (double*) malloc(N*N*sizeof(double));
    b = (double*) malloc(N*N*sizeof(double));
    c = (double*) malloc(N*N*sizeof(double));

    // Inicializace matic a, b
    for (i=0; i<N; i++)
	for (j=0; j<N; j++) {
	    a[i+j*N] = i+j;
	    b[i+j*N] = i-j;
	}

    // Vypocet soucinu c = a*b
#pragma omp parallel for private(j,k,s)
    for (i=0; i<N; i++)
	for (j=0; j<N; j++) {
	    s = 0;
	    for (k=0; k<N; k++) s += a[i+k*N]*b[k+j*N];
	    c[i+j*N] = s;
	};

    // Vypocet stopy soucinu (ma vyjit 0)
    t = 0;
    for (i=0; i<N; i++) t += c[i+i*N];

    printf("Stopa soucinu = %lf\n", t);

}
