/**********************************************************************
 * Program na reseni Laplaceovy rovnice 
 *
 *  u_xx + u_yy = f(x,y)
 *  
 * s okrajovou podminkou u(x,y) = 0 
 * na oblasti ve tvaru jednotkoveho ctverce.
 *
 * Pro reseni vyuzivam metody siti a vysledne rovnice resim
 * Gauss-Seidelvou iteracni metodou. 
 *
 **********************************************************************/

#include <stdio.h>
#include <math.h>


// Prava strana Laplaceovy rovnice
// tato prava strana odpovida reseni u = ((x-1/2)**2-1/4)*((y-1/2)**2-1/4)
// a tedy u(0.5,0.5) = -1/16 = - 0.625
double f(double x, double y) {
    return 2*(y*y-y+x*x-x);
}


int main() {

    // double *u;
    // double *b;
    static double u[201*201], b[201*201];
    double x, y, h, uc;

    int i, j, n, iter, niter;

    n = 201;

    // u = (double*) malloc(n*n*sizeof(double));
    // b = (double*) malloc(n*n*sizeof(double));

    // Pocatecni hodnota reseni
    for (i=0; i<n; i++)
	for (j=0; j<n; j++)
	    u[i+j*n] = 0;

    // Vypocet prave strany
    h = 1.0 / (n-1);
    for (i=0; i<n; i++) 
	for (j=0; j<n; j++) {
	    x = i*h;
	    y = j*h;
	    b[i+j*n] = h*h*f(x,y);
	};
    
    // 
    niter = 10000;
    for (iter=1; iter<=niter; iter++) {
#pragma omp for private(i)
        for (j=1; j<n-1; j++)
	    for (i=1; i<n-1; i++) {
	      if ((i+j) % 2 == 0)
		u[i+j*n] = (u[i-1 + j*n] + u[i + (j-1)*n] + 
			    u[i+1 + j*n] + u[i + (j+1)*n] - b[i+j*n])/4.0;
	    }

#pragma omp for private(i)
        for (j=1; j<n-1; j++)
	    for (i=1; i<n-1; i++) {
	      if ((i+j) % 2 == 1)
		u[i+j*n] = (u[i-1 + j*n] + u[i + (j-1)*n] + 
			    u[i+1 + j*n] + u[i + (j+1)*n] - b[i+j*n])/4.0;
	    }

	if (iter % 1000==0)
	    printf("%i u(0.5,0.5)= %lf\n", iter, u[(n-1)/2 + ((n-1)/2)*n]);
    }


}









