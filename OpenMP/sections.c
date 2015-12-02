/**********************************************************************
 * Tento program demonstruje pouziti konstrukce 
 * #omp parallel sections
 *
 * V kazde sekci se vykona urcity vypocet (zde soucet 0+1+2+...+random()-1)
 * a pak se na obrazovku vytiskne zprava o ukonceni sekce.
 *
 * Pri seriovem provedeni ulohy vzdy probehne nejprve vypocet prvni sekce
 * a po jeho dokonceni vypoce sekce 2. Pri paralelnim provadeni se vytiskne
 * drive zprava z te sekce, ktera byla drive dokoncena.
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    int i, j, ni, nj;
    double si, sj;

    srand((long)time(NULL));

#pragma omp parallel sections
    {
#pragma omp section
	ni = (int)(1000000.0*rand()/(RAND_MAX+1.0));
	printf("Prvni sekce, ni=%i\n",ni);
	si = 0;
	for (i=0; i<ni; i++) si += i;
	printf("Konec prvni sekce sum(0..%i)=%lf\n", ni, si);
#pragma omp section
	nj = (int)(1000000.0*rand()/(RAND_MAX+1.0));
	printf("Druha sekce, nj=%i\n",nj);
	sj = 0;
	for (j=0; j<nj; j++) sj += j;
	printf("Konec druhe sekce sum(0..%i)=%lf\n", nj, sj);
    }
    
    return 0;
}
