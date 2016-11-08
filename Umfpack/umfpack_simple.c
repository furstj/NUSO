#include <stdio.h>
#include "suitesparse/umfpack.h"

/*
 Reseni soustavy rovnic Ax=b, kde
 
     | 2  3  .  .  . |        |  8 |
     | 3  .  4  .  6 |        | 45 |
 A = | . -1 -3  2  . |,   b = | -3 |
     | .  .  1  .  . |        |  3 |
     | .  4  2  .  1 |        | 19 |
 
 Spravne reseni je x = [1, 2, 3, 4, 5]'

*/    

int n = 5;
int Ap [ ] =    { 0,     2,           5,                9, 10,     12};
int Ai [ ]    = { 0, 1,  0,   2,  4,  1,   2,  3,   4,  2,  1,  4};
double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
double b [ ] = {8., 45., -3., 3., 19.} ;
double x [5] ;

int main (void)
{
  double *null = (double *) NULL ;
  int i ;
  void *Symbolic, *Numeric ;

  umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
  umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
  umfpack_di_free_symbolic (&Symbolic) ;
  umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;
  umfpack_di_free_numeric (&Numeric) ;
  for (i = 0 ; i < n ; i++) printf ("x[%d] = %g\n", i, x [i]) ;
  return (0) ;
}

