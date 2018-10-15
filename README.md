# NUSO

This repository contains support material for my lecture "Numerical
Software" at the Czech Technical University in Prague. 


## MATMUL benchmark
* [poznámky](notebooks/matmul_benchmark.ipynb)
* matmul, matice jako ukazatel na ukazatele [dgemm_double.cpp](Benchmark/dgemm_doublepp.cpp)
* matmul, matice jako jednorozměrné pole [dgemm.cpp](Benchmark/dgemm.cpp)
* matmul, matice jako jednoduchá třída Matrix [dgemm_class.cpp](Benchmark/dgemm_class.cpp)
* matmul, násobení pomocí knihovnyBLAS [dgemm_blas.cpp](Benchmark/dgemm_blas.cpp)


## LU benchmark
* [poznámky](notebooks/lu_benchmark.ipynb)

## BLAS & LAPACK
* BLAS lvl 1 via Fortran interface [fblas1.cpp](BlasLapack/fblas1.cpp)
* BLAS lvl 1 via Cblas interface [cblas1.cpp](BlasLapack/cblas1.cpp)
* BLAS lvl 2, plná matice [cblas2.cpp](BlasLapack/cblas2.cpp)
* BLAS lvl 2, symetrická matice [cblas2_dsymv.cpp](BlasLapack/cblas2_dsymv.cpp)
* BLAS lvl 2, symetrická komprimovaná matice [cblas2_dspmv.cpp](BlasLapack/cblas2_dspmv.cpp)
* BLAS lvl 2, pásová matice [cblas2_dgbmv.cpp](BlasLapack/cblas2_dgbmv.cpp)
* BLAS lvl 3, plná matice [cblas3_dgemm.cpp](BlasLapack/cblas3_dgemm.cpp)

* LAPACKe, *A x = b*, simple driver [lapack_gesv.cpp](BlasLapack/lapack_gesv.cpp)
* LAPACKe, *A x = b* expert driver [lapack_gesvx.cpp](BlasLapack/lapack_gesvx.cpp)

* ScaLAPACK, *A x = b*, [scalapack_example.f](BlasLapack/scalapack_example.f),
  test matrix [m1000](BlasLapack/m1000)
  
## Error analysis
* [poznámky](notebooks/error_analysis.ipynb)

## MPI
* Základní příklad [basic.cpp](MPI/basic.cpp)
* Základní příklad se send/recieve [basic_send.cpp](MPI/basic_send.cpp)
* Laplaceova rovnice, statická data [laplace-mpi1.cpp](MPI/laplace-mpi1.cpp)
* Laplaceova rovnice, dynamická alokace [laplace-mpi_dynamic1.cpp](MPI/laplace-mpi_dynamic1.cpp)
* Laplaceova rovnice, dynamická alokace, vylepšeno [laplace-mpi_dynamic2.cpp](MPI/laplace-mpi_dynamic2.cpp)
* Helmholtzova rovnice, chybná implementace (deadlock)  [helmholtz-deadlock.cpp](MPI/helmholtz-deadlock.cpp)
* Helmholtzova rovnice, lepší implementace [helmholtz-serialized.cpp](MPI/helmholtz-serialized.cpp)
* Helmholtzova rovnice, správná implementace [helmholtz-fast.cpp](MPI/helmholtz-fast.cpp)

## UMFPACK
* Základní příklad [umfpack_simple.c](Umfpack/umfpack_simple.c)
* Řešení soustavy rovnic, C++,  [example.cpp](Umfpack/example.cpp)
* Řešení soustavy rovnic, velká matice,  [poisson3Da.cpp](Umfpack/poisson3Da.cpp)
* Řešení soustavy rovnic, MKL verze,  [poisson3Da_dss.cpp](Umfpack/poisson3Da_dss.cpp)
* [poznámky](notebooks/sparse_matrix_ordering.ipynb)


## [PETSc](http://www.mcs.anl.gov/petsc)

## METIS
* [poznámky](notebooks/spectral_bisection.ipynb)

## [Julia](www.julialang.org)
* Intervalová aritmetika, [skript](notebooks/intervaly.ipbn)
* Implementace neuronové sítě pomocí knihovny Knet, [skript](notebooks/neuronova_sit_1.ipbn)

## Zápočtové úlohy

* Matice *L* je pěti-diagonální matice o rozměrech $n \times n$, 
  diagonální prvky jsou *L(i,i)=4*, 
  na vedlejších diagonálách je *L(i,i±1)=-1* 
  a na vedlejších diagonálách vzdálených o
  *m=int(sqrt(n))* je *L(i,i±m}=-1*.
  
### Úloha 1 
	Určete vlastní čísla matice *L* pro *n=1000*. 
	Pro výpočet použijte knihovnu LAPACK.

### Úloha 2 
	Určete řešení soustavy *Lx=b* pro *n=1000*, *bᵢ=2*. 
	Pro výpočet použijte knihovnu LAPACK.

### Úloha 3 
	Určete řešení soustavy *Lx=b* pro *n=10000*, *bᵢ=2*. 
	Pro výpočet použijte knihovnu UMFPack.

### Úloha 4 
	Určete řešení soustavy *Lx=b* pro $n=100000$, *bᵢ=2*. 
	Pro výpočet použijte knihovnu PETSc.

### Úloha 5 
	Určete řešení soustavy *Lx=b* pro *n=1000000*, *bᵢ=2*. 
	Pro výpočet použijte knihovnu PETSc v paralelní verzi.
