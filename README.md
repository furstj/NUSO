# NUSO

This repository contains support material for my lecture "Numerical
Software" at the Czech Technical University in Prague. 


## MATMUL benchmark
* [lecture notes](notebooks/matmul_benchmark.ipynb)
* matmul with double pointer [dgemm_double.cpp](Benchmark/dgemm_doublepp.cpp)
* matmul with plain array [dgemm.cpp](Benchmark/dgemm.cpp)
* matmul with simple Matrix class [dgemm_class.cpp](Benchmark/dgemm_class.cpp)
* matmul with BLAS [dgemm_blas.cpp](Benchmark/dgemm_blas.cpp)


## LU benchmark
* [lecture notes](notebooks/lu_benchmark.ipynb)

## BLAS & LAPACK
* BLAS lvl 1 via Fortran interface [fblas1.cpp](BlasLapack/fblas1.cpp)
* BLAS lvl 1 via Cblas interface [cblas1.cpp](BlasLapack/cblas1.cpp)
* BLAS lvl 2, full matrix [cblas2.cpp](BlasLapack/cblas2.cpp)
* BLAS lvl 2, symmetric matrix [cblas2_dsymv.cpp](BlasLapack/cblas2_dsymv.cpp)
* BLAS lvl 2, symmetric packed matrix [cblas2_dspmv.cpp](BlasLapack/cblas2_dspmv.cpp)
* BLAS lvl 2, banded matrix [cblas2_dgbmv.cpp](BlasLapack/cblas2_dgbmv.cpp)
* BLAS lvl 3, full matrix [cblas3_dgemm.cpp](BlasLapack/cblas3_dgemm.cpp)

* LAPACKe, *A x = b*, simple driver [lapack_gesv.cpp](BlasLapack/lapack_gesv.cpp)
* LAPACKe, *A x = b* expert driver [lapack_gesvx.cpp](BlasLapack/lapack_gesvx.cpp)

* ScaLAPACK, *A x = b*, [scalapack_example.f](BlasLapack/scalapack_example.f),
  test matrix [m1000](BlasLapack/m1000)
  
## MPI
* Basic example [basic.cpp](MPI/basic.cpp)
* Basic example with send/recieve [basic_send.cpp](MPI/basic_send.cpp)
* Laplace equation, static data [laplace-mpi1.cpp](MPI/laplace-mpi1.cpp)
* Laplace equation, dynamic allocation [laplace-mpi_dynamic1.cpp](MPI/laplace-mpi_dynamic1.cpp)
* Laplace equation, dynamic allocation, improved [laplace-mpi_dynamic2.cpp](MPI/laplace-mpi_dynamic2.cpp)
* Helmholtz equation, wrong implementation (deadlock)  [helmholtz-deadlock.cpp](MPI/helmholtz-deadlock.cpp)
* Helmholtz equation, better implementation [helmholtz-serialized.cpp](MPI/helmholtz-serialized.cpp)
* Helmholtz equation, good implementation [helmholtz-fast.cpp](MPI/helmholtz-fast.cpp)

## UMFPACK
* Basic example [umfpack_simple.c](Umfpack/umfpack_simple.c)
* Advanced example, C++,  [example.cpp](Umfpack/example.cpp)
* Advanced example, large matrix,  [poisson3Da.cpp](Umfpack/poisson3Da.cpp)
* Advanced example, MKL version,  [poisson3Da_dss.cpp](Umfpack/poisson3Da_dss.cpp)

## Error analysis
* [lecture notes](notebooks/error_analysis.ipynb)

## [PETSc](http://www.mcs.anl.gov/petsc)
	
