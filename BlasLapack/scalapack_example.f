      program example
      
      include "mpif.h"

      parameter (N=1000)
      parameter (NL=1000)
      parameter (NB=64)

      double precision A_global(N,N), B_global(N)
      integer DESCA_g(9), DESCB_g(9)

      double precision A(NL,NL), B(NL)
      integer DESCA(9), DESCB(9)

      integer IPIV(N), ictx
      integer nprocs, myproc, i, j
 
      character*256 line
      double precision t1, t2, second
      external second

C     Inicializace procesorove site, procesory jsou usporadany do
C     jednoho sloupce
      call blacs_pinfo(iam, nprocs)

      call blacs_get(0, 0, ictx)
      call blacs_gridinit(ictx, "Row", nprocs, 1)
      call blacs_gridinfo(ictx, nprow, npcol, myprow, mypcol)
      if (iam.eq.0) then
         print *, "Sit:", nprow, " x", npcol
      endif

C     Nacteni matice a prave strany ze souboru
      if (iam.eq.0) then
         print *, "Nacitam matici ze souboru..."
         open (20, file="m1000", status="old")
         read (20,*) line
         read (20,*) line
         read (20,*) line
         read (20,*) line
         read (20,*) line
         read (20,*) ((A_global(i,j), j=1, N), i=1, N)
         read (20,*) line
         read (20,*) line
         read (20,*) line
         read (20,*) line
         read (20,*) (B_global(j), j=1, N)
         close (20)
      end if

C     Definice distribuce matice A a B
      call descinit(DESCA,   N, N, NB, NB, 0, 0, ictx, NL, info)
      if (iam.eq.0) print *, "DESCA  = ", DESCA
      call descinit(DESCB,   N, 1, NB,  1, 0, 0, ictx, N, info)
      if (iam.eq.0) print *, "DESCB  = ", DESCB

C     Definice distribuce matice A_global a B_global
      call descinit(DESCA_g, N, N,  N,  N, 0, 0, ictx, NL, info)
      if (iam.eq.0) print *, "DESCAg = ", DESCA_g
      call descinit(DESCB_g, N, 1,  N,  1, 0, 0, ictx, N, info)
      if (iam.eq.0) print *, "DESCBg = ", DESCB_g

      if (iam.eq.0) print *, "Distribuuji data na procesory..."
      call pdlacpy("All", N,1,B_global,1,1,DESCB_g, B,1,1,DESCB)
      call pdlacpy("All", N,N,A_global,1,1,DESCA_g, A,1,1,DESCA)

      if (iam.eq.0) print *, "Resim rovnici ..."
      t1 = MPI_WTIME()
      call pdgesv(N, 1, A, 1, 1, DESCA, ipiv, B, 1, 1, DESCB, info)
      t2 = MPI_WTIME()

      if (iam.eq.0) print *, "time = ", t2-t1
      if (iam.eq.0) print *, "info=0"

      call blacs_gridexit(ictx)
      call blacs_exit(0)

      end










