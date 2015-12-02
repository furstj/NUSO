program Laplace

  integer, parameter :: dbl = selected_real_kind(15)

  integer, parameter :: N = 2000

  real(dbl) :: U(N,N)
  real(dbl) :: V(N,N)
  real(dbl) :: F(N,N)

  integer :: ITER, TICKSTART, TICKSTOP, TICKRATE

  F = (1.0/N)**2
  U = 0

  call system_clock (TICKSTART, TICKRATE)

  do iter = 1, 100

     !$omp parallel do private(i) collapse(2)
     do j = 2, N-1
        do i = 2, N-1
           V(i,j) = (U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)) / 4 + F(i,j)
        end do
     end do
     
     !$omp parallel do private(i) collapse(2)
     do j = 2, N-1
        do i = 2, N-1
           U(i,j) = V(i,j)
        end do
     end do


  end do
  
  call system_clock (TICKSTOP)
  
  TICKSTOP = TICKSTOP - TICKSTART
  TIME     = float(TICKSTOP)/float(TICKRATE)
  
  print *, "U(0.5,0.5) = ", U(N/2,N/2)
  print *, "TIME = ", TIME

end program Laplace
