program Laplace

  use openacc

  implicit none
  
  integer, parameter :: n = 4096
  integer, parameter :: flt = selected_real_kind(5)
  integer, parameter :: i32 = selected_int_kind(5)

  real(flt) :: u(n,n), v(n,n), f(n,n)
  real(flt) :: h = 1.0 / (n-1)

  integer(i32) :: i, j
  integer      :: iclock_start, iclock_stop, iclock_rate, iter

  u(:,:)    = 0.0
  v(:,:) = u(:,:)
  f(:,:)    = 1.0

  call system_clock(COUNT=iclock_start, COUNT_RATE=iclock_rate)

  !$acc data copy(u, f) create(v)
  do iter = 1, 100

     !$acc kernels loop 
     do j = 2, n-1
        do i = 2, n-1
           v(i,j) = (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) + h**2*f(i,j)) / 4.0
        end do
     end do

     !$acc kernels loop 
     do j = 2, n-1
        do i = 2, n-1
           u(i,j) = (v(i-1,j) + v(i+1,j) + v(i,j-1) + v(i,j+1) + h**2*f(i,j)) / 4.0
        end do
     end do

  end do
  !$acc end data

  call system_clock(COUNT=iclock_stop)

  print *, "u_center = ", u(n/2,n/2)
  print *, "Time elapsed: ", real(iclock_stop - iclock_start) / iclock_rate

end program Laplace
