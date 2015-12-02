program MM

  integer, parameter :: dbl = kind(0.d0)
  integer, parameter :: N = 1000

  real(dbl) :: A(N,N), B(N,N), C(N,N)
  real(dbl) :: time1, time2, tr, s
  integer   :: dt1(8), dt2(8)

  do j = 1, N
     do i = 1, N
        A(i,j) = i+j
        B(i,j) = i-j
        C(i,j) = 0
     end do
  end do
  
  call cpu_time(time1)
  call date_and_time(values=dt1)
  
  do j = 1, N
     do i = 1, N
        s = 0
        !$omp parallel do shared(A,B,i,j) reduction(+:s)
        do k = 1, N
           s = s + A(i,k)*B(k,j)
        end do
        C(i,j) = s
     end do
  end do
  
  call cpu_time(time2)
  call date_and_time(values=dt2)
  
  print *, "Total CPU time :", time2 - time1
  print *, "Wall time      :", (dt2(8)-dt1(8))/1000.0 + (dt2(7)-dt1(7)) + &
       (dt2(6)-dt1(6))*60 + (dt2(5)-dt1(5))*60*24

  tr = 0
  do i = 1, N
     tr = tr + C(i,i)
  end do
  print *, "Trace(C) = ", tr

end program MM
