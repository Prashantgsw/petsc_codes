! Verification of compact finite-difference schemes
! using periodic boundaries.

! Test function:
!   f(x) = sin(kx)

! Exact derivative:
!   df/dx = k cos(kx)

! Computes L2 error and observed order
! for LELE_CD6, LELE_CD8 and LELE_CD10 schemes.
!build explicitly-- make verify_convergence


program verify_convergence !build explicitly-- make verify_convergence
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  integer :: Np_list(6) = (/20, 40, 80, 160, 320, 640/)
  integer :: idx_np, Np
  real(dp) :: err6, err8, err10
  integer :: iounit

  open(newunit=iounit, file='convergence_results.csv', status='replace', action='write')
  write(iounit,'(A)') 'Np,LELE_CD6,LELE_CD8,LELE_CD10'

  do idx_np = 1, size(Np_list)
    Np = Np_list(idx_np)
    call compute_error(Np, 'LELE_CD6',  err6)
    call compute_error(Np, 'LELE_CD8',  err8)
    call compute_error(Np, 'LELE_CD10', err10)
    write(iounit,'(I6,A,E16.8,A,E16.8,A,E16.8)') Np, ',', err6, ',', err8, ',', err10
    print '(I6,3X,E16.8,3X,E16.8,3X,E16.8)', Np, err6, err8, err10
  end do

  close(iounit)
  print *, 'Wrote convergence_results.csv'

contains

  ! Wraps any integer index i (even negative, even far out of range) into 1..N
  integer function wrap_idx(i, N) result(w)
    integer, intent(in) :: i, N
    w = mod(i - 1 + N, N) + 1
  end function wrap_idx

  subroutine compute_error(Np, scheme, l2err)
    integer, intent(in) :: Np
    character(len=*), intent(in) :: scheme
    real(dp), intent(out) :: l2err
    real(dp), allocatable :: A(:,:), rhs(:), f(:), df_exact(:)
    real(dp) :: dx, idx, x, kwave
    real(dp) :: a_m1, a_m2, a_0, a_p1, a_p2
    real(dp) :: b_m3, b_m2, b_m1
    integer :: i, im1, im2, im3, ip1, ip2, ip3
    integer, allocatable :: ipiv(:)
    integer :: info

    allocate(A(Np,Np), rhs(Np), f(Np), df_exact(Np), ipiv(Np))
    A = 0.0_dp

    dx  = 1.0_dp / real(Np, dp)
    idx = 1.0_dp/dx
    kwave = 2.0_dp*pi*3.0_dp   ! 3 full periods over the domain [0,1)

    do i = 1, Np
      x = real(i-1, dp)*dx
      f(i) = sin(kwave*x)
      df_exact(i) = kwave*cos(kwave*x)
    end do

    select case (trim(scheme))
    case ('LELE_CD6')
      a_m1 = 1.0_dp/3.0_dp;  a_0 = 1.0_dp;  a_p1 = a_m1
      b_m2 = -1.0_dp/36.0_dp;  b_m1 = -14.0_dp/18.0_dp
      do i = 1, Np
        im1 = wrap_idx(i-1, Np); ip1 = wrap_idx(i+1, Np)
        im2 = wrap_idx(i-2, Np); ip2 = wrap_idx(i+2, Np)
        A(i,im1) = a_m1; A(i,i) = a_0; A(i,ip1) = a_p1
        rhs(i) = (b_m2*f(im2) + b_m1*f(im1) - b_m1*f(ip1) - b_m2*f(ip2))*idx
      end do

    case ('LELE_CD8')
      a_m1 = 4.0_dp/9.0_dp;  a_m2 = 1.0_dp/36.0_dp;  a_0 = 1.0_dp
      a_p1 = a_m1;  a_p2 = a_m2
      b_m2 = -25.0_dp/216.0_dp;  b_m1 = -20.0_dp/27.0_dp
      do i = 1, Np
        im1 = wrap_idx(i-1, Np); ip1 = wrap_idx(i+1, Np)
        im2 = wrap_idx(i-2, Np); ip2 = wrap_idx(i+2, Np)
        A(i,im2) = a_m2; A(i,im1) = a_m1; A(i,i) = a_0; A(i,ip1) = a_p1; A(i,ip2) = a_p2
        rhs(i) = (b_m2*f(im2) + b_m1*f(im1) - b_m1*f(ip1) - b_m2*f(ip2))*idx
      end do

    case ('LELE_CD10')
      a_m1 = 0.5_dp;  a_m2 = 1.0_dp/20.0_dp;  a_0 = 1.0_dp
      a_p1 = a_m1;  a_p2 = a_m2
      b_m3 = -1.0_dp/600.0_dp;  b_m2 = -101.0_dp/600.0_dp;  b_m1 = -17.0_dp/24.0_dp 
      do i = 1, Np
        im1 = wrap_idx(i-1, Np); ip1 = wrap_idx(i+1, Np)
        im2 = wrap_idx(i-2, Np); ip2 = wrap_idx(i+2, Np)
        im3 = wrap_idx(i-3, Np); ip3 = wrap_idx(i+3, Np)
        A(i,im2) = a_m2; A(i,im1) = a_m1; A(i,i) = a_0; A(i,ip1) = a_p1; A(i,ip2) = a_p2
        rhs(i) = (b_m3*f(im3) + b_m2*f(im2) + b_m1*f(im1) &
                 - b_m1*f(ip1) - b_m2*f(ip2) - b_m3*f(ip3))*idx
      end do
    end select

    call dgesv(Np, 1, A, Np, ipiv, rhs, Np, info)
    if (info /= 0) print *, 'DGESV failed, info=', info

    l2err = sqrt(sum((rhs - df_exact)**2) / real(Np, dp))

    deallocate(A, rhs, f, df_exact, ipiv)
  end subroutine compute_error

end program verify_convergence  
!gfortran -o verify_convergence verify_convergence.f90 -llapack -lblas