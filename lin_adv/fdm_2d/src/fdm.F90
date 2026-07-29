module fdm
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none
contains

subroutine finite_diffence_method(ua, res)
  implicit none
  PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
  ! Local variables
  
  res = 0.d0
  
  select case(trim(space_disc))
  case('CD2')
    call CD2(ua, res)
  case('CD4')
    call CD4(ua, res)
  case('CD6')
    call CD6(ua, res)  
  case('LELE')
    call LELE(ua, res)
  case default
    print*, 'please select a finite difference method in space'
    print*, 'Available options: CD2, CD4, CD6, LELE'
    stop
  end select
end subroutine finite_diffence_method

subroutine CD2(ua, res)
  implicit none
  PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
  ! Local variables
  integer        :: i, j
  real(dp) :: idchi, ideta
  
  idchi = 0.5d0/dchi; ideta = 0.5d0/deta
  do j = jst, jen
  do i = ist, ien
     res(i,j) = - speed_chi * (ua(i+1,j)-ua(i-1,j)) * idchi &
                - speed_eta * (ua(i,j+1)-ua(i,j-1)) * ideta
  enddo
  enddo
end subroutine CD2

subroutine CD4(ua, res)
  implicit none
  PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
  ! Local variables
  integer  :: i, j
  real(dp) :: idchi, ideta, a_m2, a_m1, a_p1, a_p2
  
  a_m2 =  1.d0/12.d0
  a_m1 = -2.d0/3.d0
  a_p1 = - a_m1
  a_p2 = - a_m2
  idchi = 1.d0/dchi; ideta = 1.d0/deta
  do j = jst, jen
  do i = ist, ien
     res(i,j) = - speed_chi * (a_m2*ua(i-2,j)+a_m1*ua(i-1,j)+a_p1*ua(i+1,j)+a_p2*ua(i+2,j)) * idchi &
                - speed_eta * (a_m2*ua(i,j-2)+a_m1*ua(i,j-1)+a_p1*ua(i,j+1)+a_p2*ua(i,j+2)) * ideta
  enddo
  enddo
end subroutine CD4

subroutine CD6(ua, res)
  implicit none
  PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
  ! Local variables
  integer  :: i, j
  real(dp) :: idchi, ideta, a_m3, a_m2, a_m1, a_p1, a_p2, a_p3
  
  a_m3 = - 1.d0/60.d0
  a_m2 =   3.d0/20.d0
  a_m1 = - 3.d0/4.d0
  a_p1 = - a_m1
  a_p2 = - a_m2
  a_p3 = - a_m3
  idchi = 1.d0/dchi; ideta = 1.d0/deta
  do j = jst, jen
  do i = ist, ien
     res(i,j) = - speed_chi * (a_m3*ua(i-3,j)+a_m2*ua(i-2,j)+a_m1*ua(i-1,j)+a_p1*ua(i+1,j)+a_p2*ua(i+2,j)+a_p3*ua(i+3,j)) * idchi &
                - speed_eta * (a_m3*ua(i,j-3)+a_m2*ua(i,j-2)+a_m1*ua(i,j-1)+a_p1*ua(i,j+1)+a_p2*ua(i,j+2)+a_p3*ua(i,j+3)) * ideta
  enddo
  enddo
end subroutine CD6

subroutine LELE(ua, res)
  implicit none
  PetscScalar :: ua(gist:gien,gjst:gjen), res(ist:ien,jst:jen)
  ! Local variables
  integer  :: i, j
  real(dp), dimension(gist:gien) :: uf_x, ud_x
  real(dp), dimension(gjst:gjen) :: uf_y, ud_y
  real(dp) :: a_m1, a_0, a_p1, idchi, ideta
  real(dp) :: b_m2, b_m1, b_0, b_p1, b_p2
  
  !interior nodes
  a_m1 = 1.d0/3.d0; a_0  = 1.d0; a_p1 = a_m1 
  b_m2 = -1.d0/36.d0; b_m1 = -14.d0/18.d0; b_0  = 0.d0; b_p1 = - b_m1; b_p2 = - b_m2
  idchi = 1.d0/dchi; ideta = 1.d0/deta
  
  do j = jst, jen
     uf_x(:)   = ua(:, j)
     call apply_compact_scheme(gist, gien, uf_x, ud_x, idchi, a_m1, a_0, a_p1, b_m2, b_m1, b_0, b_p1, b_p2)
     res(ist:ien, j) = res(ist:ien, j) + speed_chi * ud_x(ist:ien)
  enddo
  do i = ist, ien
     uf_y(:)   = ua(i, :)
     call apply_compact_scheme(gjst, gjen, uf_y, ud_y, ideta, a_m1, a_0, a_p1, b_m2, b_m1, b_0, b_p1, b_p2)
     res(i, jst:jen) = res(i, jst:jen) + speed_eta * ud_y(jst:jen)
  enddo
end subroutine LELE

! applies compact finite difference scheme in a chosen direction
subroutine apply_compact_scheme(st, en, uf, ud, idh, a_m1, a_0, a_p1, b_m2, b_m1, b_0, b_p1, b_p2)
  implicit none
  integer, intent(in) :: st, en
  PetscScalar, intent(in)  :: uf(st:en)
  PetscScalar, intent(out) :: ud(st:en)
  real(dp), intent(in)  :: idh, a_m1, a_0, a_p1, b_m2, b_m1, b_0, b_p1, b_p2
  ! Local variables
  integer        :: i
  real(dp), dimension(st:en) :: udia, dia, ldia, rhs
  
  udia = 0.d0 ; dia = 0.d0 ; ldia = 0.d0; rhs = 0.d0
  
  ldia(st) = 0.d0 ; dia(st) = 1.d0 ; udia(st) = 0.d0 
  rhs(st) = - (uf(st+1) - uf(st)) * idh

  ldia(st+1) = 0.d0 ; dia(st+1)  = 1.d0 ; udia(st+1) = 0.d0
  rhs(st+1) = - 0.5d0 * (uf(st+2) - uf(st)) * idh

  do i = st+2, en-2
     ldia(i) = a_m1; dia(i) = a_0; udia(i) = a_p1
     rhs(i) = - ( b_m2 * uf(i-2) + b_m1 * uf(i-1) + b_0 * uf(i) + &
                  b_p1 * uf(i+1) + b_p2 * uf(i+2) ) * idh
  enddo

  ldia(en-1) = 0.d0 ; dia(en-1) = 1.d0 ; udia(en-1) = 0.d0
  rhs(en-1) = - 0.5d0 * (uf(en) - uf(en-2)) * idh

  ldia(en) = 0.d0 ; dia(en) = 1.d0 ; udia(en) = 0.d0
  rhs(en)  = - (uf(en) - uf(en-1)) * idh
  
  call tdma(ldia(st+1:en), dia(st:en), udia(st:en-1), rhs(st:en), ud(st:en), en-st+1)
end subroutine apply_compact_scheme

! Thomas algorithm to solve for tridiagonal matrix system
subroutine tdma(a3, b3, c3, d3, xa, N)
  implicit none
  integer, intent(in):: N
  real(dp), intent(in), dimension(2:N)  :: a3
  real(dp), intent(in), dimension(1:N-1):: c3
  real(dp), intent(in), dimension(1:N)  :: b3, d3
  PetscScalar, intent(out), dimension(1:N) :: xa
  real(dp), allocatable, dimension(:) :: beta1, gamma1
  integer :: i
  
  allocate(beta1(1:N),gamma1(1:N))
  beta1(1)=b3(1)
  gamma1(1)=d3(1)/b3(1)
  do i = 2, N
    beta1(i)  =   b3(i) - a3(i) * c3(i-1)/beta1(i-1)
    gamma1(i) = ( d3(i) - a3(i)*gamma1(i-1) ) / beta1(i)
  enddo
  xa(N) = gamma1(N)
  do i = N-1, 1, -1
    xa(i) = gamma1(i) - c3(i) * xa(i+1)/beta1(i)
  enddo
  deallocate(beta1, gamma1)
end subroutine TDMA

end module fdm
