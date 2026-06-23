module fdm
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none

contains

subroutine finite_diffence_method(ua, res, ctx)
implicit none
type(tsdata)         :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i

res = 0.d0

select case(trim(space_disc))
case('CD2')
  call CD2(ua, res, ctx)
case('CD4')  
  call CD4(ua, res, ctx)
case('CD6')
  call CD6(ua, res, ctx)
case('CD8')  
  call CD8(ua, res, ctx)
case('LELE_CD4')
  call LELE_CD4(ua, res, ctx)
case('LELE_CD6')
  call LELE_CD6(ua, res, ctx)
case('LELE_CD8')
  call LELE_CD8(ua, res, ctx)
case('LELE_CD10')
  call LELE_CD10(ua, res, ctx)
case default
  print*, 'please select a finite difference method in space'
  print*, 'Available options: CD2, CD4, CD6, CD8, LELE_CD4, LELE_CD6, LELE_CD8, LELE_CD10'
  stop
end select

end subroutine finite_diffence_method


subroutine CD2(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: idx

idx = 0.5d0/ctx%g%dx
do i = ist, ien
   res(i) = - speed * (ua(i+1)-ua(i-1))
enddo
res = res * idx
end subroutine CD2

subroutine CD4(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m2, a_m1, a_p1, a_p2

a_m2 =  1.d0/12.d0
a_m1 = -2.d0/3.d0
a_p1 = - a_m1
a_p2 = - a_m2
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = - speed * ( a_m2*ua(i-2) + a_m1*ua(i-1) + &
                        a_p1*ua(i+1) + a_p2*ua(i+2) )
enddo
res = res * idx
end subroutine CD4

subroutine CD6(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m3, a_m2, a_m1, a_p1, a_p2, a_p3

a_m3 = - 1.d0/60.d0
a_m2 =   3.d0/20.d0
a_m1 = - 3.d0/4.d0
a_p1 = - a_m1
a_p2 = - a_m2
a_p3 = - a_m3
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = - speed * ( a_m3*ua(i-3) + a_m2*ua(i-2) + a_m1*ua(i-1) + &
                        a_p1*ua(i+1) + a_p2*ua(i+2) + a_p3*ua(i+3) )
enddo
res = res * idx
end subroutine CD6

subroutine CD8(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer        :: i
real(dp) :: idx, a_m4, a_m3, a_m2, a_m1, a_p1, a_p2, a_p3, a_p4

a_m4 =   1.d0/280.d0
a_m3 = - 4.d0/105.d0
a_m2 =   1.d0/5.d0
a_m1 = - 4.d0/5.d0
a_p1 = - a_m1
a_p2 = - a_m2
a_p3 = - a_m3
a_p4 = - a_m4
idx  = 1.d0/ctx%g%dx
do i = ist, ien
   res(i) = - speed * ( a_m4*ua(i-4) + a_m3*ua(i-3) + a_m2*ua(i-2) + a_m1*ua(i-1) + &
                        a_p1*ua(i+1) + a_p2*ua(i+2) + a_p3*ua(i+3) + a_p4*ua(i+4) )
enddo
res = res * idx
end subroutine CD8

subroutine LELE_CD4(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: udia(gist:gien), dia(gist:gien)
real(dp)     :: ldia(gist:gien), rhs(gist:gien)
real(dp)     :: res_tmp(gist:gien)
real(dp) :: a_m1, a_m2, a_0, a_p1, idx
real(dp) :: b_m1, b_0, b_p1, b_m2, b_p2

udia = 0.0 ; dia = 0.0 ; ldia = 0.0; rhs = 0.0; res_tmp = 0.0

!interior nodes
a_m1 = 1.0/4.0
a_0  = 1.0
a_p1 = a_m1
b_m1 = -3.0/4.0
b_0  = 0.0
b_p1 = - b_m1
b_m2 = 0.0
b_p2 = 0.0

idx = 1.0/ctx%g%dx
ldia(gist) = 0.0 ; dia(gist) = 1.0 ; udia(gist) = 0.0 
rhs(gist) = - speed * (ua(gist+1) - ua(gist)) * idx
ldia(gist+1) = 0.0 ; dia(gist+1)  = 1.0 ; udia(gist+1) = 0.0
rhs(gist+1) = - speed * (0.5*ua(gist+2) - 0.5* ua(gist)) * idx
do i = gist+2, gien-2
   ldia(i) = a_m1; dia(i) = a_0; udia(i) = a_p1
   rhs(i) = - speed * ( b_m2 * ua(i-2) + b_m1 * ua(i-1) + b_0 * ua(i) + &
                        b_p1 * ua(i+1) + b_p2 * ua(i+2) ) * idx
enddo
ldia(gien-1) = 0.0 ; dia(gien-1) = 1.0 ; udia(gien-1) = 0.0
rhs(gien-1) = - speed * (0.5*ua(gien) - 0.5*ua(gien-2)) * idx
ldia(gien) = 0.0 ; dia(gien)  = 1.0 ; udia(gien) = 0.0
rhs(gien)  = - speed * (ua(gien) - ua(gien-1)) * idx

call tdma(ldia(gist+1:gien), dia(gist:gien), udia(gist:gien-1), & 
           rhs(gist:gien), res_tmp(gist:gien), gien-gist+1)

res(ist:ien) = res_tmp(ist:ien)

end subroutine LELE_CD4


subroutine LELE_CD6(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: udia(gist:gien), dia(gist:gien)
real(dp)     :: ldia(gist:gien), rhs(gist:gien)
real(dp)     :: res_tmp(gist:gien)
real(dp) :: a_m1, a_0, a_p1, idx
real(dp) :: b_m2, b_m1, b_0, b_p1, b_p2

udia = 0.0 ; dia = 0.0 ; ldia = 0.0; rhs = 0.0; res_tmp = 0.0

!interior nodes
a_m1 = 1.0/3.0
a_0  = 1.0
a_p1 = a_m1
b_m2 = -1.0/36.0
b_m1 = -14.0/18.0
b_0  = 0.0
b_p1 = - b_m1
b_p2 = - b_m2

idx = 1.0/ctx%g%dx
ldia(gist) = 0.0 ; dia(gist) = 1.0 ; udia(gist) = 0.0 
rhs(gist) = - speed * (ua(gist+1) - ua(gist)) * idx
ldia(gist+1) = 0.0 ; dia(gist+1)  = 1.0 ; udia(gist+1) = 0.0
rhs(gist+1) = - speed * (0.5*ua(gist+2) - 0.5* ua(gist)) * idx
do i = gist+2, gien-2
   ldia(i) = a_m1; dia(i) = a_0; udia(i) = a_p1
   rhs(i) = - speed * ( b_m2 * ua(i-2) + b_m1 * ua(i-1) + b_0 * ua(i) + &
                        b_p1 * ua(i+1) + b_p2 * ua(i+2) ) * idx
enddo
ldia(gien-1) = 0.0 ; dia(gien-1) = 1.0 ; udia(gien-1) = 0.0
rhs(gien-1) = - speed * (0.5*ua(gien) - 0.5*ua(gien-2)) * idx
ldia(gien) = 0.0 ; dia(gien)  = 1.0 ; udia(gien) = 0.0
rhs(gien)  = - speed * (ua(gien) - ua(gien-1)) * idx

call tdma(ldia(gist+1:gien), dia(gist:gien), udia(gist:gien-1), & 
           rhs(gist:gien), res_tmp(gist:gien), gien-gist+1)

res(ist:ien) = res_tmp(ist:ien)

end subroutine LELE_CD6

subroutine LELE_CD8(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: udia1(gist:gien), dia(gist:gien)
real(dp)     :: ldia1(gist:gien), rhs(gist:gien)
real(dp)     :: ldia2(gist:gien), udia2(gist:gien)
real(dp)     :: res_tmp(gist:gien)
real(dp) :: a_m1, a_m2, a_0, a_p1, a_p2, idx
real(dp) :: b_m2, b_m1, b_0, b_p1, b_p2

udia1 = 0.0 ; udia2 = 0.0; dia = 0.0 ; ldia1 = 0.0; ldia2 = 0.0; rhs = 0.0; res_tmp = 0.0

!interior nodes
a_m2 = 1.0/36.0
a_m1 = 4.0/9.0
a_0  = 1.0
a_p2 = a_m2
a_p1 = a_m1
b_m2 = -25.0/216.0
b_m1 = -20.0/27.0
b_0  = 0.0
b_p1 = - b_m1
b_p2 = - b_m2

idx = 1.0/ctx%g%dx
ldia2(gist) = 0.0 ; ldia1(gist) = 0.0 ; dia(gist) = 1.0 ; udia1(gist) = 0.0 ; udia2(gist) = 0.0
rhs(gist) = - speed * (ua(gist+1) - ua(gist)) * idx
ldia2(gist+1) = 0.0 ; ldia1(gist+1) = 0.0 ; dia(gist+1)  = 1.0 ; udia1(gist+1) = 0.0 ; udia2(gist+1) = 0.0
rhs(gist+1) = - speed * (0.5*ua(gist+2) - 0.5* ua(gist)) * idx
do i = gist+2, gien-2
   ldia2(i) = a_m2; ldia1(i) = a_m1; dia(i) = a_0; udia1(i) = a_p1; udia2(i) = a_p2
   rhs(i) = - speed * ( b_m2 * ua(i-2) + b_m1 * ua(i-1) + b_0 * ua(i) + &
                        b_p1 * ua(i+1) + b_p2 * ua(i+2) ) * idx
enddo
ldia2(gien-1) = 0.0 ; ldia1(gien-1) = 0.0 ; dia(gien-1) = 1.0 ; udia1(gien-1) = 0.0 ; udia2(gien-1) = 0.0
rhs(gien-1) = - speed * (0.5*ua(gien) - 0.5*ua(gien-2)) * idx
ldia2(gien) = 0.0 ; ldia1(gien) = 0.0 ; dia(gien)  = 1.0 ; udia1(gien) = 0.0 ; udia2(gien) = 0.0
rhs(gien)  = - speed * (ua(gien) - ua(gien-1)) * idx

call pdma(ldia2(gist:gien), ldia1(gist:gien), dia(gist:gien), &
            udia1(gist:gien), udia2(gist:gien), & 
           rhs(gist:gien), res_tmp(gist:gien), gien-gist+1)

res(ist:ien) = res_tmp(ist:ien)

end subroutine LELE_CD8


subroutine LELE_CD10(ua, res, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: ua(gist:gien), res(ist:ien)
! Local variables
integer      :: i
real(dp)     :: udia1(gist:gien), dia(gist:gien)
real(dp)     :: ldia1(gist:gien), rhs(gist:gien)
real(dp)     :: ldia2(gist:gien), udia2(gist:gien)
real(dp)     :: res_tmp(gist:gien)
real(dp) :: a_m1, a_m2, a_0, a_p1, a_p2, idx
real(dp) :: b_m3, b_m2, b_m1, b_0, b_p1, b_p2, b_p3

udia1 = 0.0 ; udia2 = 0.0; dia = 0.0 ; ldia1 = 0.0; ldia2 = 0.0; rhs = 0.0; res_tmp = 0.0

!interior nodes
a_m2 = 1.0/20.0
a_m1 = 1.0/2.0
a_0  = 1.0
a_p1 = a_m1
a_p2 = a_m2
b_m3 = -1.0/600.0
b_m2 = -101.0/600.0
b_m1 = -17.0/24.0
b_0  = 0.0
b_p1 = - b_m1
b_p2 = - b_m2
b_p3 = - b_m3

idx = 1.0/ctx%g%dx
ldia2(gist) = 0.0 ; ldia1(gist) = 0.0 ; dia(gist) = 1.0 ; udia1(gist) = 0.0 ; udia2(gist) = 0.0
rhs(gist) = - speed * (ua(gist+1) - ua(gist)) * idx
ldia2(gist+1) = 0.0 ; ldia1(gist+1) = 0.0 ; dia(gist+1)  = 1.0 ; udia1(gist+1) = 0.0 ; udia2(gist+1) = 0.0
rhs(gist+1) = - speed * (0.5*ua(gist+2) - 0.5* ua(gist)) * idx
ldia1(gist+2) = 0.0 ; ldia2(gist+2) = 0.0 ; dia(gist+2) = 1.0 ; udia1(gist+2) = 0.0 ; udia2(gist+2) = 0.0
rhs(gist+2) = - speed * (0.5*ua(gist+3) - 0.5* ua(gist+1)) * idx
do i = gist+3, gien-3 
   ldia2(i) = a_m2; ldia1(i) = a_m1; dia(i) = a_0; udia1(i) = a_p1; udia2(i) = a_p2
   rhs(i) = - speed * ( b_m3 * ua(i-3) + b_m2 * ua(i-2) + b_m1 * ua(i-1) + b_0 * ua(i) + &
                        b_p1 * ua(i+1) + b_p2 * ua(i+2) + b_p3 * ua(i+3)) * idx
enddo
ldia2(gien-2) = 0.0 ; ldia1(gien-2) = 0.0 ; dia(gien-2) = 1.0 ; udia2(gien-2) = 0.0 ; udia1(gien-2) = 0.0
rhs(gien-2) = - speed * (0.5*ua(gien-1) - 0.5*ua(gien-3)) * idx
ldia2(gien-1) = 0.0 ; ldia1(gien-1) = 0.0 ; dia(gien-1) = 1.0 ; udia1(gien-1) = 0.0 ; udia2(gien-1) = 0.0
rhs(gien-1) = - speed * (0.5*ua(gien) - 0.5*ua(gien-2)) * idx
ldia2(gien) = 0.0 ; ldia1(gien) = 0.0 ; dia(gien)  = 1.0 ; udia1(gien) = 0.0 ; udia2(gien) = 0.0
rhs(gien)  = - speed * (ua(gien) - ua(gien-1)) * idx


call pdma(ldia2(gist:gien), ldia1(gist:gien), dia(gist:gien), &
          udia1(gist:gien), udia2(gist:gien), &
          rhs(gist:gien), res_tmp(gist:gien), gien-gist+1)

res(ist:ien) = res_tmp(ist:ien)

end subroutine LELE_CD10


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

subroutine pdma(e3, a3, b3, c3, f3, d3, xa, N)
   use double
   implicit none

   integer,     intent(in)  :: N
   real(dp),    intent(in)  :: e3(1:N)   
   real(dp),    intent(in)  :: a3(1:N)  
   real(dp),    intent(in)  :: b3(1:N)  
   real(dp),    intent(in)  :: c3(1:N) 
   real(dp),    intent(in)  :: f3(1:N)   
   real(dp),    intent(in)  :: d3(1:N)   
   PetscScalar, intent(out) :: xa(1:N)   

   !-------------------------------------------
   ! Local working arrays
   !-------------------------------------------
   real(dp), allocatable :: beta(:), gamma(:)
   real(dp) :: c_mod(1:N)
   real(dp) :: factor1, factor2, a3_mod
   integer  :: i

   allocate(beta(1:N), gamma(1:N))

  
   c_mod(1:N) = c3(1:N)
   beta(1)  = b3(1)
   gamma(1) = d3(1) / beta(1)

   factor1  = a3(2) / beta(1)
   beta(2)  = b3(2)  - factor1 * c_mod(1)
   c_mod(2) = c_mod(2) - factor1 * f3(1)
   gamma(2) = (d3(2) - a3(2) * gamma(1)) / beta(2)

   
   do i = 3, N

      !eliminate 2nd subdiagonal using row i-2
      factor2  = e3(i) / beta(i-2)
      a3_mod   = a3(i) - factor2 * c_mod(i-2)

      !eliminate modified 1st subdiagonal using row i-1
      factor1  = a3_mod / beta(i-1)
      beta(i)  = b3(i) - factor2 * f3(i-2) - factor1 * c_mod(i-1)
      if (i <= N-1) c_mod(i) = c_mod(i) - factor1 * f3(i-1)

      gamma(i) = (d3(i) - e3(i) * gamma(i-2) - a3_mod * gamma(i-1)) / beta(i)

   enddo

   ! Back substitution
   xa(N)   = gamma(N)
   xa(N-1) = gamma(N-1) - c_mod(N-1) * xa(N) / beta(N-1)

   do i = N-2, 1, -1
      xa(i) = gamma(i) - (c_mod(i) * xa(i+1) + f3(i) * xa(i+2)) / beta(i)
   enddo

   deallocate(beta, gamma)

end subroutine pdma



end module fdm
