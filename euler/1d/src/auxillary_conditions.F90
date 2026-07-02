module auxillary_conditions
#include <petsc/finclude/petscsys.h>
use petscsys
use double
use variables
implicit none
contains

subroutine initial_condition(xp, fun)
implicit none
PetscScalar  :: xp, fun(3)
! PLACEHOLDER: same Gaussian in all 3 components, just to verify dof=3 plumbing.
! Real Sod IC (rho, rho*u, E piecewise) comes in step 3.
fun(1) = amplitude * exp(-alpha*(xp - x_0)**2)
fun(2) = fun(1)
fun(3) = fun(1)
end subroutine initial_condition

subroutine exact_solution(xp, fun, ctx)
implicit none
PetscScalar  :: xp, fun(3)
type(tsdata) :: ctx
fun(1) = amplitude * exp(-alpha*((xp-speed*ctx%g%final_time) - x_0)**2)
fun(2) = fun(1)
fun(3) = fun(1)
end subroutine exact_solution

! TEMPORARY ghost-cell fill (constant extrapolation), NOT the final Dirichlet BC.
! Real Dirichlet BC (fixed left/right conservative states) comes in step 5.
subroutine ApplyPhysicalBC(ua_bc, ng)
implicit none
integer     :: ng
PetscScalar :: ua_bc(3, gist:gien)
integer     :: k
do k = 1, ng
   ua_bc(:, ist-k) = ua_bc(:, ist)
   ua_bc(:, ien+k) = ua_bc(:, ien)
enddo
end subroutine ApplyPhysicalBC

end module auxillary_conditions