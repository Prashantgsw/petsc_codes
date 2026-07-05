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
real(dp)     :: rho, u, p

if (xp < x_disc) then
   rho = rho_L; u = u_L; p = p_L
else
   rho = rho_R; u = u_R; p = p_R
endif

fun(1) = rho
fun(2) = rho * u
fun(3) = p/(gamma_gas - 1.0d0) + 0.5d0*rho*u*u

end subroutine initial_condition

subroutine exact_solution(xp, fun, ctx)
implicit none
PetscScalar  :: xp, fun(3)
type(tsdata) :: ctx
fun(1) = amplitude * exp(-alpha*((xp-speed*ctx%g%final_time) - x_0)**2)
fun(2) = fun(1)
fun(3) = fun(1)
end subroutine exact_solution

! Dirichlet BC: fix ghost cells to the true left/right conservative states.
subroutine ApplyPhysicalBC(ua_bc, ng)
implicit none
integer     :: ng
PetscScalar :: ua_bc(3, gist:gien)
integer     :: k
real(dp)    :: rhoE_L, rhoE_R

rhoE_L = p_L/(gamma_gas - 1.0d0) + 0.5d0*rho_L*u_L*u_L
rhoE_R = p_R/(gamma_gas - 1.0d0) + 0.5d0*rho_R*u_R*u_R

 do k = 1, ng
   ua_bc(1, ist-k) = rho_L
   ua_bc(2, ist-k) = rho_L*u_L
   ua_bc(3, ist-k) = rhoE_L

   ua_bc(1, ien+k) = rho_R
   ua_bc(2, ien+k) = rho_R*u_R
   ua_bc(3, ien+k) = rhoE_R
 enddo

end subroutine ApplyPhysicalBC

!Conservative to primitive variable to flux vector
subroutine compute_euler_flux(u, flux, ctx)
implicit none
type(tsdata) :: ctx
PetscScalar  :: u(3,gist:gien), flux(3,gist:gien)
integer      :: i
real(dp)     :: rho, mom, rhoE, vel, p

do i = gist, gien
   rho  = u(1,i)
   mom  = u(2,i)
   rhoE = u(3,i)
   vel  = mom / rho
   p    = (gamma_gas - 1.0d0) * (rhoE - 0.5d0*rho*vel*vel)

   flux(1,i) = mom
   flux(2,i) = mom*vel + p
   flux(3,i) = vel*(rhoE + p)
enddo

end subroutine compute_euler_flux

end module auxillary_conditions