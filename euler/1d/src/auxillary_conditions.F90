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

!Exact Sod shock tube solution

real(dp) function f_K(p, rhoK, pK, cK)
implicit none
real(dp) :: p, rhoK, pK, cK, AK, BK
if (p > pK) then
   AK = 2.0d0/((gamma_gas+1.0d0)*rhoK)
   BK = (gamma_gas-1.0d0)/(gamma_gas+1.0d0)*pK
   f_K = (p-pK)*sqrt(AK/(p+BK))
else
   f_K = (2.0d0*cK/(gamma_gas-1.0d0)) * ( (p/pK)**((gamma_gas-1.0d0)/(2.0d0*gamma_gas)) - 1.0d0 )
endif
end function f_K

real(dp) function f_K_deriv(p, rhoK, pK, cK)
implicit none
real(dp) :: p, rhoK, pK, cK, AK, BK
if (p > pK) then
   AK = 2.0d0/((gamma_gas+1.0d0)*rhoK)
   BK = (gamma_gas-1.0d0)/(gamma_gas+1.0d0)*pK
   f_K_deriv = sqrt(AK/(p+BK)) * (1.0d0 - (p-pK)/(2.0d0*(p+BK)))
else
   f_K_deriv = (1.0d0/(rhoK*cK)) * (p/pK)**(-(gamma_gas+1.0d0)/(2.0d0*gamma_gas))
endif
end function f_K_deriv

subroutine solve_star_region(p_star, u_star)
implicit none
real(dp), intent(out) :: p_star, u_star
real(dp) :: c_L, c_R, p_old, f_total, f_deriv_total
integer  :: iter
real(dp), parameter :: tol = 1.d-10

c_L = sqrt(gamma_gas*p_L/rho_L)
c_R = sqrt(gamma_gas*p_R/rho_R)

! Initial guess (linearized, two-shock-like)
p_old = max(tol, 0.5d0*(p_L+p_R) - 0.125d0*(u_R-u_L)*(rho_L+rho_R)*(c_L+c_R))

do iter = 1, 50
   f_total       = f_K(p_old, rho_L, p_L, c_L) + f_K(p_old, rho_R, p_R, c_R) + (u_R-u_L)
   f_deriv_total = f_K_deriv(p_old, rho_L, p_L, c_L) + f_K_deriv(p_old, rho_R, p_R, c_R)
   p_star = p_old - f_total/f_deriv_total
   if (p_star < tol) p_star = tol
   if (abs(p_star-p_old) < tol) exit
   p_old = p_star
enddo

u_star = 0.5d0*(u_L+u_R) + 0.5d0*( f_K(p_star, rho_R, p_R, c_R) - f_K(p_star, rho_L, p_L, c_L) )

end subroutine solve_star_region

subroutine sample_sod(S, p_star, u_star, rho, u, p)
implicit none
real(dp), intent(in)  :: S, p_star, u_star
real(dp), intent(out) :: rho, u, p
real(dp) :: c_L, c_R, rho_star, S_wave, S_head, S_tail, c_star

c_L = sqrt(gamma_gas*p_L/rho_L)
c_R = sqrt(gamma_gas*p_R/rho_R)

if (S <= u_star) then
   ! left of contact: sample left wave (shock or rarefaction)
   if (p_star > p_L) then
      ! left shock
      rho_star = rho_L * ( (p_star/p_L + (gamma_gas-1.0d0)/(gamma_gas+1.0d0)) / &
                            ((gamma_gas-1.0d0)/(gamma_gas+1.0d0)*p_star/p_L + 1.0d0) )
      S_wave = u_L - c_L*sqrt((gamma_gas+1.0d0)/(2.0d0*gamma_gas)*p_star/p_L + (gamma_gas-1.0d0)/(2.0d0*gamma_gas))
      if (S < S_wave) then
         rho = rho_L; u = u_L; p = p_L
      else
         rho = rho_star; u = u_star; p = p_star
      endif
   else
      ! left rarefaction
      c_star = c_L*(p_star/p_L)**((gamma_gas-1.0d0)/(2.0d0*gamma_gas))
      S_head = u_L - c_L
      S_tail = u_star - c_star
      if (S < S_head) then
         rho = rho_L; u = u_L; p = p_L
      elseif (S > S_tail) then
         rho = rho_L*(p_star/p_L)**(1.0d0/gamma_gas); u = u_star; p = p_star
      else
         ! inside the fan
         u   = 2.0d0/(gamma_gas+1.0d0) * ( c_L + (gamma_gas-1.0d0)/2.0d0*u_L + S )
         c_star = 2.0d0/(gamma_gas+1.0d0) * ( c_L + (gamma_gas-1.0d0)/2.0d0*(u_L - S) )
         rho = rho_L*(c_star/c_L)**(2.0d0/(gamma_gas-1.0d0))
         p   = p_L*(c_star/c_L)**(2.0d0*gamma_gas/(gamma_gas-1.0d0))
      endif
   endif
else
   ! right of contact: sample right wave
   if (p_star > p_R) then
      ! right shock
      rho_star = rho_R * ( (p_star/p_R + (gamma_gas-1.0d0)/(gamma_gas+1.0d0)) / &
                            ((gamma_gas-1.0d0)/(gamma_gas+1.0d0)*p_star/p_R + 1.0d0) )
      S_wave = u_R + c_R*sqrt((gamma_gas+1.0d0)/(2.0d0*gamma_gas)*p_star/p_R + (gamma_gas-1.0d0)/(2.0d0*gamma_gas))
      if (S > S_wave) then
         rho = rho_R; u = u_R; p = p_R
      else
         rho = rho_star; u = u_star; p = p_star
      endif
   else
      ! right rarefaction
      c_star = c_R*(p_star/p_R)**((gamma_gas-1.0d0)/(2.0d0*gamma_gas))
      S_head = u_R + c_R
      S_tail = u_star + c_star
      if (S > S_head) then
         rho = rho_R; u = u_R; p = p_R
      elseif (S < S_tail) then
         rho = rho_R*(p_star/p_R)**(1.0d0/gamma_gas); u = u_star; p = p_star
      else
         u   = 2.0d0/(gamma_gas+1.0d0) * ( -c_R + (gamma_gas-1.0d0)/2.0d0*u_R + S )
         c_star = 2.0d0/(gamma_gas+1.0d0) * ( c_R - (gamma_gas-1.0d0)/2.0d0*(u_R - S) )
         rho = rho_R*(c_star/c_R)**(2.0d0/(gamma_gas-1.0d0))
         p   = p_R*(c_star/c_R)**(2.0d0*gamma_gas/(gamma_gas-1.0d0))
      endif
   endif
endif

end subroutine sample_sod

subroutine exact_solution(xp, fun, ctx)
implicit none
PetscScalar  :: xp, fun(3)
type(tsdata) :: ctx
real(dp)     :: p_star, u_star, S, rho, u, p

call solve_star_region(p_star, u_star)

if (ctx%g%time > 0.0d0) then
   S = (xp - x_disc) / ctx%g%time
else
   S = 0.0d0   ! at t=0 just return the initial discontinuity itself
   if (xp < x_disc) then
      rho = rho_L; u = u_L; p = p_L
   else
      rho = rho_R; u = u_R; p = p_R
   endif
endif

if (ctx%g%time > 0.0d0) call sample_sod(S, p_star, u_star, rho, u, p)

fun(1) = rho
fun(2) = rho*u
fun(3) = p/(gamma_gas-1.0d0) + 0.5d0*rho*u*u

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