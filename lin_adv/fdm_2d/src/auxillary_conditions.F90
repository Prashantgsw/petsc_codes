module auxillary_conditions
#include <petsc/finclude/petscsys.h>
use petscsys
use double
use variables
implicit none
contains

subroutine initial_condition(chi_p, eta_p, fun)
  implicit none
  integer :: i
  PetscScalar  :: chi_p, eta_p, fun
  fun = amplitude * exp(-alpha*((chi_p - chi_0)**2 + (eta_p - eta_0)**2))
end subroutine initial_condition

subroutine exact_solution(chi_p, eta_p, fun)
  implicit none
  integer :: i
  PetscScalar  :: chi_p, eta_p, fun
  fun = amplitude * exp(-alpha*(  ((chi_p-speed_chi*final_time) - chi_0)**2 &
                                + ((eta_p-speed_eta*final_time) - eta_0)**2 ) )
end subroutine exact_solution

end module auxillary_conditions
