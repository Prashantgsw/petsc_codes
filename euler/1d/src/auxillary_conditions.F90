module auxillary_conditions
#include <petsc/finclude/petscsys.h>
use petscsys
use double
use variables
implicit none
contains

subroutine initial_condition(xp, fun)
implicit none
integer :: i
PetscScalar  :: xp, fun
fun = amplitude * exp(-alpha*(xp - x_0)**2)
end subroutine initial_condition

subroutine exact_solution(xp, fun, ctx)
implicit none
integer :: i
PetscScalar  :: xp, fun
type(tsdata) :: ctx
fun = amplitude * exp(-alpha*((xp-speed*ctx%g%final_time) - x_0)**2)
end subroutine exact_solution

end module auxillary_conditions
