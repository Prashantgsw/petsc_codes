module solve
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscts.h>
#include <petsc/finclude/petscdmda.h>
use petscsys
use petscts
use petscdmda
use double
use variables
use fdm
use read_write
implicit none
contains

! Compute rhs function in y' = r(t,y)
subroutine RHSFunction(ts, time, yg, r, ctx, ierr)
#include <petsc/finclude/petscts.h>
#include <petsc/finclude/petscdmda.h>
use petscts
use petscdmda
implicit none
TS             :: ts
DM             :: da
PetscReal      :: time
Vec            :: yg
Vec            :: r
type(tsdata)   :: ctx
Vec            :: localU
PetscScalar, pointer :: u(:), res(:)
PetscErrorCode :: ierr
PetscOffset :: index_u, index_res

call TSGetDM(ts, da, ierr)
CHKERRQ(ierr)
call DMGetLocalVector(da, localU, ierr)
CHKERRQ(ierr)
call DMGlobalToLocalBegin(da, yg, INSERT_VALUES, localU, ierr)
CHKERRQ(ierr)
call DMGlobalToLocalEnd(da, yg, INSERT_VALUES, localU, ierr)
CHKERRQ(ierr)
call DMDAVecGetArrayReadF90(da, localU, u, ierr)
CHKERRQ(ierr)
call DMDAVecGetArrayF90(da, r, res, ierr)
CHKERRQ(ierr)
call finite_diffence_method(u, res, ctx)
call DMDAVecRestoreArrayReadF90(da, localU, u, ierr)
CHKERRQ(ierr)
call DMDAVecRestoreArrayF90(da, r, res, ierr)
CHKERRQ(ierr)
call DMRestoreLocalVector(da, localU, ierr)
CHKERRQ(ierr)

end subroutine RHSFunction

! This subroutine is called after every time step
! Note: This function seems to be called before time stepping starts.
subroutine Monitor(ts, step, time, u, ctx, ierr)
use petscts
implicit none
TS             :: ts
PetscInt       :: step
PetscReal      :: time, u_l2
Vec            :: u
type(tsdata)   :: ctx
PetscErrorCode :: ierr

call VecNorm(u, NORM_2, u_l2, ierr); CHKERRQ(ierr)

if(rank == 0)then
  write(*,'(a, i4, 3F8.4)') 'iter, dt, time, l2_norm = ', step, ctx%g%dt, time, u_l2
endif

if( mod(step, ctx%g%itsave) == 0 ) call save_solution(step, ctx)

end subroutine Monitor

end module solve
