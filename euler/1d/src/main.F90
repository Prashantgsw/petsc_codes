! solves linear 1D convection equation: u_t + c * u_x = 0 for periodic boundary conditions:
! using finite difference methods, PETSc DMDA and time stepping.
! written by: Ashish Bhole.
program main
! need to include petsc files as well as use modules
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscts.h>
use petscsys
use petscdm
use petscdmda
use petscts
! add modules here
use double
use variables
use auxillary_conditions
use read_write
use solve
use fdm
implicit none
! petsc datatypes
PetscInt           :: loc
PetscScalar        :: xp, fun
! datatypes
real(dp)           :: runtime
character(len=128) :: fmt1, fmt2, fmt3, fmt4
integer            :: values(8)
integer            :: i, iter
type(grid)         :: g
type(tsdata)       :: ctx

call PetscInitialize(PETSC_NULL_CHARACTER, ierr) 
CHKERRQ(ierr)
if(ierr /= 0) stop "PETSc not intialized"

call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
CHKERRQ(ierr)
call MPI_Comm_size(PETSC_COMM_WORLD, nproc, ierr)
CHKERRQ(ierr)
runtime = MPI_Wtime()
fmt1 = "(1x, 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', &
       & i2.2, ':', i2.2, ':', i2.2 )"

if(rank == 0)then
   call date_and_time(VALUES=values)
   write(*,fmt1) values(3), values(2), values(1), values(5:7)
   print*,'Running git version: ', VERSION
   print*,'Number of processes: ', nproc
endif

call set_and_braodcast_parameters(ctx)

! Creates an object that will manage the communication of one-dimensional regular array data
! that is distributed across some processors. 
call DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, ctx%g%Np, 1, stencil_width, &
                  PETSC_NULL_INTEGER, da, ierr) 
CHKERRQ(ierr)
call DMSetFromOptions(da, ierr) 
CHKERRQ(ierr)
call DMSetUp(da, ierr); CHKERRQ(ierr)
call DMCreateGlobalVector(da, ug, ierr)
CHKERRQ(ierr)
! Returns the global (x,y,z) indices of the lower left corner and size of the local region, excluding ghost points.
call DMDAGetCorners(da, ctx%g%ibeg, PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, ctx%g%nloc, &
                    PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
CHKERRQ(ierr)
! This is used to control no of grid points from command line
call DMDAGetInfo(da, PETSC_NULL_INTEGER, ctx%g%Np, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
        PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
        PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
        PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
CHKERRQ(ierr)
call DMDAGetGhostCorners(da, ctx%g%ibeg_ghosted, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                             ctx%g%nloc_ghosted, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
CHKERRQ(ierr)
ist = ctx%g%ibeg+1
gist = ctx%g%ibeg_ghosted+1
ien = ctx%g%ibeg+ctx%g%nloc ! -1
gien = ctx%g%ibeg_ghosted+ctx%g%nloc_ghosted ! -1
! setting equidistant grid and initial condition on it.
ctx%g%dx = (ctx%g%xmax - ctx%g%xmin) / dble(ctx%g%Np)
do i = ist, ien
   xp = (i-1)*ctx%g%dx ; loc = i-1
   call initial_condition(xp, fun)
   call VecSetValues(ug, one, loc, fun, INSERT_VALUES, ierr)
   CHKERRQ(ierr)
enddo
call VecAssemblyBegin(ug, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(ug, ierr); CHKERRQ(ierr)
! settin time stepping
ctx%g%time = 0.d0
ctx%g%dt   = ctx%g%cfl * ctx%g%dx / (speed + 1.d-13)
ctx%g%iter = 0

if(rank==0) call log_parameters(ctx)

if(petsc_ts)then ! solve using PETSc time stepping

  call TSCreate(PETSC_COMM_WORLD, ts, ierr); CHKERRQ(ierr)
  call TSSetDM(ts,da,ierr); CHKERRQ(ierr)
  ! by default: nonlinear problem
  call TSSetProblemType(ts, TS_NONLINEAR, ierr); CHKERRQ(ierr)
  ! in RHSfunction finite differencing in space takes place
  call TSSetRHSFunction(ts, PETSC_NULL_VEC, RHSFunction, ctx, ierr)
  CHKERRQ(ierr)
  call TSSetTime(ts, 0.d0, ierr); CHKERRQ(ierr)
  call TSSetTimeStep(ts, ctx%g%dt, ierr); CHKERRQ(ierr)
  call TSSetType(ts, TSSSP, ierr); CHKERRQ(ierr);
  call TSSetMaxTime(ts, ctx%g%final_time, ierr); CHKERRQ(ierr);
  call TSSetMaxSteps(ts, ctx%g%itmax, ierr); CHKERRQ(ierr);
  call TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP, ierr)
  CHKERRQ(ierr)
  ! in Monitor function, time stepping may be computed, solution maybe monitored for norms
  ! and solution may be written in files
  call TSSetSolution(ts, ug, ierr); CHKERRQ(ierr)
  call TSMonitorSet(ts, Monitor, ctx, PETSC_NULL_FUNCTION, ierr)
  CHKERRQ(ierr)
  ! enables the use of different PETSc implemented time integration methods from command line
  call TSGetType(ts, time_scheme, ierr); CHKERRQ(ierr)
  call TSSetFromOptions(ts, ierr); CHKERRQ(ierr)
  call TSSetUp(ts, ierr); CHKERRQ(ierr)
  call TSSolve(ts, ug, ierr); CHKERRQ(ierr)

else ! solve by implementing time integration methods

  call save_solution(ctx%g%iter, ctx)      
  do while (ctx%g%time .lt. ctx%g%final_time)
    if (ctx%g%time + ctx%g%dt .gt. ctx%g%final_time) then
    endif
    !call solve_num_meth()
    ctx%g%time = ctx%g%time + ctx%g%dt
    ctx%g%iter = ctx%g%iter + 1
    call save_solution(ctx%g%iter, ctx)
    print*, 'time = ', ctx%g%time
  enddo
  call save_solution(ctx%g%iter, ctx)

endif

! Print elapsed time
runtime = MPI_Wtime() - runtime
fmt4 = "(1x,'Iter runtime min =',e20.12)"
if(rank == 0)then
   write(*,fmt4) runtime/60.0
endif

call VecDuplicate(ug, ue, ierr); CHKERRQ(ierr)
do i = ist, ien
   xp = (i-1)*ctx%g%dx ; loc = i-1
   call exact_solution(xp, fun, ctx)
   call VecSetValues(ue, one, loc, fun, INSERT_VALUES, ierr)
   CHKERRQ(ierr)
enddo
call VecAssemblyBegin(ue, ierr); CHKERRQ(ierr)
call VecAssemblyEnd(ue, ierr); CHKERRQ(ierr)
call VecAXPY(ue, -1.d0, ug, ierr); CHKERRQ(ierr)
call VecNorm(ue, NORM_2, err_l2, ierr); CHKERRQ(ierr)
print*, ctx%g%Np, err_l2 * dsqrt(1.d0/ctx%g%Np)

call VecDestroy(ue, ierr); CHKERRQ(ierr)
call VecDestroy(ug, ierr); CHKERRQ(ierr)
call DMDestroy(da, ierr); CHKERRQ(ierr)
if(petsc_ts)then
  call TSDestroy(ts, ierr); CHKERRQ(ierr)
endif
call PetscFinalize(ierr); CHKERRQ(ierr)

end program main
