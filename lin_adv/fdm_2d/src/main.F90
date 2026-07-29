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
  PetscScalar        :: chi_p, eta_p, fun
  ! datatypes
  real(dp)           :: runtime
  character(len=128) :: fmt1, fmt2, fmt3, fmt4
  integer            :: values(8)
  integer            :: i, j
  PetscScalar, pointer :: u(:,:)
  type(tsdata)       :: ctx
  
  ! initialize petsc and read input file
  call PetscInitialize('param.in', ierr); CHKERRQ(ierr)
  if(ierr /= 0) stop "PETSc not intialized"
  
  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr); CHKERRQ(ierr)
  call MPI_Comm_size(PETSC_COMM_WORLD, nproc, ierr); CHKERRQ(ierr)
  runtime = MPI_Wtime()
  fmt1 = "(1x, 'Date ', i2.2, '/', i2.2, '/', i4.4, '; time ', &
         & i2.2, ':', i2.2, ':', i2.2 )"
  
  if(rank == 0)then
     call date_and_time(VALUES=values)
     write(*,fmt1) values(3), values(2), values(1), values(5:7)
     print*,'Running git version: ', VERSION
     print*,'Number of processes: ', nproc
  endif
  
  call get_and_set_parameters()
  
  ! Creates an object that will manage the communication of 2D regular array that is distributed across some processors. 
  if(stencil_type == 1)then
    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR, N_chi, N_eta, PETSC_DECIDE, PETSC_DECIDE, 1, stencil_width, &
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, da, ierr); CHKERRQ(ierr)
  elseif(stencil_type == 2)then
    call DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, N_chi, N_eta, PETSC_DECIDE, PETSC_DECIDE, 1, stencil_width, &
                      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, da, ierr); CHKERRQ(ierr)        
  else
    write(*,*) 'please select a stencil type: 1 or 2'
    stop
  endif
  call DMSetFromOptions(da, ierr); CHKERRQ(ierr)
  call DMSetUp(da, ierr); CHKERRQ(ierr)
  call DMCreateGlobalVector(da, ug, ierr); CHKERRQ(ierr)
  ! Returns the global (x,y,z) indices of the lower left corner and size of the local region, excluding ghost points.
  call DMDAGetCorners(da, ibeg, jbeg, PETSC_NULL_INTEGER, &
                          N_chi_loc, N_eta_loc, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
  call DMDASetUniformCoordinates(da, chi_min, chi_max, eta_min, eta_max, 0.d0, 0.d0, ierr); CHKERRQ(ierr)
  ! This is used to control no of grid points from command line
  call DMDAGetInfo(da, PETSC_NULL_INTEGER, N_chi, N_eta, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                   PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
  call DMDAGetGhostCorners(da, ibeg_ghosted, jbeg_ghosted, PETSC_NULL_INTEGER, N_chi_loc_ghosted, N_eta_loc_ghosted, PETSC_NULL_INTEGER, ierr); CHKERRQ(ierr)
  
  ist = ibeg ; ien = ibeg+N_chi_loc-1
  jst = jbeg ; jen = jbeg+N_eta_loc-1
  gist = ibeg_ghosted ; gien = ibeg_ghosted+N_chi_loc_ghosted-1
  gjst = jbeg_ghosted ; gjen = jbeg_ghosted+N_eta_loc_ghosted-1
  
  ! setting equidistant grid and initial condition on it.
  call DMDAVecGetArrayF90(da, ug, u, ierr); CHKERRQ(ierr)
  dchi = (chi_max - chi_min) / dble(N_chi) !deta is dy and dchi is dx
  deta = (eta_max - eta_min) / dble(N_eta)
  do j = jst, jen
  do i = ist, ien
     chi_p = chi_min+(i-1)*dchi ; eta_p = eta_min+(j-1)*deta
     call initial_condition(chi_p, eta_p, fun)
     u(i,j) = fun
  enddo
  enddo
  call DMDAVecRestoreArrayF90(da, ug, u, ierr); CHKERRQ(ierr)
  
  ! settin time stepping
  time = 0.d0
  iter = 0
  dt = cfl * dsqrt(dchi*dchi+deta*deta) / (dsqrt(speed_chi**2+speed_eta**2) + 1.d-13)  
  
  if(petsc_ts)then ! solve using PETSc time stepping
  
    call TSCreate(PETSC_COMM_WORLD, ts, ierr); CHKERRQ(ierr)
    call TSSetDM(ts,da,ierr); CHKERRQ(ierr)
    ! by default: nonlinear problem
    call TSSetProblemType(ts, TS_NONLINEAR, ierr); CHKERRQ(ierr)
    ! in RHSfunction finite difference methods are applied
    call TSSetRHSFunction(ts, PETSC_NULL_VEC, RHSFunction, ctx, ierr); CHKERRQ(ierr)
    call TSSetTime(ts, 0.d0, ierr); CHKERRQ(ierr)
    call TSSetTimeStep(ts, dt, ierr); CHKERRQ(ierr)
    call TSSetType(ts, TSSSP, ierr); CHKERRQ(ierr)
    call TSSetMaxTime(ts, final_time, ierr); CHKERRQ(ierr)
    call TSSetMaxSteps(ts, itmax, ierr); CHKERRQ(ierr)
    call TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP, ierr); CHKERRQ(ierr)
    call TSSetSolution(ts, ug, ierr); CHKERRQ(ierr)
    ! in Monitor function, time stepping may be computed, solution maybe monitored for norms and solution may be written in files
    call TSMonitorSet(ts, Monitor, ctx, PETSC_NULL_FUNCTION, ierr); CHKERRQ(ierr)
    ! sets the PETSc-implemented time integration methods from command line
    call TSGetType(ts, time_scheme, ierr); CHKERRQ(ierr)
    call TSSetFromOptions(ts, ierr); CHKERRQ(ierr)
    call TSSetUp(ts, ierr); CHKERRQ(ierr)
    call TSSolve(ts, ug, ierr); CHKERRQ(ierr)
  
  else ! solve by implementing time integration methods

    print*, 'Time stepping yet to be implemented'
    stop

  endif
  
  ! Print elapsed time
  runtime = MPI_Wtime() - runtime
  fmt4 = "(1x,'Iter runtime min =',e20.12)"
  if(rank == 0)then
     write(*,fmt4) runtime/60.0
  endif
  
  !compute l2 norm of error
  call VecDuplicate(ug, ue, ierr); CHKERRQ(ierr)
  call DMDAVecGetArrayF90(da, ue, u, ierr); CHKERRQ(ierr)
  do j = jst, jen
  do i = ist, ien
     chi_p = chi_min+(i-1)*dchi ; eta_p = eta_min+(j-1)*deta
     call exact_solution(chi_p, eta_p, fun)
     u(i,j) = fun
  enddo
  enddo
  call DMDAVecRestoreArrayF90(da, ue, u, ierr); CHKERRQ(ierr)
  call VecAXPY(ue, -1.d0, ug, ierr); CHKERRQ(ierr)
  call VecNorm(ue, NORM_2, err_l2, ierr); CHKERRQ(ierr)
  print*, sqrt(dchi*deta), err_l2 * dsqrt(1.d0/(N_chi*N_eta))
  
  ! destroy petsc objects
  call VecDestroy(ue, ierr); CHKERRQ(ierr)
  call VecDestroy(ug, ierr); CHKERRQ(ierr)
  call DMDestroy(da, ierr); CHKERRQ(ierr)
  if(petsc_ts)then
    call TSDestroy(ts, ierr); CHKERRQ(ierr)
  endif
  call PetscFinalize(ierr); CHKERRQ(ierr)
end program main
