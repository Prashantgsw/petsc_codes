module read_write
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
use petscsys
use petscvec
use double
use variables
implicit none
contains

subroutine save_solution(iter, ctx)
implicit none
PetscInt       :: iter
type(tsdata)   :: ctx
VecScatter     :: vsc
Vec            :: uall
character(30) :: filename
PetscScalar, pointer :: ua(:)
integer       :: i
PetscOffset   :: index_u

call VecScatterCreateToZero(ug, vsc, uall, ierr)
CHKERRQ(ierr)
call VecScatterBegin(vsc, ug, uall, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterEnd(vsc, ug, uall, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterDestroy(vsc, ierr); CHKERRQ(ierr)

if(rank==0)then
  call VecGetArrayF90(uall, ua, ierr); CHKERRQ(ierr)

  write(filename, '(a,i7.7,a)') 'solution_', iter, '.dat'
  open(10,file=trim(filename))
  write(10,*) '# x, rho, rho*u, rhoE'
  do i = 1, ctx%g%Np
     write(10, *) ctx%g%xmin+(i-1)*ctx%g%dx, &
                  ua(dof*(i-1)+1), ua(dof*(i-1)+2), ua(dof*(i-1)+3)
  enddo  
  close(10)
  call VecRestoreArrayF90(uall, ua, ierr); CHKERRQ(ierr)
endif
call VecDestroy(uall,ierr); CHKERRQ(ierr)

end subroutine save_solution

subroutine set_and_braodcast_parameters(ctx)
implicit none
type(tsdata)   :: ctx
 
if(rank == 0)then

   print*, 'please select a finite difference method in space'
   print*, 'Available options: CD2, CD4, CD6, CD8, LELE_CD4, LELE_CD6, LELE_CD8, LELE_CD10'
   read(*,*) space_disc

   print*, "Enter stencil width"
   read(*,*) stencil_width
   print*, "Enter artificial dissipation coefficient (nu_coeff)"
   read(*,*) nu_coeff

endif

call MPI_Barrier(PETSC_COMM_WORLD, ierr)

call MPI_Bcast(stencil_width, 1, MPI_int, 0, PETSC_COMM_WORLD, ierr)
call MPI_Bcast(space_disc,64,MPI_CHARACTER,0,PETSC_COMM_WORLD,ierr)
call MPI_Bcast(nu_coeff, 1, MPI_DOUBLE_PRECISION, 0, PETSC_COMM_WORLD, ierr)

end subroutine set_and_braodcast_parameters

subroutine log_parameters(ctx)
implicit none
type(tsdata)   :: ctx

! Check that floating points have correct precision
if (precision(pi) .ne. precision(pp)) then
  print*, 'Mismatch in float precision'
  stop
endif

print*,'*************About PETSc*************************'
      write(*,25) PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, &
                  PETSC_VERSION_SUBMINOR ! , PETSC_VERSION_PATCH
25    format(' PETSc Release Version: ',i1,'.',i2,'.',i1) ! ,' Patch:',i2)
print*,'****************************************************'

print*,'*************Discretization*************************'
print*, 'No of points = ', ctx%g%Np
print*, 'grid size = ', ctx%g%dx
print*, 'CFL number = ', ctx%g%cfl
print*, 'time step = ', ctx%g%dt
print*,'****************************************************'

print*,'*************Initial condition**********************'
print*, 'Amplitude of Gaussian Function= ', amplitude
print*, 'exponent = ', alpha
print*, 'centered around = ', x_0
print*,'****************************************************'

print*,'*************Finite difference methods**************'
print*, 'Use PETSc TS = ', petsc_ts 
print*, 'FDM in space = ', space_disc
if(.not. petsc_ts) print*, 'FDM in time = ', time_scheme
print*,'****************************************************'

print*,'*************Parallel computations******************'
print*, 'number of processes = ', nproc
print*, 'stencil width = ', stencil_width
print*,'****************************************************'

end subroutine log_parameters

end module read_write
