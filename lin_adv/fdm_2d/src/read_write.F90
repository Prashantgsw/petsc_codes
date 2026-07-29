module read_write
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdmda.h>
use petscsys
use petscvec
use petscdmda
use double
use variables
implicit none
contains

! save grid and solution in each rank in a seaparate file        
subroutine save_solution(iter, ug)
  implicit none
  PetscInt       :: iter
  Vec            :: ug
  DM             :: da
  Vec            :: ul
  character(30) :: filename
  PetscScalar, pointer :: ua(:,:)
  integer       :: i, j
  PetscViewer   :: viewer
  
  call TSGetDM(ts, da, ierr); CHKERRQ(ierr)
  call DMGetLocalVector(da, ul, ierr); CHKERRQ(ierr)
  call DMGlobalToLocalBegin(da, ug, INSERT_VALUES, ul, ierr); CHKERRQ(ierr)
  call DMGlobalToLocalEnd(da, ug, INSERT_VALUES, ul, ierr); CHKERRQ(ierr)
  call DMDAVecGetArrayF90(da, ul, ua, ierr); CHKERRQ(ierr)
  write(filename, '(a,i3.3,a,i3.3,a)') 'sol-', rank, '-', iter, '.plt'
  open(10,file=trim(filename))
  write(10,*) "TITLE = solution"
  write(10,*) "VARIABLES = x, y, sol"
  write(10,*) "ZONE STRANDID=1, SOLUTIONTIME=", time ,", I=", ien-ist+1, &
              ", J=", jen-jst+1, ", DATAPACKING=POINT"
  do j = jst, jen
  do i = ist, ien
     write(10, *) chi_min+(i-1)*dchi, eta_min+(j-1)*deta, ua(i,j) !deta = dy and dchi =dx
  enddo
  enddo
  close(10)
  call DMDAVecRestoreArrayF90(da, ul, ua, ierr); CHKERRQ(ierr)
  call DMRestoreLocalVector(da, ul, ierr); CHKERRQ(ierr)
  call VecDestroy(ul,ierr); CHKERRQ(ierr)
end subroutine save_solution

! dumps solution in a hdf5 file. it seems there is no way to dump grid also in the same file.
subroutine save_solution_hdf5(iter, ug)
  implicit none
  PetscInt       :: iter
  Vec            :: ug
  PetscViewer    :: viewer
  character(30) :: filename
  write(filename, '(a,i3.3,a)') 'sol-', iter, '.h5'
  call PetscViewerHDF5Open(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, viewer, ierr); CHKERRQ(ierr)
  call PetscObjectSetName(ug, "SOL", ierr); CHKERRQ(ierr)
  call VecView(ug, viewer, ierr); CHKERRQ(ierr)
  call PetscViewerDestroy(viewer, ierr); CHKERRQ(ierr)
end subroutine save_solution_hdf5

! get parameters from command line and set or braodcast them on all processors
subroutine get_and_set_parameters()
  implicit none
  PetscBool :: set
  character(len=20)  :: fmt1, fmt2, fmt3, fmt4
  
  if(rank==0) print*,'Reading parameters from param.in'

  petsc_ts = .true. ! Defualt: use petsc time stepping  
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, '-petscts', petsc_ts,set,ierr); CHKERRQ(ierr)  

  cfl = 0.5d0 ! Default cfl number
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-cfl',cfl,set,ierr); CHKERRQ(ierr)
  
  final_time = 1.d0 ! Default final time
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-tfinal',final_time,set,ierr); CHKERRQ(ierr)
  
   space_disc = 'CD4' ! Default: CD4. Options: CD2, CD4, CD6, LELE
   call PetscOptionsGetString(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-space_disc',space_disc,set,ierr); CHKERRQ(ierr)
  
  itmax = 10000000 ! Default: large value
  call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-itmax',itmax,set,ierr); CHKERRQ(ierr)
  
  itsave = 10000000 ! Default: large value
  call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-itsave',itsave,set,ierr); CHKERRQ(ierr)
 
  igridformat = 1 ! Default: 1 ascii, 2 hdf5
  call PetscOptionsGetInt (PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-igridformat',igridformat,set,ierr); CHKERRQ(ierr)

  select case(trim(space_disc))
  case('CD2')
    stencil_width = 1; stencil_type = 1
  case('CD4')
    stencil_width = 2; stencil_type = 1
  case('CD6')
    stencil_width = 3; stencil_type = 1
  case('LELE')
    stencil_width = 10; stencil_type = 1
  case default
    print*, 'please select a finite difference method in space'
    print*, 'Available options: CD2, CD4, CD6'
    stop
  end select
  
  speed_chi = 1.d0
  speed_eta = 1.d0 
  
  ! is broadcasting of variables needed?
  !call MPI_Bcast(stencil_width, 1, MPI_int, 0, PETSC_COMM_WORLD, ierr)
  
  ! log paramaters to screen
  if(rank == 0) then
  
    fmt1 = "(5x,a15,i14)"
    fmt2 = "(5x,a15,e14.4)"
    fmt3 = "(5x,a15,L)"
    fmt4 = "(5x,a15,a10)"
    write(*,fmt2) 'tfinal   =     ', final_time
    write(*,fmt3) 'petsc_ts =     ', petsc_ts
    write(*,fmt2) 'cfl      =     ', cfl
    write(*,fmt4) 'space_disc =   ', trim(space_disc)
    write(*,fmt1) 'itmax    =     ', itmax
    write(*,fmt1) 'itsave   =     ', itsave
    
    if(cfl < zero .and. dt < zero)then
       write(*,*)'Both cfl and dt cannot be zero'
    endif
    
    if(cfl > zero .and. dt > zero)then
       write(*,*)'Both cfl and dt cannot be specified'
       write(*,*)'Give one and set other to zero'
    endif
  
  endif
end subroutine get_and_set_parameters

end module read_write
