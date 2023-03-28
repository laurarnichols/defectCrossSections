module shifterMod
  
  use constants, only: dp
  use mpi

  implicit none

  ! Parameters:
  integer, parameter :: root = 0
    !! ID of the root node

  ! Global variables not passed as arguments:
  integer :: ierr
    !! Error returned by MPI
  integer :: myid
    !! ID of this process
  integer :: nProcs
    !! Number of processes
  integer :: worldComm
    !! World communicator

  logical :: ionode
    !! If this node is the root node

  ! Variables that should be passed as arguments:
  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors

  character(len=300) :: phononFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: poscarFName
    !! File name for POSCAR

  contains

!----------------------------------------------------------------------------
  subroutine mpiInitialization()
    !! Generate MPI processes and communicators 
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Output variables:
    !integer, intent(out) :: myid
      ! ID of this process
    !integer, intent(out) :: nProcs
      ! Number of processes
    !integer, intent(out) :: worldComm
      ! World communicator

    !logical, intent(out) :: ionode
      ! If this node is the root node


    call MPI_Init(ierr)
    if (ierr /= 0) call mpiExitError( 8001 )

    worldComm = MPI_COMM_WORLD

    call MPI_COMM_RANK(worldComm, myid, ierr)
    if (ierr /= 0) call mpiExitError( 8002 )
      !! * Determine the rank or ID of the calling process

    call MPI_COMM_SIZE(worldComm, nProcs, ierr)
    if (ierr /= 0) call mpiExitError( 8003 )
      !! * Determine the size of the MPI pool (i.e., the number of processes)

    ionode = (myid == root)
      ! Set a boolean for if this is the root process

    return
  end subroutine mpiInitialization

!----------------------------------------------------------------------------
  subroutine initialize(poscarFName, phononFName, shift)
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    ! Input variables:
    !integer, intent(in) :: nProcs
      ! Number of processes


    ! Output variables:
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: poscarFName
      !! File name for POSCAR


    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time

    poscarFName = 'POSCAR'
    phononFName = 'mesh.yaml'
    shift = 0.01

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(*, '(/5X,"Phonon post-processing: POSCAR shifter starts on ",A9," at ",A9)') &
             cdate, ctime

      write(*, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs


    endif

  end subroutine initialize

end module shifterMod
