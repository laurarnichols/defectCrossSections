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

end module shifterMod
