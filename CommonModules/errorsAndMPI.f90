module errorsAndMPI
  
  use constants, only: dp
  use mpi

  implicit none

  integer, parameter :: root = 0
    !! ID of the root node

  integer :: ierr
    !! Error handler
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

!----------------------------------------------------------------------------
  subroutine mpiSumIntV(msg, comm)
    !! Perform `MPI_ALLREDUCE` sum for an integer vector
    !! using a max buffer size
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input/output variables:
    integer, intent(in) :: comm
      !! MPI communicator
    integer, intent(inout) :: msg(:)
      !! Message to be sent


    ! Local variables:
    integer, parameter :: maxb = 100000
      !! Max buffer size

    integer :: ib
      !! Loop index
    integer :: buff(maxb)
      !! Buffer
    integer :: msglen
      !! Length of message to be sent
    integer :: nbuf
      !! Number of buffers

    msglen = size(msg)

    nbuf = msglen/maxb
      !! * Get the number of buffers of size `maxb` needed
  
    do ib = 1, nbuf
      !! * Send message in buffers of size `maxb`
     
        call MPI_ALLREDUCE(msg(1+(ib-1)*maxb), buff, maxb, MPI_INTEGER, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumIntV', 'error in mpi_allreduce 1', ierr)

        msg((1+(ib-1)*maxb):(ib*maxb)) = buff(1:maxb)

    enddo

    if((msglen - nbuf*maxb) > 0 ) then
      !! * Send any data left of size less than `maxb`

        call MPI_ALLREDUCE(msg(1+nbuf*maxb), buff, (msglen-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumIntV', 'error in mpi_allreduce 2', ierr)

        msg((1+nbuf*maxb):msglen) = buff(1:(msglen-nbuf*maxb))
    endif

    return
  end subroutine mpiSumIntV

!----------------------------------------------------------------------------
  subroutine mpiSumDoubleV(msg, comm)
    !! Perform `MPI_ALLREDUCE` sum for a double-precision vector
    !! using a max buffer size
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input/output variables:
    integer, intent(in) :: comm
      !! MPI communicator

    real(kind=dp), intent(inout) :: msg(:)
      !! Message to be sent


    ! Local variables:
    integer, parameter :: maxb = 20000
      !! Max buffer size

    integer :: ib
      !! Loop index
    integer :: msglen
      !! Length of message to be sent
    integer :: nbuf
      !! Number of buffers

    real(kind=dp) :: buff(maxb)
      !! Buffer

    msglen = size(msg)

    nbuf = msglen/maxb
      !! * Get the number of buffers of size `maxb` needed

    do ib = 1, nbuf
      !! * Send message in buffers of size `maxb`
     
        call MPI_ALLREDUCE(msg(1+(ib-1)*maxb), buff, maxb, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumDoubleV', 'error in mpi_allreduce 1', ierr)

        msg((1+(ib-1)*maxb):(ib*maxb)) = buff(1:maxb)

    enddo

    if((msglen - nbuf*maxb) > 0 ) then
      !! * Send any data left of size less than `maxb`

        call MPI_ALLREDUCE(msg(1+nbuf*maxb), buff, (msglen-nbuf*maxb), MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumDoubleV', 'error in mpi_allreduce 2', ierr)

        msg((1+nbuf*maxb):msglen) = buff(1:(msglen-nbuf*maxb))
    endif

    return
  end subroutine mpiSumDoubleV

!----------------------------------------------------------------------------
  subroutine mpiSumComplexV(msg, comm)
    !! Perform `MPI_ALLREDUCE` sum for a complex vector
    !! using a max buffer size
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input/output variables:
    integer, intent(in) :: comm
      !! MPI communicator

    complex(kind=dp), intent(inout) :: msg(:)
      !! Message to be sent


    ! Local variables:
    integer, parameter :: maxb = 10000
      !! Max buffer size

    integer :: ib
      !! Loop index
    integer :: msglen
      !! Length of message to be sent
    integer :: nbuf
      !! Number of buffers
    integer :: commSize

    complex(kind=dp) :: buff(maxb)
      !! Buffer


    msglen = size(msg)

    nbuf = msglen/maxb
      !! * Get the number of buffers of size `maxb` needed
  
    do ib = 1, nbuf
      !! * Send message in buffers of size `maxb`
     
        call MPI_ALLREDUCE(msg(1+(ib-1)*maxb), buff, maxb, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumComplexV', 'error in mpi_allreduce 1', ierr)

        msg((1+(ib-1)*maxb):(ib*maxb)) = buff(1:maxb)

    enddo

    if((msglen - nbuf*maxb) > 0 ) then
      !! * Send any data left of size less than `maxb`

        call MPI_ALLREDUCE(msg(1+nbuf*maxb), buff, (msglen-nbuf*maxb), MPI_DOUBLE_COMPLEX, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumComplexV', 'error in mpi_allreduce 2', ierr)

        msg((1+nbuf*maxb):msglen) = buff(1:(msglen-nbuf*maxb))
    endif

    return
  end subroutine mpiSumComplexV

!----------------------------------------------------------------------------
  subroutine mpiExitError(code)
    !! Exit on error with MPI communication

    implicit none
    
    integer, intent(in) :: code

    write(*, '( "*** MPI error ***")')
    write(*, '( "*** error code: ",I5, " ***")') code

    call MPI_ABORT(worldComm,code,ierr)
    
    stop

    return
  end subroutine mpiExitError

!----------------------------------------------------------------------------
  subroutine exitError(calledFrom, message, ierror)
    !! Output error message and abort if ierr > 0
    !!
    !! Can ensure that error will cause abort by
    !! passing abs(ierror)
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    integer, intent(in) :: ierror
      !! Error

    character(len=*), intent(in) :: calledFrom
      !! Place where this subroutine was called from
    character(len=*), intent(in) :: message
      !! Error message

    integer :: id
      !! ID of this process
    integer :: mpierr
      !! Error output from MPI

    character(len=6) :: cerr
      !! String version of error


    if ( ierror <= 0 ) return
      !! * Do nothing if the error is less than or equal to zero

    write( cerr, fmt = '(I6)' ) ierror
      !! * Write ierr to a string
    write(unit=*, fmt = '(/,1X,78("%"))' )
      !! * Output a dividing line
    write(unit=*, fmt = '(5X,"Error in ",A," (",A,"):")' ) trim(calledFrom), trim(adjustl(cerr))
      !! * Output where the error occurred and the error
    write(unit=*, fmt = '(5X,A)' ) TRIM(message)
      !! * Output the error message
    write(unit=*, fmt = '(1X,78("%"),/)' )
      !! * Output a dividing line

    write( *, '("     stopping ...")' )
  
    id = 0
  
    !> * For MPI, get the id of this process and abort
    call MPI_COMM_RANK( worldComm, id, mpierr )
    call MPI_ABORT( worldComm, mpierr, ierr )
    call MPI_FINALIZE( mpierr )

    stop 2

    return

  end subroutine exitError

end module errorsAndMPI
