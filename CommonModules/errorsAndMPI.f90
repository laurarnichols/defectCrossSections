module errorsAndMPI
  
  use constants, only: dp
  use mpi

  implicit none

  integer, parameter :: root = 0
    !! ID of the root node

  integer :: ierr
    !! Error handler
  integer :: interBgrpComm = 0
    !! Inter-band-group communicator
  integer :: intraBgrpComm = 0
    !! Intra-band-group communicator
  integer :: interPoolComm = 0
    !! Inter-pool communicator
  integer :: intraPoolComm = 0
    !! Intra-pool communicator
  integer :: indexInBgrp
    !! Process index within band group
  integer :: indexInPool
    !! Process index within pool
  integer :: myid
    !! ID of this process
  integer :: myBgrpId
    !! Band-group index for this process
  integer :: myPoolId
    !! Pool index for this process
  integer :: nBandGroups = 1
    !! Number of band groups for parallelization
  integer :: nPools = 1
    !! Number of pools for k-point parallelization
  integer :: nProcs
    !! Number of processes
  integer :: nProcPerBgrp
    !! Number of processes per band group
  integer :: nProcPerPool
    !! Number of processes per pool
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
  subroutine getCommandLineArguments()
    !! Get the command line arguments. This currently
    !! only processes the number of pools
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Output variables:
    !integer, intent(out) :: nBandGroups
      ! Number of band groups for parallelization
    !integer, intent(out) :: nPools
      ! Number of pools for k-point parallelization


    ! Local variables:
    integer :: narg = 0
      !! Arguments processed
    integer :: nargs
      !! Total number of command line arguments
    integer :: nBandGroups_ = 1
      !! Number of band groups for parallelization
    integer :: nPools_ = 1
      !! Number of k point pools for parallelization

    character(len=256) :: arg = ' '
      !! Command line argument
    character(len=256) :: command_line = ' '
      !! Command line arguments that were not processed


    nargs = command_argument_count()
      !! * Get the number of arguments input at command line

    call MPI_BCAST(nargs, 1, MPI_INTEGER, root, worldComm, ierr)

    if(ionode) then

      call get_command_argument(narg, arg)
        !! Ignore executable
      narg = narg + 1

      do while (narg <= nargs)
        call get_command_argument(narg, arg)
          !! * Get the flag
          !! @note
          !!  This program only currently processes the number of pools,
          !!  represented by `-nk`/`-nPools`/`-nPoolss`. All other flags 
          !!  will be ignored.
          !! @endnote

        narg = narg + 1

        !> * Process the flag and store the following value
        select case (trim(arg))
          case('-nk', '-nPools', '-nPoolss') 
            call get_command_argument(narg, arg)
            read(arg, *) nPools_
            narg = narg + 1
          case('-nb', '-nband', '-nbgrp', '-nband_group') 
            call get_command_argument(narg, arg)
            read(arg, *) nBandGroups_
            narg = narg + 1
          case default
            command_line = trim(command_line) // ' ' // trim(arg)
        end select
      enddo

      !> Write out unprocessed command line arguments, if there are any
      if(len_trim(command_line) /= 0) then
        write(*,*) 'Unprocessed command line arguments: ' // trim(command_line)
      endif

    endif

    call MPI_BCAST(nPools_, 1, MPI_INTEGER, root, worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8005)

    call MPI_BCAST(nBandGroups_, 1, MPI_INTEGER, root, worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8006)

    nPools = nPools_
    nBandGroups = nBandGroups_

    return
  end subroutine getCommandLineArguments

!----------------------------------------------------------------------------
  subroutine setUpPools()
    !! Split up processors between pools and generate MPI
    !! communicators for pools
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    !integer, intent(in) :: myid
      ! ID of this process
    !integer, intent(in) :: nBandGroups
      ! Number of band groups for parallelization
    !integer, intent(in) :: nPools
      ! Number of pools for k-point parallelization
    !integer, intent(in) :: nProcs
      ! Number of processes


    ! Output variables:
    !integer, intent(out) :: intraBgrpComm = 0
      ! Intra-band-group communicator
    !integer, intent(out) :: intraPoolComm = 0
      ! Intra-pool communicator
    !integer, intent(out) :: indexInBgrp
      ! Process index within band group
    !integer, intent(out) :: indexInPool
      ! Process index within pool
    !integer, intent(out) :: myBgrpId
      ! Band-group index for this process
    !integer, intent(out) :: myPoolId
      ! Pool index for this process
    !integer, intent(out) :: nProcPerBgrp
      ! Number of processes per band group
    !integer, intent(out) :: nProcPerPool
      ! Number of processes per pool


    if(nPools < 1 .or. nPools > nProcs) call exitError('mpiInitialization', &
      'invalid number of pools, out of range', 1)
      !! * Verify that the number of pools is between 1 and the number of processes

    if(mod(nProcs, nPools) /= 0) call exitError('mpiInitialization', &
      'invalid number of pools, mod(nProcs,nPools) /=0 ', 1)
      !! * Verify that the number of processes is evenly divisible by the number of pools

    nProcPerPool = nProcs / nPools
      !! * Calculate how many processes there are per pool

    myPoolId = myid / nProcPerPool
      !! * Get the pool index for this process

    indexInPool = mod(myid, nProcPerPool)
      !! * Get the index of the process within the pool

    call MPI_BARRIER(worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8007)

    call MPI_COMM_SPLIT(worldComm, myPoolId, myid, intraPoolComm, ierr)
    if(ierr /= 0) call mpiExitError(8008)
      !! * Create intra-pool communicator

    call MPI_COMM_SPLIT(worldComm, indexInPool, myid, interPoolComm, ierr)
    if(ierr /= 0) call mpiExitError(8010)
      !! * Create inter-pool communicator




    if(nBandGroups < 1 .or. nBandGroups > nProcPerPool) call exitError('mpiInitialization', &
      'invalid number of band groups, out of range', 1)
      !! * Verify that the number of band groups is between 1 and the number of processes per pool

    if(mod(nProcPerPool, nBandGroups) /= 0) call exitError('mpiInitialization', &
      'invalid number of band groups, mod(nProcPerPool,nBandGroups) /=0 ', 1)
      !! * Verify that the number of processes per pool is evenly divisible by the number of band groups

    nProcPerBgrp = nProcPerPool / nBandGroups
      !! * Calculate how many processes there are per band group

    myBgrpId = indexInPool / nProcPerBgrp
      !! * Get the band-group index for this process

    indexInBgrp = mod(indexInPool, nProcPerBgrp)
      !! * Get the index of the process within the band group

    call MPI_BARRIER(worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8009)

    call MPI_COMM_SPLIT(intraPoolComm, myBgrpId, indexInPool, intraBgrpComm, ierr)
    if(ierr /= 0) call mpiExitError(8010)
      !! * Create intra-band group communicator

    call MPI_BARRIER(worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8011)

    call MPI_COMM_SPLIT(intraPoolComm, indexInBgrp, indexInPool, interBgrpComm, ierr)
    if(ierr /= 0) call mpiExitError(8012)
      !! * Create inter-band-group communicator

    return
  end subroutine setUpPools

!----------------------------------------------------------------------------
  function checkDirInitialization(variableName, variableValue, fileInDir) result(abortExecution)

    implicit none

    ! Input variables:
    character(*), intent(in) :: fileInDir
      !! File to look for in directory 
      !! (for compiler compatibility)
    character(*), intent(in) :: variableName
      !! Name of variable to be tested
    character(*), intent(in) :: variableValue
      !! Value of variable to be tested

    ! Output variables:
    logical :: abortExecution
      !! If execution should stop
    logical :: fileExists
      !! If file exists


    abortExecution = .false.

    if(trim(variableValue) == '' ) then

      write(*,*)
      write(*,'(" Variable ",a," is not defined!")') trim(variableName)
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    else
      
      inquire(file=trim(variableValue)//'/'//trim(fileInDir), exist=fileExists)
      
      if(fileExists .eqv. .false.) then

        write(*,'(" File ", a, " , does not exist!")') trim(variableValue)//'/'//trim(fileInDir)
        write(*,'(" This variable is mandatory and thus the program will not be executed!")')

        abortExecution = .true.

      endif
    endif


    write(*, '(a," = ''", a, "''")') trim(variableName), trim(variableValue)

  end function checkDirInitialization

!----------------------------------------------------------------------------
  function checkFileInitialization(variableName, variableValue) result(abortExecution)

    implicit none

    ! Input variables:
    character(*), intent(in) :: variableName
      !! Name of variable to be tested
    character(*), intent(in) :: variableValue
      !! Value of variable to be tested

    ! Output variables:
    logical :: abortExecution
      !! If execution should stop
    logical :: fileExists
      !! If file exists


    abortExecution = .false.

    if(trim(variableValue) == '' ) then

      write(*,*)
      write(*,'(" Variable ",a," is not defined!")') trim(variableName)
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    else
      
      inquire(file=trim(variableValue), exist=fileExists)
      
      if(fileExists .eqv. .false.) then

        write(*,'(" File ", a, " , does not exist!")') trim(variableValue)
        write(*,'(" This variable is mandatory and thus the program will not be executed!")')

        abortExecution = .true.

      endif
    endif


    write(*, '(a," = ''", a, "''")') trim(variableName), trim(variableValue)

  end function checkFileInitialization

!----------------------------------------------------------------------------
  function checkStringInitialization(variableName, variableValue) result(abortExecution)

    implicit none

    ! Input variables:
    character(*), intent(in) :: variableName
      !! Name of variable to be tested
    character(*), intent(in) :: variableValue
      !! Value of variable to be tested

    ! Output variables:
    logical :: abortExecution
      !! If execution should stop


    abortExecution = .false.

    if(trim(variableValue) == '' ) then

      write(*,*)
      write(*,'(" Variable ",a," is not defined!")') trim(variableName)
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    endif


    write(*, '(a," = ''", a, "''")') trim(variableName), trim(variableValue)

  end function checkStringInitialization

!----------------------------------------------------------------------------
  function checkIntInitialization(variableName, variableValue, minValue, maxValue) result(abortExecution)

    implicit none

    ! Input variables:
    integer, intent(in) :: minValue, maxValue
      !! Max and min value of variable to test
    integer, intent(in) :: variableValue
      !! Value of variable to be tested

    character(*), intent(in) :: variableName
      !! Name of variable to be tested

    ! Output variables:
    logical :: abortExecution
      !! If execution should stop


    abortExecution = .false.

    if(variableValue < minValue .or. variableValue > maxValue) then

      write(*,*)
      write(*,'(" Variable ",a," is not defined or is out of range!")') trim(variableName)
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    endif


    write(*, '(a," = ", i7)') trim(variableName), variableValue

  end function checkIntInitialization

!----------------------------------------------------------------------------
  function checkDoubleInitialization(variableName, variableValue, minValue, maxValue) result(abortExecution)

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: minValue, maxValue
      !! Max and min value of variable to test
    real(kind=dp), intent(in) :: variableValue
      !! Value of variable to be tested

    character(*), intent(in) :: variableName
      !! Name of variable to be tested

    ! Output variables:
    logical :: abortExecution
      !! If execution should stop


    abortExecution = .false.

    if(variableValue < minValue .or. variableValue > maxValue) then

      write(*,*)
      write(*,'(" Variable ",a," is not defined or is out of range!")') trim(variableName)
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    endif


    write(*, '(a," = ", ES10.3E1)') trim(variableName), variableValue

  end function checkDoubleInitialization

!----------------------------------------------------------------------------
  subroutine distributeItemsInSubgroups(mySubgroupId, nItemsToDistribute, nProcPerLargerGroup, nProcPerSubgroup, nSubgroups, &
        iItemStart_subgroup, iItemEnd_subgroup, nItemsPerSubgroup)
    !! Distribute items across a subgroup
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: mySubgroupId
      !! Process ID in subgroup
    integer, intent(in) :: nItemsToDistribute
      !! Total number of items to distribute
    integer, intent(in) :: nProcPerLargerGroup
      !! Number of processes in group next
      !! larger than this subgroup
    integer, intent(in) :: nProcPerSubgroup
      !! Number of processes per subgroup
    integer, intent(in) :: nSubgroups
      !! Number of subgroups


    ! Output variables:
    integer, intent(out) :: iItemStart_subgroup
      !! Starting index for items in single subgroup
    integer, intent(out) :: iItemEnd_subgroup
      !! Ending index for items in single subgroup
    integer, intent(out) :: nItemsPerSubgroup
      !! Number of items in each subgroup


    ! Local variables:
    integer :: nr
      !! Number of items left over after evenly divided across subgroups


    if(nItemsToDistribute > 0) then

      if(nProcPerSubgroup > nProcPerLargerGroup .or. mod(nProcPerLargerGroup, nProcPerSubgroup) /= 0) &
        call exitError('distributeItemsInSubgroups','nProcPerSubgroup', 1)

      nItemsPerSubgroup = nItemsToDistribute / nSubgroups
        !!  * Calculate items per subgroup

      nr = nItemsToDistribute - nItemsPerSubgroup * nSubgroups
        !! * Calculate the remainder 

      IF( mySubgroupId < nr ) nItemsPerSubgroup = nItemsPerSubgroup + 1
        !! * Assign the remainder to the first `nr` subgroups

      !>  * Calculate the index of the first item in this subgroup
      iItemStart_subgroup = nItemsPerSubgroup * mySubgroupId + 1
      IF( mySubgroupId >= nr ) iItemStart_subgroup = iItemStart_subgroup + nr

      iItemEnd_subgroup = iItemStart_subgroup + nItemsPerSubgroup - 1
        !!  * Calculate the index of the last k-point in this pool

    endif

    return
  end subroutine distributeItemsInSubgroups

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
