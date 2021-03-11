module wfcExportVASPMod
  
  use constants, only: dp, iostd, angToBohr, eVToRy, ryToHartree, pi
  use mpi

  implicit none

  ! Parameters:
  integer, parameter :: root = 0
    !! ID of the root node
  integer, parameter :: rootInPool = 0
    !! Index of the root process within each pool
  integer, parameter :: mainOutFileUnit = 50
    !! Main output file unit
  integer, parameter :: potcarUnit = 86
    !! POTCAR unit for I/O
  integer, parameter :: wavecarUnit = 86
    !! WAVECAR unit for I/O

  real(kind = dp), parameter :: twoPiSquared = (2.0_dp*pi)**2
    !! This is used in place of \(2\pi/a\) which assumes that \(a=1\)


  ! Global variables not passed as arguments:
  integer :: ierr
    !! Error returned by MPI
  integer :: ikEnd
    !! Ending index for k-points in single pool 
  integer :: ikStart
    !! Starting index for k-points in single pool 
  integer :: ios
    !! Error for input/output
  integer :: indexInPool
    !! Process index within pool
  integer :: interPoolComm = 0
    !! Inter-pool communicator
  integer :: intraPoolComm = 0
    !! Intra-pool communicator
  integer :: myid
    !! ID of this process
  integer :: myPoolId
    !! Pool index for this process
  integer :: nkPerPool
    !! Number of k-points in each pool
  integer :: nPools = 1
    !! Number of pools for k-point parallelization
  integer :: nProcs
    !! Number of processes
  integer :: nProcPerPool
    !! Number of processes per pool
  integer :: worldComm
    !! World communicator

  logical :: ionode
    !! If this node is the root node


  ! Variables that should be passed as arguments:
  real(kind=dp) :: realSpaceLatticeVectors(3,3)
    !! Real space lattice vectors
  real(kind=dp) :: recipSpaceLatticeVectors(3,3)
    !! Reciprocal lattice vectors
  real(kind=dp) :: wfcECut
    !! Plane wave energy cutoff in Ry
  real(kind=dp) :: eFermi
    !! Fermi energy
  real(kind=dp), allocatable :: gVecInCart(:,:)
    !! G-vectors in Cartesian coordinates
  real(kind=dp), allocatable :: bandOccupation(:,:)
    !! Occupation of band
  real(kind=dp) :: omega
    !! Volume of unit cell
  real(kind=dp), allocatable :: atomPositions(:,:)
    !! Atom positions
  real(kind=dp) :: tStart
    !! Start time
  real(kind=dp) :: wfcVecCut
    !! Energy cutoff converted to vector cutoff
  real(kind=dp), allocatable :: kPosition(:,:)
    !! Position of k-points in reciprocal space
  real(kind=dp), allocatable :: kWeight(:)
    !! Weight of k-points

  complex*16, allocatable :: eigenE(:,:,:)
    !! Band eigenvalues
  
  integer, allocatable :: gIndexLocalToGlobal(:)
    !! Converts local index `ig` to global index
  integer, allocatable :: gKIndexLocalToGlobal(:,:)
    !! Local to global indices for \(G+k\) vectors 
    !! ordered by magnitude at a given k-point
  integer, allocatable :: gToGkIndexMap(:,:)
    !! Index map from \(G\) to \(G+k\);
    !! indexed up to `nGVecsLocal` which
    !! is greater than `maxNumPWsPool` and
    !! stored for each k-point
  integer, allocatable :: gKIndexGlobal(:,:)
    !! Indices of \(G+k\) vectors for each k-point
    !! and all processors
  integer, allocatable :: iType(:)
    !! Atom type index
  integer, allocatable :: gVecMillerIndicesGlobal(:,:)
    !! Integer coefficients for G-vectors on all processors
  integer :: nb1max, nb2max, nb3max
    !! Not sure what this is??
  integer :: nAtoms
    !! Number of atoms
  integer :: nBands
    !! Total number of bands
  integer, allocatable :: nGkLessECutGlobal(:)
    !! Global number of \(G+k\) vectors with energy
    !! less than `wfcECut` for each k-point
  integer, allocatable :: nGkLessECutLocal(:)
    !! Number of \(G+k\) vectors with energy
    !! less than `wfcECut` for each
    !! k-point, on this processor
  integer :: maxGkNum
    !! Maximum number of \(G+k\) combinations
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors on this processor
  integer :: nKPoints
    !! Total number of k-points
  integer, allocatable :: nAtomsEachType(:)
    !! Number of atoms of each type
  integer, allocatable :: nPWs1kGlobal(:)
    !! Input number of plane waves for a single k-point for all processors
  integer :: maxGIndexGlobal
    !! Maximum G-vector index among all \(G+k\)
    !! and processors
  integer :: maxNumPWsGlobal
    !! Max number of \(G+k\) vectors with energy
    !! less than `wfcECut` among all k-points
  integer :: maxNumPWsPool
    !! Maximum number of \(G+k\) vectors
    !! across all k-points for just this
    !! ppool
  integer :: nAtomTypes
    !! Number of types of atoms
  integer :: nSpins
    !! Number of spins
  
  character(len=256) :: exportDir
    !! Directory to be used for export
  character(len=256) :: mainOutputFile
    !! Main output file
  character(len=256) :: VASPDir
    !! Directory with VASP files

  type pseudo
    integer :: angmom(16) = 0
      !! Angular momentum of projectors
    integer :: iRAugMax
      !! Max index of augmentation sphere
    integer :: lmmax
      !! Total number of nlm channels
    integer :: nChannels
      !! Number of l channels;
      !! also number of projectors
    integer :: nmax
      !! Number of radial grid points

    real(kind=dp), allocatable :: dRadGrid(:)
      !! Derivative of radial grid
    real(kind=dp) :: psRMax
      !! Max r for non-local contribution
    real(kind=dp), allocatable :: radGrid(:)
      !! Radial grid points
    real(kind=dp) :: rAugMax
      !! Maximum radius of augmentation sphere
    real(kind=dp) :: recipProj(16,100)
      !! Reciprocal-space projectors
    real(kind=dp) :: realProj(16,100)
      !! Real-space projectors
    real(kind=dp), allocatable :: wae(:,:)
      !! AE wavefunction
    real(kind=dp), allocatable :: wps(:,:)
      !! PS wavefunction

    character(len=2) :: element
  end type pseudo

  type (pseudo), allocatable :: ps(:)

  namelist /inputParams/ VASPDir, exportDir


  contains

!----------------------------------------------------------------------------
  subroutine mpiInitialization()
    !! Generate MPI processes and communicators 
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Output variables:
    !logical, intent(out) :: ionode
      ! If this node is the root node
    !integer, intent(out) :: interPoolComm = 0
      ! Inter-pool communicator
    !integer, intent(out) :: intraPoolComm = 0
      ! Intra-pool communicator
    !integer, intent(out) :: indexInPool
      ! Process index within pool
    !integer, intent(out) :: myid
      ! ID of this process
    !integer, intent(out) :: myPoolId
      ! Pool index for this process
    !integer, intent(out) :: nPools
      ! Number of pools for k-point parallelization
    !integer, intent(out) :: nProcs
      ! Number of processes
    !integer, intent(out) :: nProcPerPool
      ! Number of processes per pool
    !integer, intent(out) :: worldComm
      ! World communicator


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

    call getCommandLineArguments()
      !! * Get the number of pools from the command line

    call setUpPools()
      !! * Split up processors between pools and generate MPI
      !!   communicators for pools

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
    !integer, intent(out) :: nPools
      ! Number of pools for k-point parallelization


    ! Local variables:
    integer :: narg = 0
      !! Arguments processed
    integer :: nargs
      !! Total number of command line arguments
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
          case default
            command_line = trim(command_line) // ' ' // trim(arg)
        end select
      enddo

      write(iostd,*) 'Unprocessed command line arguments: ' // trim(command_line)
    endif

    call MPI_BCAST(nPools_, 1, MPI_INTEGER, root, worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8005)

    nPools = nPools_

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
    !integer, intent(in) :: nPools
      ! Number of pools for k-point parallelization
    !integer, intent(in) :: nProcs
      ! Number of processes


    ! Output variables:
    !integer, intent(out) :: interPoolComm = 0
      ! Inter-pool communicator
    !integer, intent(out) :: intraPoolComm = 0
      ! Intra-pool communicator
    !integer, intent(out) :: indexInPool
      ! Process index within pool
    !integer, intent(out) :: myPoolId
      ! Pool index for this process
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

    call MPI_BARRIER(worldComm, ierr)
    if(ierr /= 0) call mpiExitError(8009)

    call MPI_COMM_SPLIT(worldComm, indexInPool, myid, interPoolComm, ierr)
    if(ierr /= 0) call mpiExitError(8010)
      !! * Create inter-pool communicator

    return
  end subroutine setUpPools

!----------------------------------------------------------------------------
  subroutine initialize(exportDir, VASPDir)
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    ! Input variables:
    !integer, intent(in) :: nPools
      ! Number of pools for k-point parallelization
    !integer, intent(in) :: nProcs
      ! Number of processes


    ! Output variables:
    character(len=256), intent(out) :: exportDir
      !! Directory to be used for export
    character(len=256), intent(out) :: VASPDir
      !! Directory with VASP files


    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time


    VASPDir = './'
    exportDir = './Export'

    call cpu_time(tStart)

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(iostd, '(/5X,"VASP wavefunction export program starts on ",A9," at ",A9)') &
             cdate, ctime

      write(iostd, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs

      if(nPools > 1) write(iostd, '(5X,"K-points division:     nPools     = ",I7)') nPools

    else

      open(unit = iostd, file='/dev/null', status='unknown')
        ! Make the iostd unit point to null for non-root processors
        ! to avoid tons of duplicate output

    endif

  end subroutine initialize

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
  subroutine mpiExitError(code)
    !! Exit on error with MPI communication

    implicit none
    
    integer, intent(in) :: code

    write( iostd, '( "*** MPI error ***")' )
    write( iostd, '( "*** error code: ",I5, " ***")' ) code

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


    if ( ierr <= 0 ) return
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
  
    call flush( iostd )
  
    id = 0
  
    !> * For MPI, get the id of this process and abort
    call MPI_COMM_RANK( worldComm, id, mpierr )
    call MPI_ABORT( worldComm, mpierr, ierr )
    call MPI_FINALIZE( mpierr )

    stop 2

    return

  end subroutine exitError

!----------------------------------------------------------------------------
  subroutine readWAVECAR(VASPDir, realSpaceLatticeVectors, recipSpaceLatticeVectors, wfcECut, bandOccupation, omega, wfcVecCut, &
        kPosition, nb1max, nb2max, nb3max, nBands, maxGkNum, nKPoints, nPWs1kGlobal, nSpins, &
        eigenE)
    !! Read cell and wavefunction data from the WAVECAR file
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    
    ! Output variables:
    real(kind=dp), intent(out) :: realSpaceLatticeVectors(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: recipSpaceLatticeVectors(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(out) :: wfcECut
      !! Plane wave energy cutoff in Ry
    real(kind=dp), allocatable, intent(out) :: bandOccupation(:,:)
      !! Occupation of band
    real(kind=dp), intent(out) :: omega
      !! Volume of unit cell
    real(kind=dp), intent(out) :: wfcVecCut
      !! Energy cutoff converted to vector cutoff
    real(kind=dp), allocatable, intent(out) :: kPosition(:,:)
      !! Position of k-points in reciprocal space

    integer, intent(out) :: nb1max, nb2max, nb3max
      !! Not sure what this is??
    integer, intent(out) :: nBands
      !! Total number of bands
    integer, intent(out) :: maxGkNum
      !! Maximum number of \(G+k\) combinations
    integer, intent(out) :: nKPoints
      !! Total number of k-points
    integer, allocatable, intent(out) :: nPWs1kGlobal(:)
      !! Input number of plane waves for a single k-point 
      !! for all processors
    integer, intent(out) :: nSpins
      !! Number of spins

    complex*16, allocatable, intent(out) :: eigenE(:,:,:)
      !! Band eigenvalues


    ! Local variables:
    real(kind=dp) :: c = 0.26246582250210965422
      !! \(2m/\hbar^2\) converted from J\(^{-1}\)m\(^{-2}\)
      !! to eV\(^{-1}\)A\(^{-2}\)
    real(kind=dp) :: nRecords_real, nspin_real, prec_real, nkstot_real 
      !! Real version of integers for reading from file
    real(kind=dp) :: nbnd_real
      !! Real version of integers for reading from file

    integer :: j
      !! Index used for reading lattice vectors
    integer :: nRecords
      !! Number of records in WAVECAR file
    integer :: prec
      !! Precision of plane wave coefficients

    character(len=256) :: fileName
      !! Full WAVECAR file name including path


    if(ionode) then
      fileName = trim(VASPDir)//'/WAVECAR'

      nRecords = 24
        ! Set a starting value for the number of records

      open(unit=wavecarUnit, file=fileName, access='direct', recl=nRecords, iostat=ierr, status='old')
      if (ierr .ne. 0) write(iostd,*) 'open error - iostat =', ierr
        !! * If root node, open the `WAVECAR` file

      read(unit=wavecarUnit,rec=1) nRecords_real, nspin_real, prec_real
        !! @note Must read in as real first then convert to integer @endnote

      close(unit=wavecarUnit)

      nRecords = nint(nRecords_real)
      nSpins = nint(nspin_real)
      prec = nint(prec_real)
        ! Convert input variables to integers

      !if(prec .eq. 45210) call exitError('readWAVECAR', 'WAVECAR_double requires complex*16', 1)

      open(unit=wavecarUnit, file=fileName, access='direct', recl=nRecords, iostat=ierr, status='old')
      if (ierr .ne. 0) write(iostd,*) 'open error - iostat =', ierr

      read(unit=wavecarUnit,rec=2) nkstot_real, nbnd_real, wfcECut,(realSpaceLatticeVectors(j,1),j=1,3),&
          (realSpaceLatticeVectors(j,2),j=1,3), (realSpaceLatticeVectors(j,3),j=1,3)
        !! * Read total number of k-points, plane wave cutoff energy, and real
        !!   space lattice vectors

      wfcECut = wfcECut*eVToRy
        !! * Convert energy from VASP from eV to Rydberg to match QE/TME expectation

      wfcVecCut = sqrt(wfcECut/evToRy*c)/angToBohr
        !! * Calculate vector cutoff from energy cutoff

      realSpaceLatticeVectors = realSpaceLatticeVectors*angToBohr

      nKPoints = nint(nkstot_real)
      nBands = nint(nbnd_real)
        ! Convert input variables to integers

      call calculateOmega(realSpaceLatticeVectors, omega)
        !! * Calculate the cell volume as \(a_1\cdot a_2\times a_3\)

      call getReciprocalVectors(realSpaceLatticeVectors, omega, recipSpaceLatticeVectors)
        !! * Calculate the reciprocal lattice vectors from the real-space
        !!   lattice vectors and the cell volume

      call estimateMaxNumPlanewaves(recipSpaceLatticeVectors, wfcECut, nb1max, nb2max, nb3max, maxGkNum)
        !! * Get the maximum number of plane waves

      !> * Write out total number of k-points, number of bands, 
      !>   the energy cutoff, the real-space-lattice vectors,
      !>   the cell volume, and the reciprocal lattice vectors
      write(iostd,*) 'no. k points =', nKPoints
      write(iostd,*) 'no. bands =', nBands
      write(iostd,*) 'max. energy (eV) =', sngl(wfcECut/eVToRy)
        !! @note 
        !!  The energy cutoff is currently output to the `iostd` file
        !!  in eV to compare with output from WaveTrans.
        !! @endnote
      write(iostd,*) 'real space lattice vectors:'
      write(iostd,*) 'a1 =', (sngl(realSpaceLatticeVectors(j,1)),j=1,3)
      write(iostd,*) 'a2 =', (sngl(realSpaceLatticeVectors(j,2)),j=1,3)
      write(iostd,*) 'a3 =', (sngl(realSpaceLatticeVectors(j,3)),j=1,3)
      write(iostd,*) 
      write(iostd,*) 'volume unit cell =', sngl(omega)
      write(iostd,*) 
      write(iostd,*) 'reciprocal lattice vectors:'
      write(iostd,*) 'b1 =', (sngl(recipSpaceLatticeVectors(j,1)),j=1,3)
      write(iostd,*) 'b2 =', (sngl(recipSpaceLatticeVectors(j,2)),j=1,3)
      write(iostd,*) 'b3 =', (sngl(recipSpaceLatticeVectors(j,3)),j=1,3)
      write(iostd,*) 
        !! @note
        !!  I made an intentional choice to stick with the unscaled lattice
        !!  vectors until I see if it will be convenient to scale them down.
        !!  QE uses the `alat` and `tpiba` scaling quite a bit though, so I
        !!  will have to be careful with the scaling/units.
        !! @endnote

      write(mainOutFileUnit, '("# Cell volume (a.u.)^3. Format: ''(ES24.15E3)''")')
      write(mainOutFileUnit, '(ES24.15E3)' ) omega
      flush(mainOutFileUnit)

    endif

    call MPI_BCAST(nb1max, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nb2max, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nb3max, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(maxGkNum, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nSpins, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nBands, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(wfcECut, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(wfcVecCut, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(realSpaceLatticeVectors, size(realSpaceLatticeVectors), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(recipSpaceLatticeVectors, size(recipSpaceLatticeVectors), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if (ionode) then
      write(iostd,*) 
      write(iostd,*) '******'
      write(iostd,*) 'Starting to read wavefunction'
    endif

    call readWavefunction(nBands, maxGkNum, nKPoints, nSpins, bandOccupation, kPosition, nPWs1kGlobal, eigenE)
      !! * Get the position of each k-point in reciprocal space 
      !!   and the number of \(G+k) vectors below the cutoff 
      !!   energy for each k-point

    if(ionode) close(wavecarUnit)

    if (ionode) then
      write(iostd,*) 'Finished reading wavefunction'
      write(iostd,*) '******'
      write(iostd,*) 
    endif

    return
  end subroutine readWAVECAR

!----------------------------------------------------------------------------
  subroutine calculateOmega(realSpaceLatticeVectors, omega)
    !! Calculate the cell volume as \(a_1\cdot a_2\times a_3\)

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: realSpaceLatticeVectors(3,3)
      !! Real space lattice vectors


    ! Output variables:
    real(kind=dp), intent(out) :: omega
      !! Volume of unit cell


    ! Local variables:
    real(kind=dp) :: vtmp(3)
      !! \(a_2\times a_3\)


    call vcross(realSpaceLatticeVectors(:,2), realSpaceLatticeVectors(:,3), vtmp)

    omega = sum(realSpaceLatticeVectors(:,1)*vtmp(:))

    return
  end subroutine calculateOmega

!----------------------------------------------------------------------------
  subroutine getReciprocalVectors(realSpaceLatticeVectors, omega, recipSpaceLatticeVectors)
    !! Calculate the reciprocal lattice vectors from the real-space
    !! lattice vectors and the cell volume

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: realSpaceLatticeVectors(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: omega
      !! Volume of unit cell


    ! Output variables:
    real(kind=dp), intent(out) :: recipSpaceLatticeVectors(3,3)
      !! Reciprocal lattice vectors


    ! Local variables:
    integer :: i
      !! Loop index
    

    call vcross(2.0d0*pi*realSpaceLatticeVectors(:,2)/omega, realSpaceLatticeVectors(:,3), recipSpaceLatticeVectors(:,1))
      ! \(b_1 = 2\pi/\Omega a_2\times a_3\)
    call vcross(2.0d0*pi*realSpaceLatticeVectors(:,3)/omega, realSpaceLatticeVectors(:,1), recipSpaceLatticeVectors(:,2))
      ! \(b_2 = 2\pi/\Omega a_3\times a_1\)
    call vcross(2.0d0*pi*realSpaceLatticeVectors(:,1)/omega, realSpaceLatticeVectors(:,2), recipSpaceLatticeVectors(:,3))
      ! \(b_3 = 2\pi/\Omega a_1\times a_2\)


    return
  end subroutine getReciprocalVectors

!----------------------------------------------------------------------------
  subroutine vcross(vec1, vec2, crossProd)
    !! Calculate the cross product `crossProd` of 
    !! two vectors `vec1` and `vec2`

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: vec1(3), vec2(3)
      !! Input vectors


    ! Output variables:
    real(kind=dp), intent(out) :: crossProd(3)
      !! Cross product of input vectors


    crossProd(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    crossProd(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    crossProd(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

    return
  end subroutine vcross

!----------------------------------------------------------------------------
  subroutine estimateMaxNumPlanewaves(recipSpaceLatticeVectors, wfcECut, nb1max, nb2max, nb3max, maxGkNum)
    !! Get the maximum number of plane waves. I'm not sure how 
    !! this is done completely. It seems to be just basic vector
    !! stuff, but I haven't been able to make sense of it.
    !! 
    !! @todo Figure out how `estimateMaxNumPlanewaves` works @endtodo
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: recipSpaceLatticeVectors(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: wfcECut
      !! Plane wave energy cutoff in Ry


    ! Output variables:
    integer, intent(out) :: nb1max, nb2max, nb3max
      !! Not sure what this is??
    integer, intent(out) :: maxGkNum
      !! Maximum number of \(G+k\) combinations


    ! Local variables:
    real(kind=dp) :: b1mag, b2mag, b3mag
      !! Reciprocal vector magnitudes
    real(kind=dp) :: c = 0.26246582250210965422
      !! \(2m/\hbar^2\) converted from J\(^{-1}\)m\(^{-2}\)
      !! to eV\(^{-1}\)A\(^{-2}\)
    real(kind=dp) :: phi12, phi13, phi23
      !! Angle between vectors
    real(kind=dp) :: sinphi123
      !! \(\sin\phi_{123}\)
    real(kind=dp) :: vmag
      !! Magnitude of temporary vector
    real(kind=dp) :: vtmp(3)
      !! Temporary vector for calculating angles

    integer :: nb1maxA, nb2maxA, nb3maxA
      !! Not sure what this is??
    integer :: nb1maxB, nb2maxB, nb3maxB
      !! Not sure what this is??
    integer :: nb1maxC, nb2maxC, nb3maxC
      !! Not sure what this is??
    integer :: npmaxA, npmaxB, npmaxC
      !! Not sure what this is??


    b1mag = sqrt(sum(recipSpaceLatticeVectors(:,1)**2))
    b2mag = sqrt(sum(recipSpaceLatticeVectors(:,2)**2))
    b3mag = sqrt(sum(recipSpaceLatticeVectors(:,3)**2))

    write(iostd,*) 'reciprocal lattice vector magnitudes:'
    write(iostd,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)
      !! * Calculate and output reciprocal vector magnitudes


    phi12 = acos(sum(recipSpaceLatticeVectors(:,1)*recipSpaceLatticeVectors(:,2))/(b1mag*b2mag))
      !! * Calculate angle between \(b_1\) and \(b_2\)

    call vcross(recipSpaceLatticeVectors(:,1), recipSpaceLatticeVectors(:,2), vtmp)
    vmag = sqrt(sum(vtmp(:)**2))
    sinphi123 = sum(recipSpaceLatticeVectors(:,3)*vtmp(:))/(vmag*b3mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxA = (dsqrt(wfcECut/eVToRy*c)/(b1mag*abs(sin(phi12)))) + 1
    nb2maxA = (dsqrt(wfcECut/eVToRy*c)/(b2mag*abs(sin(phi12)))) + 1
    nb3maxA = (dsqrt(wfcECut/eVToRy*c)/(b3mag*abs(sinphi123))) + 1
    npmaxA = nint(4.0*pi*nb1maxA*nb2maxA*nb3maxA/3.0)
      !! * Get first set of max values


    phi13 = acos(sum(recipSpaceLatticeVectors(:,1)*recipSpaceLatticeVectors(:,3))/(b1mag*b3mag))
      !! * Calculate angle between \(b_1\) and \(b_3\)

    call vcross(recipSpaceLatticeVectors(:,1), recipSpaceLatticeVectors(:,3), vtmp)
    vmag = sqrt(sum(vtmp(:)**2))
    sinphi123 = sum(recipSpaceLatticeVectors(:,2)*vtmp(:))/(vmag*b2mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxB = (dsqrt(wfcECut/eVToRy*c)/(b1mag*abs(sin(phi13)))) + 1
    nb2maxB = (dsqrt(wfcECut/eVToRy*c)/(b2mag*abs(sinphi123))) + 1
    nb3maxB = (dsqrt(wfcECut/eVToRy*c)/(b3mag*abs(sin(phi13)))) + 1
    npmaxB = nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      !! * Get first set of max values


    phi23 = acos(sum(recipSpaceLatticeVectors(:,2)*recipSpaceLatticeVectors(:,3))/(b2mag*b3mag))
      !! * Calculate angle between \(b_2\) and \(b_3\)

    call vcross(recipSpaceLatticeVectors(:,2), recipSpaceLatticeVectors(:,3), vtmp)
    vmag = sqrt(sum(vtmp(:)**2))
    sinphi123 = sum(recipSpaceLatticeVectors(:,1)*vtmp(:))/(vmag*b1mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxC = (dsqrt(wfcECut/eVToRy*c)/(b1mag*abs(sinphi123))) + 1
    nb2maxC = (dsqrt(wfcECut/eVToRy*c)/(b2mag*abs(sin(phi23)))) + 1
    nb3maxC = (dsqrt(wfcECut/eVToRy*c)/(b3mag*abs(sin(phi23)))) + 1
    npmaxC = nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)
      !! * Get first set of max values


    nb1max = max0(nb1maxA,nb1maxB,nb1maxC)
    nb2max = max0(nb2maxA,nb2maxB,nb2maxC)
    nb3max = max0(nb3maxA,nb3maxB,nb3maxC)
    maxGkNum = min0(npmaxA,npmaxB,npmaxC)

    write(iostd,*) 'max. no. G values; 1,2,3 =', nb1max, nb2max, nb3max
    write(iostd,*) ' '

    write(iostd,*) 'estimated max. no. plane waves =', maxGkNum

    return
  end subroutine estimateMaxNumPlanewaves

!----------------------------------------------------------------------------
  subroutine readWavefunction(nBands, maxGkNum, nKPoints, nSpins, bandOccupation, kPosition, nPWs1kGlobal, eigenE)
    !! For each spin and k-point, read the number of
    !! \(G+k\) vectors below the energy cutoff, the
    !! position of the k-point in reciprocal space, 
    !! and the eigenvalue, occupation, and plane
    !! wave coefficients for each band
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: maxGkNum
      !! Maximum number of \(G+k\) combinations
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nSpins
      !! Number of spins


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: bandOccupation(:,:)
      !! Occupation of band
    real(kind=dp), allocatable, intent(out) :: kPosition(:,:)
      !! Position of k-points in reciprocal space

    integer, allocatable, intent(out) :: nPWs1kGlobal(:)
      !! Input number of plane waves for a single k-point 
      !! for all processors

    complex*16, allocatable, intent(out) :: eigenE(:,:,:)
      !! Band eigenvalues


    ! Local variables:
    real(kind=dp) :: nPWs1kGlobal_real
      !! Real version of integers for reading from file

    complex*8, allocatable :: coeff(:,:)
      !! Plane wave coefficients

    integer :: irec, isp, ik, i, iband, iplane
      !! Loop indices


    allocate(bandOccupation(nBands, nKPoints))
    allocate(kPosition(3,nKPoints))
    allocate(nPWs1kGlobal(nKPoints))
    allocate(eigenE(nSpins,nKPoints,nBands))

    if(ionode) then

      allocate(coeff(maxGkNum,nBands))
    
      irec=2

      do isp = 1, nSpins
        !! * For each spin:
        !!    * Go through each k-point
        !!       * Read in the number of \(G+k\) plane wave
        !!         vectors below the energy cutoff
        !!       * Read the position of the k-point in 
        !!         reciprocal space
        !!       * Read in the eigenvalue and occupation for
        !!         each band
        !!       * Read in the plane wave coefficients for
        !!         each band

        write(iostd,*) '  Reading spin ', isp

        do ik = 1, nKPoints

          write(iostd,*) '    Reading k-point ', ik

          irec = irec + 1
       
          read(unit=wavecarUnit,rec=irec) nPWs1kGlobal_real, (kPosition(i,ik),i=1,3), &
                 (eigenE(isp,ik,iband), bandOccupation(iband, ik), iband=1,nBands)
            ! Read in the number of \(G+k\) plane wave vectors below the energy
            ! cutoff, the position of the k-point in reciprocal space, and
            ! the eigenvalue and occupation for each band

          nPWs1kGlobal(ik) = nint(nPWs1kGlobal_real)

          do iband = 1, nBands

            irec = irec + 1

            read(unit=wavecarUnit,rec=irec) (coeff(iplane,iband), iplane=1,nPWs1kGlobal(ik))
              ! Read in the plane wave coefficients for each band

          enddo
        enddo
      enddo

      eigenE(:,:,:) = eigenE(:,:,:)*eVToRy

      deallocate(coeff)
        !! @note 
        !!  The plane wave coefficients are not currently used anywhere.
        !! @endnote

    endif

    call MPI_BCAST(kPosition, size(kPosition), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(bandOccupation, size(bandOccupation), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eigenE, size(eigenE), MPI_COMPLEX, root, worldComm, ierr)
    call MPI_BCAST(nPWs1kGlobal, size(nPWs1kGlobal), MPI_INTEGER, root, worldComm, ierr)

    return
  end subroutine readWavefunction

!----------------------------------------------------------------------------
  subroutine distributeKpointsInPools(nKPoints)
    !! Figure out how many k-points there should be per pool
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nKPoints
      !! Total number of k-points


    ! Output variables:
    !integer, intent(out) :: ikEnd
      ! Ending index for k-points in single pool 
    !integer, intent(out) :: ikStart
      ! Starting index for k-points in single pool 
    !integer, intent(out) :: nkPerPool
      ! Number of k-points in each pool


    ! Local variables:
    integer :: nkr
      !! Number of k-points left over after evenly divided across pools


    if( nKPoints > 0 ) then

      IF( ( nProcPerPool > nProcs ) .or. ( mod( nProcs, nProcPerPool ) /= 0 ) ) &
        CALL exitError( 'distributeKpointsInPools','nProcPerPool', 1 )

      nkPerPool = nKPoints / nPools
        !!  * Calculate k-points per pool

      nkr = nKPoints - nkPerPool * nPools 
        !! * Calculate the remainder `nkr`

      IF( myPoolId < nkr ) nkPerPool = nkPerPool + 1
        !! * Assign the remainder to the first `nkr` pools

      !>  * Calculate the index of the first k-point in this pool
      ikStart = nkPerPool * myPoolId + 1
      IF( myPoolId >= nkr ) ikStart = ikStart + nkr

      ikEnd = ikStart + nkPerPool - 1
        !!  * Calculate the index of the last k-point in this pool

    endif

    return
  end subroutine distributeKpointsInPools

!----------------------------------------------------------------------------
  subroutine calculateGvecs(nb1max, nb2max, nb3max, recipSpaceLatticeVectors, gVecInCart, gIndexLocalToGlobal, gVecMillerIndicesGlobal, &
      nGVecsGlobal, nGVecsLocal)
    !! Calculate Miller indices and G-vectors and split
    !! over processors
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nb1max, nb2max, nb3max
      !! Not sure what this is??

    real(kind=dp), intent(in) :: recipSpaceLatticeVectors(3,3)
      !! Reciprocal lattice vectors


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: gVecInCart(:,:)
      !! G-vectors in Cartesian coordinates

    integer, allocatable, intent(out) :: gIndexLocalToGlobal(:)
      !! Converts local index `ig` to global index
    integer, allocatable, intent(out) :: gVecMillerIndicesGlobal(:,:)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors
    integer, intent(out) :: nGVecsLocal
      !! Local number of G-vectors on this processor


    ! Local variables:
    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp), allocatable :: millSum(:)
      !! Sum of integer coefficients for G-vectors

    integer :: ig1, ig2, ig3, ig1p, ig2p, ig3p, ig, ix
      !! Loop indices
    integer, allocatable :: iMill(:)
      !! Indices of miller indices after sorting
    integer, allocatable :: gVecMillerIndicesGlobal_tmp(:,:)
      !! Integer coefficients for G-vectors on all processors
    integer, allocatable :: mill_local(:,:)
      !! Integer coefficients for G-vectors
    integer :: npmax
      !! Max number of plane waves


    npmax = (2*nb1max+1)*(2*nb2max+1)*(2*nb3max+1) 
    allocate(gVecMillerIndicesGlobal(3,npmax))
    allocate(millSum(npmax))

    if(ionode) then

      allocate(gVecMillerIndicesGlobal_tmp(3,npmax))

      write(iostd,*)
      write(iostd,*) "***************"
      write(iostd,*) "Calculating miller indices"

      nGVecsGlobal = 0
      gVecMillerIndicesGlobal_tmp = 0

      do ig3 = 0, 2*nb3max

        ig3p = ig3

        if (ig3 .gt. nb3max) ig3p = ig3 - 2*nb3max - 1

        do ig2 = 0, 2*nb2max

          ig2p = ig2

          if (ig2 .gt. nb2max) ig2p = ig2 - 2*nb2max - 1

          do ig1 = 0, 2*nb1max

            ig1p = ig1

            if (ig1 .gt. nb1max) ig1p = ig1 - 2*nb1max - 1

            nGVecsGlobal = nGVecsGlobal + 1

            gVecMillerIndicesGlobal_tmp(1,nGVecsGlobal) = ig1p
            gVecMillerIndicesGlobal_tmp(2,nGVecsGlobal) = ig2p
            gVecMillerIndicesGlobal_tmp(3,nGVecsGlobal) = ig3p
              !! * Calculate Miller indices

            millSum(nGVecsGlobal) = sqrt(real(ig1p**2 + ig2p**2 + ig3p**2))
              !! * Calculate the sum of the Miller indices
              !!   for sorting

          enddo
        enddo
      enddo

      if (nGVecsGlobal .ne. npmax) call exitError('calculateGvecs', & 
        '*** error - computed no. of G-vectors != estimated number of plane waves', 1)
        !! * Check that number of G-vectors are the same as the number of plane waves

      write(iostd,*) "Sorting miller indices"

      allocate(iMill(nGVecsGlobal))

      do ig = 1, nGVecsGlobal
        !! * Initialize the index array that will track elements
        !!   after sorting

        iMill(ig) = ig

      enddo

      call hpsort_eps(nGVecsGlobal, millSum, iMill, eps8)
        !! * Order vector `millSum` keeping initial position in `iMill`

      do ig = 1, nGVecsGlobal
        !! * Rearrange the miller indices to match order of `millSum`

        gVecMillerIndicesGlobal(:,ig) = gVecMillerIndicesGlobal_tmp(:,iMill(ig))

      enddo

      deallocate(iMill)
      deallocate(gVecMillerIndicesGlobal_tmp)

      write(iostd,*) "Done calculating and sorting miller indices"
      write(iostd,*) "***************"
      write(iostd,*)
    endif

    call MPI_BCAST(nGVecsGlobal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(gVecMillerIndicesGlobal, size(gVecMillerIndicesGlobal), MPI_INTEGER, root, worldComm, ierr)

    if (ionode) then
      write(iostd,*)
      write(iostd,*) "***************"
      write(iostd,*) "Distributing G-vecs over processors"
    endif

    call distributeGvecsOverProcessors(nGVecsGlobal, gVecMillerIndicesGlobal, gIndexLocalToGlobal, mill_local, nGVecsLocal)
      !! * Split up the G-vectors and Miller indices over processors 

    if (ionode) write(iostd,*) "Calculating G-vectors"

    allocate(gVecInCart(3,nGVecsLocal))

    do ig = 1, nGVecsLocal

      do ix = 1, 3
        !! * Calculate \(G = m_1b_1 + m_2b_2 + m_3b_3\)

        gVecInCart(ix,ig) = sum(mill_local(:,ig)*recipSpaceLatticeVectors(ix,:))

      enddo
      
    enddo

    if (ionode) then
      write(iostd,*) "***************"
      write(iostd,*)
    endif

    deallocate(mill_local)

    return
  end subroutine calculateGvecs

!----------------------------------------------------------------------------
  subroutine distributeGvecsOverProcessors(nGVecsGlobal, gVecMillerIndicesGlobal, gIndexLocalToGlobal, mill_local, nGVecsLocal)
    !! Figure out how many G-vectors there should be per processor
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors
      
    integer, intent(in) :: gVecMillerIndicesGlobal(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors

    
    ! Output variables:
    integer, allocatable, intent(out) :: gIndexLocalToGlobal(:)
      ! Converts local index `ig` to global index
    integer, allocatable, intent(out) :: mill_local(:,:)
      !! Integer coefficients for G-vectors
    integer, intent(out) :: nGVecsLocal
      !! Local number of G-vectors on this processor


    ! Local variables:
    integer :: ig_l, ig_g
      !! Loop indices
    integer :: ngr
      !! Number of G-vectors left over after evenly divided across processors


    if( nGVecsGlobal > 0 ) then
      nGVecsLocal = nGVecsGlobal/nProcs
        !!  * Calculate number of G-vectors per processor

      ngr = nGVecsGlobal - nGVecsLocal*nProcs 
        !! * Calculate the remainder

      if( myid < ngr ) nGVecsLocal = nGVecsLocal + 1
        !! * Assign the remainder to the first `ngr` processors

      !> * Generate an array to map a local index
      !>   (`ig` passed to `gIndexLocalToGlobal`) to a global
      !>   index (the value stored at `gIndexLocalToGlobal(ig)`)
      !>   and get local miller indices
      allocate(gIndexLocalToGlobal(nGVecsLocal))
      allocate(mill_local(3,nGVecsLocal))

      ig_l = 0
      do ig_g = 1, nGVecsGlobal

        if (myid == mod(ig_g-1,nProcs)) then
        
          ig_l = ig_l + 1
          gIndexLocalToGlobal(ig_l) = ig_g
          mill_local(:,ig_l) = gVecMillerIndicesGlobal(:,ig_g)

        endif

      enddo

      if (ig_l /= nGVecsLocal) call exitError('distributeGvecsOverProcessors', 'unexpected number of G-vecs for this processor', 1)

    endif

    return
  end subroutine distributeGvecsOverProcessors

!----------------------------------------------------------------------------
  subroutine reconstructFFTGrid(nGVecsLocal, gIndexLocalToGlobal, maxGkNum, nKPoints, nPWs1kGlobal, recipSpaceLatticeVectors, gVecInCart, &
      wfcVecCut, kPosition, gKIndexLocalToGlobal, gToGkIndexMap, nGkLessECutLocal, nGkLessECutGlobal, maxGIndexGlobal, maxNumPWsGlobal, maxNumPWsPool)
    !! Determine which G-vectors result in \(G+k\)
    !! below the energy cutoff for each k-point and
    !! sort the indices based on \(|G+k|^2\)
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsLocal
      !! Number of G-vectors on this processor

    integer, intent(in) :: gIndexLocalToGlobal(nGVecsLocal)
      ! Converts local index `ig` to global index
    integer, intent(in) :: maxGkNum
      !! Maximum number of \(G+k\) combinations
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nPWs1kGlobal(nKPoints)
      !! Input number of plane waves for a single k-point

    real(kind=dp), intent(in) :: recipSpaceLatticeVectors(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: gVecInCart(3,nGVecsLocal)
      !! G-vectors in Cartesian coordinates
    real(kind=dp), intent(in) :: wfcVecCut
      !! Energy cutoff converted to vector cutoff
    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space


    ! Output variables:
    integer, allocatable, intent(out) :: gKIndexLocalToGlobal(:,:)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point
    integer, allocatable, intent(out) :: gToGkIndexMap(:,:)
      !! Index map from \(G\) to \(G+k\);
      !! indexed up to `nGVecsLocal` which
      !! is greater than `maxNumPWsPool` and
      !! stored for each k-point
    integer, allocatable, intent(out) :: nGkLessECutLocal(:)
      !! Number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each
      !! k-point, on this processor
    integer, allocatable, intent(out) :: nGkLessECutGlobal(:)
      !! Global number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each k-point
    integer, intent(out) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(out) :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with energy
      !! less than `wfcECut` among all k-points
    integer, intent(out) :: maxNumPWsPool
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! pool


    ! Local variables:
    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp) :: gkMod(nkPerPool,nGVecsLocal)
      !! \(|G+k|^2\);
      !! only stored if less than `wfcVecCut`
    real(kind=dp) :: q
      !! \(|q|^2\) where \(q = G+k\)
    real(kind=dp) :: xkCart(3)
      !! Cartesian coordinates for given k-point

    integer :: ik, ig, ix
      !! Loop indices
    integer, allocatable :: igk(:)
      !! Index map from \(G\) to \(G+k\)
      !! indexed up to `maxNumPWsPool`
    integer :: ngk_tmp
      !! Temporary variable to hold `nGkLessECutLocal`
      !! value so that don't have to keep accessing
      !! array
    integer :: maxGIndexLocal
      !! Maximum G-vector index among all \(G+k\)
      !! for just this processor
    integer :: maxNumPWsLocal
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    allocate(nGkLessECutLocal(nkPerPool))
    allocate(gToGkIndexMap(nkPerPool,nGVecsLocal))
    
    maxNumPWsLocal = 0
    nGkLessECutLocal(:) = 0
    gToGkIndexMap(:,:) = 0

    if (ionode) then
      write(iostd,*)
      write(iostd,*) "***************"
      write(iostd,*) "Determining G+k combinations less than energy cutoff"
    endif

    do ik = 1, nkPerPool
      !! * For each \(G+k\) combination, calculate the 
      !!   magnitude and, if it is less than the energy
      !!   cutoff, store the G index and magnitude and 
      !!   increment the number of \(G+k\) vectors at
      !!   this k-point. Also, keep track of the maximum 
      !!   number of \(G+k\) vectors among all k-points
      !!
      !! @note
      !!  All of the above calculations are local to a single
      !!  processor.
      !! @endnote

      if (ionode) write(iostd,*) "Processing k-point ", ik

      do ix = 1, 3
        xkCart(ix) = sum(kPosition(:,ik+ikStart-1)*recipSpaceLatticeVectors(ix,:))
      enddo

      ngk_tmp = 0

      do ig = 1, nGVecsLocal

        q = sqrt(sum((xkCart(:) + gVecInCart(:,ig))**2))
          ! Calculate \(|G+k|\)

        if (q <= eps8) q = 0.d0

        !if (ionode) write(89,*) ik+ikStart-1, ig, q <= wfcVecCut

        if (q <= wfcVecCut) then

          ngk_tmp = ngk_tmp + 1
            ! If \(|G+k| \leq \) `vcut` increment the count for
            ! this k-point

          gkMod(ik,ngk_tmp) = q
            ! Store the modulus for sorting

          gToGkIndexMap(ik,ngk_tmp) = ig
            ! Store the index for this G-vector

        !else

          !if (sqrt(sum(gVecInCart(:, ig)**2)) .gt. &
            !sqrt(sum(kPosition(:,ik+ikStart-1)**2) + sqrt(wfcVecCut))) goto 100
            ! if |G| > |k| + sqrt(Ecut)  stop search
            !! @todo Figure out if there is valid exit check for `ig` loop @endtodo

        endif
      enddo

      if (ngk_tmp == 0) call exitError('reconstructFFTGrid', 'no G+k vectors on this processor', 1) 

100   maxNumPWsLocal = max(maxNumPWsLocal, ngk_tmp)
        ! Track the maximum number of \(G+k\)
        ! vectors among all k-points

      nGkLessECutLocal(ik) = ngk_tmp
        ! Store the total number of \(G+k\)
        ! vectors for this k-point

    enddo

    allocate(nGkLessECutGlobal(nKPoints))
    nGkLessECutGlobal = 0
    nGkLessECutGlobal(ikStart:ikEnd) = nGkLessECutLocal(1:nkPerPool)
    CALL mpiSumIntV(nGkLessECutGlobal, worldComm)
      !! * Calculate the global number of \(G+k\) 
      !!   vectors for each k-point
      
    if (ionode) then

      do ik = 1, nKPoints

        if (nGkLessECutGlobal(ik) .ne. nPWs1kGlobal(ik)) call exitError('reconstructFFTGrid', &
          'computed no. of G-vectors != input no. of plane waves', 1)
          !! * Make sure that number of G-vectors isn't higher than the calculated maximum

        if (nGkLessECutGlobal(ik) .gt. maxGkNum) call exitError('reconstructFFTGrid', &
          'G-vector count exceeds estimate', 1)
          !! * Make sure that number of G-vectors isn't higher than the calculated maximum

      enddo
    endif

    if (ionode) then
      write(iostd,*) "Done determining G+k combinations less than energy cutoff"
      write(iostd,*) "***************"
      write(iostd,*)
    endif

    if (maxNumPWsLocal <= 0) call exitError('reconstructFFTGrid', &
                'No plane waves found: running on too many processors?', 1)
      !! * Make sure that each processor gets some \(G+k\) vectors. If not,
      !!   should rerun with fewer processors.

    call MPI_ALLREDUCE(maxNumPWsLocal, maxNumPWsPool, 1, MPI_INTEGER, MPI_MAX, interPoolComm, ierr)
    if(ierr /= 0) call exitError('reconstructFFTGrid', 'error in mpi_allreduce 1', ierr)
      !! * When using pools, set `maxNumPWsPool` to the maximum value of `maxNumPWsLocal` 
      !!   in the pool 


    allocate(gKIndexLocalToGlobal(maxNumPWsPool,nkPerPool))
    allocate(igk(maxNumPWsPool))

    gKIndexLocalToGlobal = 0
    igk = 0

    if (ionode) then
      write(iostd,*)
      write(iostd,*) "***************"
      write(iostd,*) "Sorting G+k combinations by magnitude"
    endif

    do ik = 1, nkPerPool
      !! * Reorder the indices of the G-vectors so that
      !!   they are sorted by \(|G+k|^2\) for each k-point

      ngk_tmp = nGkLessECutLocal(ik)

      igk(1:ngk_tmp) = gToGkIndexMap(ik,1:ngk_tmp)

      call hpsort_eps(ngk_tmp, gkMod(ik,:), igk, eps8)
        ! Order vector `gkMod` keeping initial position in `igk`

      do ig = 1, ngk_tmp
        
        gKIndexLocalToGlobal(ig,ik) = gIndexLocalToGlobal(igk(ig))
        
      enddo
     
      gKIndexLocalToGlobal(ngk_tmp+1:maxNumPWsPool, ik) = 0

    enddo

    if (ionode) then
      write(iostd,*) "Done sorting G+k combinations by magnitude"
      write(iostd,*) "***************"
      write(iostd,*)
    endif

    deallocate(igk)

    maxGIndexLocal = maxval(gKIndexLocalToGlobal(:,:))
    call MPI_ALLREDUCE(maxGIndexLocal, maxGIndexGlobal, 1, MPI_INTEGER, MPI_MAX, worldComm, ierr)
    if(ierr /= 0) call exitError('reconstructFFTGrid', 'error in mpi_allreduce 2', ierr)
      !! * Calculate the maximum G-vector index 
      !!   among all \(G+k\) and processors

    maxNumPWsGlobal = maxval(nGkLessECutGlobal(1:nKPoints))
      !! * Calculate the maximum number of G-vectors 
      !!   among all k-points

    return
  end subroutine reconstructFFTGrid

!----------------------------------------------------------------------------
  subroutine hpsort_eps(n, ra, ind, eps)
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! considering two elements equal if their difference is less than `eps`
    !!
    !! n is input, ra is replaced on output by its sorted rearrangement.
    !! Create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (ra).
    !! In case of equal values in the data array (ra) the values in the
    !! index array (ind) are used to order the entries.
    !!
    !! if on input ind(1)  = 0 then indices are initialized in the routine,
    !! if on input ind(1) != 0 then indices are assumed to have been
    !!                initialized before entering the routine and these
    !!                indices are carried around during the sorting process
    !!
    !!
    !! From QE code, adapted from Numerical Recipes pg. 329 (new edition)
    !!

    implicit none

    ! Input/Output variables:
    real(kind=dp), intent(in) :: eps
    integer, intent(in) :: n

    integer, intent(inout) :: ind(:)
    real(kind=dp), intent(inout) :: ra (:)


    ! Local variables
    integer :: i, ir, j, l, iind
    real(kind=dp) :: rra

    ! initialize index array
    if (ind (1) .eq.0) then
      do i = 1, n
        ind (i) = i
      enddo
    endif
    ! nothing to order
    if (n.lt.2) return
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1

    ir = n

    sorting: do

      ! still in hiring phase
      if ( l .gt. 1 ) then
        l    = l - 1
        rra  = ra (l)
        iind = ind (l)
        ! in retirement-promotion phase.
      else
        ! clear a space at the end of the array
        rra  = ra (ir)
        !
        iind = ind (ir)
        ! retire the top of the heap into it
        ra (ir) = ra (1)
        !
        ind (ir) = ind (1)
        ! decrease the size of the corporation
        ir = ir - 1
        ! done with the last promotion
        if ( ir .eq. 1 ) then
          ! the least competent worker at all !
          ra (1)  = rra
          !
          ind (1) = iind
          exit sorting
        endif
      endif
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      do while ( j .le. ir )
        if ( j .lt. ir ) then
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then
            if (ra(j).lt.ra(j+1)) j = j + 1
          else
            ! this means ra(j) == ra(j+1) within tolerance
            if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
        endif
        ! demote rra
        if ( abs(rra - ra(j)).ge.eps ) then
          if (rra.lt.ra(j)) then
            ra (i) = ra (j)
            ind (i) = ind (j)
            i = j
            j = j + j
          else
            ! set j to terminate do-while loop
            j = ir + 1
          end if
        else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind (j) ) then
            ra (i) = ra (j)
            ind (i) = ind (j)
            i = j
            j = j + j
          else
            ! set j to terminate do-while loop
            j = ir + 1
          endif
        end if
      enddo
      ra (i) = rra
      ind (i) = iind

    end do sorting

    return 
  end subroutine hpsort_eps

!----------------------------------------------------------------------------
  subroutine read_vasprun_xml(realSpaceLatticeVectors, nKPoints, VASPDir, eFermi, kWeight, iType, nAtoms, nAtomTypes)
    !! Read the k-point weights and cell info from the `vasprun.xml` file
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: realSpaceLatticeVectors(3,3)
      !! Real space lattice vectors

    integer, intent(in) :: nKPoints
      !! Total number of k-points

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files


    ! Output variables:
    real(kind=dp), intent(out) :: eFermi
      !! Fermi energy
    real(kind=dp), allocatable, intent(out) :: kWeight(:)
      !! K-point weights

    integer, allocatable, intent(out) :: iType(:)
      !! Atom type index
    integer, intent(out) :: nAtoms
      !! Number of atoms
    integer, intent(out) :: nAtomTypes
      !! Number of types of atoms


    ! Local variables:
    real(kind=dp) :: dir(3)
      !! Direct coordinates read from file

    integer :: ik, ia, i, ix
      !! Loop indices

    character(len=256) :: cDum
      !! Dummy variable to ignore input
    character(len=256) :: fileName
      !! `vasprun.xml` with path
    character(len=256) :: line
      !! Line read from file

    logical :: fileExists
      !! If the `vasprun.xml` file exists
    logical :: found
      !! If the required tag was found

    allocate(kWeight(nKPoints))

    if (ionode) then

      fileName = trim(VASPDir)//'/vasprun.xml'

      inquire(file = fileName, exist = fileExists)

      if (.not. fileExists) call exitError('getKPointWeights', 'Required file vasprun.xml does not exist', 1)

      open(57, file=fileName)
        !! * If root node, open `vasprun.xml`


      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'weights'`, indicating the
        !!   tag surrounding the k-point weights
        
        read(57, '(A)') line

        if (index(line,'weights') /= 0) found = .true.
        
      enddo

      do ik = 1, nKPoints
        !! * Read in the weight for each k-point

        read(57,*) cDum, kWeight(ik), cDum

      enddo


      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'atominfo'`, indicating the
        !!   tag surrounding the cell info
        
        read(57, '(A)') line

        if (index(line,'atominfo') /= 0) found = .true.
        
      enddo

      read(57,*) cDum, nAtoms, cDum
      read(57,*) cDum, nAtomTypes, cDum
      read(57,*) 
      read(57,*) 
      read(57,*) 
      read(57,*) 
      read(57,*) 

      allocate(iType(nAtoms))

      do ia = 1, nAtoms
        !! * Read in the atom type index for each atom

        read(57,'(a21,i3,a9)') cDum, iType(ia), cDum

      enddo

      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'efermi'`, indicating the
        !!   tag with the Fermi energy
        
        read(57, '(A)') line

        if (index(line,'efermi') /= 0) found = .true.
        
      enddo

      read(line,*) cDum, cDum, eFermi, cDum
      eFermi = eFermi*eVToRy

      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'finalpos'`, indicating the
        !!   tag surrounding the final cell parameters
        !!   and positions
        
        read(57, '(A)') line

        if (index(line,'finalpos') /= 0) found = .true.
        
      enddo

      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'positions'`, indicating the
        !!   tag surrounding the final positions
        
        read(57, '(A)') line

        if (index(line,'positions') /= 0) found = .true.
        
      enddo

      allocate(atomPositions(3,nAtoms))

      do ia = 1, nAtoms
        !! * Read in the final position for each atom

        read(57,*) cDum, (dir(i),i=1,3), cDum
          !! @note
          !!  I assume that the coordinates are always direct
          !!  in the `vasprun.xml` file and that the scaling
          !!  factor is already included as I cannot find it 
          !!  listed anywhere in that file. Extensive testing
          !!  needs to be done to confirm this assumption.
          !! @endnote

        do ix = 1, 3

          atomPositions(ix,ia) = sum(dir(:)*realSpaceLatticeVectors(ix,:))
            !! @todo Test logic of direct to cartesian coordinates with scaling factor @endtodo

        enddo

      enddo

    endif

    call MPI_BCAST(kWeight, size(kWeight), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nAtomTypes, 1, MPI_INTEGER, root, worldComm, ierr)

    if (.not. ionode) allocate(iType(nAtoms))
    if (.not. ionode) allocate(atomPositions(3,nAtoms))

    call MPI_BCAST(iType, size(iType), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(atomPositions, size(atomPositions), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return
  end subroutine read_vasprun_xml

!----------------------------------------------------------------------------
  subroutine readPOTCAR(nAtomTypes, VASPDir, ps)
    !! Read PAW pseudopotential information from POTCAR
    !! file
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    
    ! Output variables:
    type (pseudo) :: ps(nAtomTypes)
      !! Holds all information needed from pseudopotential


    ! Local variables:
    real(kind=dp) :: dummyD(1000)
      !! Dummy variable to ignore input
    real(kind=dp), allocatable :: dummyDA1(:), dummyDA2(:,:)
      !! Allocatable dummy variable to ignore input
    real(kind=dp) :: H
      !! Factor for generating derivative of 
      !! radial grid

    integer :: angMom
      !! Angular momentum of projectors
    integer :: iT, i, j, ip, ir
      !! Loop indices
    integer :: nProj
      !! Number of projectors with given angular momentum

    character(len=1) :: charSwitch
      !! Switch to determine what section reading
    character(len=256) :: dummyC
      !! Dummy character to ignore input
    character(len=256) :: fileName
      !! Full WAVECAR file name including path

    logical :: found
      !! If the required tag was found


    if(ionode) then
      fileName = trim(VASPDir)//'/POTCAR'

      open(unit=potcarUnit, file=fileName, iostat=ierr, status='old')
      if (ierr .ne. 0) write(iostd,*) 'open error - iostat =', ierr
        !! * If root node, open the `POTCAR` file

      do iT = 1, nAtomTypes

        ps(iT)%nChannels = 0
        ps(iT)%lmmax = 0

        read(potcarUnit,*) dummyC, ps(iT)%element, dummyC
          !! * Read in the header
        read(potcarUnit,*)
          !! * Ignore the valence line
        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read in character switch to determine if there is a 
          !!   PSCRT section (switch not used)
          !! @note
          !!  Some of the switches do not actually seem to be used
          !!  as a switch because the following code does not include
          !!  any logic to process the switch. If the POTCAR files 
          !!  ever have a different form than assumed here, the logic
          !!  will need to be updated.
          !! @endnote

        found = .false.
        do while (.not. found)
          !! * Ignore all lines until you get to the `END` of
          !!   the PSCRT section
        
          read(potcarUnit, '(A)') dummyC

          if (dummyC(1:3) == 'END') found = .true.
        
        enddo

        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch (switch not used)
        read(potcarUnit,*)
          !! * Ignore the max G for local potential
        read(potcarUnit,*) (dummyD(i), i=1,1000)
          !! * Ignore the local pseudopotential in reciprocal
          !!   space
        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch

        if (charSwitch == 'g') then
          !! * Ignore gradient correction

          read(potcarUnit,*)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        if (charSwitch == 'c') then
          !! * Ignore core charge density

          read(potcarUnit,*) (dummyD(i), i=1,1000)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        if (charSwitch == 'k') then
          !! * Ignore partial kinetic energy density

          read(potcarUnit,*) (dummyD(i), i=1,1000)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        if (charSwitch == 'K') then
          !! * Ignore kinetic energy density

          read(potcarUnit,*) (dummyD(i), i=1,1000)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        read(potcarUnit,*) (dummyD(i), i=1,1000)
          !! * Ignore the atomic pseudo charge density

        read(potcarUnit,*) 
          !! * Ignore the max G for non-local potential and 
          !!   unused boolean (`LDUM` in VASP)
        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch

        do while (charSwitch /= 'D' .and. charSwitch /= 'A' .and. charSwitch /= 'P' &
          .and. charSwitch /= 'E')
            !! * Until you have read in all of the momentum channels
            !!   (i.e. you get to a character switch that is not `'N'`)
            !!     * Read in the angular momentum, the number of 
            !!       projectors at this angular momentum, and the max
            !!       r for the non-local contribution
            !!     * Increment the number of nlm channels
            !!     * Ignore non-local strength multipliers
            !!     * Read in the reciprocal-space and real-space
            !!       projectors
            !!     * Increment the number of l channels
            !!     * Read the next character switch

          read(potcarUnit,*) angMom, nProj, ps(iT)%psRMax
            ! Read in angular momentum, the number of projectors
            ! at this angular momentum, and the max r for the 
            ! non-local contribution

          ps(iT)%lmmax = ps(iT)%lmmax + (2*angMom+1)*nProj
            ! Increment the number of nlm channels

          allocate(dummyDA2(nProj,nProj))

          read(potcarUnit,*) dummyDA2(:,:)
            ! Ignore non-local strength multipliers

          do ip = 1, nProj
            ! Read in the reciprocal-space and real-space
            ! projectors

            ps(iT)%angmom(ps(iT)%nChannels+ip) = angmom

            read(potcarUnit,*) 
            read(potcarUnit,*) (ps(iT)%recipProj(ps(iT)%nChannels+ip,i), i=1,100)
            read(potcarUnit,*) 
            read(potcarUnit,*) (ps(iT)%realProj(ps(iT)%nChannels+ip,i), i=1,100)
              !! @todo Figure out if you actually need to read these projectors @endtodo

          enddo

          ps(iT)%nChannels = ps(iT)%nChannels + nProj
            ! Increment the number of l channels

          deallocate(dummyDA2)

          read(potcarUnit,'(1X,A1)') charSwitch
            ! Read character switch

        enddo

        if (charSwitch /= 'P') then
          !! * Ignore depletion charges

          read(potcarUnit,*)
          read(potcarUnit,*)
          
        else

          read(potcarUnit,*) ps(iT)%nmax, ps(iT)%rAugMax  
            !! * Read the number of mesh grid points and
            !!   the maximum radius in the augmentation sphere
          read(potcarUnit,*)
            !! * Ignore format specifier
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch (not used)

          allocate(ps(iT)%radGrid(ps(iT)%nmax))
          allocate(ps(iT)%dRadGrid(ps(iT)%nmax))
          allocate(ps(iT)%wps(ps(iT)%nChannels,ps(iT)%nmax))
          allocate(ps(iT)%wae(ps(iT)%nChannels,ps(iT)%nmax))
          allocate(dummyDA2(ps(iT)%nChannels, ps(iT)%nChannels))

          read(potcarUnit,*) dummyDA2(:,:)
            !! * Ignore augmentation charges
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch

          if (charSwitch == 't') then
            !! * Ignore total charge in each channel 

            read(potcarUnit,*) dummyDA2(:,:)
            read(potcarUnit,*) 

          endif

          read(potcarUnit,*) dummyDA2
            !! * Ignore initial occupancies in atom

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'g') call exitError('readPOTCAR', 'expected grid section', 1)

          read(potcarUnit,*) (ps(iT)%radGrid(i), i=1,ps(iT)%nmax)

          ps(iT)%rAugMax = ps(iT)%rAugMax*angToBohr 
          ps(iT)%radGrid(:) = ps(iT)%radGrid(:)*angToBohr

          H = log(ps(iT)%radGrid(ps(iT)%nmax)/ps(iT)%radGrid(1))/(ps(iT)%nmax - 1)
            !! * Calculate \(H\) which is used to generate the derivative of the grid
            !!
            !! @note
            !!  The grid in VASP is defined as \(R_i = R_0e^{H(i-1)}\), so we define the
            !!  derivative as \(dR_i = R_0He^{H(i-1)}\)
            !! @endnote
          
          found = .false.
          do ir = 1, ps(iT)%nmax
            !! * Calculate the max index of the augmentation sphere and
            !!   the derivative of the radial grid

            if (.not. found .and. ps(iT)%radGrid(ir) > ps(iT)%rAugMax) then
              ps(iT)%iRAugMax = ir - 1
              found = .true.
            endif

            ps(iT)%dRadGrid(ir) = ps(iT)%radGrid(1)*H*exp(H*(ir-1))

          enddo

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'a') call exitError('readPOTCAR', 'expected aepotential section', 1)

          allocate(dummyDA1(ps(iT)%nmax))

          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore AE potential

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'c') call exitError('readPOTCAR', 'expected core charge-density section', 1)

          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore the frozen core charge
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch

          if (charSwitch == 'k') then
            !! * Ignore kinetic energy density

            read(potcarUnit,*) dummyDA1(:)
            read(potcarUnit,'(1X,A1)') charSwitch

          endif

          if (charSwitch == 'm') then
            !! * Ignore pseudo-ized kinetic energy density

            read(potcarUnit,*) dummyDA1(:)
            read(potcarUnit,'(1X,A1)') charSwitch

          endif

          if (charSwitch == 'l') then
            !! * Ignore local pseudopotential core

            read(potcarUnit,*) dummyDA1(:)
            read(potcarUnit,'(1X,A1)') charSwitch

          endif

          if (charSwitch /= 'p') call exitError('readPOTCAR', 'expected pspotential section', 1)
          
          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore PS potential

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'c') call exitError('readPOTCAR', 'expected core charge-density section', 1)
          
          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore core charge density

          do ip = 1, ps(iT)%nChannels
            !! * Read the AE and PS partial waves for each projector
            
            read(potcarUnit,'(1X,A1)') charSwitch
            if (charSwitch /= 'p') call exitError('readPOTCAR', 'expected pseudowavefunction section', 1)
            read(potcarUnit,*) (ps(iT)%wps(ip,i), i=1,ps(iT)%nmax)

            read(potcarUnit,'(1X,A1)') charSwitch
            if (charSwitch /= 'a') call exitError('readPOTCAR', 'expected aewavefunction section', 1)
            read(potcarUnit,*) (ps(iT)%wae(ip,i), i=1,ps(iT)%nmax)

          enddo

          ps(iT)%wps(:,:) = ps(iT)%wps(:,:)/sqrt(angToBohr)
          ps(iT)%wae(:,:) = ps(iT)%wae(:,:)/sqrt(angToBohr)
            !! @todo Make sure that partial wave unit conversion makes sense @endtodo

          deallocate(dummyDA1)
          deallocate(dummyDA2)

        endif

        found = .false.
        do while (.not. found)
          !! * Ignore all lines until you get to the `End of Dataset`
        
          read(potcarUnit, '(A)') dummyC

          if (index(dummyC,'End of Dataset') /= 0) found = .true.
        
        enddo

      enddo

    endif    

    return
  end subroutine readPOTCAR

!----------------------------------------------------------------------------
  subroutine writeKInfo(nKPoints, maxNumPWsPool, gKIndexLocalToGlobal, nBands, nGkLessECutGlobal, nGkLessECutLocal, &
      maxGIndexGlobal, maxNumPWsGlobal, bandOccupation, kWeight, kPosition, gKIndexGlobal)
    !! Calculate the highest occupied band for each k-point,
    !! gather the \(G+k\) vector indices in single, global 
    !! array, and write out k-point information
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: maxNumPWsPool
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    integer, intent(in) :: gKIndexLocalToGlobal(maxNumPWsPool, nkPerPool)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point;
      !! the first index goes up to `maxNumPWsPool`,
      !! but only valid values are up to `nGkLessECutLocal`
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nGkLessECutGlobal(nKPoints)
      !! Global number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each k-point
    integer, intent(in) :: nGkLessECutLocal(nkPerPool)
      !! Number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each
      !! k-point, on this processor
    integer, intent(in) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(in) :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with energy
      !! less than `wfcECut` among all k-points

    real(kind=dp), intent(in) :: bandOccupation(nBands, nKPoints)
      !! Occupation of band
    real(kind=dp), intent(in) :: kWeight(nKPoints)
      !! K-point weights
    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space


    ! Output variables:
    integer, allocatable, intent(out) :: gKIndexGlobal(:,:)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors


    ! Local variables:
    integer, allocatable :: groundState(:)
      !! Holds the highest occupied band
      !! for each k-point
    integer :: ik, ig
      !! Loop indices


    if(ionode) then

      write(iostd,*)
      write(iostd,*) "***************"
      write(iostd,*) "Getting ground state bands"
    
      write(mainOutFileUnit, '("# Number of K-points. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nKPoints
      write(mainOutFileUnit, '("# ik, groundState, nGkLessECutGlobal(ik), wk(ik), xk(1:3,ik). Format: ''(3i10,4ES24.15E3)''")')
      flush(mainOutFileUnit)
    
      allocate(groundState(nKPoints))

      call getGroundState(nBands, nKPoints, bandOccupation, groundState)
        !! * For each k-point, find the index of the 
        !!   highest occupied band
        !!
        !! @note
        !!  Although `groundState` is written out in `Export`,
        !!  it is not currently used by the `TME` program.
        !! @endtodo

      write(iostd,*) "Done getting ground state bands"
      write(iostd,*) "***************"
      write(iostd,*)

    endif

    if(ionode) then

      write(iostd,*)
      write(iostd,*) "***************"
      write(iostd,*) "Getting global G+k indices"

    endif
  
    allocate(gKIndexGlobal(maxNumPWsGlobal, nKPoints))
  
    gKIndexGlobal(:,:) = 0
    do ik = 1, nKPoints

      if (ionode) write(iostd,*) "Processing k-point ", ik

      call getGlobalGkIndices(nKPoints, maxNumPWsPool, gKIndexLocalToGlobal, ik, nGkLessECutGlobal, nGkLessECutLocal, maxGIndexGlobal, &
          maxNumPWsGlobal, gKIndexGlobal)
        !! * For each k-point, gather all of the \(G+k\) indices
        !!   among all processors in a single global array
    
      if (ionode) write(mainOutFileUnit, '(3i10,4ES24.15E3)') ik, groundState(ik), nGkLessECutGlobal(ik), kWeight(ik), kPosition(1:3,ik)
      if (ionode) flush(mainOutFileUnit)
        !! * Write the k-point index, the ground state band, and
        !!   the number of G-vectors, weight, and position for this 
        !!   k-point
    
    enddo

    if(ionode) then

      write(iostd,*) "Done getting global G+k indices"
      write(iostd,*) "***************"
      write(iostd,*)
      flush(iostd)

    endif

    if (ionode) deallocate(groundState)

    return
  end subroutine writeKInfo

!----------------------------------------------------------------------------
  subroutine getGroundState(nBands, nKPoints, bandOccupation, groundState)
    !! * For each k-point, find the index of the 
    !!   highest occupied band

    implicit none

    ! Input variables:
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nKPoints
      !! Total number of k-points

    real(kind=dp), intent(in) :: bandOccupation(nBands, nKPoints)
      !! Occupation of band

    
    ! Output variables:
    integer, intent(out) :: groundState(nKPoints)
      !! Holds the highest occupied band
      !! for each k-point


    ! Local variables:
    integer :: ik, ibnd
      !! Loop indices


    groundState(:) = 0
    do ik = 1, nKPoints

      do ibnd = 1, nBands

        if (bandOccupation(ibnd,ik) < 0.5_dp) then
          !! @todo Figure out if boundary for "occupied" should be 0.5 or less @endtodo
        !if (et(ibnd,ik) > ef) then

          groundState(ik) = ibnd - 1
          goto 10

        endif
      enddo

10    continue

    enddo

    return
  end subroutine getGroundState

!----------------------------------------------------------------------------
  subroutine getGlobalGkIndices(nKPoints, maxNumPWsPool, gKIndexLocalToGlobal, ik, nGkLessECutGlobal, nGkLessECutLocal, maxGIndexGlobal, &
      maxNumPWsGlobal, gKIndexGlobal)
    !! Gather the \(G+k\) vector indices in single, global 
    !! array
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: maxNumPWsPool
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    integer, intent(in) :: gKIndexLocalToGlobal(maxNumPWsPool, nkPerPool)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point;
      !! the first index goes up to `maxNumPWsPool`,
      !! but only valid values are up to `nGkLessECutLocal`
    integer, intent(in) :: ik
      !! Index of current k-point
    integer, intent(in) :: nGkLessECutGlobal(nKPoints)
      !! Global number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each k-point
    integer, intent(in) :: nGkLessECutLocal(nkPerPool)
      !! Number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each
      !! k-point, on this processor
    integer, intent(in) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(in) :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with energy
      !! less than `wfcECut` among all k-points


    ! Output variables:
    integer, intent(out) :: gKIndexGlobal(maxNumPWsGlobal, nKPoints)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors


    ! Local variables:
    integer :: ig
    integer, allocatable :: itmp1(:)
      !! Global \(G+k\) indices for single
      !! k-point with zeros for G-vector indices
      !! where \(G+k\) was greater than the cutoff
    integer :: ngg 
      !! Counter for \(G+k\) vectors for given
      !! k-point; should equal `nGkLessECutGlobal`

    
    allocate(itmp1(maxGIndexGlobal), stat=ierr)
    if (ierr/= 0) call exitError('getGlobalGkIndices','allocating itmp1', abs(ierr))

    itmp1 = 0
    if(ik >= ikStart .and. ik <= ikEnd) then

      do ig = 1, nGkLessECutLocal(ik-ikStart+1)

        itmp1(gKIndexLocalToGlobal(ig, ik-ikStart+1)) = gKIndexLocalToGlobal(ig, ik-ikStart+1)
          !! * For each k-point and \(G+k\) vector for this processor,
          !!   store the local to global indices (`gKIndexLocalToGlobal`) in an
          !!   array that will later be combined globally
          !!
          !! @note
          !!  This will leave zeros in spots where the \(G+k\) 
          !!  combination for this k-point was greater than the energy 
          !!  cutoff.
          !! @endnote

      enddo
    endif

    call mpiSumIntV(itmp1, worldComm)

    ngg = 0
    do  ig = 1, maxGIndexGlobal

      if(itmp1(ig) == ig) then
        !! * Go through and find all of the non-zero
        !!   indices in the now-global `itmp1` array,
        !!   and store them in a new array that won't
        !!   have the extra zeros

        ngg = ngg + 1

        gKIndexGlobal(ngg, ik) = ig

      endif
    enddo


    if(ionode .and. ngg /= nGkLessECutGlobal(ik)) call exitError('writeKInfo', 'Unexpected number of G+k vectors', 1)
      !! * Make sure that the total number of non-zero
      !!   indices matches the global number of \(G+k\)
      !!   vectors for this k-point
    
    deallocate( itmp1 )

    return
  end subroutine getGlobalGkIndices

!----------------------------------------------------------------------------
  subroutine writeGridInfo(nGVecsGlobal, nKPoints, maxNumPWsGlobal, gKIndexGlobal, gVecMillerIndicesGlobal, nGkLessECutGlobal, maxGIndexGlobal, exportDir)
    !! Write out grid boundaries and miller indices
    !! for just \(G+k\) combinations below cutoff energy
    !! in one file and all miller indices in another 
    !! file
    !!
    !! <h2>Walkthrough</h2>
    !!

    use miscUtilities, only: int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with energy
      !! less than `wfcECut` among all k-points

    integer, intent(in) :: gKIndexGlobal(maxNumPWsGlobal, nKPoints)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors
    integer, intent(in) :: gVecMillerIndicesGlobal(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(in) :: nGkLessECutGlobal(nKPoints)
      !! Global number of \(G+k\) vectors with energy
      !! less than `wfcECut` for each k-point
    integer, intent(in) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors

    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export


    ! Output variables:


    ! Local variables:
    integer :: ik, ig, igk
      !! Loop indices

    character(len=300) :: indexC
      !! Character index


    if (ionode) then
    
      !> * Write the global number of G-vectors, the maximum
      !>   G-vector index, and the max/min miller indices
      write(mainOutFileUnit, '("# Number of G-vectors. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nGVecsGlobal
    
      write(mainOutFileUnit, '("# Number of PW-vectors. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') maxGIndexGlobal
    
      write(mainOutFileUnit, '("# Number of min - max values of fft grid in x, y and z axis. Format: ''(6i10)''")')
      write(mainOutFileUnit, '(6i10)') minval(gVecMillerIndicesGlobal(1,1:nGVecsGlobal)), maxval(gVecMillerIndicesGlobal(1,1:nGVecsGlobal)), &
                          minval(gVecMillerIndicesGlobal(2,1:nGVecsGlobal)), maxval(gVecMillerIndicesGlobal(2,1:nGVecsGlobal)), &
                          minval(gVecMillerIndicesGlobal(3,1:nGVecsGlobal)), maxval(gVecMillerIndicesGlobal(3,1:nGVecsGlobal))
      flush(mainOutFileUnit)
    
      do ik = 1, nKPoints
        !! * For each k-point, write out the miller indices
        !!   resulting in \(G+k\) vectors less than the energy
        !!   cutoff in a `grid.ik` file
      
        call int2str(ik, indexC)
        open(72, file=trim(exportDir)//"/grid."//trim(indexC))
        write(72, '("# Wave function G-vectors grid")')
        write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
      
        do igk = 1, nGkLessECutGlobal(ik)
          write(72, '(4i10)') gKIndexGlobal(igk,ik), gVecMillerIndicesGlobal(1:3,gKIndexGlobal(igk,ik))
          flush(72)
        enddo
      
        close(72)
      
      enddo

      !> * Output all miller indices in `mgrid` file
      open(72, file=trim(exportDir)//"/mgrid")
      write(72, '("# Full G-vectors grid")')
      write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
    
      do ig = 1, nGVecsGlobal
        write(72, '(4i10)') ig, gVecMillerIndicesGlobal(1:3,ig)
        flush(72)
      enddo
    
      close(72)

    endif

    return
  end subroutine writeGridInfo


!----------------------------------------------------------------------------
  subroutine writeCellInfo(iType, nAtoms, nBands, nAtomTypes, nSpins, realSpaceLatticeVectors, recipSpaceLatticeVectors, atomPositions, nAtomsEachType)
    !! Write out the real- and reciprocal-space lattice vectors, 
    !! the number of atoms, the number of types of atoms, the
    !! final atom positions, number of bands, and number of spins,
    !! then calculate the number of atoms of each type

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms
    integer, intent(in) :: nSpins
      !! Number of spins

    real(kind=dp), intent(in) :: realSpaceLatticeVectors(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: recipSpaceLatticeVectors(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: atomPositions(3,nAtoms)
      !! Atom positions


    ! Output variables:
    integer, allocatable, intent(out) :: nAtomsEachType(:)
      !! Number of atoms of each type


    ! Local variables:
    integer :: i
      !! Loop index


    if (ionode) then
    
      write(mainOutFileUnit, '("# Cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
      write(mainOutFileUnit, '("# a1 ",3ES24.15E3)') realSpaceLatticeVectors(:,1)
      write(mainOutFileUnit, '("# a2 ",3ES24.15E3)') realSpaceLatticeVectors(:,2)
      write(mainOutFileUnit, '("# a3 ",3ES24.15E3)') realSpaceLatticeVectors(:,3)
    
      write(mainOutFileUnit, '("# Reciprocal cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
      write(mainOutFileUnit, '("# b1 ",3ES24.15E3)') recipSpaceLatticeVectors(:,1)
      write(mainOutFileUnit, '("# b2 ",3ES24.15E3)') recipSpaceLatticeVectors(:,2)
      write(mainOutFileUnit, '("# b3 ",3ES24.15E3)') recipSpaceLatticeVectors(:,3)
    
      write(mainOutFileUnit, '("# Number of Atoms. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nAtoms
    
      write(mainOutFileUnit, '("# Number of Types. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nAtomTypes
    
      write(mainOutFileUnit, '("# Atoms type, position(1:3) (a.u.). Format: ''(i10,3ES24.15E3)''")')
      do i = 1, nAtoms
        write(mainOutFileUnit,'(i10,3ES24.15E3)') iType(i), atomPositions(:,i)
      enddo
    
      write(mainOutFileUnit, '("# Number of Bands. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nBands

      write(mainOutFileUnit, '("# Spin. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nSpins
    
      allocate( nAtomsEachType(nAtomTypes) )
      nAtomsEachType = 0
      do i = 1, nAtoms
        nAtomsEachType(iType(i)) = nAtomsEachType(iType(i)) + 1
      enddo

    endif

    return
  end subroutine writeCellInfo

!----------------------------------------------------------------------------
  subroutine writePseudoInfo(nAtomTypes, nAtomsEachType, ps)
    !! For each atom type, write out the element name,
    !! number of atoms of this type, projector info,
    !! radial grid info, and partial waves

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms

    integer, intent(in) :: nAtomsEachType(nAtomTypes)
      !! Number of atoms of each type

    type (pseudo) :: ps(nAtomTypes)
      !! Holds all information needed from pseudopotential


    ! Output variables:


    ! Local variables:
    integer :: iT, ip, ir
      !! Loop index

  
    if (ionode) then

      do iT = 1, nAtomTypes
        
        write(mainOutFileUnit, '("# Element")')
        write(mainOutFileUnit, *) trim(ps(iT)%element)
        write(mainOutFileUnit, '("# Number of Atoms of this type. Format: ''(i10)''")')
        write(mainOutFileUnit, '(i10)') nAtomsEachType(iT)
        write(mainOutFileUnit, '("# Number of projectors. Format: ''(i10)''")')
        write(mainOutFileUnit, '(i10)') ps(iT)%nChannels
        
        write(mainOutFileUnit, '("# Angular momentum, index of the projectors. Format: ''(2i10)''")')
        do ip = 1, ps(iT)%nChannels

          write(mainOutFileUnit, '(2i10)') ps(iT)%angmom(ip), ip

        enddo
        
        write(mainOutFileUnit, '("# Number of channels. Format: ''(i10)''")')
        write(mainOutFileUnit, '(i10)') ps(iT)%lmmax
        
        write(mainOutFileUnit, '("# Number of radial mesh points. Format: ''(2i10)''")')
        write(mainOutFileUnit, '(2i10)') ps(iT)%nmax, ps(iT)%iRAugMax
          ! Number of points in the radial mesh, number of points inside the aug sphere
        
        write(mainOutFileUnit, '("# Radial grid, Integratable grid. Format: ''(2ES24.15E3)''")')
        do ir = 1, ps(iT)%nmax
          write(mainOutFileUnit, '(2ES24.15E3)') ps(iT)%radGrid(ir), ps(iT)%dRadGrid(ir) 
            ! Radial grid, derivative of radial grid
        enddo
        
        write(mainOutFileUnit, '("# AE, PS radial wfc for each beta function. Format: ''(2ES24.15E3)''")')
        do ip = 1, ps(iT)%nChannels
          do ir = 1, ps(iT)%nmax
            write(mainOutFileUnit, '(2ES24.15E3)') ps(iT)%wae(ip,ir), ps(iT)%wps(ip,ir)
          enddo
        enddo
      
      enddo
    
    endif

    return
  end subroutine writePseudoInfo

!----------------------------------------------------------------------------
  subroutine writeEigenvalues(nBands, nKPoints, nSpins, eFermi, bandOccupation, eigenE)
    !! Write Fermi energy and eigenvalues and occupations for each band

    use miscUtilities

    implicit none

    ! Input variables:
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nSpins
      !! Number of spins
      
    real(kind=dp), intent(in) :: eFermi
      !! Fermi energy
    real(kind=dp), intent(in) :: bandOccupation(nBands,nKPoints)
      !! Occupation of band

    complex*16, intent(in) :: eigenE(nSpins,nKPoints,nBands)
      !! Band eigenvalues


    ! Output variables:


    ! Local variables:
    integer :: ik, ib
      !! Loop indices
    integer :: ispin
      !! Spin index

    character(len=300) :: indexC
      !! Character index


    if (ionode ) then
    
      write(mainOutFileUnit, '("# Fermi Energy (Hartree). Format: ''(ES24.15E3)''")')
      write(mainOutFileUnit, '(ES24.15E3)') eFermi*ryToHartree
      flush(mainOutFileUnit)
    
      do ik = 1, nKPoints
      
        !ispin = isk(ik)
        ispin = 1
          !! @todo Figure out if spin needs to be incorporated for eigenvalues @endtodo
      
        call int2str(ik, indexC)
        open(72, file=trim(exportDir)//"/eigenvalues."//trim(indexC))
      
        write(72, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
        write(72, '("# Eigenvalues (Hartree), band occupation number. Format: ''(2ES24.15E3)''")')
      
        do ib = 1, nBands

          write(72, '(2ES24.15E3)') real(eigenE(ispin,ik,ib))*ryToHartree, bandOccupation(ib,ik)
          flush(72)

        enddo
      
        close(72)
      
      enddo
    
    endif

    return
  end subroutine writeEigenvalues

!----------------------------------------------------------------------------
!  subroutine writeProjections()
!
!    implicit none
!
!    IF ( nkb > 0 ) THEN
!      !! @todo Add calculation of total number of projectors #thistask @endtodo
!      !! @todo Add `nkb` as a local variable #thistask @endtodo
!    
!      CALL init_us_1
!        !   This routine performs the following tasks:
!        !   a) For each non vanderbilt pseudopotential it computes the D and
!        !      the betar in the same form of the Vanderbilt pseudopotential.
!        !   b) It computes the indices indv which establish the correspondence
!        !      nh <-> beta in the atom
!        !   c) It computes the indices nhtol which establish the correspondence
!        !      nh <-> angular momentum of the beta function
!        !   d) It computes the indices nhtolm which establish the correspondence
!        !      nh <-> combined (l,m) index for the beta function.
!        !   e) It computes the coefficients c_{LM}^{nm} which relates the
!        !      spherical harmonics in the Q expansion
!        !   f) It computes the radial fourier transform of the Q function on
!        !      all the g vectors
!        !   g) It computes the q terms which define the S matrix.
!        !   h) It fills the interpolation table for the beta functions
!        !
!        ! `nh` is found in the `uspp_param` module and is the number of beta
!        ! functions per atom type which corresponds to `ps(iT)%nChannels` I think
!
!      CALL init_at_1
!        ! This routine computes a table with the radial Fourier transform of the 
!        ! atomic wavefunctions
!    
!      CALL allocate_bec_type (nkb,nBands, becp)
!        ! `nkb` may be the total number of projectors over all atom types 
!        !
!        ! This routine allocates space for the projections \(<\beta|\psi>\),
!        ! and since we are ignoring noncollinear calculations for now, we
!        ! will either have `becp%r(nkb,nbnd)` (real) or `becp%k(nkb,nbnd)`
!        ! (complex)
!
!      allocate(igk(maxNumPWsPool))
!        !! @todo Add `maxNumPWsPool` as an input variable #thistask @endtodo
!      
!      igk = 0
!    
!      DO ik = 1, nKPoints
!        !! @todo Add `ik` as a local variable #thistask @endtodo
!        !! @todo Add `nKPoints` as an input variable #thistask @endtodo
!      
!        local_pw = 0
!          !! @todo Figure out the meaning of `local_pw` #thistask @endtodo
!        IF ( (ik >= ikStart) .and. (ik <= ikEnd) ) THEN
!          CALL davcio (evc, nwordwfc, iunwfc, (ik-ikStart+1), - 1)
!            ! Read the wavefunction for this k-point
!            !! @todo Figure out what the "wavefunction" is here #thistask @endtodo
!
!          igk(1:nGkLessECutLocal(ik-ikStart+1)) = gToGkIndexMap(ik-ikStart+1,1:nGkLessECutLocal(ik-ikStart+1))
!            !! @todo Add `igk` as a local variable #thistask @endtodo
!            !! @todo Add `gToGkIndexMap` as an input variable #thistask @endtodo
!
!          CALL init_us_2(npw, igk, kPosition(1,ik), vkb)
!          local_pw = ngk(ik-ikStart+1)
!
!          IF ( gamma_only ) THEN
!            !! @todo Figure out if need to have `gamma_only` variable here and how to set @endtodo
!            CALL calbec ( nGkLessECutGlobal(ik), vkb, evc, becp )
!            WRITE(0,*) 'Gamma only PW_EXPORT not yet tested'
!              !! @todo Figure out what about `gamma_only` makes the rest of the program break @endtodo
!          ELSE
!            CALL calbec ( npw, vkb, evc, becp )
!            if ( ionode ) then
!
!              WRITE(iostd,*) "Writing projectors of kpt", ik
!
!              file_exists = .false.
!              inquire(file =trim(exportDir)//"/projections"//iotk_index(ik), exist = file_exists)
!              if ( .not. file_exists ) then
!                open(72, file=trim(exportDir)//"/projections"//iotk_index(ik))
!                write(72, '("# Complex projections <beta|psi>. Format: ''(2ES24.15E3)''")')
!                do j = 1,  nBands ! number of bands
!                  !! @todo Verify that `nBands` is equal to `becp%nbnd` #thistask @endtodo
!                  do i = 1, nkb      ! number of projections
!                    write(72,'(2ES24.15E3)') becp%k(i,j)
!                  enddo
!                enddo
!              
!                close(72)
!              
!              endif
!            endif
!          ENDIF
!        ENDIF
!      enddo
!    endif
!
!    return
!  end subroutine writeProjections

!----------------------------------------------------------------------------
  subroutine subroutineTemplate()
    implicit none


    return
  end subroutine subroutineTemplate

end module wfcExportVASPMod
