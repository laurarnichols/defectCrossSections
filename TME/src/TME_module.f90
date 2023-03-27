module declarations
  use mpi
  
  implicit none

  ! Parameters:
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: iostd = 16
  integer, parameter :: root  = 0
    !! ID of the root node

  real(kind = dp), parameter ::          pi = 3.141592653589793_dp
  real(kind = dp), parameter ::       sq4pi = 3.544907701811032_dp
  real(kind = dp), parameter :: evToHartree = 0.03674932538878_dp
  
  complex(kind = dp), parameter :: ii = cmplx(0.0_dp, 1.0_dp, kind = dp)
  
  character(len = 6), parameter :: output = 'output'


  ! Global variables not passed as arguments:
  integer :: ierr
    !! Error returned by MPI
  integer :: iGkStart_poolPC, iGkEnd_poolPC, iGkStart_poolSD, iGkEnd_poolSD
    !! Start and end G+k vector for process in pool
  integer :: iGStart_pool, iGEnd_pool
    !! Start and end G-vector for process in pool
  integer :: ikStart_pool, ikEnd_pool
    !! Start and end k-points in pool
  integer :: indexInPool
    !! Process index within pool
  integer :: intraPoolComm = 0
    !! Intra-pool communicator
  integer :: interPoolComm = 0
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
  real(kind=dp) :: realLattVec(3,3)
    !! Real space lattice vectors
  real(kind=dp) :: recipLattVec(3,3)
    !! Reciprocal lattice vectors
  real(kind = dp) :: t1, t2
    !! For timing different processes


  integer :: maxGIndexGlobal
    !! Maximum G-vector index among all \(G+k\)
    !! and processors for PC and SD
  integer, allocatable :: mill_local(:,:)
    !! Local Miller indices
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors
  integer :: nGkVecsLocalPC, nGkVecsLocalSD
    !! Local number of G+k vectors on this processor
  integer :: nKPoints
    !! Total number of k-points
  integer :: nSpins
    !! Max number of spins for both systems
  integer :: nSpinsPC
    !! Number of spins for PC system
  integer :: nSpinsSD
    !! Number of spins for SD system

  complex(kind=dp), allocatable :: betaPC(:,:)
    !! PC projectors
  complex(kind=dp), allocatable :: betaSD(:,:)
    !! SD projectors
  complex(kind=dp), allocatable :: Ufi(:,:,:,:)
    !! All-electron overlap
  

  character(len = 200) :: exportDirSD, exportDirPC, VfisOutput
  character(len = 300) :: input, inputPC, textDum, elementsPath
  character(len = 320) :: mkdir
  
  integer :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ik, ki, kf, ig, ibi, ibf
  integer :: JMAX, iTypes, iPn
  integer :: nIonsSD, nIonsPC, nProjsPC, numOfTypesPC
  integer :: numOfTypes, nBands, nProjsSD
  integer :: numOfUsedGvecsPP, ios, npwNi, npwNf, npwMi, npwMf
  integer :: gx, gy, gz, nGvsI, nGvsF, nGi, nGf
  integer :: np, nI, nF, nPP, ind2
  integer :: i, j, n1, n2, n3, n4, n, id, npw
  
  integer, allocatable :: counts(:), displmnt(:), nPWsI(:), nPWsF(:)
  
  real(kind = dp) t0, tf
  
  real(kind = dp) :: omega, threej
  
  real(kind = dp), allocatable :: eigvI(:), eigvF(:), gvecs(:,:), posIonSD(:,:), posIonPC(:,:)
  real(kind = dp), allocatable :: wk(:), xk(:,:)
  real(kind = dp), allocatable :: DE(:,:,:), absVfi2(:,:,:)
  
  complex(kind = dp), allocatable :: wfcPC(:,:), wfcSD(:,:), paw_SDKKPC(:,:), paw_id(:,:)
  complex(kind = dp), allocatable :: pawKPC(:,:,:), pawSDK(:,:,:), pawPsiPC(:,:), pawSDPhi(:,:)
  complex(kind = dp), allocatable :: cProjPC(:,:), cProjSD(:,:)
  complex(kind = dp), allocatable :: paw_PsiPC(:,:), paw_SDPhi(:,:)
  complex(kind = dp), allocatable :: cProjBetaPCPsiSD(:,:)
  complex(kind = dp), allocatable :: cProjBetaSDPhiPC(:,:)
  
  integer, allocatable :: TYPNISD(:), TYPNIPC(:), igvs(:,:,:), pwGvecs(:,:), iqs(:)
  integer, allocatable :: npwsSD(:), pwGs(:,:), nIs(:,:), nFs(:,:), ngs(:,:)
  integer, allocatable :: npwsPC(:)
  real(kind = dp), allocatable :: wkPC(:), xkPC(:,:)
  
  type :: atom
    integer :: numOfAtoms, lMax, lmMax, nMax, iRc
    integer, allocatable :: lps(:)
    real(kind = dp), allocatable :: r(:), rab(:), wae(:,:), wps(:,:), F(:,:), F1(:,:,:), F2(:,:,:), bes_J_qr(:,:)
  end type atom
  
  TYPE(atom), allocatable :: atoms(:), atomsPC(:)
  
  type :: vec
    integer :: ind
    integer, allocatable :: igN(:), igM(:)
  end type vec
  
  TYPE(vec), allocatable :: vecs(:), newVecs(:)
  
  real(kind = dp) :: eBin
  complex(kind = dp) :: paw, pseudo1, pseudo2, paw2
  
  logical :: gamma_only, master, calculateVfis, coulomb, tmes_file_exists
  
  NAMELIST /TME_Input/ exportDirSD, exportDirPC, elementsPath, &
                       iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, &
                       ki, kf, calculateVfis, VfisOutput, eBin
  
  
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
  subroutine readInput(maxGIndexGlobal, nKPoints, nGVecsGlobal, realLattVec, recipLattVec)
    
    implicit none

    ! Output variables:
    integer, intent(out) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC and SD
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors
    integer, intent(out) :: nKPoints
      !! Total number of k-points

    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors

    ! Local variables:    
    integer :: maxGIndexGlobalPC
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC
    integer :: maxGIndexGlobalSD
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for SD

    logical :: file_exists
    
    
    if(ionode) then
      !> Check if file output exists. If it does, delete it.
      inquire(file = output, exist = file_exists)
      if ( file_exists ) then
        open (unit = 11, file = output, status = "old")
        close(unit = 11, status = "delete")
      endif
    
      open(iostd, file = output)
        !! Open new output file.

    endif

    if(ionode) then
    
      call initialize()
    
      READ(5, TME_Input, iostat = ios)
    
      call checkInitialization()

    endif

    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(calculateVfis, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(eBin, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(exportDirSD, len(exportDirSD), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirPC, len(exportDirPC), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(elementsPath, len(elementsPath), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(VfisOutput, len(exportDirSD), MPI_CHARACTER, root, worldComm, ierr)

    call readInputPC(nKPoints, maxGIndexGlobalPC)
    call readInputSD(nKPoints, maxGIndexGlobalSD, nGVecsGlobal, realLattVec, recipLattVec)
    
    if(ionode) then

      maxGIndexGlobal = max(maxGIndexGlobalPC, maxGIndexGlobalSD)

      if(maxGIndexGlobal > nGVecsGlobal) call exitError('readInput', &
          'Trying to reference G vecs outside of max grid size. Try switching which grid is read.', 1)

    endif

    call MPI_BCAST(maxGIndexGlobal, 1, MPI_INTEGER, root, worldComm, ierr)

    nSpins = max(nSpinsPC,nSpinsSD)
    
    return
    
  end subroutine readInput
  
!----------------------------------------------------------------------------
  subroutine initialize()
    
    implicit none

    ! Output variables:
    

    exportDirSD = ''
    exportDirPC = ''
    elementsPath = ''
    VfisOutput = ''
    
    ki = -1
    kf = -1
    
    eBin = -1.0_dp
    
    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    
    calculateVfis = .false.
    
    return
    
  end subroutine initialize
  
!----------------------------------------------------------------------------
  subroutine checkInitialization()
    
    implicit none
    
    logical :: file_exists, abortExecution
    
    abortExecution = .false.
    
    write(iostd, '(" Inputs : ")')
    
    if ( trim(exportDirSD) == '' ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""exportDirSD"" is not defined!")')
      write(iostd, '(" usage : exportDirSD = ''./Export/''")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    else
      input = trim(trim(exportDirSD)//'/input')
      
      inquire(file =trim(input), exist = file_exists)
      
      if ( file_exists .eqv. .false. ) then
        write(iostd, '(" File : ", a, " , does not exist!")') trim(input)
        write(iostd, '(" Please make sure that folder : ", a, " has been created successfully !")') trim(exportDirSD)
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
      endif
    endif
    
    write(iostd, '("exportDirSD = ''", a, "''")') trim(exportDirSD)
    
    if ( trim(exportDirPC) == '' ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""exportDirPC"" is not defined!")')
      write(iostd, '(" usage : exportDirPC = ''./Export/''")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    else
      inputPC = trim(trim(exportDirPC)//'/input')
      
      inquire(file =trim(inputPC), exist = file_exists)
      
      if ( file_exists .eqv. .false. ) then
        write(iostd, '(" File : ", a, " , does not exist!")') trim(inputPC)
        write(iostd, '(" Please make sure that folder : ", a, " has been created successfully !")') trim(exportDirPC)
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
      endif
    endif
    
    write(iostd, '("exportDirPC = ''", a, "''")') trim(exportDirPC)
    
    if ( trim(elementsPath) == '' ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""elementsPath"" is not defined!")')
      write(iostd, '(" usage : elementsPath = ''./''")')
      write(iostd, '(" The homepath will be used as elementsPath.")')
      elementsPath = './TMEs'
    endif
    inquire(file= trim(elementsPath), exist = file_exists)
      !! @note 
      !!   `inquire` only works for directories using the `gfortran` compiler, not
      !!   `ifort`. For portability, I updated the `inquire` check for `exportDirSD`
      !!   and `exportDirPC` to check for the `input` file in those directories. 
      !!   that test will work for both compilers. However, I can't change the `inquire`
      !!   test for `elementsPath` in the same way because `elementsPath` is expected
      !!   to be empty. With `ifort`, the `mkdir` command will always be run, but I
      !!   don't think that will cause any issues. 
      !! @endnote

    if ( .not.file_exists ) then
      write(mkDir, '("mkdir -p ", a)') trim(elementsPath) 
      call system(mkDir)
    endif
    
    write(iostd, '("elementsPath = ''", a, "''")') trim(elementsPath)
    
    if ( iBandIinit < 0 ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandIinit"" is not defined!")')
      write(iostd, '(" usage : iBandIinit = 10")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    endif
    
    write(iostd, '("iBandIinit = ", i4)') iBandIinit
    
    if ( iBandIfinal < 0 ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandIfinal"" is not defined!")')
      write(iostd, '(" usage : iBandIfinal = 20")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    endif
    
    write(iostd, '("iBandIfinal = ", i4)') iBandIfinal
    
    if ( iBandFinit < 0 ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandFinit"" is not defined!")')
      write(iostd, '(" usage : iBandFinit = 9")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    endif
    
    write(iostd, '("iBandFinit = ", i4)') iBandFinit
    
    if ( iBandFfinal < 0 ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandFfinal"" is not defined!")')
      write(iostd, '(" usage : iBandFfinal = 9")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    endif
    
    write(iostd, '("iBandFfinal = ", i4)') iBandFfinal
    
    if ( ( calculateVfis ) .and. ( iBandFinit /= iBandFfinal ) ) then
      write(iostd, *)
      write(iostd, '(" Vfis can be calculated only if the final state is one and only one!")')
      write(iostd, '(" ''iBandFInit'' = ", i10)') iBandFinit
      write(iostd, '(" ''iBandFfinal'' = ", i10)') iBandFfinal
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
    endif
    
    write(iostd, '("calculateVfis = ", l )') calculateVfis
    
    if ( trim(VfisOutput) == '' ) then
      write(iostd, *)
      write(iostd, '(" Variable : ""VfisOutput"" is not defined!")')
      write(iostd, '(" usage : VfisOutput = ''VfisVsE''")')
      write(iostd, '(" The default value ''VfisOutput'' will be used.")')
      VfisOutput = 'VfisVsE'
    endif
    
    write(iostd, '("VfisOutput = ''", a, "''")') trim(VfisOutput)
    
    if ( eBin < 0.0_dp ) then
      eBin = 0.01_dp ! eV
      write(iostd,'(" Variable : ""eBin"" is not defined!")')
      write(iostd,'(" usage : eBin = 0.01")')
      write(iostd,'(" A default value of 0.01 eV will be used !")')
    endif
    
    write(iostd, '("eBin = ", f8.4, " (eV)")') eBin
    
    eBin = eBin*evToHartree
    
    if ( abortExecution ) then
      write(iostd, '(" Program stops!")')
      stop
    endif
    
    flush(iostd)
    
    return
    
  end subroutine checkInitialization
  
!----------------------------------------------------------------------------
  subroutine readInputPC(nKPoints, maxGIndexGlobalPC)
    
    implicit none

    ! Outpu variables:
    integer, intent(out) :: nKPoints
      !! Total number of k-points

    ! Output variables:
    integer, intent(out) :: maxGIndexGlobalPC
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC
    
    ! Local variables:
    integer :: i, j, l, ind, ik, iDum, iType, ni, irc
    
    real(kind = dp) :: t1, t2 
    
    character(len = 300) :: textDum
    
    logical :: file_exists
    
    if(ionode) then
      call cpu_time(t1)
    
      write(iostd,*)
      write(iostd, '(" Reading perfect crystal inputs.")')
      write(iostd, *)
    
      inputPC = trim(trim(exportDirPC)//'/input')
    
      inquire(file =trim(inputPC), exist = file_exists)
    
      if ( file_exists .eqv. .false. ) then
        write(iostd, '(" File : ", a, " , does not exist!")') trim(inputPC)
        write(iostd, '(" Please make sure that folder : ", a, " has been created successfully !")') trim(exportDirPC)
        write(iostd, '(" Program stops!")')
        flush(iostd)
      endif
    
      open(50, file=trim(inputPC), status = 'old')
    
      read(50, '(a)') textDum
      read(50, * ) 
    
      read(50, '(a)') textDum
      read(50, '(i10)') nSpinsPC
    
      read(50, '(a)') textDum
      read(50, '(i10)') nKPoints

    endif

    call MPI_BCAST(nSpinsPC, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(npwsPC(nKPoints), wkPC(nKPoints), xkPC(3,nKPoints))
    
    if(ionode) then

      read(50, '(a)') textDum
    
      do ik = 1, nKPoints
      
        read(50, '(2i10,4ES24.15E3)') iDum, npwsPC(ik), wkPC(ik), xkPC(1:3,ik)
      
      enddo

    endif
    
    call MPI_BCAST(npwsPC, nKPoints, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(wkPC, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(xkPC, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if(ionode) then

      read(50, '(a)') textDum
      read(50, * ) ! nGVecsGlobal
    
      read(50, '(a)') textDum
      read(50, '(i10)') maxGIndexGlobalPC
    
      read(50, '(a)') textDum     
      read(50, * ) ! fftxMin, fftxMax, fftyMin, fftyMax, fftzMin, fftzMax
    
      read(50, '(a)') textDum
      read(50,  * )
      read(50,  * )
      read(50,  * )
    
      read(50, '(a)') textDum
      read(50, * )
      read(50, * )
      read(50, * )
    
      read(50, '(a)') textDum
      read(50, '(i10)') nIonsPC
    
      read(50, '(a)') textDum
      read(50, '(i10)') numOfTypesPC
    
    endif

    call MPI_BCAST(nIonsPC, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(numOfTypesPC, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(posIonPC(3,nIonsPC), TYPNIPC(nIonsPC))


    if(ionode) then
    
      read(50, '(a)') textDum
      do ni = 1, nIonsPC
        read(50,'(i10, 3ES24.15E3)') TYPNIPC(ni), (posIonPC(j,ni) , j = 1,3)
      enddo
    
      read(50, '(a)') textDum
      read(50, * )

    endif

    call MPI_BCAST(TYPNIPC,  size(TYPNIPC),  MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(posIonPC, size(posIonPC), MPI_DOUBLE_PRECISION,root,worldComm,ierr)
    
    allocate(atomsPC(numOfTypesPC))

    nProjsPC = 0
    do iType = 1, numOfTypesPC
      
      if(ionode) then

        read(50, '(a)') textDum
        read(50, *) 
      
        read(50, '(a)') textDum
        read(50, '(i10)') atomsPC(iType)%numOfAtoms

        read(50, '(a)') textDum
        read(50, '(i10)') atomsPC(iType)%lMax              ! number of projectors

      endif

      call MPI_BCAST(atomsPC(iType)%numOfAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%lMax, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(atomsPC(iType)%lps(atomsPC(iType)%lMax))
      
      if(ionode) then

        read(50, '(a)') textDum

        do i = 1, atomsPC(iType)%lMax 

          read(50, '(2i10)') l, ind
          atomsPC(iType)%lps(ind) = l

        enddo

      endif

      call MPI_BCAST(atomsPC(iType)%lps, atomsPC(iType)%lMax, MPI_INTEGER, root, worldComm, ierr)

      if(ionode) then
      
        read(50, '(a)') textDum
        read(50, '(i10)') atomsPC(iType)%lmMax
      
        read(50, '(a)') textDum
        read(50, '(2i10)') atomsPC(iType)%nMax, atomsPC(iType)%iRc

      endif
    
      call MPI_BCAST(atomsPC(iType)%lmMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%nMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%iRc, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(atomsPC(iType)%r(atomsPC(iType)%nMax))
      
      if(ionode) then
      
        allocate(atomsPC(iType)%rab(atomsPC(iType)%nMax))

        read(50, '(a)') textDum

        do i = 1, atomsPC(iType)%nMax
          read(50, '(2ES24.15E3)') atomsPC(iType)%r(i), atomsPC(iType)%rab(i)
        enddo
       
        allocate(atomsPC(iType)%wae(atomsPC(iType)%nMax, atomsPC(iType)%lMax))
        allocate(atomsPC(iType)%wps(atomsPC(iType)%nMax, atomsPC(iType)%lMax))
      
        read(50, '(a)') textDum
        do j = 1, atomsPC(iType)%lMax
          do i = 1, atomsPC(iType)%nMax
            read(50, '(2ES24.15E3)') atomsPC(iType)%wae(i, j), atomsPC(iType)%wps(i, j) 
          enddo
        enddo
        
      endif

      allocate(atomsPC(iType)%F(atomsPC(iType)%iRc, atomsPC(iType)%lMax))
      allocate(atomsPC(iType)%F1(atomsPC(iType)%iRc, atomsPC(iType)%lMax, atomsPC(iType)%lMax))
      allocate(atomsPC(iType)%F2(atomsPC(iType)%iRc, atomsPC(iType)%lMax, atomsPC(iType)%lMax))
      
      if(ionode) then

        atomsPC(iType)%F = 0.0_dp
        atomsPC(iType)%F1 = 0.0_dp
        atomsPC(iType)%F2 = 0.0_dp
      
        do j = 1, atomsPC(iType)%lMax
        
          irc = atomsPC(iType)%iRc
          atomsPC(iType)%F(1:irc,j)=(atomsPC(iType)%wae(1:irc,j)-atomsPC(iType)%wps(1:irc,j))* &
                atomsPC(iType)%r(1:irc)*atomsPC(iType)%rab(1:irc)
        
          do i = 1, atomsPC(iType)%lMax
            atomsPC(iType)%F1(1:irc,i,j) = ( atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wae(1:irc,j) - &
                                            atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wps(1:irc,j))*atomsPC(iType)%rab(1:irc)
          
            atomsPC(iType)%F2(1:irc,i,j) = ( atomsPC(iType)%wae(1:irc,i)*atomsPC(iType)%wae(1:irc,j) - &
                                             atomsPC(iType)%wae(1:irc,i)*atomsPC(iType)%wps(1:irc,j) - &
                                             atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wae(1:irc,j) + &
                                             atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wps(1:irc,j))*atomsPC(iType)%rab(1:irc)

          enddo
        enddo

        deallocate(atomsPC(iType)%wae)
        deallocate(atomsPC(iType)%wps)
        deallocate(atomsPC(iType)%rab)
      
        nProjsPC = nProjsPC + atomsPC(iType)%numOfAtoms*atomsPC(iType)%lmMax

      endif

      call MPI_BCAST(atomsPC(iType)%r, size(atomsPC(iType)%r), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%F, size(atomsPC(iType)%F), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%F1, size(atomsPC(iType)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%F2, size(atomsPC(iType)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      
    enddo
  
    call MPI_BCAST(nProjsPC, 1, MPI_INTEGER, root, worldComm, ierr)

    if(ionode) then
    
      close(50)
    
      JMAX = 0
      do iType = 1, numOfTypesPC
        do i = 1, atomsPC(iType)%lMax
          if ( atomsPC(iType)%lps(i) > JMAX ) JMAX = atomsPC(iType)%lps(i)
        enddo
      enddo
    
      JMAX = 2*JMAX + 1

    endif

    call MPI_BCAST(JMAX, 1, MPI_INTEGER, root, worldComm, ierr)
    
    do iType = 1, numOfTypesPC

      allocate(atomsPC(iType)%bes_J_qr(0:JMAX, atomsPC(iType)%iRc))
      atomsPC(iType)%bes_J_qr(:,:) = 0.0_dp
      
    enddo

    if(ionode) then
      call cpu_time(t2)
      write(iostd, '(" Reading input files done in:                ", f10.2, " secs.")') t2-t1
      write(iostd, *)
      flush(iostd)

    endif
    
    return
    
  end subroutine readInputPC
  
!----------------------------------------------------------------------------
  subroutine readInputSD(nKPoints, maxGIndexGlobalSD, nGVecsGlobal, realLattVec, recipLattVec)
    !
    implicit none
    
    ! Input variables:
    integer, intent(in) :: nKPoints
      !! Total number of k-points

    ! Output variables:
    integer, intent(out) :: maxGIndexGlobalSD
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for SD
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors

    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors
    
    ! Local variables:
    integer :: nKpts
      !! Number of k-points read from SD file

    integer :: i, j, l, ind, ik, iDum, iType, ni, irc
    
    real(kind = dp) :: t1, t2
    
    character(len = 300) :: textDum
    
    logical :: file_exists
    

    if(ionode) then
      call cpu_time(t1)
    
      write(iostd, *)
      write(iostd, '(" Reading solid defect inputs.")')
      write(iostd, *)
    
      input = trim(trim(exportDirSD)//'/input')
    
      inquire(file = trim(input), exist = file_exists)
    
      if ( file_exists .eqv. .false. ) then
        write(iostd, '(" File : ", a, " , does not exist!")') trim(input)
        write(iostd, '(" Please make sure that folder : ", a, " has been created successfully !")') trim(exportDirSD)
        write(iostd, '(" Program stops!")')
        flush(iostd)
      endif
    
      open(50, file=trim(input), status = 'old')
    
      read(50, '(a)') textDum
      read(50, '(ES24.15E3)' ) omega

    endif

    call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    
    if(ionode) then

      read(50, '(a)') textDum
      read(50, '(i10)') nSpinsSD

      read(50, '(a)') textDum
      read(50, '(i10)') nKpts

      if(nKpts /= nKPoints) call exitError('readInputsSD', 'Number of k-points in systems must match', 1)
    
      read(50, '(a)') textDum

    endif

    call MPI_BCAST(nSpinsSD, 1, MPI_INTEGER, root, worldComm, ierr)

    
    allocate(npwsSD(nKPoints), wk(nKPoints), xk(3,nKPoints))

    if(ionode) then
    
      do ik = 1, nKPoints
      
        read(50, '(2i10,4ES24.15E3)') iDum, npwsSD(ik), wk(ik), xk(1:3,ik)
      
      enddo

    endif
    
    call MPI_BCAST(npwsSD, nKPoints, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(wk, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(xk, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if(ionode) then
    
      read(50, '(a)') textDum
      read(50, '(i10)') nGVecsGlobal
    
      read(50, '(a)') textDum
      read(50, '(i10)') maxGIndexGlobalSD
    
      read(50, '(a)') textDum     
      read(50,*) 
      !read(50, '(6i10)') fftxMin, fftxMax, fftyMin, fftyMax, fftzMin, fftzMax
    
      read(50, '(a)') textDum
      read(50, '(a5, 3ES24.15E3)') textDum, realLattVec(1:3,1)
      read(50, '(a5, 3ES24.15E3)') textDum, realLattVec(1:3,2)
      read(50, '(a5, 3ES24.15E3)') textDum, realLattVec(1:3,3)
    
      read(50, '(a)') textDum
      read(50, '(a5, 3ES24.15E3)') textDum, recipLattVec(1:3,1)
      read(50, '(a5, 3ES24.15E3)') textDum, recipLattVec(1:3,2)
      read(50, '(a5, 3ES24.15E3)') textDum, recipLattVec(1:3,3)
    
      read(50, '(a)') textDum
      read(50, '(i10)') nIonsSD
    
      read(50, '(a)') textDum
      read(50, '(i10)') numOfTypes

    endif

    call MPI_BCAST(nGVecsGlobal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(recipLattVec, size(recipLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(nIonsSD, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(numOfTypes, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(posIonSD(3,nIonsSD), TYPNISD(nIonsSD))

    if(ionode) then
    
      read(50, '(a)') textDum
      do ni = 1, nIonsSD
        read(50,'(i10, 3ES24.15E3)') TYPNISD(ni), (posIonSD(j,ni), j = 1,3)
      enddo

    endif
  
    call MPI_BCAST(TYPNISD, size(TYPNISD), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(posIonSD, size(posIonSD), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if(ionode) then
    
      read(50, '(a)') textDum
      read(50, '(i10)') nBands

    endif

    call MPI_BCAST(nBands, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(atoms(numOfTypes))
    
    nProjsSD = 0
    do iType = 1, numOfTypes
      
      if(ionode) then

        read(50, '(a)') textDum
        read(50, *) 
      
        read(50, '(a)') textDum
        read(50, '(i10)') atoms(iType)%numOfAtoms
      
        read(50, '(a)') textDum
        read(50, '(i10)') atoms(iType)%lMax              ! number of projectors

      endif

      call MPI_BCAST(atoms(iType)%numOfAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%lMax, 1, MPI_INTEGER, root, worldComm, ierr)
      
      allocate(atoms(iType)%lps(atoms(iType)%lMax))
      
      if(ionode) then

        read(50, '(a)') textDum
        do i = 1, atoms(iType)%lMax 
          read(50, '(2i10)') l, ind
          atoms(iType)%lps(ind) = l
        enddo

      endif

      call MPI_BCAST(atoms(iType)%lps, size(atoms(iType)%lps), MPI_INTEGER, root, worldComm, ierr)

      if(ionode) then
      
        read(50, '(a)') textDum
        read(50, '(i10)') atoms(iType)%lmMax
      
        read(50, '(a)') textDum
        read(50, '(2i10)') atoms(iType)%nMax, atoms(iType)%iRc

      endif
    
      call MPI_BCAST(atoms(iType)%lmMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%nMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%iRc, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(atoms(iType)%r(atoms(iType)%nMax))

      if(ionode) then
      
        allocate(atoms(iType)%rab(atoms(iType)%nMax))
      
        read(50, '(a)') textDum
        do i = 1, atoms(iType)%nMax
          read(50, '(2ES24.15E3)') atoms(iType)%r(i), atoms(iType)%rab(i)
        enddo
       
        allocate(atoms(iType)%wae(atoms(iType)%nMax, atoms(iType)%lMax))
        allocate(atoms(iType)%wps(atoms(iType)%nMax, atoms(iType)%lMax))
      
        read(50, '(a)') textDum
        do j = 1, atoms(iType)%lMax
          do i = 1, atoms(iType)%nMax
            read(50, '(2ES24.15E3)') atoms(iType)%wae(i, j), atoms(iType)%wps(i, j) 
          enddo
        enddo

      endif

        
      allocate(atoms(iType)%F(atoms(iType)%iRc, atoms(iType)%lMax))
      allocate(atoms(iType)%F1(atoms(iType)%iRc, atoms(iType)%lMax, atoms(iType)%lMax))
      allocate(atoms(iType)%F2(atoms(iType)%iRc, atoms(iType)%lMax, atoms(iType)%lMax))
      
      if(ionode) then

        atoms(iType)%F = 0.0_dp
        atoms(iType)%F1 = 0.0_dp
        atoms(iType)%F2 = 0.0_dp
      
        do j = 1, atoms(iType)%lMax
          
          irc = atoms(iType)%iRc
          atoms(iType)%F(1:irc,j)=(atoms(iType)%wae(1:irc,j)-atoms(iType)%wps(1:irc,j))*atoms(iType)%r(1:irc) * &
              atoms(iType)%rab(1:irc)
          
          do i = 1, atoms(iType)%lMax
                  
            atoms(iType)%F1(1:irc,i,j) = ( atoms(iType)%wae(1:irc,i)*atoms(iType)%wps(1:irc,j) - &
                                           atoms(iType)%wps(1:irc,i)*atoms(iType)%wps(1:irc,j))*atoms(iType)%rab(1:irc)
          
            atoms(iType)%F2(1:irc,i,j) = ( atoms(iType)%wae(1:irc,i)*atoms(iType)%wae(1:irc,j) - &
                                           atoms(iType)%wae(1:irc,i)*atoms(iType)%wps(1:irc,j) - &
                                           atoms(iType)%wps(1:irc,i)*atoms(iType)%wae(1:irc,j) + &
                                           atoms(iType)%wps(1:irc,i)*atoms(iType)%wps(1:irc,j))*atoms(iType)%rab(1:irc)
          enddo
        enddo
      
        nProjsSD = nProjsSD + atoms(iType)%numOfAtoms*atoms(iType)%lmMax
      
        deallocate(atoms(iType)%wae)
        deallocate(atoms(iType)%wps)
        deallocate(atoms(iType)%rab)

      endif

      call MPI_BCAST(atoms(iType)%r, size(atoms(iType)%r), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%F, size(atoms(iType)%F), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%F1, size(atoms(iType)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%F2, size(atoms(iType)%F2), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      
    enddo

    call MPI_BCAST(nProjsSD, 1, MPI_INTEGER, root, worldComm, ierr)
    
    if(ionode) then
    
      close(50)

      JMAX = 0
      do iType = 1, numOfTypes
        do i = 1, atoms(iType)%lMax
          if ( atoms(iType)%lps(i) > JMAX ) JMAX = atoms(iType)%lps(i)
        enddo
      enddo
    
      JMAX = 2*JMAX + 1

    endif

    call MPI_BCAST(JMAX, 1, MPI_INTEGER, root, worldComm, ierr)
    
    do iType = 1, numOfTypes

      allocate(atoms(iType)%bes_J_qr( 0:JMAX, atoms(iType)%iRc))
      atoms(iType)%bes_J_qr(:,:) = 0.0_dp
      
    enddo
    
    if(ionode) then

      call cpu_time(t2)
      write(iostd, '(" Reading solid defect inputs done in:                ", f10.2, " secs.")') t2-t1
      write(iostd, *)
      flush(iostd)

    endif
    
    return
    
  end subroutine readInputSD

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
  subroutine checkIfCalculated(ikGlobal, isp, tmes_file_exists)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Output variables:
    logical, intent(out) :: tmes_file_exists
      !! If the current overlap file exists

    ! Local variables:
    character(len=300) :: ikC, ispC
      !! Character indices

    
    call int2str(ikGlobal, ikC)
    call int2str(isp, ispC)
    
    inquire(file=trim(elementsPath)//"/allElecOverlap."//trim(ispC)//"."//trim(ikC), exist = tmes_file_exists)
    
    return
    
  end subroutine checkIfCalculated

!----------------------------------------------------------------------------
  subroutine getFullPWGrid(iGStart_pool, iGEnd_pool, nGVecsLocal, nGVecsGlobal, mill_local)
    !! Read full PW grid from mgrid file
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iGStart_pool, iGEnd_pool
      !! Start and end G-vectors for each process in pool
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors

    ! Output variables:
    integer, intent(out) :: mill_local(3,nGVecsLocal)
      !! Integer coefficients for G-vectors
    
    ! Local variables:
    integer :: gVecMillerIndicesGlobal(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer :: iDum
      !! Ignore dummy integer
    integer :: ig
      !! Loop index
    

    if(ionode) then

      open(72, file=trim(exportDirSD)//"/mgrid")
        !! Read full G-vector grid from defect folder.
        !! This assumes that the grids are the same.
    
      read(72, * )
      read(72, * )
    
      do ig = 1, nGVecsGlobal

        read(72, '(4i10)') iDum, gVecMillerIndicesGlobal(1:3,ig)

      enddo
    
      close(72)

    endif

    call MPI_BCAST(gVecMillerIndicesGlobal, size(gVecMillerIndicesGlobal), MPI_INTEGER, root, worldComm, ierr)

    mill_local(:,:) = gVecMillerIndicesGlobal(:,iGStart_pool:iGEnd_pool)
    
    return
    
  end subroutine getFullPWGrid

!----------------------------------------------------------------------------
  subroutine readProjectors(crystalType, iGkStart_pool, ikGlobal, nGkVecsLocal, nProjs, npws, betaLocalPWs)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iGkStart_pool
      !! Start G+k index for this process in pool
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor
    integer, intent(in) :: nProjs
      !! Number of projectors
    integer, intent(in) :: npws
      !! Number of G+k vectors less than
      !! the cutoff at each k-point

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD) for projectors

    ! Output variables:
    complex(kind=dp), intent(out) :: betaLocalPWs(nGkVecsLocal,nProjs)
      !! Projector of `projCrystalType` with all projectors
      !! and only local PWs/G+k vectors
    
    ! Local variables:
    integer :: displacement(nProcPerPool)
      !! Offset from beginning of array for
      !! scattering coefficients to each process
    integer :: endingProj(nProcPerPool)
      !! End projector for each process in pool
    integer :: iprStart_pool, iprEnd_pool
      !! Start and end projector for process in pool
    integer :: nProjsLocal
      !! Number of projectors read by this process
    integer :: sendCount(nProcPerPool)
      !! Number of items to send to each process
      !! in the pool
    integer :: suffixLength
      !! Length of split-file suffix
    integer :: ipr, igk, iproc
      !! Loop indices

    complex(kind=dp), allocatable :: betaLocalProjs(:,:)
      !! Projector of `projCrystalType` with all PWs/G+k 
      !! vectors and only local projectors
    
    character(len=300) :: splitCommand
      !! Command for splitting `projectors.ik` file
    character(len=300) :: fNameExport, splitFilePrefix
      !! File names
    character(len=300) :: ikC, ngkC, sLC, iprC
      !! Character index
    

    call int2str(ikGlobal, ikC)
    
    if(crystalType == 'PC') then
      fNameExport = trim(exportDirPC)//"/projectors."//trim(ikC) 
    else
      fNameExport = trim(exportDirSD)//"/projectors."//trim(ikC)
    endif


    if(indexInPool == 0) then

      call int2str(npws,ngkC)

      splitCommand = 'cat '//trim(fNameExport)//' | tail -n +3 | split -l '//trim(ngkC)//' -d -a '
        !! Split the projectors file into separate files for each `nProj`
        !! so that the length of each file is equal to the number of 
        !! G+k vectors < cutoff at this k-point. Ignore the two header
        !! lines using `tail`. `-d` makes suffix numeric. `-a plus 
        !! `suffixLength` below sets the length of the suffix with leading
        !! zeros as needed.

    endif

    if(nProjs < 10) then
      suffixLength = 1
    else if(nProjs < 100) then
      suffixLength = 2
    else if(nProjs < 1000) then
      suffixLength = 3
    else if(nProjs < 10000) then
      suffixLength = 4
    else if(nProjs < 100000) then
      suffixLength = 5
    endif

    splitFilePrefix = './workingFiles/projectors.'//trim(ikC)//'.split.'

    if(indexInPool == 0) then

      call int2str(suffixLength, sLC)
      splitCommand = trim(splitCommand)//trim(sLC)//' - '//trim(splitFilePrefix)
        !! Finish the split command with the suffix length. `-` stands for
        !! `stdin` and allows the piped output from `tail` to be used as the
        !! input file to be split.

      call execute_command_line('mkdir -p workingFiles')
        !!  Make the `workingFiles` directory
      call execute_command_line(splitCommand)
        !! Split files

    endif


    call distributeItemsInSubgroups(indexInPool, nProjs, nProcPerPool, nProcPerPool, nProcPerPool, iprStart_pool, iprEnd_pool, nProjsLocal)

    allocate(betaLocalProjs(npws,iprStart_pool:iprEnd_pool))

    betaLocalProjs(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)
    betaLocalPWs(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)

    call MPI_Barrier(intraPoolComm, ierr)


    !> Have each process open the split files for just their
    !> projectors and then delete the files
    do ipr = iprStart_pool, iprEnd_pool
      
      call int2strLeadZero(ipr-1, suffixLength, iprC)

      open(72,file=trim(splitFilePrefix)//trim(iprC))

      do igk = 1, npws

        read(72, '(2ES24.15E3)') betaLocalProjs(igk,ipr)

      enddo

      close(72, status="delete")

    enddo
    
    call MPI_Barrier(intraPoolComm, ierr)

    sendCount = 0
    sendCount(indexInPool+1) = nGkVecsLocal
    call mpiSumIntV(sendCount, intraPoolComm)
      !! * Put the number of G+k vectors on each process
      !!   in a single array per pool

    displacement = 0
    displacement(indexInPool+1) = iGkStart_pool-1
    call mpiSumIntV(displacement, intraPoolComm)
      !! * Put the displacement from the beginning of the array
      !!   for each process in a single array per pool

    endingProj = 0
    endingProj(indexInPool+1) = iprEnd_pool
    call mpiSumIntV(endingProj, intraPoolComm)


    !> Distribute projectors across processors so that PWs are local
    !> instead of projectors
    iproc = 0
    do ipr = 1, nProjs

      if(ipr == endingProj(iproc+1)+1) iproc = iproc + 1

      call MPI_SCATTERV(betaLocalProjs(:,ipr), sendCount, displacement, MPI_COMPLEX, betaLocalPWs(1:nGkVecsLocal,ipr), nGkVecsLocal, &
          MPI_COMPLEX, iproc, intraPoolComm, ierr)

    enddo

    deallocate(betaLocalProjs)
    
  end subroutine readProjectors

!----------------------------------------------------------------------------
  subroutine readWfc(crystalType, iBandinit, iBandfinal, iGkStart_pool, ikGlobal, isp, nGkVecsLocal, npws, wfc)
    !! Read wave function for given `crystalType` from `iBandinit`
    !! to `iBandfinal`
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: iBandinit
      !! Starting band
    integer, intent(in) :: iBandfinal
      !! Ending band
    integer, intent(in) :: iGkStart_pool
      !! Start G+k index for this process in pool
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor
    integer, intent(in) :: npws
      !! Number of G+k vectors less than
      !! the cutoff at this k-point

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    ! Output variables
    complex(kind=dp), intent(out) :: wfc(nGkVecsLocal,iBandinit:iBandfinal)
      !! Wave function coefficients for 
      !! local G-vectors

    ! Local variables:
    integer :: displacement(nProcPerPool)
      !! Offset from beginning of array for
      !! scattering coefficients to each process
    integer :: sendCount(nProcPerPool)
      !! Number of items to send to each process
      !! in the pool
    integer :: ib, igk, iproc
      !! Loop indices

    complex(kind=dp) :: wfcAllPWs(npws)
      !! Wave function read from file
 
    character(len = 300) :: ikC, ispC
      !! Character indices


    if(indexInPool == 0) then
      !! Have the root node in the pool read the wave
      !! function coefficients. Ignore the bands before 
      !! `iBandinit`

      call int2str(ikGlobal, ikC)
      call int2str(isp, ispC)
    
      if(crystalType == 'PC') then
        open(72, file=trim(exportDirPC)//"/wfc."//trim(ispC)//"."//trim(ikC))
      else
        open(72, file=trim(exportDirSD)//"/wfc."//trim(ispC)//"."//trim(ikC))
      endif
    
      read(72, * )
      read(72, * )
        ! Ignore the header
    
      do ib = 1, iBandinit - 1
        do igk = 1, npws
          read(72, *)
        enddo
      enddo

    endif
    
    wfc(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)

    sendCount = 0
    sendCount(indexInPool+1) = nGkVecsLocal
    call mpiSumIntV(sendCount, intraPoolComm)
      !! * Put the number of G+k vectors on each process
      !!   in a single array per pool

    displacement = 0
    displacement(indexInPool+1) = iGkStart_pool-1
    call mpiSumIntV(displacement, intraPoolComm)
      !! * Put the displacement from the beginning of the array
      !!   for each process in a single array per pool
    
    !> Have the root node read the PW coefficients
    !> then broadcast to other processes
    do ib = iBandinit, iBandfinal

      if(indexInPool == 0) then

        do igk = 1, npws
          read(72,'(2ES24.15E3)') wfcAllPWs(igk)
        enddo

      endif

      call MPI_SCATTERV(wfcAllPWs(:), sendCount, displacement, MPI_COMPLEX, wfc(1:nGkVecsLocal,ib), nGkVecsLocal, &
        MPI_COMPLEX, 0, intraPoolComm, ierr)

    enddo

    if(indexInPool == 0) close(72)
    
    return
    
  end subroutine readWfc
  
!----------------------------------------------------------------------------
  subroutine calculatePWsOverlap(ikLocal,isp)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ibi, ibf

    
    Ufi(:,:,ikLocal,isp) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    do ibi = iBandIinit, iBandIfinal 
      
      do ibf = iBandFinit, iBandFfinal

        Ufi(ibf, ibi, ikLocal,isp) = sum(conjg(wfcSD(:,ibf))*wfcPC(:,ibi))
          !! Calculate local overlap

      enddo
      
    enddo
    
    return
    
  end subroutine calculatePWsOverlap
  
!----------------------------------------------------------------------------
  subroutine readProjections(crystalType, iBandinit, iBandfinal, ikGlobal, isp, nProjs, cProj)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandinit
      !! Starting band
    integer, intent(in) :: iBandfinal
      !! Ending band
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nProjs
      !! Number of projectors

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    ! Output variables:
    complex(kind=dp), intent(out) :: cProj(nProjs,iBandinit:iBandfinal)
      !! Projections <beta|wfc>

    ! Local variables:
    integer :: ipr, ib
      !! Loop indices
    
    character(len = 300) :: ikC, ispC
      !! Character indices
    

    call int2str(ikGlobal, ikC)
    call int2str(isp, ispC)
    
    cProj(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
    
    if(indexInPool == 0) then

      ! Open the projections file for the given crystal type
      if(crystalType == 'PC') then
        open(72, file=trim(exportDirPC)//"/projections."//trim(ispC)//"."//trim(ikC))
      else
        open(72, file=trim(exportDirSD)//"/projections."//trim(ispC)//"."//trim(ikC))
      endif
    
      read(72, *)
        ! Ignore the header

      ! Ignore bands before initial band
      do ib = 1, iBandinit-1
        do ipr = 1, nProjs
          read(72,*)
        enddo
      enddo
    
      ! Read the projections
      do ib = iBandinit, iBandfinal
        do ipr = 1, nProjs 

          read(72,'(2ES24.15E3)') cProj(ipr,ib)

        enddo
      enddo
    
      close(72)

    endif

    call MPI_BCAST(cProj, size(cProj), MPI_DOUBLE_COMPLEX, root, intraPoolComm, ierr)
      ! Broadcast entire array to all processes
    
    return
    
  end subroutine readProjections
  
!----------------------------------------------------------------------------
  subroutine calculateCrossProjection(iBandinit, iBandfinal, ikGlobal, nGkVecsLocal1, nGkVecsLocal2, nProjs, beta, wfc, crossProjection)
    !! Calculate the cross projection of one crystal's projectors
    !! on the other crystal's wave function coefficients, distributing
    !! the result to all processors
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandinit
      !! Starting band for crystal wfc comes from
      !! (not `projCrystalType`)
    integer, intent(in) :: iBandfinal
      !! Ending band for crystal wfc comes from
      !! (not `projCrystalType`)
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: nGkVecsLocal1
      !! Local number of G+k vectors on this processor
      !! for projectors
    integer, intent(in) :: nGkVecsLocal2
      !! Local number of G+k vectors on this processor
      !! for wave function
    integer, intent(in) :: nProjs
      !! Number of projectors

    complex(kind=dp) :: beta(nGkVecsLocal1,nProjs)
      !! Projector of one crystal type
    complex(kind=dp), intent(in) :: wfc(nGkVecsLocal2,iBandinit:iBandfinal)
      !! Wave function coefficients for local G-vectors
      !! for other crystal type

    ! Output variables:
    complex(kind=dp), intent(out) :: crossProjection(nProjs,iBandinit:iBandfinal)
      !! Projections <beta|wfc>
    
    ! Local variables:
    integer :: ib, ipr
      !! Loop indices

    complex(kind=dp) :: crossProjectionLocal
      !! Local version of cross projection to
      !! be summed across processors in pool
    
    character(len = 300) :: ikC
    

    crossProjection(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do ib = iBandinit, iBandfinal
      do ipr = 1, nProjs

        crossProjectionLocal = sum(conjg(beta(:,ipr))*wfc(:,ib))

        call MPI_ALLREDUCE(crossProjectionLocal, crossProjection(ipr,ib), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, intraPoolComm, ierr)

      enddo
    enddo
    
    return
    
  end subroutine calculateCrossProjection
  
!----------------------------------------------------------------------------
  subroutine pawCorrectionWfc(nAtoms, iType, nProjs, cProjI, cProjF, atoms, pawWfc)
    ! calculates the augmentation part of the transition matrix element
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: nProjs
      !! First index for `cProjI` and `cProjF`

    complex(kind = dp) :: cProjI(nProjs,iBandIinit:iBandIfinal)
      !! Initial-system (PC) projection
    complex(kind = dp) :: cProjF(nProjs,iBandFinit:iBandFfinal)
      !! Final-system (SD) projection

    type(atom), intent(in) :: atoms(nAtoms)
      !! Structure to hold details for each atom

    ! Output variables:
    complex(kind=dp), intent(out) :: pawWfc(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal)
      !! Augmentation part of the transition matrix element

    ! Local variables:
    integer :: ia, ispin
      !! Loop index

    complex(kind = dp) :: cProjIe
      !! Single element of initial-system (PC)
      !! projection
    complex(kind = dp) :: cProjFe
      !! Single element of final-system (SD)
      !! projection

    integer :: ibi, ibf
    integer :: LL, LLP, LMBASE, LM, LMP
    integer :: L, M, LP, MP, iT
    real(kind = dp) :: atomicOverlap
    
    
    pawWfc(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    LMBASE = 0
    
    do ia = 1, nAtoms
      
      iT = iType(ia)

      LM = 0
      DO LL = 1, atoms(iT)%lMax
        L = atoms(iT)%LPS(LL)
        DO M = -L, L
          LM = LM + 1 !1st index for CPROJ
          
          LMP = 0
          DO LLP = 1, atoms(iT)%lMax
            LP = atoms(iT)%LPS(LLP)
            DO MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) then 
                atomicOverlap = sum(atoms(iT)%F1(:,LL, LLP))
                
                do ibi = iBandIinit, iBandIfinal
                  cProjIe = cProjI(LMP + LMBASE, ibi)
                  
                  do ibf = iBandFinit, iBandFfinal
                    cProjFe = conjg(cProjF(LM + LMBASE, ibf))
                    
                    pawWfc(ibf, ibi) = pawWfc(ibf, ibi) + cProjFe*atomicOverlap*cProjIe
                    
                  enddo
                  
                enddo
                
              endif
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      LMBASE = LMBASE + atoms(iT)%lmMax
    ENDDO
    
    return
    
  end subroutine pawCorrectionWfc

!----------------------------------------------------------------------------
  subroutine pawCorrectionK(crystalType, nAtoms, iType, numOfTypes, atomPositions, atoms, atomsSD, pawK)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: numOfTypes
      !! Number of different atom types

    real(kind=dp), intent(in) :: atomPositions(3,nAtoms)
      !! Atom positions in Cartesian coordinates

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    type(atom), intent(inout) :: atoms(nAtoms)
      !! Structure to hold details for each atom
    type(atom), intent(in) :: atomsSD(nAtoms)
      !! Structure to hold details for each atom
      !! for SD (needed for grid for Bessel functions)


    ! Output variables:
    complex(kind=dp), intent(out) :: pawK(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal,nGVecsLocal)

    ! Local variables:
    real(kind=dp) :: gCart(3)
      !! G-vector in Cartesian coordinates
    integer :: ibi, ibf, ispin, ig, ix
    integer :: LL, I, NI, LMBASE, LM
    integer :: L, M, ind, iT
    real(kind = dp) :: q, qDotR, FI, t1, t2
    
    real(kind = dp) :: JL(0:JMAX), v_in(3)
    complex(kind = dp) :: Y( (JMAX+1)**2 )
    complex(kind = dp) :: VifQ_aug, ATOMIC_CENTER
    
    
    call cpu_time(t1)

    pawK(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    do ig = 1, nGVecsLocal
      
      do ix = 1, 3
        gCart(ix) = sum(mill_local(:,ig)*recipLattVec(ix,:))
      enddo

      q = sqrt(sum(gCart(:)*gCart(:)))
      
      v_in(:) = gCart
      if ( abs(q) > 1.0e-6_dp ) v_in = v_in/q ! i have to determine v_in = q
      Y = cmplx(0.0_dp, 0.0_dp, kind = dp)
      CALL ylm(v_in, JMAX, Y) ! calculates all the needed spherical harmonics once
      
      LMBASE = 0
      
      do iT = 1, numOfTypes
        
        DO I = 1, atoms(iT)%iRc
          
          JL = 0.0_dp
          CALL bessel_j(q*atomsSD(iT)%r(I), JMAX, JL) ! returns the spherical bessel at qr point
            ! Previously used SD atoms structure here for both PC and SD
          atoms(iT)%bes_J_qr(:,I) = JL(:)
          
        ENDDO
        
      enddo
      
      do ni = 1, nAtoms ! LOOP OVER THE IONS
        
        qDotR = sum(gCart(:)*atomPositions(:,ni))
        
        if(crystalType == 'PC') then
          ATOMIC_CENTER = exp( -ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        else
          ATOMIC_CENTER = exp( ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        endif
        
        iT = iType(ni)
        LM = 0
        DO LL = 1, atoms(iT)%lMax
          L = atoms(iT)%LPS(LL)
          DO M = -L, L
            LM = LM + 1 !1st index for CPROJ
            
            FI = 0.0_dp
            
            FI = sum(atoms(iT)%bes_J_qr(L,:)*atoms(iT)%F(:,LL)) ! radial part integration F contains rab
            
            ind = L*(L + 1) + M + 1 ! index for spherical harmonics

            if(crystalType == 'PC') then
              VifQ_aug = ATOMIC_CENTER*Y(ind)*(-II)**L*FI
            else
              VifQ_aug = ATOMIC_CENTER*conjg(Y(ind))*(II)**L*FI
            endif

            do ibi = iBandIinit, iBandIfinal
              
              do ibf = iBandFinit, iBandFfinal
                
                if(crystalType == 'PC') then
                  pawK(ibf, ibi, ig) = pawK(ibf, ibi, ig) + VifQ_aug*cProjPC(LM + LMBASE, ibi)
                else
                  pawK(ibf, ibi, ig) = pawK(ibf, ibi, ig) + VifQ_aug*conjg(cProjSD(LM + LMBASE, ibf))
                endif
                
              enddo
            enddo
          ENDDO
        ENDDO

        LMBASE = LMBASE + atoms(iT)%lmMax
      ENDDO
      
    enddo
    
    return
    
  end subroutine pawCorrectionK
  
!----------------------------------------------------------------------------
  subroutine writeResults(ikLocal, isp)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ikGlobal
      !! Current global k-point
 
    character(len = 300) :: ikC, ispC
      !! Character indices
    
    integer :: ibi, ibf, totalNumberOfElements
    real(kind = dp) :: t1, t2
    
    character(len = 300) :: text


    ikGlobal = ikLocal+ikStart_pool-1
    
    call cpu_time(t1)
    
    call readEigenvalues(ikGlobal, isp)
    
    write(*, '(" Writing Ufi(:,:) of k-point ", i2, " and spin ", i1, ".")') ikGlobal, isp

    call int2str(ikGlobal, ikC)
    call int2str(isp, ispC)
    
    open(17, file=trim(elementsPath)//"/allElecOverlap."//trim(ispC)//"."//trim(ikC), status='unknown')
    
    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') omega
    
    text = "# Total number of <f|U|i> elements, Initial States (bandI, bandF), Final States (bandI, bandF)"
    write(17,'(a, " Format : ''(5i10)''")') trim(text)
    
    totalNumberOfElements = (iBandIfinal - iBandIinit + 1)*(iBandFfinal - iBandFinit + 1)
    write(17,'(5i10)') totalNumberOfElements, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    
    write(17, '("# Final Band, Initial Band, Delta energy, Complex <f|U|i>, |<f|U|i>|^2 Format : ''(2i10,4ES24.15E3)''")')
    
    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal
        
        write(17, 1001) ibf, ibi, eigvI(ibi) - eigvF(ibf), Ufi(ibf,ibi,ikLocal,isp), abs(Ufi(ibf,ibi,ikLocal,isp))**2
            
      enddo
    enddo
    
    close(17)
    
    call cpu_time(t2)
    write(*, '(" Writing Ufi(:,:) of k-point ", i4, "and spin ", i1, " done in:                   ", f10.2, " secs.")') &
      ikGlobal, isp, t2-t1
    
 1001 format(2i10,4ES24.15E3)
    
    return
    
  end subroutine writeResults
  
!----------------------------------------------------------------------------
  subroutine readEigenvalues(ikGlobal,isp)
    
    implicit none
    
    ! Input variables
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ib
      !! Loop index
    
    character(len = 300) :: ikC, ispC
      !! Character indices
    
    call int2str(ikGlobal, ikC)
    call int2str(min(isp,nSpinsSD), ispC)
    
    open(72, file=trim(exportDirSD)//"/eigenvalues."//trim(ispC)//"."//trim(ikC))
    
    read(72, * )
    read(72, * )
    
    do ib = 1, iBandIinit - 1
      read(72, *)
    enddo
    
    do ib = iBandIinit, iBandIfinal
      read(72, '(ES24.15E3)') eigvI(ib)
    enddo
    
    close(72)
    
    open(72, file=trim(exportDirSD)//"/eigenvalues."//trim(ispC)//"."//trim(ikC))
    
    read(72, * )
    read(72, * ) 
    
    do ib = 1, iBandFinit - 1
      read(72, *)
    enddo
    
    do ib = iBandFinit, iBandFfinal
      read(72, '(ES24.15E3)') eigvF(ib)
    enddo
    
    close(72)
    
    return
    
  end subroutine readEigenvalues
  
!----------------------------------------------------------------------------
  subroutine readUfis(ikLocal,isp)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ikGlobal
      !! Current global k-point
    
    character(len = 300) :: ikC, ispC
      !! Character indices
    
    integer :: ibi, ibf, totalNumberOfElements, iDum, i
    real(kind = dp) :: rDum, t1, t2
    complex(kind = dp):: cUfi


    ikGlobal = ikLocal+ikStart_pool-1
    
    call cpu_time(t1)
    write(*, '(" Reading Ufi(:,:) of k-point ", i4, " and spin ", i1)') ikGlobal, isp
    
    call int2str(ikGlobal, ikC)
    call int2str(isp, ispC)
    
    open(17, file=trim(elementsPath)//"/allElecOverlap."//trim(ispC)//"."//trim(ikC), status='unknown')
    
    read(17, *) 
    read(17, *) 
    read(17,'(5i10)') totalNumberOfElements, iDum, iDum, iDum, iDum
    read(17, *) 
    
    do i = 1, totalNumberOfElements
      
      read(17, 1001) ibf, ibi, rDum, cUfi, rDum
      Ufi(ibf,ibi,ikLocal,isp) = cUfi
          
    enddo
    
    close(17)
    
    call cpu_time(t2)
    write(*, '(" Reading Ufi(:,:) of k-point ", i4, " and spin ", i1, " done in:                   ", f10.2, " secs.")') &
      ikGlobal, isp, t2-t1
    
 1001 format(2i10,4ES24.15E3)
    
    return
    
  end subroutine readUfis
  
!----------------------------------------------------------------------------
  subroutine calculateVfiElements()
    
    implicit none
    
    integer :: ikLocal, ikGlobal, ib, nOfEnergies, iE, isp
    
    real(kind = dp) :: eMin, eMax, E, av, sd, x, EiMinusEf, A, DHifMin
    
    real(kind = dp), allocatable :: sumWk(:), sAbsVfiOfE2(:), absVfiOfE2(:)
    integer, allocatable :: nKsInEbin(:)
    
    character (len = 300) :: text
    character (len = 300) :: fNameBase
    character (len = 300) :: fNameSK
    character(len = 300) :: ikC, ispC


    allocate(DE(iBandIinit:iBandIfinal, nKPerPool, nSpins))
    allocate(absVfi2(iBandIinit:iBandIfinal, nKPerPool, nSpins))
     
    DE(:,:,:) = 0.0_dp
    absVfi2(:,:,:) = 0.0_dp 
    
    do isp = 1, nSpins
  
      do ikLocal = 1, nKPerPool

        ikGlobal = ikLocal+ikStart_pool-1
      
        eigvI(:) = 0.0_dp
        eigvF(:) = 0.0_dp
      
        call readEigenvalues(ikGlobal, isp)
      
        do ib = iBandIinit, iBandIfinal

          EiMinusEf = eigvI(ib) - eigvF(iBandFinit)
          absVfi2(ib,ikLocal,isp) = EiMinusEf**2*( abs(Ufi(iBandFinit,ib,ikLocal,isp))**2 - abs(Ufi(iBandFinit,ib,ikLocal,isp))**4 )
        
          DE(ib,ikLocal,isp) = sqrt(EiMinusEf**2 - 4.0_dp*absVfi2(ib,ikLocal,isp))

        enddo

      enddo
      
    enddo
    
    eMin = minval(DE(:,:,:))
    call MPI_ALLREDUCE(MPI_IN_PLACE, eMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, interPoolComm, ierr)

    eMax = maxval(DE(:,:,:))
    call MPI_ALLREDUCE(MPI_IN_PLACE, eMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, interPoolComm, ierr)

    nOfEnergies = int((eMax-eMin)/eBin) + 1
    
    allocate(absVfiOfE2(0:nOfEnergies), nKsInEbin(0:nOfEnergies), sumWk(0:nOfEnergies))
    
    absVfiOfE2(:) = 0.0_dp
    nKsInEbin(:) = 0
    sumWk(:) = 0.0_dp
    DHifMin = 0.0_dp
    
    do isp = 1, nSpins
      
      do ikLocal = 1, nKPerPool

        ikGlobal = ikLocal+ikStart_pool-1
      
        do ib = iBandIinit, iBandIfinal
        
          if(abs(eMin - DE(ib, ikLocal, isp)) < 1.0e-3_dp) DHifMin = absVfi2(ib, ikLocal,isp)

          iE = int((DE(ib, ikLocal, isp) - eMin)/eBin)

          if(absVfi2(ib, ikLocal, isp) > 0.0_dp) then

            absVfiOfE2(iE) = absVfiOfE2(iE) + wkPC(ikGlobal)*absVfi2(ib, ikLocal, isp)
  
            sumWk(iE) = sumWk(iE) + wkPC(ikGlobal)

            nKsInEbin(iE) = nKsInEbin(iE) + 1

          else
            write(*,*) 'absVfi2', absVfi2(ib, ikLocal, isp)
          endif
        
        enddo
      
      enddo

    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE, DHifMin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, interPoolComm, ierr)
    call mpiSumDoubleV(absVfiOfE2, interPoolComm)
    call mpiSumDoubleV(sumWk, interPoolComm)
    call mpiSumIntV(nKsInEbin, interPoolComm)

    allocate(sAbsVfiOfE2(0:nOfEnergies))
    
    sAbsVfiOfE2 = 0.0_dp
    
    do isp = 1, nSpins
  
      do ikLocal = 1, nKPerPool

        ikGlobal = ikLocal+ikStart_pool-1
    
        call int2str(ikGlobal, ikC)
        call int2str(isp, ispC)

        open(11, file=trim(VfisOutput)//'ofKpt.'//trim(ispC)//'.'//trim(ikC), status='unknown')
      
        do ib = iBandIinit, iBandIfinal
        
          iE = int((DE(ib, ikLocal, isp) - eMin)/eBin)

          av = absVfiOfE2(iE)/sumWk(iE)

          x = absVfi2(ib, ikLocal, isp)

          write(11, '(2ES24.15E3,2i10)') (eMin + iE*eBin), x, isp, ikGlobal
          !write(12, '(2ES24.15E3,i10)') DE(ib,ik), absVfi2(ib, ik), ik

          sAbsVfiOfE2(iE) = sAbsVfiOfE2(iE) + wkPC(ikGlobal)*(x - av)**2/sumWk(iE)
        
        enddo

        close(11)
      
      enddo

    enddo

    call mpiSumDoubleV(sAbsVfiOfE2, interPoolComm)
    
    if(myPoolId == 0) then

      open(11, file=trim(VfisOutput)//'ofKpt', status='unknown')
    
      write(11, '("# |<f|V|i>|^2 versus energy for all the k-points.")')
      write(text, '("# Energy (shifted by eBin/2) (Hartree), |<f|V|i>|^2 (Hartree)^2,")')
      write(11, '(a, " spin index, k-point index. Format : ''(2ES24.15E3,,i2,i10)''")') trim(text)

      close(11)

      do isp = 1, nSpins

        do ikGlobal = 1, nKPoints

          call int2str(ikGlobal, ikC)
          call int2str(isp, ispC)

          fNameBase = trim(VfisOutput)//'ofKpt'
          fNameSK = trim(fNameBase)//'.'//trim(ispC)//'.'//trim(ikC)

          call execute_command_line('cat '//trim(fNameSK)//' >> '//trim(fNameBase)//'&& rm '//trim(fNameSK))

        enddo

      enddo

      open(63, file=trim(VfisOutput), status='unknown')
    
      write(63, '("# Averaged |<f|V|i>|^2 over K-points versus energy.")')
      write(63, '("#                 Cell volume : ", ES24.15E3, " (a.u.)^3,   Format : ''(ES24.15E3)''")') omega
      write(63, '("#   Minimun transition energy : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eMin
      write(63, '("# |DHif|^2 at minimum Tr. En. : ", ES24.15E3, " (Hartree^2),Format : ''(ES24.15E3)''")') DHifMin
      write(63, '("#                  Energy bin : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eBin
      write(text, '("# Energy (Hartree), averaged |<f|V|i>|^2 over K-points (Hartree)^2,")')
      write(63, '(a, " standard deviation (Hartree)^2. Format : ''(3ES24.15E3)''")') trim(text)
    
      do iE = 0, nOfEnergies
        E = iE*eBin
        av = 0.0_dp
        sd = 0.0_dp
        if (nKsInEbin(iE) > 0) then
          av = absVfiOfE2(iE)/sumWk(iE)
          sd = sqrt(sAbsVfiOfE2(iE))
        endif
        write(63,'(3ES24.15E3)') eMin + E, av, sd
      enddo
    
      close(63)
    
    endif  
    
    return
    
  end subroutine calculateVfiElements
   
!----------------------------------------------------------------------------
  subroutine bessel_j (x, lmax, jl)
    
    ! x is the argument of j, jl(0:lmax) is the output values.
    implicit none
    integer, intent(in) :: lmax
    real(kind = dp), intent(in) :: x
    real(kind = dp), intent(out) :: jl(0:lmax)
    integer :: l
    
    if (x <= 0.0_dp) then
      jl = 0.0_dp
      jl(0) = 1.0_dp
      return
    end if
    
    jl(0) = sin(x)/x
    if (lmax <= 0) return
    jl(1) = (jl(0)-cos(x))/x
    if (lmax == 1) return
    
    do l = 2, lmax
      jl(l) = dble(2*l-1)*jl(l-1)/x - jl(l-2)
    enddo
    
    return
    
  end subroutine bessel_j
  
  
!----------------------------------------------------------------------------
  subroutine ylm(v_in,lmax,y)
  !
  ! lmax   : spherical harmonics are calculated for l = 0 to lmax
  ! v      : vector, argument of the spherical harmonics (we calculate
  ! Ylm(v/norm(v))
  ! y      : array containing Ylm(v) for several l,m
  !
  ! !DESCRIPTION:
  !   1.  PURPOSE
  !        The spherical harmonics (Condon and Shortley convention)
  !          Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
  !        for vector V (given in Cartesian coordinates)
  !        are calculated. In the Condon Shortley convention the
  !        spherical harmonics are defined as
  !        $$ Y(l,m) = (-1)^m \sqrt{\frac{1}{\pi}} P_{lm}(\cos{\theta})
  !        \rm
  !        e^{\rm i m \phi} $$
  !                        
  !        where  $P_{lm}(\cos{\theta})$ is the normalized Associated
  !        Legendre
  !                  
  !        function. Thus,
  !                                             
  !                     $$  Y(l,-m) = (-1)^m Y^*(l,m) $$
  !
  !   2.  USAGE
  !        DOUBLE PRECISION V(3), Y(5*5)
  !        V(1) = ...
  !        V(2) = ...
  !        V(3) = ...
  !        CALL YLM(V,4,Y)
  !
  !       ARGUMENT-DESCRIPTION
  !          V      - DOUBLE PRECISION vector, dimension 3        (input)
  !                   Must be given in Cartesian coordinates.
  !                   Conversion of V to polar coordinates gives the
  !                   angles Theta and Phi necessary for the calculation
  !                   of the spherical harmonics.
  !          LMAX   - INTEGER value                               (input)
  !                   upper bound of L for which spherical harmonics
  !                   will be calculated
  !                   constraint:
  !                      LMAX >= 0
  !          Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
  !                   contains the calculated spherical harmonics
  !                   Y(1)                   for L .EQ. 0 (M = 0)
  !                   Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
  !                   ...
  !                   Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
  !                                          for L .EQ. LMAX
  !                                              (M = -L,...,L)
  !                   constraint:
  !                      Dimension of Y .GE. (LMAX+1)**2 (not checked)
  !        USED SUBROUTINES (DIRECTLY CALLED)
  !           none
  !
  !        INDIRECTLY CALLED SUBROUTINES
  !           none
  !
  !        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
  !           none
  !
  !        INPUT/OUTPUT (READ/WRITE)
  !           none
  !
  !        MACHINENDEPENDENT PROGRAMPARTS
  !           Type COMPLEX*16 is used which does not conform to the
  !           FORTRAN 77 standard.
  !           Also the non-standard type conversion function DCMPLX()
  !           is used which combines two double precision values into
  !           one double complex value.
  !
  !   3.     METHOD
  !           The basic algorithm used to calculate the spherical
  !           harmonics for vector V is as follows:
  !
  !           Y(0,0)
  !           Y(1,0)
  !           Y(1,1)
  !           Y(1,-1) = -Y(1,1)
  !           DO L = 2, LMAX
  !              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
  !              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
  !              DO M = L-2, 0, -1
  !                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
  !                 Y(L,-M)= (-1)**M*Y(L,M)
  !              ENDDO
  !           ENDDO
  !
  !           In the following the necessary recursion formulas and
  !           starting values are given:
  !
  !        Start:
  !%                        +------+
  !%                        |   1     
  !%           Y(0,0) =  -+ | -----  
  !%                       \| 4(Pi)  
  !%
  !%                                   +------+
  !%                                   |   3     
  !%           Y(1,0) =  cos(Theta) -+ | -----  
  !%                                  \| 4(Pi)  
  !%
  !%                                     +------+
  !%                                     |   3    i(Phi)
  !%           Y(1,1) =  - sin(Theta) -+ | ----- e
  !%                                    \| 8(Pi)  
  !%
  !%        Formula 1:
  !%
  !%           Y(l,l) =
  !%                           +--------+
  !%                           | (2l+1)   i(Phi)
  !%            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
  !%                          \|   2l  
  !%
  !%        Formula 2:
  !%                                  +---------------+  
  !%                                  |  (2l-1)(2l+1)   
  !%           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
  !%                                 \|   (l-m)(l+m)       
  !%
  !%                                    +--------------------+  
  !%                                    |(l-1+m)(l-1-m)(2l+1)
  !%                              -  -+ |-------------------- Y(l-2,m)
  !%                                   \|  (2l-3)(l-m)(l+m)                 
  !%
  !%        Formula 3: (not used in the algorithm because of the division
  !%                    by sin(Theta) which may be zero)
  !%
  !%                                    +--------------+  
  !%                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
  !%           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
  !%                      sin(Theta)   \| (l+m+1)(l-m)       
  !%
  !%                                    +--------------+  
  !%                                    |(l-m-1)(l+m+2)  -2i(Phi)
  !%                              -  -+ |-------------- e        Y(l,m+2)
  !%                                   \| (l-m)(l+m+1)                         
  !%                                  
  !%
  ! !REVISION HISTORY:
  !   26. April 1994                                   Version 1.2
  !   Taken 8 1 98 from SRC_lapw2 to SRC_telnes
  !   Updated November 2004 (Kevin Jorissen)
  !   cosmetics March 2005 (Kevin Jorissen)
  !
      implicit none
  
  !   In/Output :
  
      integer, intent(in) :: LMAX
      real(kind = dp), intent(in) :: V_in(3)
      complex(kind = dp), intent(out) :: Y(*)
  !   Local variables :
      real(kind = dp), parameter :: pi = 3.1415926535897932384626433_dp
  
      INTEGER         ::  I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      real(kind = dp) ::  A, B, C, AB, ABC, ABMAX, ABCMAX, V(3)
      real(kind = dp) ::  D4LL1C, D2L13
      real(kind = dp) ::  COSTH, SINTH, COSPH, SINPH
      real(kind = dp) ::  TEMP1, TEMP2, TEMP3
      real(kind = dp) ::  YLLR, YLL1R, YL1L1R, YLMR
      real(kind = dp) ::  YLLI, YLL1I, YL1L1I, YLMI
      
      ! Y(0,0)
      
      do INDEX = 1,3
        V(INDEX) = dble(V_in(INDEX))
      enddo
      YLLR = 1.0_dp/sqrt(4.0_dp*PI)
      YLLI = 0.0_dp
      Y(1) = CMPLX(YLLR, YLLI, kind = dp)
      
      ! continue only if spherical harmonics for (L .GT. 0) are desired
      
      IF (LMAX .LE. 0) GOTO 999
      
      ! calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
      ! Theta, Phi ... polar angles of vector V
      
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0_dp) THEN
        A = V(1)/ABMAX
        B = V(2)/ABMAX
        AB = SQRT(A*A+B*B)
        COSPH = A/AB
        SINPH = B/AB
      ELSE
        COSPH = 1.0_dp
        SINPH = 0.0_dp
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. dble(0)) THEN
        A = V(1)/ABCMAX
        B = V(2)/ABCMAX
        C = V(3)/ABCMAX
        AB = A*A + B*B
        ABC = SQRT(AB + C*C)
        COSTH = C/ABC
        SINTH = SQRT(AB)/ABC
      ELSE
        COSTH = 1.0_dp
        SINTH = 0.0_dp
      ENDIF
      
      ! Y(1,0)
      
      Y(3) = CMPLX(sqrt(3.0_dp)*YLLR*COSTH, 0.0_dp, kind = dp)
      
      ! Y(1,1) ( = -DCONJG(Y(1,-1)))
      
      TEMP1 = -SQRT(1.5_dp)*YLLR*SINTH
      Y(4) = CMPLX(TEMP1*COSPH,TEMP1*SINPH, kind = dp)
      Y(2) = -CONJG(Y(4))
      
      DO L = 2, LMAX
        INDEX  = L*L + 1
        INDEX2 = INDEX + 2*L
        MSIGN  = 1 - 2*MOD(L,2)
        
        ! YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
        
        YL1L1R = DBLE(Y(INDEX-1))
        YL1L1I = DIMAG(Y(INDEX-1))
        TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
        YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
        YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
        Y(INDEX2) = CMPLX(YLLR,YLLI, kind = dp)
        Y(INDEX)  = cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
        !Y(INDEX)  = dble(MSIGN)*CONJG(Y(INDEX2))
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        
        ! YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
        ! (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
        
        TEMP2 = SQRT(DBLE(2*L+1))*COSTH
        YLL1R = TEMP2*YL1L1R
        YLL1I = TEMP2*YL1L1I
        Y(INDEX2) = CMPLX(YLL1R,YLL1I, kind = dp)
        Y(INDEX)  = -cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
  !      Y(INDEX)  = -dble(MSIGN)*CONJG(Y(INDEX2))
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        
        I4L2 = INDEX2 - 4*L + 2
        I2L  = INDEX2 - 2*L
        D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
        D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
        
        DO M = L - 2, 0, -1
          
          ! YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
          
          TEMP1 = 1.0_dp/SQRT(DBLE((L+M)*(L-M)))
          TEMP2 = D4LL1C*TEMP1
          TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
          YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
          YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
          Y(INDEX2) = CMPLX(YLMR,YLMI, kind = dp)
          Y(INDEX)  = cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
    !      Y(INDEX)  = dble(MSIGN)*CONJG(Y(INDEX2))
          
          MSIGN  = -MSIGN
          INDEX2 = INDEX2 - 1
          INDEX  = INDEX  + 1
          I4L2   = I4L2   - 1
          I2L    = I2L    - 1
        ENDDO
      ENDDO
      
  999 RETURN
  END subroutine ylm
  
!----------------------------------------------------------------------------
  subroutine finalizeCalculation()
    
    implicit none

    integer :: iType
      !! Loop index


    deallocate(npwsPC)
    deallocate(wkPC)
    deallocate(xkPC)
    deallocate(posIonPC)
    deallocate(TYPNIPC)
    deallocate(mill_local)

    do iType = 1, numOfTypesPC
      deallocate(atomsPC(iType)%r)
      deallocate(atomsPC(iType)%lps)
      deallocate(atomsPC(iType)%F)
      deallocate(atomsPC(iType)%F1)
      deallocate(atomsPC(iType)%F2)
      deallocate(atomsPC(iType)%bes_J_qr)
    enddo

    deallocate(atomsPC)

    deallocate(npwsSD)
    deallocate(wk)
    deallocate(xk)
    deallocate(posIonSD)
    deallocate(TYPNISD)

    do iType = 1, numOfTypes
      deallocate(atoms(iType)%lps)
      deallocate(atoms(iType)%r)
      deallocate(atoms(iType)%F)
      deallocate(atoms(iType)%F1)
      deallocate(atoms(iType)%F2)
      deallocate(atoms(iType)%bes_J_qr)
    enddo

    deallocate(atoms)

    if(indexInPool == 0) deallocate(eigvI, eigvF)

    call MPI_Barrier(worldComm, ierr)
    
    if(ionode) then
    
      call execute_command_line('rm -r ./workingFiles')

      write(iostd,'("-----------------------------------------------------------------")')
    
      call cpu_time(tf)
      write(iostd, '(" Total time needed:                         ", f10.2, " secs.")') tf-t0
    
      close(iostd)

    endif

    call MPI_Barrier(worldComm, ierr)
    
    return
    
  end subroutine finalizeCalculation

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
  subroutine int2str(integ, string)
    
    implicit none
    integer :: integ
    character(len = 300) :: string
    
    if ( integ < 10 ) then
      write(string, '(i1)') integ
    else if ( integ < 100 ) then
      write(string, '(i2)') integ
    else if ( integ < 1000 ) then
      write(string, '(i3)') integ
    else if ( integ < 10000 ) then
      write(string, '(i4)') integ
    else if ( integ < 100000 ) then
      write(string, '(i5)') integ
    else if ( integ < 1000000 ) then
      write(string, '(i6)') integ
    else if ( integ < 10000000 ) then
      write(string, '(i7)') integ
    endif
    
    string = trim(string)
    
    return
    
  end subroutine int2str
  
!----------------------------------------------------------------------------
  subroutine int2strLeadZero(integ, length, string)
    
    implicit none
    integer :: integ
    integer :: length
    character(len = 300) :: string
    character(len = 300) :: lengthC
    character(len = 300) :: formatString
    
    call int2str(length, lengthC)
    formatString = '(i'//trim(lengthC)//'.'//trim(lengthC)//')'

    write(string, formatString) integ
    
    string = trim(string)
    
    return
    
  end subroutine int2strLeadZero
  
  
end module declarations
