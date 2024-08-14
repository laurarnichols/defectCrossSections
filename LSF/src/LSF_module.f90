module LSFmod
  
  use constants, only: dp, HartreeToEv, ii, hbar_atomic, time_atomicToSI
  use base, only: nKPoints, order
  use TMEmod, only: getMatrixElementFNameWPath, getMatrixElementFName, readMatrixElement
  use PhononPPMod, only: diffOmega, readSjOneFreq, readSjTwoFreq, omega, omegaPrime, Sj, SjPrime, readNj
  use miscUtilities, only: int2strLeadZero, int2str
  use errorsAndMPI

  use energyTabulatorMod, only: energyTableDir, readCaptureEnergyTable

  implicit none 

  integer :: iTime_start, iTime_end
    !! Start and end time steps for this process
  integer :: nStepsLocal
    !! Number of time steps for each
    !! process to complete


  integer, allocatable :: ibi(:)
    !! Initial-state indices
  integer :: ibf
    !! Final-state index
  integer :: iSpin
    !! Spin channel to use
  integer, allocatable :: jReSort(:)
    !! Indices to optionally resort matrix elements
  integer :: mDim
    !! Size of first dimension for matrix element
  integer :: nModes
    !! Number of phonon modes
  integer :: nTransitions
    !! Total number of transitions 
  integer :: suffixLength
    !! Length of shifted POSCAR file suffix

  real(kind=dp), allocatable :: dE(:,:,:)
    !! All energy differences from energy table
  real(kind=dp) :: dt
    !! Time step size
  real(kind=dp) :: gamma0
    !! \(\gamma\) for Lorentzian smearing
  real(kind=dp) :: hbarGamma
    !! \(\hbar\gamma\) for Lorentzian smearing
    !! to guarantee convergence
  real(kind=dp), allocatable :: matrixElement(:,:,:)
    !! Electronic matrix element
  real(kind=dp) :: maxTime
    !! Max time for integration
  real(kind=dp), allocatable :: nj(:)
    !! \(n_j\) occupation number
  real(kind=dp) :: SjThresh
    !! Threshold for Sj to determine which modes to calculate
  real(kind=dp) :: smearingExpTolerance
    !! Tolerance for the Lorentzian-smearing
    !! exponential used to calculate max time

  logical :: newEnergyTable
    !! If this code and TME are being run with a different
    !! energy table
  logical :: oldFormat
    !! If the old format of the matrix element files
    !! should be used
  logical :: rereadDq
    !! If dq should be read from matrix element file
    !! (.false.) or from the dq.txt file (.true.)
  logical :: reSortMEs
    !! If matrix elements should be resorted

  character(len=300) :: matrixElementDir
    !! Path to matrix element file `allElecOverlap.isp.ik`. 
    !! For first-order term, the path is just within each 
    !! subdirectory.
  character(len=300) :: MjBaseDir
    !! Path to the base directory for the first-order
    !! matrix element calculations
  character(len=300) :: njInput
    !! Path to nj file
  character(len=300) :: outputDir
    !! Path to output transition rates
  character(len=300) :: PhononPPDir
    !! Path to PhononPP output dir to get Sj.out
    !! and potentially optimalPairs.out
  character(len=300) :: prefix
    !! Prefix of directories for first-order matrix
    !! elements
  character(len=300) :: volumeLine
    !! Volume line from overlap file to be
    !! output exactly in transition rate file


  namelist /inputParams/ energyTableDir, matrixElementDir, MjBaseDir, PhononPPDir, njInput, hbarGamma, dt, &
                         smearingExpTolerance, outputDir, order, prefix, iSpin, diffOmega, newEnergyTable, &
                         suffixLength, reSortMEs, oldFormat, rereadDq, SjThresh

contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iSpin, order, dt, gamma0, hbarGamma, maxTime, SjThresh, smearingExpTolerance, diffOmega, &
        newEnergyTable, oldFormat, rereadDq, reSortMEs, energyTableDir, matrixElementDir, MjBaseDir, njInput, outputDir, &
        PhononPPDir, prefix)

    implicit none

    ! Output variables
    integer, intent(out) :: iSpin
      !! Spin channel to use
    integer, intent(out) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(out) :: dt
      !! Time step size
    real(kind=dp), intent(out) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(out) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: maxTime
      !! Max time for integration
    real(kind=dp), intent(out) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(out) :: newEnergyTable
      !! If this code and TME are being run with a different
      !! energy table
    logical, intent(out) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used
    logical, intent(out) :: rereadDq
      !! If dq should be read from matrix element file
      !! (.false.) or from the dq.txt file (.true.)
    logical, intent(out) :: reSortMEs
      !! If matrix elements should be resorted

    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(out) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: njInput
      !! Path to nj file
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: PhononPPDir
      !! Path to PhononPP output dir to get Sj.out
      !! and potentially optimalPairs.out
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements

  
    call initialize(iSpin, order, dt, hbarGamma, SjThresh, smearingExpTolerance, diffOmega, newEnergyTable, &
          oldFormat, rereadDq, reSortMEs, energyTableDir, matrixElementDir, MjBaseDir, njInput, outputDir, &
          PhononPPDir, prefix)

    if(ionode) then

      read(5, inputParams, iostat=ierr)
        !! * Read input variables

    
      if(ierr /= 0) call exitError('LSF module', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(iSpin, order, dt, hbarGamma, SjThresh, smearingExpTolerance, diffOmega, newEnergyTable, &
            oldFormat, rereadDq, reSortMEs, energyTableDir, matrixElementDir, MjBaseDir, njInput, outputDir, &
            PhononPPDir, prefix)

      gamma0 = hbarGamma*1e-3/HartreeToEv
        ! Input expected in meV

      maxTime = -log(smearingExpTolerance)/gamma0
      write(*,'("Max time: ", ES24.15E3)') maxTime

    endif

    call MPI_BCAST(iSpin, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)
  
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(hbarGamma, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(gamma0, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(SjThresh, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(smearingExpTolerance, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(maxTime, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(diffOmega, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(newEnergyTable, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(oldFormat, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(reSortMEs, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(rereadDq, 1, MPI_LOGICAL, root, worldComm, ierr)
  
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(matrixElementDir, len(matrixElementDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(MjBaseDir, len(MjBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(njInput, len(njInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(PhononPPDir, len(PhononPPDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)
    
    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(iSpin, order, dt, hbarGamma, SjThresh, smearingExpTolerance, diffOmega, newEnergyTable, &
        oldFormat, rereadDq, reSortMEs, energyTableDir, matrixElementDir, MjBaseDir, njInput, outputDir, &
        PhononPPDir, prefix)

    implicit none

    ! Output variables
    integer, intent(out) :: iSpin
      !! Spin channel to use
    integer, intent(out) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(out) :: dt
      !! Time step size
    real(kind=dp), intent(out) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(out) :: newEnergyTable
      !! If this code and TME are being run with a different
      !! energy table
    logical, intent(out) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used
    logical, intent(out) :: rereadDq
      !! If dq should be read from matrix element file
      !! (.false.) or from the dq.txt file (.true.)
    logical, intent(out) :: reSortMEs
      !! If matrix elements should be resorted

    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(out) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: njInput
      !! Path to nj file
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: PhononPPDir
      !! Path to PhononPP output dir to get Sj.out
      !! and potentially optimalPairs.out
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements


    iSpin = 1
    order = -1

    dt = 1d-4
    hbarGamma = 0.0_dp
    SjThresh = 0.0_dp
    smearingExpTolerance = 0.0_dp

    diffOmega = .false.
    newEnergyTable = .false.
    oldFormat = .false.
    reSortMEs = .false.
    rereadDq = .false.

    energyTableDir = ''
    matrixElementDir = ''
    MjBaseDir = ''
    njInput = ''
    outputDir = './'
    PhononPPDir = ''
    prefix = 'disp-'

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iSpin, order, dt, hbarGamma, SjThresh, smearingExpTolerance, diffOmega, newEnergyTable, &
        oldFormat, rereadDq, reSortMEs, energyTableDir, matrixElementDir, MjBaseDir, njInput, outputDir, &
        PhononPPDir, prefix)

    implicit none

    ! Input variables
    integer, intent(in) :: iSpin
      !! Spin channel to use
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dt
      !! Time step size
    real(kind=dp), intent(in) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(in) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(in) :: newEnergyTable
      !! If this code and TME are being run with a different
      !! energy table
    logical, intent(in) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used
    logical, intent(in) :: rereadDq
      !! If dq should be read from matrix element file
      !! (.false.) or from the dq.txt file (.true.)
    logical, intent(in) :: reSortMEs
      !! If matrix elements should be resorted

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(in) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(in) :: njInput
      !! Path to nj file
    character(len=300), intent(in) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(in) :: PhononPPDir
      !! Path to PhononPP output dir to get Sj.out
      !! and potentially optimalPairs.out
    character(len=300), intent(in) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements

    ! Local variables:
    character(len=300) :: fName
      !! File name for matrix element file to get
      !! nKPoints from 

    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iSpin', iSpin, 1, 2)
    abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution 

    abortExecution = checkDoubleInitialization('dt', dt, 1.0d-6, 1.0d-2) .or. abortExecution
    abortExecution = checkDoubleInitialization('hbarGamma', hbarGamma, 0.1_dp, 20.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('SjThresh', SjThresh, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('smearingExpTolerance', smearingExpTolerance, 0.0_dp, 1.0_dp) .or. abortExecution
      ! These limits are my best guess as to what is reasonable; they are not
      ! hard and fast, but you should think about the application of the theory
      ! to numbers outside these ranges.

    abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.'//trim(int2str(iSpin))//'.1') .or. abortExecution
    abortExecution = checkDirInitialization('PhononPPDir', PhononPPDir, 'Sj.out') .or. abortExecution
    abortExecution = checkFileInitialization('njInput', njInput) .or. abortExecution

    write(*,'("diffOmega = ",L)') diffOmega
    write(*,'("oldFormat = ",L)') oldFormat
    write(*,'("newEnergyTable = ",L)') newEnergyTable

    if(order == 0) then 
      abortExecution = checkDirInitialization('matrixElementDir', matrixElementDir, getMatrixElementFName(1,iSpin)) .or. abortExecution

      fName = getMatrixElementFNameWPath(1,iSpin,matrixElementDir)

    else if(order == 1) then

      write(*,'("reSortMEs = ",L)') reSortMEs
      if(reSortMEs) abortExecution = checkFileInitialization('optimalPairsFile', trim(PhononPPDir)//'/optimalPairs.out') .or. abortExecution

      write(*,'("rereadDq = ",L)') rereadDq
      if(rereadDq) abortExecution = checkFileInitialization('dqFile', trim(PhononPPDir)//'/dq.txt') .or. abortExecution

      abortExecution = checkIntInitialization('suffixLength', suffixLength, 1, 5) .or. abortExecution 
      abortExecution = checkDirInitialization('MjBaseDir', MjBaseDir, &
            '/'//trim(prefix)//trim(int2strLeadZero(1,suffixLength))//'/'//trim(getMatrixElementFNameWPath(1,iSpin,matrixElementDir))) .or. abortExecution
      write(*,'("prefix = ''",a,"''")') trim(prefix)
      write(*,'("matrixElementDir = ''",a,"''")') trim(matrixElementDir)

      fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(1,suffixLength))//'/'//trim(getMatrixElementFNameWPath(1,iSpin,matrixElementDir))

    endif

    open(unit=12,file=trim(fName))
    read(12,*)
    read(12,*) nKPoints
    close(12)

    write(*,'("nKPoints = ", i10)') nKPoints

    call system('mkdir -p '//trim(outputDir))


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif

    return 

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readEnergyTable(iSpin, energyTableDir, nTransitions, ibi, ibf, dE)

    implicit none

    ! Input variables:
    integer, intent(in) :: iSpin
      !! Spin channel to use

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to read

    ! Output variables:
    integer, intent(out) :: nTransitions
      !! Total number of transitions 
    integer, allocatable, intent(out) :: ibi(:)
      !! Initial-state indices
    integer, intent(out) :: ibf
      !! Final-state index

    real(kind=dp), allocatable, intent(out) :: dE(:,:,:)
      !! All energy differences from energy table

    ! Local variables:
    integer, allocatable :: iDum1D(:)
      !! Integer to ignore input
    integer :: iDum1, iDum2
      !! Integers to ignore input
    integer :: ikLocal, ikGlobal
      !! Loop indices

    real(kind=dp), allocatable :: dE2D(:,:)
      !! 2D array to read energy at single k-point


    ! Get the number of transitions and the state indices
    if(ionode) then
      call readCaptureEnergyTable(1, iSpin, energyTableDir, ibi, ibf, nTransitions, dE2D)
       ! Assume that band bounds and number of transitions do not depend on k-points or spin
       ! We ignore the energy here because we are only reading ibi, ibf, and nTransitions

      deallocate(dE2D)
    endif

    call MPI_BCAST(nTransitions, 1, MPI_INTEGER, root, worldComm, ierr)
    if(.not. ionode) allocate(ibi(nTransitions))
    call MPI_BCAST(ibi, nTransitions, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ibf, 1, MPI_INTEGER, root, worldComm, ierr)


    allocate(dE(3,nTransitions,nkPerPool))


    if(indexInPool == 0) then

      dE = 0.0_dp

      do ikLocal = 1, nkPerPool
    
        ! Get the global `ik` index from the local one
        ikGlobal = ikLocal+ikStart_pool-1


        call readCaptureEnergyTable(ikGlobal, iSpin, energyTableDir, iDum1D, iDum1, iDum2, dE2D)
          ! Assume that band bounds and number of transitions do not depend on k-points or spin

        dE(:,:,ikLocal) = dE2D
        deallocate(dE2D)

      enddo

    endif


    call MPI_BCAST(dE, size(dE), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)


    return

  end subroutine readEnergyTable

!----------------------------------------------------------------------------
  subroutine readSj(diffOmega, PhononPPDir, nModes, omega, omegaPrime, Sj, SjPrime)

    implicit none

    ! Input variables:
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(in) :: PhononPPDir
      !! Path to PhononPP output dir to get Sj.out
      !! and potentially optimalPairs.out

    ! Output variables:
    integer, intent(out) :: nModes
      !! Number of phonon modes

    real(kind=dp),allocatable, intent(out) :: omega(:), omegaPrime(:)
      !! Frequency for each mode
    real(kind=dp), allocatable, intent(out) :: Sj(:), SjPrime(:)
      !! Huang-Rhys factor for each mode

    ! Local variables:
    !real(kind=dp), allocatable :: randVal(:)
      ! Random adjustment to be made to frequencies
      ! to test sensitivity

    character(len=300) :: fName
      !! File name to read


    fName = trim(PhononPPDir)//'/Sj.out' 
    if(diffOmega) then
      call readSjTwoFreq(fName, nModes, omega, omegaPrime, Sj, SjPrime)

      ! This is the code that I used to test what difference different
      ! frequencies would have. I input the same two frequencies twice
      ! and randomly adjusted the omegaPrime. 
      !
      ! If doing a random adjustment, must do with a single process then
      ! broadcast, otherwise the random variables will be different across
      ! processes.
      !allocate(randVal(nModes))
      !if(ionode) then
        !call random_seed()
        !call random_number(randVal)
        !randVal = (randVal*2.0_dp - 1.0_dp)*0.50_dp  ! the number multiplied here is the % adjustment
        !SjPrime(:) = SjPrime(:)/omegaPrime(:)
        !omegaPrime(:) = omegaPrime(:)*(1.0_dp + randVal) ! this applies the random adjustment
        !omegaPrime(:) = omegaPrime(:)*(1.0_dp + 0.5_dp) ! this applies a uniform adjustment
        !SjPrime(:) = SjPrime(:)*omegaPrime(:)
      !endif
      !deallocate(randVal)

      ! Need to rebroadcast omegaPrime if we adjust it
      !call MPI_BCAST(omegaPrime, size(omegaPrime), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      !call MPI_BCAST(SjPrime, size(SjPrime), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    else
      allocate(omegaPrime(nModes), SjPrime(nModes))
        ! Need to allocate to avoid issues with passing variables
        ! and deallocating

      call readSjOneFreq(fName, nModes, omega, Sj)
    endif


    return

  end subroutine readSj

!----------------------------------------------------------------------------
  subroutine getjReSort(nModes, PhononPPDir, jReSort)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of phonon modes

    character(len=300), intent(in) :: PhononPPDir
      !! Path to PhononPP output dir to get Sj.out
      !! and potentially optimalPairs.out

    ! Output variables:
    integer, intent(out) :: jReSort(nModes)
      !! Indices to optionally resort matrix elements

    ! Local variables:
    integer :: j
      !! Loop index
    integer :: jCurrent
      !! Current-ordered index
    integer :: jNew
      !! New-ordered index


    open(unit=32,file=trim(PhononPPDir)//'/optimalPairs.out')

    read(32,*)

    do j = 1, nModes
      read(32,*) jCurrent, jNew
      jReSort(jCurrent) = jNew
    enddo

    close(32)

    return

  end subroutine

!----------------------------------------------------------------------------
  subroutine readAllMatrixElements(iSpin, nTransitions, ibi, nModes, jReSort, order, suffixLength, dE, newEnergyTable, &
          oldFormat, rereadDq, reSortMEs, matrixElementDir, MjBaseDir, PhononPPDir, prefix, mDim, matrixElement, volumeLine)

    implicit none

    ! Input variables:
    integer, intent(in) :: iSpin
      !! Spin channel to use
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jReSort(nModes)
      !! Indices to optionally resort matrix elements
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)
    integer, intent(in) :: suffixLength
      !! Length of shifted POSCAR file suffix

    real(kind=dp), intent(in) :: dE(3,nTransitions,nkPerPool)
      !! All energy differences from energy table

    logical, intent(in) :: newEnergyTable
      !! If this code and TME are being run with a different
      !! energy table
    logical, intent(in) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used
    logical, intent(in) :: rereadDq
      !! If dq should be read from matrix element file
      !! (.false.) or from the dq.txt file (.true.)
    logical, intent(in) :: reSortMEs
      !! If matrix elements should be resorted

    character(len=300), intent(in) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(in) :: PhononPPDir
      !! Path to PhononPP output dir to get Sj.out
      !! and potentially optimalPairs.out
    character(len=300), intent(in) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements

    ! Output variables:
    integer, intent(out) :: mDim
      !! Size of first dimension for matrix element

    real(kind=dp), allocatable, intent(out) :: matrixElement(:,:,:)
      !! Electronic matrix element

    character(len=300), intent(out) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: ikLocal, ikGlobal, j, jStore
      !! Loop indices

    real(kind=dp), allocatable :: dENew(:)
      !! New energy to update matrix element; needed
      !! not to create a temporary array when passing
      !! a slice of dE
    real(kind=dp), allocatable :: ME_tmp(:)
      !! Temporary storage of matrix element

    character(len=300) :: fName
      !! File name to read


    allocate(dENew(nTransitions))

    if(order == 0) then
      mDim = 1
    else if(order == 1) then
      mDim = nModes
    endif

    allocate(matrixElement(mDim,nTransitions,nkPerPool))
    allocate(ME_tmp(nTransitions))


    if(indexInPool == 0) then

      matrixElement = 0.0_dp

      do ikLocal = 1, nkPerPool
    
        ikGlobal = ikLocal+ikStart_pool-1
          !! Get the global `ik` index from the local one

        ! For both the zeroth-order and first-order matrix elements, there
        ! is the option to specify `newEnergyTable`. If this variable is true,
        ! that means that the matrix elements (TME) and LSF code were not run
        ! with the same energy table as input. This might be beneficial if tweaks
        ! were made to the energies used or you only want to use a subset of 
        ! the transitions for the LSF code. 
        !
        ! With `newEnergyTable = .true.`, the states to be read from the matrix
        ! element file(s) come from the energy table input to the LSF code. All
        ! of the states in the energy table input to the LSF code must be present
        ! in the matrix element file. 
        !
        ! Both terms also allow you to specify the old format where both the
        ! initial and final band indices were given for capture and the headers
        ! were different. I needed this for the Si calculations and left it here
        ! just in case.
        !..........................................................................
        ! The zeroth-order term reads only one matrix element, with or without the
        ! new energy table
        if(order == 0) then
          fName = getMatrixElementFNameWPath(ikGlobal, iSpin, matrixElementDir)

          dENew = dE(2,:,ikLocal)
            ! The second index holds the zeroth-order energy
          if(newEnergyTable) then
            call readMatrixElement(minval(ibi), maxval(ibi), nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
          else
            call readMatrixElement(-1, -1, nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
              ! dENew will be ignored here
          endif

          matrixElement(1,:,ikLocal) = ME_tmp

        !..........................................................................
        ! The first-order term reads a matrix element for each mode. In addition to
        ! the options `newEnergyTable` and `oldFormat` in the zeroth-order term, 
        ! the first-order term also has options `reSortMEs` and `rereadDq`.
        !
        ! `reSortMEs` assumes that the modes were labeled differently when the
        ! matrix elements were calculated vs what is to be used for the LSF run.
        ! This can be the case if you want to compare the difference between
        ! allowing the frequencies to change for different states (e.g., before and
        ! after capture). 
        !
        ! `rereadDq` allows you to read a new delta q for each mode than was used
        ! in the calculation of the matrix elements. You shouldn't need this in 
        ! general, but there was a bug in the code that only messed up the dq.txt
        ! file previously, so it was useful for me to be able to fix the dq here 
        ! rather than having to redo all of the matrix element calculations. 
        !
        ! The volume line will get overwritten each time through the loop, but 
        ! that's okay because the volume doesn't change between the files. 
        else if(order == 1) then

          dENew = dE(3,:,ikLocal)
            ! The third index holds the first-order energy
    
          do j = 1, nModes

            fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(j,suffixLength))//'/'&
                    //trim(getMatrixElementFNameWPath(ikGlobal,iSpin,matrixElementDir))

            ! If resorting the matrix element files based on a different PhononPP output
            ! order, make sure to pass the resorted mode index if re-reading dq's
            if(reSortMEs) then
              jStore = jReSort(j)
            else 
              jStore = j
            endif


            ! For new energy table and new dqs, must pass the band bounds, mode index
            ! to read (jStore) and PhononPPDir
            if(newEnergyTable .and. rereadDq) then
              call readMatrixElement(minval(ibi), maxval(ibi), nTransitions, order, dENew, newEnergyTable, oldFormat, &
                    fName, ME_tmp, volumeLine, jStore, PhononPPDir)

            ! For just new energy table, only pass band bounds
            else if(newEnergyTable .and. (.not. rereadDq)) then
              call readMatrixElement(minval(ibi), maxval(ibi), nTransitions, order, dENew, newEnergyTable, oldFormat, &
                    fName, ME_tmp, volumeLine)

            ! For just new dq, only pass jStore and PhononPPDir
            else if((.not. newEnergyTable) .and. rereadDq) then
              call readMatrixElement(-1, -1, nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine, &
                jStore, PhononPPDir)
                  ! dENew will be ignored here

            ! For neither, don't pass anything and just read matrix elements as-is
            else if((.not. newEnergyTable) .and. (.not. rereadDq)) then
              call readMatrixElement(-1, -1, nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
                ! dENew will be ignored here
            endif

            ! Store the matrix element for this mode
            matrixElement(jStore,:,ikLocal) = ME_tmp

          enddo
        endif
      enddo
    endif

    deallocate(dENew)
    deallocate(ME_tmp)

    call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)

    return

  end subroutine readAllMatrixElements

!----------------------------------------------------------------------------
  subroutine getAndWriteTransitionRate(nTransitions, ibi, iSpin, mDim, nModes, order, dE, dt, gamma0, & 
          matrixElement, nj, omega, omegaPrime, Sj, SjPrime, SjThresh, diffOmega, volumeLine)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: iSpin
      !! Spin index
    integer, intent(in) :: mDim
      !! Size of first dimension for matrix element
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dE(3,nTransitions,nkPerPool)
      !! All energy differences from energy table
    real(kind=dp), intent(in) :: dt
      !! Time step size
    real(kind=dp), intent(in) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(in) :: matrixElement(mDim,nTransitions,nkPerPool)
      !! Electronic matrix element
    real(kind=dp), intent(in) :: nj(nModes)
      !! \(n_j\) occupation number
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: Sj(nModes), SjPrime(nModes)
      !! Huang-Rhys factor for each mode
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(in) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: countTrue
      !! Number of modes where Sj > SjThresh
    integer :: iTime, iE, ikLocal, ikGlobal
      !! Loop indices
    integer, allocatable :: jTrue(:)
      !! Mode indices where Sj > SjThresh
    integer :: updateFrequency
      !! Frequency of steps to write status update

    real(kind=dp) :: Eif
      !! Local storage of energy for delta function
    real(kind=dp) :: multFact
      !! Multiplication factor for each term per
      !! Simpson's integration method
    real(kind=dp) :: t0
      !! Initial time for this process
    real(kind=dp) :: time
      !! Time for each time step
    real(kind=dp) :: timer1, timer2
      !! Timers
    real(kind=dp) :: transitionRate(nTransitions,nkPerPool)
      !! \(Gamma_i\) transition rate

    complex(kind=dp) :: Dj0_t(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
    complex(kind=dp) :: Dj1_t(nModes)
      !! Factor multiplying |M_j|^2*G_j^0(t) in G_j^1(t)
    complex(kind=dp) :: expArg_base
      !! Base exponential argument for each time step
    complex(kind=dp) :: expArg
      !! Exponential argument for each time step
    complex(kind=dp) :: expPrefactor
      !! Prefactor for the exponential inside the integral.
      !! For zeroth-order, this is just the matrix element,
      !! but for first-order this is \(\sum_j M_j A_j\).

    logical :: mask(nModes)
      !! Select the modes to calculate

      
    updateFrequency = ceiling(nStepsLocal/10.0)
    call cpu_time(timer1)


    ! Create a mask to determine which modes to calculate
    ! based on an optional threshold given (default 0.0).
    mask = Sj(:) >= SjThresh
    call getTrueIndices(nModes, mask, countTrue, jTrue)


    t0 = indexInPool*float(nStepsLocal-1)*dt

    transitionRate(:,:) = 0.0_dp
    do iTime = 0, nStepsLocal-1

      if(ionode .and. mod(iTime,updateFrequency) == 0) then

        call cpu_time(timer2)
        write(*,'("    ", i2,"% complete with transition-rate loop. Time in loop: ",f10.2," secs")') int((iTime*100.0)/nStepsLocal), timer2-timer1

      endif

      time = t0 + float(iTime)*dt
        ! Must do this arithmetic with floats to avoid
        ! integer overflow

      call setupTimeTables(countTrue, nModes, jTrue, nj, omega, omegaPrime, Sj, SjPrime, time, diffOmega, Dj0_t, Dj1_t)

      expArg_base = ii*sum(Dj0_t(jTrue)) - gamma0*time

      if(iTime == 0 .or. iTime == nStepsLocal-1) then
        multFact = 1.0_dp
      else if(mod(iTime,2) == 0) then
        multFact = 2.0_dp
      else 
        multFact = 4.0_dp
      endif


      do ikLocal = 1, nkPerPool
        do iE = 1, nTransitions

          if(order == 0) then
            expPrefactor = matrixElement(1,iE,ikLocal)
          else if(order == 1) then
            expPrefactor = sum(matrixElement(jTrue,iE,ikLocal)*Dj1_t(jTrue))
          endif

          Eif = dE(1,iE,ikLocal)

          expArg = expArg_base + ii*Eif/hbar_atomic*time

          transitionRate(iE,ikLocal) = transitionRate(iE,ikLocal) + Real(multFact*expPrefactor*exp(expArg))
            ! We are doing multiple sums (integral and sum over final states), 
            ! but they are all commutative. Here we add in the contribution 
            ! to the integral at this time step from a given final state. The 
            ! loop over final states adds in the contributions from all final 
            ! states. 

        enddo
      enddo
    enddo

    
    ! Combine results from all processes in pool
    do ikLocal = 1, nkPerPool

      call mpiSumDoubleV(transitionRate(:,ikLocal), intraPoolComm)

    enddo
    

    if(indexInPool == 0) then

      ! Multiply by prefactor for Simpson's integration method 
      ! and prefactor for time-domain integral
      transitionRate(:,:) = transitionRate(:,:)*(dt/3.0_dp)*(2.0_dp/(hbar_atomic*hbar_atomic))/time_atomicToSI

      do ikLocal = 1, nkPerPool

        ikGlobal = ikLocal+ikStart_pool-1

        open(unit=37, file=trim(outputDir)//'transitionRate.'//trim(int2str(iSpin))//"."//trim(int2str(ikGlobal)))

        write(37,'(a)') trim(volumeLine)

        write(37,'("# Total number of transitions, Initial states (bandI, bandF) Format : ''(3i10)''")')
        write(37,'(3i10)') nTransitions, ibi(1), ibi(nTransitions)

        write(37,'("# Initial state, Transition rate Format : ''(i10, f10.5, ES24.15E3)''")')

        do iE = 1, nTransitions
          write(37,'(i10, ES24.14E3)') ibi(iE), transitionRate(iE,ikLocal)
            ! Plotting energy doesn't depend on final band, so
            ! just pick one
        enddo

        write(*, '("  Transition rate of k-point ", i4, " and spin ", i1, " written.")') ikGlobal, iSpin

      enddo
    endif

    return

  end subroutine getAndWriteTransitionRate

!----------------------------------------------------------------------------
  subroutine getTrueIndices(arrSize, mask, countTrue, trueIndices)
    ! Get an array of only the indices where the input mask is true

    implicit none

    ! Input variables:
    integer, intent(in) :: arrSize
      !! Size of the mask array

    logical, intent(in) :: mask(arrSize)
      !! Mask used to determine which indices to use

    ! Output variables:
    integer, intent(out) :: countTrue
      !! Number of true values in mask
    integer, allocatable, intent(out) :: trueIndices(:)
      !! Indices where mask is true

    ! Local variables:
    integer :: j, jTrue
      !! Loop indices

    countTrue = count(mask)

    allocate(trueIndices(countTrue))

    jTrue = 0
    do j = 1, arrSize
      if(mask(j)) then
        jTrue = jTrue + 1
        trueIndices(jTrue) = j
      endif
    enddo

    return

  end subroutine

!----------------------------------------------------------------------------
  subroutine setupTimeTables(countTrue, nModes, jTrue, nj, omega, omegaPrime, Sj, SjPrime, time, diffOmega, Dj0_t, Dj1_t)

    implicit none

    ! Input variables:
    integer, intent(in) :: countTrue
      !! Number of indices where Sj > SjThresh
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jTrue(countTrue)
      !! Mode indices where Sj > SjThresh

    real(kind=dp), intent(in) :: nj(nModes)
      !! \(n_j\) occupation number
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: Sj(nModes), SjPrime(nModes)
      !! Huang-Rhys factor for each mode
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    ! Output variables:
    complex(kind=dp), intent(out) :: Dj0_t(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
    complex(kind=dp), intent(out) :: Dj1_t(nModes)
      !! Factor multiplying |M_j|^2*G_j^0(t) in G_j^1(t)

    ! Local variables:
    real(kind=dp) :: cosOmegaPrime(nModes)
      !! cos(omega' t/2)
    real(kind=dp) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)

    complex(kind=dp) :: Aj_t(nModes)
      !! Aj from text. See equation below.
    complex(kind=dp) :: Dj0OverSinOmegaPrime_t(nModes)
      !! Needed to keep from getting NaNs from cot()
    complex(kind=dp) :: expTimesNBarPlus1(nModes)
      !! Local storage of \(e^{i\omega_j t}(\bar{n}_j + 1)\) 
    complex(kind=dp) :: njOverPosExp_t(nModes)
      !! n_j*e^{-i\omega_j t}
    complex(kind=dp) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed


    posExp_t(jTrue) = exp(ii*omega(jTrue)*time)


    expTimesNBarPlus1(jTrue) = posExp_t(jTrue)*(nj(jTrue)+1.0_dp)

    if(diffOmega) then
      Aj_t(jTrue) = (expTimesNBarPlus1(jTrue) + nj(jTrue))/(expTimesNBarPlus1(jTrue) - nj(jTrue))

      sinOmegaPrime(jTrue) = sin(omegaPrime(jTrue)*time/2.0_dp)
      cosOmegaPrime(jTrue) = cos(omegaPrime(jTrue)*time/2.0_dp)

      Dj0OverSinOmegaPrime_t(jTrue) = -2.0_dp/(ii*Aj_t*sinOmegaPrime(jTrue)/Sj(jTrue) - cosOmegaPrime(jTrue)/SjPrime(jTrue))

      Dj0_t(jTrue) = sinOmegaPrime(jTrue)*Dj0OverSinOmegaPrime_t(jTrue)
        ! Need to factor out the sin() to avoid getting NaNs
        ! when calculating cot() here and in Dj1_t

      if(order == 1) then

        Dj1_t(jTrue) = -(hbar_atomic/(2.0_dp*omega(jTrue)*Sj(jTrue)))* &
                        Dj0OverSinOmegaPrime_t(jTrue)*(sinOmegaPrime(jTrue)*Dj0_t(jTrue)*Aj_t(jTrue)**2 - &
                        0.5_dp*omega(jTrue)/omegaPrime(jTrue)*(Aj_t(jTrue)*cosOmegaPrime(jTrue) - sinOmegaPrime(jTrue)* &
                          (omega(jTrue)*cosOmegaPrime(jTrue) - ii*omegaPrime(jTrue)*Aj_t(jTrue)*sinOmegaPrime(jTrue))/ &
                          (omega(jTrue)*Aj_t(jTrue)*sinOmegaPrime(jTrue) + ii*omegaPrime(jTrue)*cosOmegaPrime(jTrue))))
          ! I don't have access to (Delta q_j) here, and I don't want to 
          ! get another variable to deal with. Instead, I rearranged this
          ! to not be in terms of (Delta q_j). Also have to rearrange to 
          ! get rid of cot()

      endif

    else
      njOverPosExp_t(jTrue) = nj(jTrue)/posExp_t(jTrue)

      Dj0_t(jTrue) = Sj(jTrue)/ii*(expTimesNBarPlus1(jTrue) + njOverPosExp_t(jTrue) - (2.0_dp*nj(jTrue) + 1.0_dp))

      if(order == 1) then

        Dj1_t(jTrue) = (hbar_atomic/omega(jTrue))/2.0_dp*(njOverPosExp_t(jTrue) + expTimesNBarPlus1(jTrue) + &
            Sj(jTrue)*(1 + njOverPosExp_t(jTrue) - expTimesNBarPlus1(jTrue))**2)

      endif
    endif

    return

  end subroutine
  
!----------------------------------------------------------------------------
  function transitionRateFileExists(ikGlobal, isp) result(fileExists)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Output variables:
    logical :: fileExists
      !! If the overlap file exists for the given 
      !! k-point and spin channel


    inquire(file=trim(outputDir)//'transitionRate.'//trim(int2str(isp))//"."//trim(int2str(ikGlobal)), exist=fileExists)
    
  end function transitionRateFileExists

end module LSFmod
