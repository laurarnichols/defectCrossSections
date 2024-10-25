module LSFmod
  
  use constants, only: dp, HartreeToEv, ii, hbar_atomic, time_atomicToSI, BohrToMeter
  use base, only: nKPoints, order
  use TMEmod, only: getMatrixElementFNameWPath, getMatrixElementFName, readSingleKMatrixElements
  use PhononPPMod, only: diffOmega, readSjOneFreq, readSjTwoFreq, omega, omegaPrime, readNj
  use energyTabulatorMod, only: energyTableDir, readCaptureEnergyTable, readScatterEnergyTable, readDEPlot
  use miscUtilities, only: int2strLeadZero, int2str
  use errorsAndMPI


  implicit none 

  integer :: iTime_start, iTime_end
    !! Start and end time steps for this process
    !! for the transition-rate integration
  integer :: nStepsLocal_transRate
    !! Number of time steps for each
    !! process to complete


  integer, allocatable :: ibi(:), ibf(:), iki(:), ikf(:)
    !! State indices
  integer :: iSpin
    !! Spin channel to use
  integer, allocatable :: jReSort(:)
    !! Indices to optionally resort matrix elements
  integer :: mDim
    !! Size of first dimension for matrix element
  integer :: nModes
    !! Number of phonon modes
  integer :: nRealTimeSteps
    !! Number of real-time steps for updating occupations
    !! if applicable
  integer :: nTransitions
    !! Total number of transitions 
  integer :: suffixLength
    !! Length of shifted POSCAR file suffix

  real(kind=dp), allocatable :: dE(:,:,:)
    !! All energy differences from energy table
  real(kind=dp), allocatable :: deltaNjInitApproach(:,:)
    !! Optional delta nj from adjustment to 
    !! carrier approach in initial state
  real(kind=dp) :: dt
    !! Optional real time step for generating 
    !! new phonon occupations
  real(kind=dp) :: dtau
    !! Time step size for time-domain integration
  real(kind=dp) :: energyAvgWindow
    !! Size of window to average over for carrier 
    !! density evaluation
  real(kind=dp) :: gamma0
    !! \(\gamma\) for Lorentzian smearing
  real(kind=dp) :: hbarGamma
    !! \(\hbar\gamma\) for Lorentzian smearing
    !! to guarantee convergence
  real(kind=dp), allocatable :: matrixElement(:,:,:)
    !! Electronic matrix element
  real(kind=dp) :: maxTime_transRate
    !! Max time for integration
  real(kind=dp), allocatable :: njBase(:)
    !! Base \(n_j\) occupation number for all states
  real(kind=dp), allocatable :: Sj(:,:), SjPrime(:,:)
    !! Huang-Rhys factor for each mode (and transition
    !! for scattering)
  real(kind=dp) :: SjThresh
    !! Threshold for Sj to determine which modes to calculate
  real(kind=dp) :: smearingExpTolerance
    !! Tolerance for the Lorentzian-smearing
    !! exponential used to calculate max time
  real(kind=dp), allocatable :: totalDeltaNj(:,:)
    !! Optional total change in occupation numbers
    !! for each mode and transition
  real(kind=dp), allocatable :: transitionRate(:,:)
    !! \(Gamma_i\) transition rate

  logical :: addDeltaNj
    !! Add change in occupations for different scattering states
  logical :: captured
    !! If carrier is captured as opposed to scattered
  logical :: generateNewOccupations
    !! If new occupation numbers should be calculated based on
    !! summing over initial and final scattering states
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
  logical :: thermalize
    !! If should convert energy transfer to local temp. or
    !! channel directly into modes
  logical :: writeEiRate
    !! If the energy transfer rate should be output as a 
    !! function of energy

  character(len=300) :: carrierDensityInput
    !! Path to carrier density if generating 
    !! new phonon occupations
  character(len=300) :: deltaNjBaseDir
    !! Path to base directory for deltaNj files
  character(len=300) :: dqInput
    !! Input file for dq.txt if rereading
  character(len=300) :: EiRateOutDir
    !! Output directory for energy transfer rate by 
    !! initial-state energy if writeEiRate == .true.
  character(len=300) :: matrixElementDir
    !! Path to matrix element file `allElecOverlap.isp.ik`. 
    !! For first-order term, the path is just within each 
    !! subdirectory.
  character(len=300) :: MjBaseDir
    !! Path to the base directory for the first-order
    !! matrix element calculations
  character(len=300) :: njBaseInput
    !! Path to base nj file
  character(len=300) :: njNewOutDir
    !! Path to output new occupations if applicable
  character(len=300) :: optimalPairsInput
    !! Path to get optimalPairs.out
  character(len=300) :: prefix
    !! Prefix of directories for first-order matrix
    !! elements
  character(len=300) :: SjBaseDir
    !! Path to directory holding Sj.out file(s)
  character(len=300) :: transRateOutDir
    !! Path to output transition rates
  character(len=300) :: volumeLine
    !! Volume line from overlap file to be
    !! output exactly in transition rate file


  namelist /inputParams/ energyTableDir, matrixElementDir, MjBaseDir, SjBaseDir, njBaseInput, hbarGamma, dtau, &
                         smearingExpTolerance, transRateOutDir, order, prefix, iSpin, diffOmega, newEnergyTable, &
                         suffixLength, reSortMEs, oldFormat, rereadDq, SjThresh, captured, addDeltaNj, &
                         optimalPairsInput, dqInput, deltaNjBaseDir, generateNewOccupations, dt, carrierDensityInput, &
                         energyAvgWindow, njNewOutDir, nRealTimeSteps, thermalize, writeEiRate, EiRateOutDir

contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iSpin, nRealTimeSteps, order, dt, dtau, energyAvgWindow, gamma0, hbarGamma, maxTime_transRate, &
        SjThresh, smearingExpTolerance, addDeltaNj, captured, diffOmega, generateNewOccupations, newEnergyTable, oldFormat, &
        rereadDq, reSortMEs, thermalize, writeEiRate, carrierDensityInput, deltaNjBaseDir, dqInput, energyTableDir, &
        EiRateOutDir, matrixElementDir, MjBaseDir, njBaseInput, njNewOutDir, optimalPairsInput, prefix, SjBaseDir, transRateOutDir)

    implicit none

    ! Output variables
    integer, intent(out) :: iSpin
      !! Spin channel to use
    integer, intent(out) :: nRealTimeSteps
      !! Number of real-time steps for updating occupations
      !! if applicable
    integer, intent(out) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(out) :: dt
      !! Optional real time step for generating 
      !! new phonon occupations
    real(kind=dp), intent(out) :: dtau
      !! Time step size for time-domain integration
    real(kind=dp), intent(out) :: energyAvgWindow
      !! Size of window to average over for carrier 
      !! density evaluation
    real(kind=dp), intent(out) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(out) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: maxTime_transRate
      !! Max time for integration
    real(kind=dp), intent(out) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    
    logical, intent(out) :: addDeltaNj
      !! Add change in occupations for different scattering states
    logical, intent(out) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(out) :: generateNewOccupations
      !! If new occupation numbers should be calculated based on
      !! summing over initial and final scattering states
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
    logical, intent(out) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes
    logical, intent(out) :: writeEiRate
      !! If the energy transfer rate should be output as a 
      !! function of energy

    character(len=300), intent(out) :: carrierDensityInput
      !! Path to carrier density if generating 
      !! new phonon occupations
    character(len=300), intent(out) :: deltaNjBaseDir
      !! Path to base directory for deltaNj files
    character(len=300), intent(out) :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(out) :: EiRateOutDir
      !! Output directory for energy transfer rate by 
      !! initial-state energy if writeEiRate == .true.
    character(len=300), intent(out) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: njBaseInput
      !! Path to base nj file
    character(len=300), intent(out) :: njNewOutDir
      !! Path to output new occupations if applicable
    character(len=300), intent(out) :: optimalPairsInput
      !! Path to get optimalPairs.out
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(out) :: SjBaseDir
      !! Path to directory holding Sj.out file(s)
    character(len=300), intent(out) :: transRateOutDir
      !! Path to store transition rates

  
    call initialize(iSpin, nRealTimeSteps, order, dt, dtau, energyAvgWindow, hbarGamma, SjThresh, smearingExpTolerance, &
          addDeltaNj, captured, diffOmega, generateNewOccupations, newEnergyTable, oldFormat, rereadDq, reSortMEs, &
          thermalize, writeEiRate, carrierDensityInput, deltaNjBaseDir, dqInput, energyTableDir, EiRateOutDir, &
          matrixElementDir, MjBaseDir, njBaseInput, njNewOutDir, optimalPairsInput, prefix, SjBaseDir, transRateOutDir)

    if(ionode) then

      read(5, inputParams, iostat=ierr)
        !! * Read input variables

    
      if(ierr /= 0) call exitError('LSF module', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(iSpin, nRealTimeSteps, order, dt, dtau, energyAvgWindow, hbarGamma, SjThresh, smearingExpTolerance, &
            addDeltaNj, captured, diffOmega, generateNewOccupations, newEnergyTable, oldFormat, rereadDq, reSortMEs, &
            thermalize, writeEiRate, carrierDensityInput, deltaNjBaseDir, dqInput, energyTableDir, EiRateOutDir, &
            matrixElementDir, MjBaseDir, njBaseInput, njNewOutDir, optimalPairsInput, prefix, SjBaseDir, transRateOutDir)

      gamma0 = hbarGamma*1e-3/HartreeToEv
        ! Input expected in meV

      maxTime_transRate = -log(smearingExpTolerance)/gamma0
      write(*,'("Max time for transition-rate integration: ", ES24.15E3)') maxTime_transRate

    endif

    call MPI_BCAST(iSpin, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
      ! nKPoints is not meaningful for scattering
    call MPI_BCAST(nRealTimeSteps, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)
  
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(dtau, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(energyAvgWindow, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(hbarGamma, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(gamma0, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(SjThresh, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(smearingExpTolerance, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(maxTime_transRate, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(addDeltaNj, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(captured, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(diffOmega, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(generateNewOccupations, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(newEnergyTable, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(oldFormat, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(reSortMEs, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(rereadDq, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(thermalize, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(writeEiRate, 1, MPI_LOGICAL, root, worldComm, ierr)
  
    call MPI_BCAST(carrierDensityInput, len(carrierDensityInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(deltaNjBaseDir, len(deltaNjBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(dqInput, len(dqInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(EiRateOutDir, len(EiRateOutDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(matrixElementDir, len(matrixElementDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(MjBaseDir, len(MjBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(njBaseInput, len(njBaseInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(njNewOutDir, len(njNewOutDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(optimalPairsInput, len(optimalPairsInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(SjBaseDir, len(SjBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(transRateOutDir, len(transRateOutDir), MPI_CHARACTER, root, worldComm, ierr)
    
    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(iSpin, nRealTimeSteps, order, dt, dtau, energyAvgWindow, hbarGamma, SjThresh, smearingExpTolerance, &
        addDeltaNj, captured, diffOmega, generateNewOccupations, newEnergyTable, oldFormat, rereadDq, reSortMEs, &
        thermalize, writeEiRate, carrierDensityInput, deltaNjBaseDir, dqInput, energyTableDir, EiRateOutDir, &
        matrixElementDir, MjBaseDir, njBaseInput, njNewOutDir, optimalPairsInput, prefix, SjBaseDir, transRateOutDir)

    implicit none

    ! Output variables
    integer, intent(out) :: iSpin
      !! Spin channel to use
    integer, intent(out) :: nRealTimeSteps
      !! Number of real-time steps for updating occupations
      !! if applicable
    integer, intent(out) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(out) :: dt
      !! Optional real time step for generating 
      !! new phonon occupations
    real(kind=dp), intent(out) :: dtau
      !! Time step size for time-domain integration
    real(kind=dp), intent(out) :: energyAvgWindow
      !! Size of window to average over for carrier 
      !! density evaluation
    real(kind=dp), intent(out) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    
    logical, intent(out) :: addDeltaNj
      !! Add change in occupations for different scattering states
    logical, intent(out) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(out) :: generateNewOccupations
      !! If new occupation numbers should be calculated based on
      !! summing over initial and final scattering states
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
    logical, intent(out) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes
    logical, intent(out) :: writeEiRate
      !! If the energy transfer rate should be output as a 
      !! function of energy

    character(len=300), intent(out) :: carrierDensityInput
      !! Path to carrier density if generating 
      !! new phonon occupations
    character(len=300), intent(out) :: deltaNjBaseDir
      !! Path to base directory for deltaNj files
    character(len=300), intent(out) :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(out) :: EiRateOutDir
      !! Output directory for energy transfer rate by 
      !! initial-state energy if writeEiRate == .true.
    character(len=300), intent(out) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: njBaseInput
      !! Path to base nj file
    character(len=300), intent(out) :: njNewOutDir
      !! Path to output new occupations if applicable
    character(len=300), intent(out) :: optimalPairsInput
      !! Path to get optimalPairs.out
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(out) :: SjBaseDir
      !! Path to directory holding Sj.out file(s)
    character(len=300), intent(out) :: transRateOutDir
      !! Path to store transition rates


    iSpin = 1
    nRealTimeSteps = 1
    order = -1

    dt = 1d-4
    dtau = 1d-4
    energyAvgWindow = 1e-2
    hbarGamma = 0.0_dp
    SjThresh = 0.0_dp
    smearingExpTolerance = 0.0_dp

    addDeltaNj = .false.
    captured = .true.
    diffOmega = .false.
    generateNewOccupations = .false.
    newEnergyTable = .false.
    oldFormat = .false.
    reSortMEs = .false.
    rereadDq = .false.
    thermalize = .false.
    writeEiRate = .false.

    carrierDensityInput = ''
    deltaNjBaseDir = ''
    dqInput = ''
    energyTableDir = ''
    EiRateOutDir = './EiRates'
    matrixElementDir = ''
    MjBaseDir = ''
    njBaseInput = ''
    njNewOutDir = './njNew'
    optimalPairsInput = ''
    prefix = 'disp-'
    SjBaseDir = ''
    transRateOutDir = './transitionRates'

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iSpin, nRealTimeSteps, order, dt, dtau, energyAvgWindow, hbarGamma, SjThresh, &
        smearingExpTolerance, addDeltaNj, captured, diffOmega, generateNewOccupations, newEnergyTable, oldFormat, &
        rereadDq, reSortMEs, thermalize, writeEiRate, carrierDensityInput, deltaNjBaseDir, dqInput, energyTableDir, &
        EiRateOutDir, matrixElementDir, MjBaseDir, njBaseInput, njNewOutDir, optimalPairsInput, prefix, SjBaseDir, &
        transRateOutDir)

    implicit none

    ! Input variables
    integer, intent(in) :: iSpin
      !! Spin channel to use
    integer, intent(in) :: nRealTimeSteps
      !! Number of real-time steps for updating occupations
      !! if applicable
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dt
      !! Optional real time step for generating 
      !! new phonon occupations
    real(kind=dp), intent(in) :: dtau
      !! Time step size for time-domain integration
    real(kind=dp), intent(in) :: energyAvgWindow
      !! Size of window to average over for carrier 
      !! density evaluation
    real(kind=dp), intent(in) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(in) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    
    logical, intent(in) :: addDeltaNj
      !! Add change in occupations for different scattering states
    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(in) :: generateNewOccupations
      !! If new occupation numbers should be calculated based on
      !! summing over initial and final scattering states
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
    logical, intent(in) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes
    logical, intent(in) :: writeEiRate
      !! If the energy transfer rate should be output as a 
      !! function of energy

    character(len=300), intent(in) :: carrierDensityInput
      !! Path to carrier density if generating 
      !! new phonon occupations
    character(len=300), intent(in) :: deltaNjBaseDir
      !! Path to base directory for deltaNj files
    character(len=300), intent(in) :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(in) :: EiRateOutDir
      !! Output directory for energy transfer rate by 
      !! initial-state energy if writeEiRate == .true.
    character(len=300), intent(in) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(in) :: njBaseInput
      !! Path to base nj file
    character(len=300), intent(in) :: njNewOutDir
      !! Path to output new occupations if applicable
    character(len=300), intent(in) :: optimalPairsInput
      !! Path to get optimalPairs.out
    character(len=300), intent(in) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(in) :: SjBaseDir
      !! Path to directory holding Sj.out file(s)
    character(len=300), intent(in) :: transRateOutDir
      !! Path to store transition rates

    ! Local variables:
    integer :: ikTest
      !! Index to pass to getMatrixElementFName
    character(len=300) :: fName
      !! File name for matrix element file to get
      !! nKPoints from 

    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iSpin', iSpin, 1, 2)
    abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution 

    abortExecution = checkDoubleInitialization('dtau', dtau, 1.0d-6, 1.0d-2) .or. abortExecution
    abortExecution = checkDoubleInitialization('hbarGamma', hbarGamma, 0.1_dp, 20.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('SjThresh', SjThresh, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('smearingExpTolerance', smearingExpTolerance, 0.0_dp, 1.0_dp) .or. abortExecution
      ! These limits are my best guess as to what is reasonable; they are not
      ! hard and fast, but you should think about the application of the theory
      ! to numbers outside these ranges.

    abortExecution = checkFileInitialization('njBaseInput', njBaseInput) .or. abortExecution

    
    write(*,'("captured = ",L)') captured
    write(*,'("addDeltaNj = ",L)') addDeltaNj
    write(*,'("generateNewOccupations = ",L)') generateNewOccupations

    if(addDeltaNj .or. generateNewOccupations) then
      if(captured) &
        call exitError('checkInitialization','Can only treat change in occupations for different states with scattering!',1)

      write(*,'("deltaNjBaseDir = ''",a,"''")') trim(deltaNjBaseDir)
      if(trim(deltaNjBaseDir) == '') then
        abortExecution = .true.
        write(*,'("Must have deltaNjBaseDir for addDeltaNj or generateNewOccupations!!")')
      endif

      if(generateNewOccupations) then
        abortExecution = checkDoubleInitialization('dt', dt, 0.0_dp, 10.0_dp) .or. abortExecution
        abortExecution = checkIntInitialization('nRealTimeSteps', nRealTimeSteps, 1, int(1e9))
        abortExecution = checkFileInitialization('carrierDensityInput', carrierDensityInput) .or. abortExecution
        abortExecution = checkDoubleInitialization('energyAvgWindow', energyAvgWindow, 0.0_dp, 1.0_dp) .or. abortExecution
        abortExecution = checkFileInitialization('njBaseInput', njBaseInput) .or. abortExecution

        write(*,'("njNewOutDir = ''",a,"''")') trim(njNewOutDir)
        call system('mkdir -p '//trim(njNewOutDir))

        write(*,'("thermalize = ",L)') thermalize

        write(*,'("EiRateOutDir = ''",a,"''")') trim(EiRateOutDir)
        call system('mkdir -p '//trim(EiRateOutDir))
      endif

    endif


    ! We don't know the band indices here to check for the Sj files for
    ! scattering, but we can check for Sj.analysis.out that should always 
    ! be there when calculating Sj's as well.
    if(captured) then
      abortExecution = checkDirInitialization('SjBaseDir', SjBaseDir, 'Sj.out') .or. abortExecution
    else
      abortExecution = checkDirInitialization('SjBaseDir', SjBaseDir, 'Sj.analysis.out') .or. abortExecution
    endif


    ! Check for the energy table file for the first k-point for
    ! capture and without the k-point index for scattering
    if(captured) then
      abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.'//trim(int2str(iSpin))//'.1') .or. abortExecution
    else
      abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.'//trim(int2str(iSpin))) .or. abortExecution
    endif


    write(*,'("diffOmega = ",L)') diffOmega
    write(*,'("oldFormat = ",L)') oldFormat
    write(*,'("newEnergyTable = ",L)') newEnergyTable


    ! Below, we want to test the presence of the required matrix
    ! element files. The only difference between the capture and
    ! scattering case is that the capture case has a k-point index.
    ! For capture, we look for the presence of the k-point 1 file,
    ! and for scattering we pass -1 to getMatrixElementFName to 
    ! remove the k-point index from the file name.
    if(captured) then
      ikTest = 1
    else
      ikTest = -1
    endif


    if(order == 0) then 
      abortExecution = checkDirInitialization('matrixElementDir', matrixElementDir, getMatrixElementFName(ikTest,iSpin)) .or. abortExecution

    else if(order == 1) then

      write(*,'("reSortMEs = ",L)') reSortMEs
      if(reSortMEs) abortExecution = checkFileInitialization('optimalPairsInput', optimalPairsInput) .or. abortExecution

      write(*,'("rereadDq = ",L)') rereadDq
      if(rereadDq) abortExecution = checkFileInitialization('dqInput', dqInput) .or. abortExecution

      abortExecution = checkIntInitialization('suffixLength', suffixLength, 1, 5) .or. abortExecution 
      abortExecution = checkDirInitialization('MjBaseDir', MjBaseDir, &
            '/'//trim(prefix)//trim(int2strLeadZero(1,suffixLength))//'/'//trim(getMatrixElementFNameWPath(ikTest,iSpin,matrixElementDir))) .or. abortExecution
      write(*,'("prefix = ''",a,"''")') trim(prefix)
      write(*,'("matrixElementDir = ''",a,"''")') trim(matrixElementDir)

    endif


    ! For capture, we need to know nKPoints, so read this from the
    ! matrix element file
    if(captured) then
      if(order == 0) then
        fName = getMatrixElementFNameWPath(1,iSpin,matrixElementDir)
      else if(order == 1) then
        fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(1,suffixLength))//'/'//&
                trim(getMatrixElementFNameWPath(1,iSpin,matrixElementDir))
      endif

      open(unit=12,file=trim(fName))
      read(12,*)
      read(12,*) nKPoints
      close(12)

      write(*,'("nKPoints = ", i10)') nKPoints
    endif


    write(*,'("transRateOutDir = ''",a,"''")') trim(transRateOutDir)
    call system('mkdir -p '//trim(transRateOutDir))


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif
    write(*,*) 
    write(*,*) 
    write(*,*) 
    write(*,*) 

    return 

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readEnergyTable(iSpin, captured, energyTableDir, nTransitions, ibi, ibf, iki, ikf, dE)

    implicit none

    ! Input variables:
    integer, intent(in) :: iSpin
      !! Spin channel to use

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to read

    ! Output variables:
    integer, intent(out) :: nTransitions
      !! Total number of transitions 
    integer, allocatable, intent(out) :: ibi(:), ibf(:), iki(:), ikf(:)
      !! State indices

    real(kind=dp), allocatable, intent(out) :: dE(:,:,:)
      !! All energy differences from energy table

    ! Local variables:
    integer :: ibf1
      !! Scalar integer for reading single final state
      !! from capture energy table
    integer, allocatable :: iDum1D(:)
      !! Integer to ignore input
    integer :: iDum1, iDum2
      !! Integers to ignore input
    integer :: ikLocal, ikGlobal
      !! Loop indices

    real(kind=dp), allocatable :: dE2D(:,:)
      !! 2D array to read energy at single k-point


    ! For capture, assume that the band bounds and number of transitions do not depend
    ! on k-points or spin, so read those from the first k-point and ignore those inputs
    ! from different files.
    if(captured) then

      allocate(ibf(1),iki(1),ikf(1))

      ! Get the number of transitions and the state indices
      if(ionode) then
        call readCaptureEnergyTable(1, iSpin, energyTableDir, ibi, ibf1, nTransitions, dE2D)
         ! Assume that band bounds and number of transitions do not depend on k-points or spin
         ! We ignore the energy here because we are only reading ibi, ibf, and nTransitions

        deallocate(dE2D)

        ibf(1) = ibf1
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

    ! Scattering should only have one input file because the file does not depend on 
    ! the k-point and the user selects a single spin, so read from the scatter energy
    ! table and get all of the information (nTransitions, state indices, energy) in a
    ! single shot.
    else
      if(ionode) &
        call readScatterEnergyTable(iSpin, .false., energyTableDir, ibi, ibf, iki, ikf, nTransitions, dE2D)
          ! dE2D will get allocated here with (5,nTransitions) and will read all
          ! energies from the table


      call MPI_BCAST(nTransitions, 1, MPI_INTEGER, root, worldComm, ierr)

      ! For scattering, all state indices are indexed with nTransitions
      if(.not. ionode) allocate(ibi(nTransitions), ibf(nTransitions),iki(nTransitions),ikf(nTransitions))
      call MPI_BCAST(ibi, nTransitions, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(ibf, nTransitions, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(iki, nTransitions, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(ikf, nTransitions, MPI_INTEGER, root, worldComm, ierr)


      allocate(dE(5,nTransitions,1))
        ! dE is expected to be 3D here and in later subroutines to 
        ! be consistent with capture. The last index is the number
        ! of k-points per pool, but scattering currently assumes
        ! nkPerPool = 1.

      if(ionode) then
        dE(:,:,1) = dE2D
        deallocate(dE2D)
      endif

      call MPI_BCAST(dE, size(dE), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    endif

    return

  end subroutine readEnergyTable

!----------------------------------------------------------------------------
  subroutine readSj(ibi, ibf, iki, ikf, nTransitions, captured, diffOmega, SjBaseDir, nModes, omega, &
          omegaPrime, Sj, SjPrime)

    implicit none

    ! Input variables:
    integer, intent(in) :: ibi(:), ibf(:), iki(:), ikf(:)
      !! State indices
    integer, intent(in) :: nTransitions
      !! Total number of transitions 

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(in) :: SjBaseDir
      !! Path to directory holding Sj.out file(s)

    ! Output variables:
    integer, intent(out) :: nModes
      !! Number of phonon modes

    real(kind=dp),allocatable, intent(out) :: omega(:), omegaPrime(:)
      !! Frequency for each mode
    real(kind=dp), allocatable, intent(out) :: Sj(:,:), SjPrime(:,:)
      !! Huang-Rhys factor for each mode (and transition
      !! for scattering)

    ! Local variables:
    integer :: iDum
      !! Dummy integer to ignore input
    integer :: iE
      !! Loop index

    !real(kind=dp), allocatable :: randVal(:)
      ! Random adjustment to be made to frequencies
      ! to test sensitivity
    real(kind=dp), allocatable :: rDum1D_1(:), rDum1D_2(:)
      !! Dummy variables to ignore input
    real(kind=dp), allocatable :: Sj1D(:), SjPrime1D(:)
      !! 1D arrays to pass to readSjOneFreq and readSjTwoFreq

    character(len=300) :: fName
      !! File name to read


    ! Check that the arrays have the expected size
    if(captured) then
      if(.not. (size(iki) == 1 .and. size(ibf) == 1 .and. size(ikf) == 1)) &
          call exitError('readSj','For capture, the iki, ibf, and ikf arrays should have length 1!', 1)
    endif

    ! Go through reading one file for capture and scattering in the
    ! same way but with different file names. This will get the
    ! Sj/Sj' and omega/omega'. For scattering, we only save the 
    ! frequencies and number of modes from the first file because 
    ! they do not depend on the transition
    if(captured) then
      fName = trim(SjBaseDir)//'/Sj.out' 
    else
      fName = trim(SjBaseDir)//'/Sj.k'//trim(int2str(iki(1)))//'_b'//trim(int2str(ibi(1)))//'.k'&
                               //trim(int2str(ikf(1)))//'_b'//trim(int2str(ibf(1)))//'.out'
    endif

    if(diffOmega) then
      call readSjTwoFreq(fName, nModes, omega, omegaPrime, Sj1D, SjPrime1D)

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
      allocate(omegaPrime(nModes),SjPrime1D(nModes))
        ! Need to allocate to avoid issues with passing variables
        ! and deallocating

      call readSjOneFreq(fName, nModes, omega, Sj1D)
    endif


    ! Store in full Sj/Sj' arrays with second index for optional
    ! transitions
    if(captured) then
      allocate(Sj(nModes,1),SjPrime(nModes,1))
    else
      allocate(Sj(nModes,nTransitions),SjPrime(nModes,nTransitions))
    endif

    Sj(:,1) = Sj1D
    if(diffOmega) SjPrime(:,1) = SjPrime1D
    deallocate(Sj1D,SjPrime1D)


    ! For scattering, read the Sj for the rest of the transitions
    if(.not. captured) then
      do iE = 2, nTransitions
        fName = trim(SjBaseDir)//'/Sj.k'//trim(int2str(iki(iE)))//'_b'//trim(int2str(ibi(iE)))//'.k'&
                                 //trim(int2str(ikf(iE)))//'_b'//trim(int2str(ibf(iE)))//'.out'

        if(diffOmega) then
          call readSjTwoFreq(fName, iDum, rDum1D_1, rDum1D_2, Sj1D, SjPrime1D)
          deallocate(rDum1D_1,rDum1D_2)
        else
          allocate(SjPrime1D(nModes))
          call readSjOneFreq(fName, iDum, rDum1D_1, Sj1D)
          deallocate(rDum1D_1)
        endif

        Sj(:,iE) = Sj1D
        if(diffOmega) SjPrime(:,iE) = SjPrime1D
        deallocate(Sj1D,SjPrime1D)

      enddo
    endif

    return

  end subroutine readSj

!----------------------------------------------------------------------------
  subroutine getjReSort(nModes, optimalPairsInput, jReSort)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of phonon modes

    character(len=300), intent(in) :: optimalPairsInput
      !! Path to get optimalPairs.out

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


    open(unit=32,file=trim(optimalPairsInput))

    read(32,*)

    do j = 1, nModes
      read(32,*) jCurrent, jNew
      jReSort(jCurrent) = jNew
    enddo

    close(32)

    return

  end subroutine

!----------------------------------------------------------------------------
  subroutine readAllMatrixElements(iSpin, nTransitions, ibi, nModes, jReSort, order, suffixLength, dE, captured, newEnergyTable, &
          oldFormat, rereadDq, reSortMEs, dqInput, matrixElementDir, MjBaseDir, prefix, mDim, matrixElement, volumeLine)
    ! For both the zeroth-order and first-order matrix elements, there is the option 
    ! to specify `newEnergyTable`. If this variable is true, that means that the 
    ! matrix elements (TME) and LSF code were not run with the same energy table as 
    ! input. This might be beneficial if tweaks were made to the energies used or 
    ! you only want to use a subset of the transitions for the LSF code. 
    !
    ! With `newEnergyTable = .true.`, the states to be read from the matrix element 
    ! file(s) come from the energy table input to the LSF code. All of the states in 
    ! the energy table input to the LSF code must be present in the matrix element file. 
    !
    ! Both terms also allow you to specify the old format where both the initial and 
    ! final band indices were given for capture and the headers were different. I 
    ! needed this for the Si calculations and left it here just in case.
    !
    ! In addition to the options `newEnergyTable` and `oldFormat` in the zeroth-order 
    ! term, the first-order term also has options `reSortMEs` and `rereadDq`.
    !
    ! `reSortMEs` assumes that the modes were labeled differently when the matrix 
    ! elements were calculated vs what is to be used for the LSF run. This can be the 
    ! case if you want to compare the difference between allowing the frequencies to 
    ! change for different states (e.g., before and after capture) vs not. 
    !
    ! `rereadDq` allows you to read a new delta q for each mode than was used in the 
    ! calculation of the matrix elements. You shouldn't need this in general, but 
    ! there was a bug in the code that only messed up the dq.txt file previously, so 
    ! it was useful for me to be able to fix the dq here rather than having to redo 
    ! all of the matrix element calculations. 
    !
    ! I have many layers of encapsulation here because otherwise I would have been 
    ! repeating a lot of code. It was also really hard for me to understand without
    ! it, even as the person who wrote it. I know the many layers can also be 
    ! confusing, but I just had to make a choice. Hopefully with the comments you can
    ! follow the layers.

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

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
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

    character(len=300), intent(in) :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(in) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
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
    integer :: ikLocal, ikGlobal
      !! Loop indices


    if(order == 0) then
      mDim = 1
    else if(order == 1) then
      mDim = nModes
    endif

    allocate(matrixElement(mDim,nTransitions,nkPerPool))


    if(captured) then

      if(indexInPool == 0) then

        matrixElement = 0.0_dp

        do ikLocal = 1, nkPerPool
    
          ikGlobal = ikLocal+ikStart_pool-1
            !! Get the global `ik` index from the local one

           call readSingleKMatrixElements(ikGlobal, iSpin, nTransitions, ibi, nModes, jReSort, mDim, order, suffixLength, &
                dE(:,:,ikLocal), captured, newEnergyTable, oldFormat, rereadDq, reSortMEs, dqInput, matrixElementDir, &
                MjBaseDir, prefix, matrixElement(:,:,ikLocal), volumeLine)

        enddo
      endif

      call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)

    else
      if(ionode) then
        matrixElement = 0.0_dp

        call readSingleKMatrixElements(-1, iSpin, nTransitions, ibi, nModes, jReSort, mDim, order, suffixLength, &
              dE(:,:,1), captured, newEnergyTable, oldFormat, rereadDq, reSortMEs, dqInput, matrixElementDir, MjBaseDir, &
              prefix, matrixElement(:,:,1), volumeLine)

      endif

      call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    endif


    return

  end subroutine readAllMatrixElements

!----------------------------------------------------------------------------
  subroutine getAndWriteTransitionRate(nTransitions, ibi, ibf, iki, ikf, iRt, iSpin, mDim, nModes, order, dE, deltaNjInitApproach, &
          dtau, gamma0, matrixElement, njBase, omega, omegaPrime, Sj, SjPrime, SjThresh, addDeltaNj, &
          captured, diffOmega, generateNewOccupations, transRateOutDir, volumeLine, transitionRate)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions), ibf(:), iki(:), ikf(:)
      !! State indices
    integer, intent(in) :: iRt
      !! Loop index
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
    real(kind=dp), intent(in) :: deltaNjInitApproach(:,:)
      !! Optional delta nj from adjustment to 
      !! carrier approach in initial state
    real(kind=dp), intent(in) :: dtau
      !! Time step size for time-domain integration
    real(kind=dp), intent(in) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(in) :: matrixElement(mDim,nTransitions,nkPerPool)
      !! Electronic matrix element
    real(kind=dp), intent(in) :: njBase(nModes)
      !! Base \(n_j\) occupation number for all states
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: Sj(:,:), SjPrime(:,:)
      !! Huang-Rhys factor for each mode (and 
      !! transition for scattering)
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate

    logical, intent(in) :: addDeltaNj
      !! Add change in occupations for different scattering states
    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(in) :: generateNewOccupations
      !! If new occupation numbers should be calculated based on
      !! summing over initial and final scattering states

    character(len=300), intent(in) :: transRateOutDir
      !! Path to store transition rates
    character(len=300), intent(in) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Output variables:
    real(kind=dp), intent(out) :: transitionRate(nTransitions,nkPerPool)
      !! \(Gamma_i\) transition rate

    ! Local variables:
    integer :: countTrue
      !! Number of modes where Sj > SjThresh
    integer :: iTime, iE, ikLocal, ikGlobal
      !! Loop indices
    integer, allocatable :: jTrue(:)
      !! Mode indices where Sj > SjThresh
    integer :: updateFrequency
      !! Frequency of steps to write status update

    real(kind=dp) :: cosOmegaPrime(nModes)
      !! cos(omega' t/2)
    real(kind=dp) :: Eif
      !! Local storage of energy for delta function
    real(kind=dp) :: multFact
      !! Multiplication factor for each term per
      !! Simpson's integration method
    real(kind=dp) :: njPlusDelta1D(nModes)
      !! Optional nj plus delta nj from adjustment to 
      !! carrier approach in initial state. Indexed only
      !! over modes and not stored for every transition.
    real(kind=dp) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)
    real(kind=dp) :: Sj1D(nModes), SjPrime1D(nModes)
      !! 1D arrays to pass to setUpTimeTables
    real(kind=dp) :: t0
      !! Initial time for this process
    real(kind=dp) :: time
      !! Time for each time step
    real(kind=dp) :: timer1, timer2
      !! Timers

    complex(kind=dp) :: Aj_t(nModes)
      !! Aj from text. See equation below.
    complex(kind=dp) :: Dj0_tOverSj(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
      !! divided by Sj
    complex(kind=dp) :: Dj0_t(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
    complex(kind=dp) :: Dj1_t(nModes)
      !! Factor multiplying |M_j|^2*G_j^0(t) in G_j^1(t)
    complex(kind=dp) :: Dj1_termsOneAndTwo(nModes)
      !! First and second terms in the parentheses of D_j^1(t)
    complex(kind=dp) :: Dj1_termThree(nModes)
      !! Third term in the parentheses of D_j^1(t) (without Sj)
    complex(kind=dp) :: expArg_base
      !! Base exponential argument for each time step
    complex(kind=dp) :: expArg
      !! Exponential argument for each time step
    complex(kind=dp) :: expPrefactor
      !! Prefactor for the exponential inside the integral.
      !! For zeroth-order, this is just the matrix element,
      !! but for first-order this is \(\sum_j M_j A_j\).
    complex(kind=dp) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed

    logical :: mask(nModes)
      !! Select the modes to calculate

    character(len=300) :: text
      !! Text for long header

      
    updateFrequency = ceiling(nStepsLocal_transRate/10.0)
    call cpu_time(timer1)


    ! Create a mask to determine which modes to calculate
    ! based on an optional threshold given (default 0.0).
    mask = Sj(:,1) >= SjThresh
    call getTrueIndices(nModes, mask, countTrue, jTrue)


    t0 = indexInPool*float(nStepsLocal_transRate-1)*dtau

    
    transitionRate(:,:) = 0.0_dp
    do iTime = 0, nStepsLocal_transRate-1

      if(ionode .and. mod(iTime,updateFrequency) == 0) then

        call cpu_time(timer2)
        write(*,'("    ", i2,"% complete with transition-rate loop. Time in loop: ",f10.2," secs")') &
              int((iTime*100.0)/nStepsLocal_transRate), timer2-timer1

      endif

      time = t0 + float(iTime)*dtau
        ! Must do this arithmetic with floats to avoid
        ! integer overflow

      if(iTime == 0 .or. iTime == nStepsLocal_transRate-1) then
        multFact = 1.0_dp
      else if(mod(iTime,2) == 0) then
        multFact = 2.0_dp
      else 
        multFact = 4.0_dp
      endif

      if(captured) then
        Sj1D = Sj(:,1)
        SjPrime1D = SjPrime(:,1)
        call setupAllTimeTables(countTrue, nModes, jTrue, njBase, omega, omegaPrime, Sj1D, SjPrime1D, time, diffOmega, Dj0_t, Dj1_t)

        expArg_base = ii*sum(Dj0_t(jTrue)) - gamma0*time
      else 
        ! I had a choice here to have one subroutine called and pass addDeltaNj to 
        ! switch what is calculated, but that would significantly affect the input
        ! and output variables. I don't, in general, like to have different copies
        ! of subroutines, but I felt having two separate subroutines would be
        ! clearer and safer. Same below for setupStateDepTimeTables
        if(addDeltaNj) then
          call setupStateIndTimeTablesDeltaNj(countTrue, nModes, jTrue, omega, omegaPrime, time, diffOmega, cosOmegaPrime, &
                    sinOmegaPrime, posExp_t)
        else
          call setupStateIndTimeTablesNoDeltaNj(countTrue, nModes, jTrue, njBase, omega, omegaPrime, time, diffOmega, &
                    cosOmegaPrime, sinOmegaPrime, Aj_t, Dj0_tOverSj, Dj1_termsOneAndTwo, Dj1_termThree)
        endif

      endif



      ! For scattering, this loop is basically non-functional because the number
      ! of pools is 1 and nkPerPool is also set to 1. That means that anywhere
      ! ikLocal is seen below, it can be replaced with 1 for scattering.
      do ikLocal = 1, nkPerPool

        do iE = 1, nTransitions

          if(.not. captured) then
            Sj1D = Sj(:,iE)
            SjPrime1D = SjPrime(:,iE)

            if(addDeltaNj) then
              njPlusDelta1D = njBase(:) + deltaNjInitApproach(:,iE)
              call setupStateDepTimeTablesDeltaNj(countTrue, nModes, jTrue, cosOmegaPrime, njPlusDelta1D, omega, omegaPrime, &
                    sinOmegaPrime, Sj1D, SjPrime1D, posExp_t, diffOmega, Dj0_t, Dj1_t)
            else
              call setupStateDepTimeTablesNoDeltaNj(countTrue, nModes, jTrue, cosOmegaPrime, omega, omegaPrime, sinOmegaPrime, Sj1D, &
                    SjPrime1D, Aj_t, Dj0_tOverSj, Dj1_termsOneAndTwo, Dj1_termThree, diffOmega, Dj0_t, Dj1_t)
            endif

            expArg_base = ii*sum(Dj0_t(jTrue)) - gamma0*time
          endif

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
      transitionRate(:,:) = transitionRate(:,:)*(dtau/3.0_dp)*(2.0_dp/(hbar_atomic*hbar_atomic))/time_atomicToSI

      do ikLocal = 1, nkPerPool

        if(captured) then
          ikGlobal = ikLocal+ikStart_pool-1
          open(unit=37, file=trim(transRateOutDir)//'/transitionRate.'//trim(int2str(iSpin))//"."//trim(int2str(ikGlobal)))
        else
          if(generateNewOccupations) then
            open(unit=37, file=trim(transRateOutDir)//'/transitionRate.'//trim(int2str(iSpin))//'.'//trim(int2str(iRt)))
          else
            open(unit=37, file=trim(transRateOutDir)//'/transitionRate.'//trim(int2str(iSpin)))
          endif
        endif

        write(37,'(a)') trim(volumeLine)

        if(captured) then
          write(37,'("# Total number of transitions, Initial states (bandI, bandF) Format : ''(3i10)''")')
          write(37,'(3i10)') nTransitions, ibi(1), ibi(nTransitions)

          write(37,'("# Initial state, Transition rate Format : ''(i10, ES24.15E3)''")')
        else
          text = "# Total number of transitions, Initial States (kI, kF, bandI, bandF), Final States (kI, kF, bandI, bandF)"
          write(37,'(a, " Format : ''(9i10)''")') trim(text)   
  
          write(37,'(9i10)') nTransitions, iki(1), iki(nTransitions), ibi(1), ibi(nTransitions), &
                                           ikf(1), ikf(nTransitions), ibf(1), ibf(nTransitions)

          write(37,'("# iki, ibi, ikf, ibf, Transition rate Format : ''(4i10, ES24.15E3)''")')
        endif 


        do iE = 1, nTransitions
          if(captured) then
            write(37,'(i10, ES24.14E3)') ibi(iE), transitionRate(iE,ikLocal)
          else
            write(37,'(4i10, ES24.14E3)') iki(iE), ibi(iE), ikf(iE), ibf(iE), transitionRate(iE,ikLocal)
          endif
        enddo

        if(captured) then
          write(*, '("  Transition rate of k-point ", i4, " and spin ", i1, " written.")') ikGlobal, iSpin
        else
          write(*, '("  Transition rate written.")')
        endif

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
  subroutine setupAllTimeTables(countTrue, nModes, jTrue, nj, omega, omegaPrime, Sj, SjPrime, time, diffOmega, Dj0_t, Dj1_t)

    implicit none

    ! Input variables:
    integer, intent(in) :: countTrue
      !! Number of indices where Sj > SjThresh
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jTrue(countTrue)
      !! Mode indices where Sj > SjThresh

    real(kind=dp), intent(in) :: nj(nModes)
      !! Base \(n_j\) occupation number for all states
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

  end subroutine setupAllTimeTables

!----------------------------------------------------------------------------
  subroutine setupStateIndTimeTablesNoDeltaNj(countTrue, nModes, jTrue, nj, omega, omegaPrime, time, diffOmega, &
            cosOmegaPrime, sinOmegaPrime, Aj_t, Dj0_tOverSj, Dj1_termsOneAndTwo, Dj1_termThree)

    implicit none

    ! Input variables:
    integer, intent(in) :: countTrue
      !! Number of indices where Sj > SjThresh
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jTrue(countTrue)
      !! Mode indices where Sj > SjThresh

    real(kind=dp), intent(in) :: nj(nModes)
      !! Base \(n_j\) occupation number for all states
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    ! Output variables:
    real(kind=dp), intent(out) :: cosOmegaPrime(nModes)
      !! cos(omega' t/2)
    real(kind=dp), intent(out) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)

    complex(kind=dp), intent(out) :: Aj_t(nModes)
      !! Aj from text. See equation below.
    complex(kind=dp), intent(out) :: Dj0_tOverSj(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
      !! divided by Sj
    complex(kind=dp), intent(out) :: Dj1_termsOneAndTwo(nModes)
      !! First and second terms in the parentheses of D_j^1(t)
    complex(kind=dp), intent(out) :: Dj1_termThree(nModes)
      !! Third term in the parentheses of D_j^1(t) (without Sj)

    ! Local variables:
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

    else
      njOverPosExp_t(jTrue) = nj(jTrue)/posExp_t(jTrue)

      Dj0_tOverSj(jTrue) = 1.0_dp/ii*(expTimesNBarPlus1(jTrue) + njOverPosExp_t(jTrue) - (2.0_dp*nj(jTrue) + 1.0_dp))

      if(order == 1) then

        Dj1_termsOneAndTwo(jTrue) = njOverPosExp_t(jTrue) + expTimesNBarPlus1(jTrue)
        Dj1_termThree(jTrue) = (1 + njOverPosExp_t(jTrue) - expTimesNBarPlus1(jTrue))**2

      endif
    endif

    return

  end subroutine setupStateIndTimeTablesNoDeltaNj

!----------------------------------------------------------------------------
  subroutine setupStateDepTimeTablesNoDeltaNj(countTrue, nModes, jTrue, cosOmegaPrime, omega, omegaPrime, sinOmegaPrime, Sj, &
            SjPrime, Aj_t, Dj0_tOverSj, Dj1_termsOneAndTwo, Dj1_termThree, diffOmega, Dj0_t, Dj1_t)

    implicit none

    ! Input variables:
    integer, intent(in) :: countTrue
      !! Number of indices where Sj > SjThresh
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jTrue(countTrue)
      !! Mode indices where Sj > SjThresh

    real(kind=dp), intent(in) :: cosOmegaPrime(nModes)
      !! cos(omega' t/2)
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)
    real(kind=dp), intent(in) :: Sj(nModes), SjPrime(nModes)
      !! Huang-Rhys factor for each mode

    complex(kind=dp), intent(in) :: Aj_t(nModes)
      !! Aj from text. See equation below.
    complex(kind=dp), intent(in) :: Dj0_tOverSj(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
      !! divided by Sj
    complex(kind=dp), intent(in) :: Dj1_termsOneAndTwo(nModes)
      !! First and second terms in the parentheses of D_j^1(t)
    complex(kind=dp), intent(in) :: Dj1_termThree(nModes)
      !! Third term in the parentheses of D_j^1(t) (without Sj)

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    ! Output variables:
    complex(kind=dp), intent(out) :: Dj0_t(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
    complex(kind=dp), intent(out) :: Dj1_t(nModes)
      !! Factor multiplying |M_j|^2*G_j^0(t) in G_j^1(t)

    ! Local variables:
    complex(kind=dp) :: Dj0OverSinOmegaPrime_t(nModes)
      !! Needed to keep from getting NaNs from cot()


    if(diffOmega) then

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

      Dj0_t(jTrue) = Sj(jTrue)*Dj0_tOverSj(jTrue)

      if(order == 1) then

        Dj1_t(jTrue) = (hbar_atomic/omega(jTrue))/2.0_dp*(Dj1_termsOneAndTwo(jTrue) + &
            Sj(jTrue)*Dj1_termThree(jTrue))

      endif
    endif

    return

  end subroutine setupStateDepTimeTablesNoDeltaNj

!----------------------------------------------------------------------------
  subroutine setupStateIndTimeTablesDeltaNj(countTrue, nModes, jTrue, omega, omegaPrime, time, diffOmega, cosOmegaPrime, &
            sinOmegaPrime, posExp_t)

    implicit none

    ! Input variables:
    integer, intent(in) :: countTrue
      !! Number of indices where Sj > SjThresh
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jTrue(countTrue)
      !! Mode indices where Sj > SjThresh

    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    ! Output variables:
    real(kind=dp), intent(out) :: cosOmegaPrime(nModes)
      !! cos(omega' t/2)
    real(kind=dp), intent(out) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)

    complex(kind=dp), intent(out) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed


    posExp_t(jTrue) = exp(ii*omega(jTrue)*time)

    if(diffOmega) then
      sinOmegaPrime(jTrue) = sin(omegaPrime(jTrue)*time/2.0_dp)
      cosOmegaPrime(jTrue) = cos(omegaPrime(jTrue)*time/2.0_dp)
    endif

    return

  end subroutine setupStateIndTimeTablesDeltaNj

!----------------------------------------------------------------------------
  subroutine setupStateDepTimeTablesDeltaNj(countTrue, nModes, jTrue, cosOmegaPrime, nj, omega, omegaPrime, sinOmegaPrime, Sj, &
            SjPrime, posExp_t, diffOmega, Dj0_t, Dj1_t)

    implicit none

    ! Input variables:
    integer, intent(in) :: countTrue
      !! Number of indices where Sj > SjThresh
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jTrue(countTrue)
      !! Mode indices where Sj > SjThresh

    real(kind=dp), intent(in) :: cosOmegaPrime(nModes)
      !! cos(omega' t/2)
    real(kind=dp), intent(in) :: nj(nModes)
      !! Base \(n_j\) occupation number for all states
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)
    real(kind=dp), intent(in) :: Sj(nModes), SjPrime(nModes)
      !! Huang-Rhys factor for each mode

    complex(kind=dp), intent(in) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    ! Output variables:
    complex(kind=dp), intent(out) :: Dj0_t(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
    complex(kind=dp), intent(out) :: Dj1_t(nModes)
      !! Factor multiplying |M_j|^2*G_j^0(t) in G_j^1(t)

    ! Local variables:
    complex(kind=dp) :: Aj_t(nModes)
      !! Aj from text. See equation below.
    complex(kind=dp) :: Dj0_tOverSj(nModes)
      !! Argument of exponential in G_j^0(t) = e^{i*D_j^0(t)}
      !! divided by Sj
    complex(kind=dp) :: Dj0OverSinOmegaPrime_t(nModes)
      !! Needed to keep from getting NaNs from cot()
    complex(kind=dp) :: Dj1_termsOneAndTwo(nModes)
      !! First and second terms in the parentheses of D_j^1(t)
    complex(kind=dp) :: Dj1_termThree(nModes)
      !! Third term in the parentheses of D_j^1(t) (without Sj)
    complex(kind=dp) :: expTimesNBarPlus1(nModes)
      !! Local storage of \(e^{i\omega_j t}(\bar{n}_j + 1)\) 
    complex(kind=dp) :: njOverPosExp_t(nModes)
      !! n_j*e^{-i\omega_j t}


    expTimesNBarPlus1(jTrue) = posExp_t(jTrue)*(nj(jTrue)+1.0_dp)

    if(diffOmega) then
      Aj_t(jTrue) = (expTimesNBarPlus1(jTrue) + nj(jTrue))/(expTimesNBarPlus1(jTrue) - nj(jTrue))

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

      Dj0_tOverSj(jTrue) = 1.0_dp/ii*(expTimesNBarPlus1(jTrue) + njOverPosExp_t(jTrue) - (2.0_dp*nj(jTrue) + 1.0_dp))

      if(order == 1) then

        Dj1_termsOneAndTwo(jTrue) = njOverPosExp_t(jTrue) + expTimesNBarPlus1(jTrue)
        Dj1_termThree(jTrue) = (1 + njOverPosExp_t(jTrue) - expTimesNBarPlus1(jTrue))**2

      endif

      Dj0_t(jTrue) = Sj(jTrue)*Dj0_tOverSj(jTrue)

      if(order == 1) then

        Dj1_t(jTrue) = (hbar_atomic/omega(jTrue))/2.0_dp*(Dj1_termsOneAndTwo(jTrue) + &
            Sj(jTrue)*Dj1_termThree(jTrue))

      endif
    endif

    return

  end subroutine setupStateDepTimeTablesDeltaNj

!----------------------------------------------------------------------------
  subroutine realTimeIntegration(mDim, nModes, nRealTimeSteps, nTransitions, order, ibi, ibf, iki, ikf, iSpin, &
          dE, deltaNjInitApproach, dt, dtau, energyAvgWindow, gamma0, matrixElement, njBase, omega, omegaPrime, Sj, SjPrime, &
          SjThresh, totalDeltaNj, transitionRate, addDeltaNj, captured, diffOmega, thermalize, writeEiRate, carrierDensityInput, &
          EiRateOutDir, energyTableDir, njNewOutDir, transRateOutDir, volumeLine)

    implicit none

    ! Input variables:
    integer, intent(in) :: mDim
      !! Size of first dimension for matrix element
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: nRealTimeSteps
      !! Number of real-time steps for updating occupations
      !! if applicable
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)
    integer, intent(in) :: ibi(nTransitions), ibf(nTransitions), iki(nTransitions), ikf(nTransitions)
      !! State indices
    integer, intent(in) :: iSpin
      !! Spin channel to use

    real(kind=dp), intent(in) :: dE(5,nTransitions,nkPerPool)
      !! All energy differences from energy table
    real(kind=dp), intent(in) :: deltaNjInitApproach(:,:)
      !! Optional delta nj from adjustment to 
      !! carrier approach in initial state
    real(kind=dp), intent(in) :: dt
      !! Optional real time step for generating 
      !! new phonon occupations
    real(kind=dp), intent(in) :: dtau
      !! Time step size for time-domain integration
    real(kind=dp), intent(in) :: energyAvgWindow
      !! Size of window to average over for carrier 
      !! density evaluation
    real(kind=dp), intent(in) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(in) :: matrixElement(mDim,nTransitions,nkPerPool)
      !! Electronic matrix element
    real(kind=dp), intent(inout) :: njBase(nModes)
      !! Base \(n_j\) occupation number for all states
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: Sj(:,:), SjPrime(:,:)
      !! Huang-Rhys factor for each mode (and 
      !! transition for scattering)
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold for Sj to determine which modes to calculate
    real(kind=dp), intent(in) :: totalDeltaNj(nModes,nTransitions)
      !! Optional total change in occupation numbers
      !! for each mode and transition
    real(kind=dp), intent(inout) :: transitionRate(nTransitions)
      !! \(Gamma_i\) transition rate

    logical, intent(in) :: addDeltaNj
      !! Add change in occupations for different scattering states
    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(in) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes
    logical, intent(in) :: writeEiRate
      !! If the energy transfer rate should be output as a 
      !! function of energy

    character(len=300), intent(in) :: carrierDensityInput
      !! Path to carrier density if generating 
      !! new phonon occupations
    character(len=300), intent(in) :: EiRateOutDir
      !! Output directory for energy transfer rate by 
      !! initial-state energy if writeEiRate == .true.
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(in) :: njNewOutDir
      !! Path to output new occupations if applicable
    character(len=300), intent(in) :: transRateOutDir
      !! Path to store transition rates
    character(len=300), intent(in) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: iRt
      !! Loop index
    integer :: nUniqueInitStates
      !! Number of unique initial states, defined by
      !! iki and ibi pairs
    integer, allocatable :: uniqueInitStates_ib(:)
      !! Band indices for unique initial states
    integer, allocatable :: uniqueInitStates_ik(:)
      !! k-point indices for unique initial states

    real(kind=dp), allocatable :: carrierDensity(:)
      !! Carrier density for each of the initial states, averaged
      !! over a given window
    real(kind=dp), allocatable :: dEEigInit(:)
      !! Eigenvalue difference of initial states
      !! relative to band edge
    real(kind=dp) :: energyTransferRate
      !! Total energy transfer combining all states
    real(kind=dp), allocatable :: energyTransferRate_i(:)
      !! Average energy transfer rate from each initial state
    real(kind=dp) :: k1(nModes), k2(nModes), k3(nModes), k4(nModes)
      !! Intermediate slopes for RK4 integration
    real(kind=dp) :: njNextEst(nModes)
      !! Next estimate for njBase
    real(kind=dp) :: njRateOfChange(nModes)
      !! Rate of change of occupations due to average
      !! effect of all transitions


    call getUniqueInitialStates(nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, uniqueInitStates_ik)

    allocate(dEEigInit(nUniqueInitStates), carrierDensity(nUniqueInitStates), energyTransferRate_i(nUniqueInitStates))

    call readDEPlot(iSpin, nUniqueInitStates, energyTableDir, dEEigInit)

    call readCarrierDensity(nUniqueInitStates, dEEigInit, energyAvgWindow, carrierDensityInput, volumeLine, carrierDensity)


    ! Get initial rate of change of occupations
    call getExcitationRate(nModes, nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, &
            uniqueInitStates_ik, carrierDensity, dE, dEEigInit, omega, totalDeltaNj, transitionRate, thermalize, &
            writeEiRate, njRateOfChange, energyTransferRate, energyTransferRate_i)

    ! Write occupations file with rate of change for first point
    call writeNewOccupations(1, nModes, njBase, njRateOfChange, dt, njNewOutDir)


    if(ionode) write(*, '("--------------------Beginning real-time integration ")')

    do iRt = 2, nRealTimeSteps
      !-----------------------------------------------------------
      ! Based on initial derivative estimate
      if(ionode) then
        ! First estimated slope is at the current point
        k1(:) = njRateOfChange(:)
        ! Update estimate of njBase after half time step
        njNextEst(:) = njBase(:) + k1(:)*dt/2.0d0
      endif

      ! Broadcast to all processes and get new transition rate and nj rate of change
      call MPI_BCAST(njNextEst, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

      if(ionode) write(*, '("Beginning transition-rate calculation for part 1 of RK4 time-integration step ",i5)') iRt
      call getAndWriteTransitionRate(nTransitions, ibi, ibf, iki, ikf, iRt, iSpin, mDim, nModes, order, dE(1:3,:,:), deltaNjInitApproach, &
            dtau, gamma0, matrixElement, njNextEst, omega, omegaPrime, Sj, SjPrime, SjThresh, addDeltaNj, &
            captured, diffOmega, generateNewOccupations, transRateOutDir, volumeLine, transitionRate)

      call getExcitationRate(nModes, nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, &
            uniqueInitStates_ik, carrierDensity, dE, dEEigInit, omega, totalDeltaNj, transitionRate, thermalize, &
            writeEiRate, njRateOfChange, energyTransferRate, energyTransferRate_i)

      !-----------------------------------------------------------
      ! Based on first midpoint estimate
      if(ionode) then
        ! Second estimated slope is from the midpoint
        k2(:) = njRateOfChange(:)
        ! Update estimate of njBase after half time step
        njNextEst(:) = njBase(:) + k2(:)*dt/2.0d0
      endif

      ! Broadcast to all processes and get new transition rate
      call MPI_BCAST(njNextEst, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

      if(ionode) write(*, '("Beginning transition-rate calculation for part 2 of RK4 time-integration step ",i5)') iRt
      call getAndWriteTransitionRate(nTransitions, ibi, ibf, iki, ikf, iRt, iSpin, mDim, nModes, order, dE(1:3,:,:), deltaNjInitApproach, &
            dtau, gamma0, matrixElement, njNextEst, omega, omegaPrime, Sj, SjPrime, SjThresh, addDeltaNj, &
            captured, diffOmega, generateNewOccupations, transRateOutDir, volumeLine, transitionRate)

      call getExcitationRate(nModes, nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, &
            uniqueInitStates_ik, carrierDensity, dE, dEEigInit, omega, totalDeltaNj, transitionRate, thermalize, &
            writeEiRate, njRateOfChange, energyTransferRate, energyTransferRate_i)

      !-----------------------------------------------------------
      ! Based on second midpoint estimate
      if(ionode) then
        ! Third estimated slope is from the midpoint
        k3(:) = njRateOfChange(:)
        ! Update estimate of njBase after full step
        njNextEst(:) = njBase(:) + k3(:)*dt
      endif

      ! Broadcast to all processes and get new transition rate
      call MPI_BCAST(njNextEst, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

      if(ionode) write(*, '("Beginning transition-rate calculation for part 3 of RK4 time-integration step ",i5)') iRt
      call getAndWriteTransitionRate(nTransitions, ibi, ibf, iki, ikf, iRt, iSpin, mDim, nModes, order, dE(1:3,:,:), deltaNjInitApproach, &
            dtau, gamma0, matrixElement, njNextEst, omega, omegaPrime, Sj, SjPrime, SjThresh, addDeltaNj, &
            captured, diffOmega, generateNewOccupations, transRateOutDir, volumeLine, transitionRate)

      call getExcitationRate(nModes, nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, &
            uniqueInitStates_ik, carrierDensity, dE, dEEigInit, omega, totalDeltaNj, transitionRate, thermalize, &
            writeEiRate, njRateOfChange, energyTransferRate, energyTransferRate_i)

      !-----------------------------------------------------------
      ! Based on endpoint estimate
      if(ionode) then
        ! Fourth  estimated slope is at the next time point
        k4(:) = njRateOfChange(:)

        ! Calculate the total estimated rate of change based on a weighted
        ! average of the four estimated slopes
        njRateOfChange(:) = (k1(:) + 2.0d0*k2(:) + 2.0d0*k3(:) + k4(:))/6.0d0 

        ! Update estimate of njBase with final estimated rate of change
        njBase(:) = njBase(:) + njRateOfChange(:)*dt
      endif

      ! Broadcast to all processes and write
      call MPI_BCAST(njBase, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(njRateOfChange, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

      if(ionode) write(*, '("Writing occupations for time-integration step ",i5)') iRt

      call writeNewOccupations(iRt, nModes, njBase, njRateOfChange, dt, njNewOutDir)

    enddo


    deallocate(dEEigInit)
    deallocate(carrierDensity)
    deallocate(energyTransferRate_i)

    return

  end subroutine realTimeIntegration

!----------------------------------------------------------------------------
  subroutine getUniqueInitialStates(nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, uniqueInitStates_ik)

    use miscUtilities, only: getUniqueInts

    implicit none

    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions), iki(nTransitions)
      !! Initial-state indices

    ! Output variables:
    integer, intent(out) :: nUniqueInitStates
      !! Number of unique initial states, defined by
      !! iki and ibi pairs
    integer, allocatable, intent(out) :: uniqueInitStates_ib(:)
      !! Band indices for unique initial states
    integer, allocatable, intent(out) :: uniqueInitStates_ik(:)
      !! k-point indices for unique initial states

    ! Local variables:
    integer, allocatable :: ibi_thisUniqueK(:)
      !! All values in ibi corresponding to the indices
      !! in iki that equal the current unique k index
    integer, allocatable :: ibiUnique(:)
      !! Unique initial bands for each unique k-point
    integer, allocatable :: ikiUnique(:)
      !! Unique initial k-points
    integer :: iU_iki, iU_ibi
      !! Loop indices
    integer :: nUnique_iki, nUnique_ibi
      !! Number of unique initial k-points and bands 
    integer :: uniqueInitStates_ik_large(nTransitions), uniqueInitStates_ib_large(nTransitions)
      !! Arrays of max possible size. Will only pass out
      !! arrays of size nUniqueInitStates


    ! First, get the unique initial k-points for each system
    if(ionode) then
      call getUniqueInts(nTransitions, iki, nUnique_iki, ikiUnique)

      nUniqueInitStates = 0
      do iU_iki = 1, nUnique_iki

        ibi_thisUniqueK = PACK(ibi, iki == ikiUnique(iU_iki))
          ! Includes allocate(ibi_thisUniqueK)

        call getUniqueInts(SIZE(ibi_thisUniqueK), ibi_thisUniqueK, nUnique_ibi, ibiUnique)
          ! Includes allocate(ibiUnique)

        do iU_ibi = 1, nUnique_ibi

          nUniqueInitStates = nUniqueInitStates + 1

          ! Store unique state-index pairs in an array with the maximum
          ! possible size of 
          uniqueInitStates_ik_large(nUniqueInitStates) = ikiUnique(iU_iki)
          uniqueInitStates_ib_large(nUniqueInitStates) = ibiUnique(iU_ibi)

        enddo


        ! Make sure to deallocate these variables allocated withing
        ! the loop before repeating the loop
        deallocate(ibi_thisUniqueK)
        deallocate(ibiUnique)
      enddo

      ! Copy larger arrays into arrays of size nUniqueInitStates to return
      allocate(uniqueInitStates_ik(nUniqueInitStates), uniqueInitStates_ib(nUniqueInitStates))
      uniqueInitStates_ik(:) = uniqueInitStates_ik_large(1:nUniqueInitStates)
      uniqueInitStates_ib(:) = uniqueInitStates_ib_large(1:nUniqueInitStates)
    endif

    call MPI_BCAST(nUniqueInitStates, 1, MPI_INTEGER, root, worldComm, ierr)
    if(.not. ionode) &
      allocate(uniqueInitStates_ik(nUniqueInitStates), uniqueInitStates_ib(nUniqueInitStates))
    call MPI_BCAST(uniqueInitStates_ik, nUniqueInitStates, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(uniqueInitStates_ib, nUniqueInitStates, MPI_INTEGER, root, worldComm, ierr)


    return

  end subroutine getUniqueInitialStates

!----------------------------------------------------------------------------
  subroutine readCarrierDensity(nUniqueInitStates, dEEigInit, energyAvgWindow, carrierDensityInput, volumeLine, carrierDensity)

    use cell, only: volume

    implicit none

    ! Input variables:
    integer, intent(in) :: nUniqueInitStates
      !! Number of unique initial states, defined by
      !! iki and ibi pairs

    real(kind=dp), intent(in) :: dEEigInit(nUniqueInitStates)
      !! Eigenvalue difference of initial states
      !! relative to band edge
    real(kind=dp), intent(in) :: energyAvgWindow
      !! Size of window to average over for carrier 
      !! density evaluation

    character(len=300), intent(in) :: carrierDensityInput
      !! Path to carrier density if generating 
      !! new phonon occupations
    character(len=300), intent(in) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Output variables:
    real(kind=dp), intent(out) :: carrierDensity(nUniqueInitStates)
      !! Carrier density for each of the initial states, averaged
      !! over a given window

    ! Local variables:
    integer :: iS, iUInit
      !! Loop index
    integer :: nSamplesEachInitState(nUniqueInitStates)
      !! Number of energy samples within the averaging window
      !! for each initial state
    integer :: nSampleEnergies
      !! Number of energy samples in carrier density file

    real(kind=dp) :: carrierDensity_iS
      !! Carrier density of this energy sample
    real(kind=dp) :: dEEigInit_iS
      !! Eigenvalue difference of this energy sample
      !! relative to band edge
    real(kind=dp) :: rDum
      !! Dummy real to ignore input

    character(len=300) :: textDum
      !! Dummy text to ignore input


    if(ionode) then
      open(unit=37, file=trim(carrierDensityInput))

      ! Ignore 3 header lines
      read(37,*)
      read(37,*)
      read(37,*)

      ! Read in the number of points in the carrier density file
      read(37,*) nSampleEnergies

      ! Ignore final header line
      read(37,*)


      ! Initialize the carrier density for each initial state to zero
      carrierDensity = 0.0_dp
      nSamplesEachInitState = 0

      ! Note that the arrays in our system are index in order of increasing
      ! energy (i.e., band index), but the carrier-density file is in order
      ! of increasing energy from the band edge, so for holes they will be
      ! in the opposite order.
      do iS = 1, nSampleEnergies
        read(37,*) dEEigInit_iS, rDum, carrierDensity_iS

        do iUInit = 1, nUniqueInitStates
          ! Have to use the negative sign because energies are positive in carrier
          ! density file
          if(abs(-dEEigInit_iS - dEEigInit(iUInit)) <= energyAvgWindow/2.0_dp) then
            nSamplesEachInitState(iUInit) = nSamplesEachInitState(iUInit) + 1
            carrierDensity(iUInit) = carrierDensity(iUInit) + carrierDensity_iS
          endif

          ! Use this to investigate your energyAvgWindow
          !write(*,'(2i5,4ES12.3E3, L5, i5, ES12.3E3)') iS, iUInit, -dEEigInit_iS, dEEigInit(iUInit), &
          !                           abs(-dEEigInit_iS - dEEigInit(iUInit)), energyAvgWindow, &
          !                           abs(-dEEigInit_iS - dEEigInit(iUInit)) <= energyAvgWindow/2.0_dp, &
          !                           nSamplesEachInitState(iUInit), carrierDensity(iUInit) 

        enddo

      enddo

      close(37)


      read(volumeLine,'(a51, ES24.15E3)') textDum, volume

      where(nSamplesEachInitState /= 0) carrierDensity = carrierDensity/nSamplesEachInitState*volume*(BohrToMeter*1.0d2)**3

    endif

    call MPI_BCAST(carrierDensity, nUniqueInitStates, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine readCarrierDensity

!----------------------------------------------------------------------------
  subroutine getExcitationRate(nModes, nTransitions, ibi, iki, nUniqueInitStates, uniqueInitStates_ib, &
          uniqueInitStates_ik, carrierDensity, dE, dEEigInit, omega, totalDeltaNj, transitionRate, thermalize, &
          writeEiRate, njRateOfChange, energyTransferRate, energyTransferRate_i)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions), iki(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: nUniqueInitStates
      !! Number of unique initial states, defined by
      !! iki and ibi pairs
    integer, intent(in) :: uniqueInitStates_ib(nUniqueInitStates)
      !! Band indices for unique initial states
    integer, intent(in) :: uniqueInitStates_ik(nUniqueInitStates)
      !! k-point indices for unique initial states

    real(kind=dp), intent(in) :: carrierDensity(nUniqueInitStates)
      !! Carrier density for each of the initial states, averaged
      !! over a given window
    real(kind=dp), intent(in) :: dE(5,nTransitions,nkPerPool)
      !! All energy differences from energy table
    real(kind=dp), intent(in) :: dEEigInit(nUniqueInitStates)
      !! Eigenvalue difference of initial states
      !! relative to band edge
    real(kind=dp), intent(in) :: omega(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: totalDeltaNj(nModes,nTransitions)
      !! Optional total change in occupation numbers
      !! for each mode and transition
    real(kind=dp), intent(in) :: transitionRate(nTransitions)
      !! \(Gamma_i\) transition rate

    logical, intent(in) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes
    logical, intent(in) :: writeEiRate
      !! If the energy transfer rate should be output as a 
      !! function of energy

    ! Output variables:
    real(kind=dp), intent(out) :: njRateOfChange(nModes)
      !! Rate of change of occupations due to average
      !! effect of all transitions
    real(kind=dp), intent(out) :: energyTransferRate
      !! Total energy transfer combining all states
    real(kind=dp), intent(out) :: energyTransferRate_i(nUniqueInitStates)
      !! Average energy transfer rate from each initial state

    ! Local variables:
    real(kind=dp), allocatable :: njRateOfChange_i(:,:)
      !! Rate of change of occupations due to transitions
      !! from each initial state, i


    call sumOverFinalStates(nModes, nTransitions, nUniqueInitStates, ibi, iki, uniqueInitStates_ib, uniqueInitStates_ik, &
          dE, energyTransferRate_i, totalDeltaNj, transitionRate, thermalize, writeEiRate, njRateOfChange_i)

    call integrateOverInitialStates(nModes, nUniqueInitStates, carrierDensity, dEEigInit, energyTransferRate_i, &
          njRateOfChange_i, omega, thermalize, energyTransferRate, njRateOfChange)

    deallocate(njRateOfChange_i)

    return

  end subroutine getExcitationRate

!----------------------------------------------------------------------------
  subroutine sumOverFinalStates(nModes, nTransitions, nUniqueInitStates, ibi, iki, uniqueInitStates_ib, uniqueInitStates_ik, &
        dE, energyTransferRate_i, totalDeltaNj, transitionRate, thermalize, writeEiRate, njRateOfChange_i)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: nUniqueInitStates
      !! Number of unique initial states, defined by
      !! iki and ibi pairs
    integer, intent(in) :: ibi(nTransitions), iki(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: uniqueInitStates_ib(nUniqueInitStates)
      !! Band indices for unique initial states
    integer, intent(in) :: uniqueInitStates_ik(nUniqueInitStates)
      !! k-point indices for unique initial states

    real(kind=dp), intent(in) :: dE(5,nTransitions,nkPerPool)
      !! All energy differences from energy table
    real(kind=dp), intent(out) :: energyTransferRate_i(nUniqueInitStates)
      !! Average energy transfer rate from each initial state
    real(kind=dp), intent(in) :: totalDeltaNj(nModes,nTransitions)
      !! Optional total change in occupation numbers
      !! for each mode and transition
    real(kind=dp), intent(in) :: transitionRate(nTransitions)
      !! \(Gamma_i\) transition rate

    logical, intent(in) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes
    logical, intent(in) :: writeEiRate
      !! If the energy transfer rate should be output as a 
      !! function of energy

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: njRateOfChange_i(:,:)
      !! Rate of change of occupations due to transitions
      !! from each initial state, i

    ! Local variables:
    integer :: iUInit, iE
      !! Loop indices


    allocate(njRateOfChange_i(nModes,nUniqueInitStates))
    njRateOfChange_i = 0.0_dp

    if(ionode) then
      do iUInit = 1, nUniqueInitStates
        do iE = 1, nTransitions
        
          ! Pick out the elements where the initial state (iki,ibi) corresponds to
          ! a given unique initial state
          if(iki(iE) == uniqueInitStates_ik(iUInit) .and. ibi(iE) == uniqueInitStates_ib(iUInit)) then

            ! Then sum over all of the final states for each unique initial state.
            ! For thermalized energy distribution, calculate the energy rate of change.
            ! Otherwise, channel the energy into the modes.
            if(thermalize .or. writeEiRate) then
              energyTransferRate_i(iUInit) = energyTransferRate_i(iUInit) + (dE(1,iE,1) + dE(4,iE,1) + dE(5,iE,1))*transitionRate(iE)
            else
              ! The (:) here indcates the phonon modes. They are all independent.
              njRateOfChange_i(:,iUInit) = njRateOfChange_i(:,iUInit) + totalDeltaNj(:,iE)*transitionRate(iE)
            endif

          endif
        enddo
      enddo
    endif

    return

  end subroutine sumOverFinalStates

!----------------------------------------------------------------------------
  subroutine integrateOverInitialStates(nModes, nUniqueInitStates, carrierDensity, dEEigInit, energyTransferRate_i, &
      njRateOfChange_i, omega, thermalize, energyTransferRate, njRateOfChange)

    use generalComputations, only: trapezoidIntegrationVariableDx 

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: nUniqueInitStates
      !! Number of unique initial states, defined by
      !! iki and ibi pairs

    real(kind=dp), intent(in) :: carrierDensity(nUniqueInitStates)
      !! Carrier density for each of the initial states, averaged
      !! over a given window
    real(kind=dp), intent(in) :: dEEigInit(nUniqueInitStates)
      !! Eigenvalue difference of initial states
      !! relative to band edge
    real(kind=dp), intent(in) :: energyTransferRate_i(nUniqueInitStates)
      !! Average energy transfer rate from each initial state
    real(kind=dp), intent(in) :: njRateOfChange_i(nModes,nUniqueInitStates)
      !! Rate of change of occupations due to transitions
      !! from each initial state, i
    real(kind=dp), intent(in) :: omega(nModes)
      !! Frequency for each mode

    logical, intent(in) :: thermalize
      !! If should convert energy transfer to local temp. or
      !! channel directly into modes

    ! Output variables:
    real(kind=dp), intent(out) :: energyTransferRate
      !! Total energy transfer combining all states
    real(kind=dp), intent(out) :: njRateOfChange(nModes)
      !! Rate of change of occupations due to average
      !! effect of all transitions

    ! Local variables:
    integer :: j
      !! Loop index

    real(kind=dp) :: integral
      !! Result of integration for each mode
    real(kind=dp) :: integrand(nUniqueInitStates)
      !! Integrand to pass to integration subroutine


    if(ionode) then
      if(thermalize) then
        integrand(:) = carrierDensity(:)*energyTransferRate_i(:)
        call trapezoidIntegrationVariableDx(nUniqueInitStates, integrand, dEEigInit, integral)
        energyTransferRate = integral
      else
        do j = 1, nModes
          integrand(:) = carrierDensity(:)*njRateOfChange_i(j,:)
          call trapezoidIntegrationVariableDx(nUniqueInitStates, integrand, dEEigInit, integral)
          njRateOfChange(j) = integral
        enddo

        energyTransferRate = sum(njRateOfChange(:)*omega(:))
      endif
    endif

    call MPI_BCAST(njRateOfChange, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(energyTransferRate, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine integrateOverInitialStates
  
!----------------------------------------------------------------------------
  subroutine writeNewOccupations(iRt, nModes, njBase, njRateOfChange, dt, njNewOutDir)

    implicit none

    ! Input variables:
    integer, intent(in) :: iRt
      !! Real-time integration step index
    integer, intent(in) :: nModes
      !! Number of phonon modes

    real(kind=dp), intent(in) :: njBase(nModes)
      !! Base \(n_j\) occupation number for all states
    real(kind=dp), intent(in) :: njRateOfChange(nModes)
      !! Rate of change of occupations due to average
      !! effect of all transitions
    real(kind=dp), intent(in) :: dt
      !! Optional real time step for generating 
      !! new phonon occupations

    character(len=300), intent(in) :: njNewOutDir
      !! Path to output new occupations if applicable

    ! Local variables:
    integer :: j
      !! Loop index

    if(ionode) then
      open(unit=37, file=trim(njNewOutDir)//'/nj.'//trim(int2str(iRt))//'.out')
      write(37,'("# dt = ",ES24.15E3)') dt
      write(37,'("# j, New nj, Rate of change (1/s)")')

      do j = 1, nModes
        write(37,'(i7,3ES24.15E3)') j, njBase(j), njRateOfChange(j)
      enddo

      close(37)
    endif

    return

  end subroutine writeNewOccupations
  
!----------------------------------------------------------------------------
  function transitionRateFileExists(ikGlobal, isp, transRateOutDir) result(fileExists)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    character(len=300), intent(in) :: transRateOutDir
      !! Path to store transition rates

    ! Output variables:
    logical :: fileExists
      !! If the overlap file exists for the given 
      !! k-point and spin channel


    inquire(file=trim(transRateOutDir)//'transitionRate.'//trim(int2str(isp))//"."//trim(int2str(ikGlobal)), exist=fileExists)
    
  end function transitionRateFileExists

end module LSFmod
