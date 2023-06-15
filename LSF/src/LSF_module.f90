module LSF0mod
  
  use constants, only: dp, HartreeToJ, HartreeToEv, eVToJ, ii, hbar
  use errorsAndMPI

  implicit none 
  real(kind=dp),parameter :: Kb =  1.38064852d-23
  real(kind=dp),parameter :: pi= 3.14159265358979
  real(kind=dp),parameter :: tpi = 6.2831853071795864769 
  real(kind=dp),parameter :: Thz = 1.0d12


  integer :: iTime_start, iTime_end
    !! Start and end time steps for this process
  integer :: nStepsLocal
    !! Number of time steps for each
    !! process to complete


  integer :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !! Energy band bounds for initial and final state
  integer :: nModes
    !! Number of phonon modes
  integer :: order
    !! Order to calculate (0 or 1)

  real(kind=dp) :: beta
    !! 1/kb*T
  real(kind=dp), allocatable :: dEDelta(:,:)
    !! Energy for delta function
  real(kind=dp), allocatable :: dEPlot(:)
    !! Energy for plotting
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
  real(kind=dp),allocatable :: modeFreq(:)
    !! Frequency for each mode
  real(kind=dp), allocatable :: Sj(:)
    !! Huang-Rhys factor for each mode
  real(kind=dp) :: smearingExpTolerance
    !! Tolerance for the Lorentzian-smearing
    !! exponential used to calculate max time
  real(kind=dp) :: temperature

  character(len=300) :: EInput
    !! Path to energy table to read
  character(len=300) :: MifInput
    !! Path to matrix element file `allElecOverlap.isp.ik`. 
    !! For first-order term, the path is just within each 
    !! subdirectory.
  character(len=300) :: MjDir
    !! Path to the base directory for the first-order
    !! matrix element calculations
  character(len=300) :: outputDir
    !! Path to output transition rates
  character(len=300) :: prefix
    !! Prefix of directories for first-order matrix
    !! elements
  character(len=300) :: SjInput
    !! Path to Sj.out file


  namelist /inputParams/ iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, EInput, MifInput, MjDir, SjInput, &
                        temperature, hbarGamma, dt, smearingExpTolerance, outputDir, order, prefix

contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, beta, dt, gamma0, hbarGamma, maxTime, &
        smearingExpTolerance, temperature, EInput, MifInput, MjDir, outputDir, prefix, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(out) :: beta
      !! 1/kb*T
    real(kind=dp), intent(out) :: dt
      !! Time step size
    real(kind=dp), intent(out) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(out) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: maxTime
      !! Max time for integration
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    real(kind=dp), intent(out) :: temperature

    character(len=300), intent(out) :: EInput
      !! Path to energy table to read
    character(len=300), intent(out) :: MifInput
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(out) :: SjInput
      !! Path to Sj.out file

  
    call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, dt, hbarGamma, smearingExpTolerance, temperature, EInput, &
          MifInput, MjDir, outputDir, prefix, SjInput)

    if(ionode) then

      read(5, inputParams, iostat=ierr)
        !! * Read input variables

    
      if(ierr /= 0) call exitError('LSF0 main', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, dt, hbarGamma, smearingExpTolerance, temperature, EInput, &
            MifInput, MjDir, outputDir, prefix, SjInput)

      dt = dt/Thz

      gamma0 = hbarGamma*1e-3*eVToJ/hbar
        ! Input expected in meV

      beta = 1.0d0/Kb/temperature

      maxTime = -log(smearingExpTolerance)/gamma0
      write(*,'("Max time: ", ES24.15E3)') maxTime

    endif

    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)
  
    call MPI_BCAST(beta, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(hbarGamma, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(gamma0, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(smearingExpTolerance, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(temperature, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(maxTime, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  
    call MPI_BCAST(EInput, len(EInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(MifInput, len(MifInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(MjDir, len(MjDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(SjInput, len(SjInput), MPI_CHARACTER, root, worldComm, ierr)
    
    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, dt, hbarGamma, smearingExpTolerance, temperature, EInput, &
        MifInput, MjDir, outputDir, prefix, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(out) :: dt
      !! Time step size
    real(kind=dp), intent(out) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    real(kind=dp), intent(out) :: temperature

    character(len=300), intent(out) :: EInput
      !! Path to energy table to read
    character(len=300), intent(out) :: MifInput
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(out) :: SjInput
      !! Path to Sj.out file

    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time


    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    order = -1

    dt = 1d-6
    hbarGamma = 0.0_dp
    smearingExpTolerance = 0.0_dp
    temperature = 0.0_dp

    EInput = ''
    MifInput = ''
    MjDir = ''
    SjInput = ''
    outputDir = './'
    prefix = 'disp-'

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(*, '(/5X,"LSF0 starts on ",A9," at ",A9)') &
             cdate, ctime

      write(*, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs


    endif

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, dt, hbarGamma, smearingExpTolerance, temperature, &
        EInput, MifInput, MjDir, outputDir, prefix, SjInput)

    implicit none

    ! Input variables
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dt
      !! Time step size
    real(kind=dp), intent(in) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(in) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    real(kind=dp), intent(in) :: temperature

    character(len=300), intent(in) :: EInput
      !! Path to energy table to read
    character(len=300), intent(in) :: MifInput
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(in) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(in) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(in) :: SjInput
      !! Path to Sj.out file

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9))
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 
    abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution 

    abortExecution = checkDoubleInitialization('dt', dt, 1.0d-10, 1.0d-4) .or. abortExecution
    abortExecution = checkDoubleInitialization('hbarGamma', hbarGamma, 0.1_dp, 20.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('smearingExpTolerance', smearingExpTolerance, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('temperature', temperature, 0.0_dp, 1500.0_dp) .or. abortExecution
      ! These limits are my best guess as to what is reasonable; they are not
      ! hard and fast, but you should think about the application of the theory
      ! to numbers outside these ranges.

    abortExecution = checkFileInitialization('EInput', EInput) .or. abortExecution
    abortExecution = checkFileInitialization('SjInput', SjInput) .or. abortExecution

    if(order == 0) then 
      abortExecution = checkFileInitialization('MifInput', MifInput) .or. abortExecution
    else if(order == 1) then
      abortExecution = checkDirInitialization('MjDir', MjDir, '/'//trim(prefix)//'0001/'//trim(MifInput)) .or. abortExecution
      write(*,'("prefix = ''",a,"''")') trim(prefix)
      write(*,'("MifInput = ''",a,"''")') trim(MifInput)
    endif

    call system('mkdir -p '//trim(outputDir))


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif

    return 

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readSj(SjInput, nModes, modeFreq, Sj)

    implicit none

    ! Input variables
    character(len=300), intent(in) :: SjInput
      !! Path to Sj.out file

    ! Output variables:
    integer, intent(out) :: nModes
      !! Number of phonon modes

    real(kind=dp),allocatable, intent(out) :: modeFreq(:)
      !! Frequency for each mode
    real(kind=dp), allocatable, intent(out) :: Sj(:)
      !! Huang-Rhys factor for each mode

    ! Local variables:
    integer :: iDum
      !! Dummy integer
    integer :: j
      !! Loop index
  
  
    if(ionode) then

      open(12,file=trim(SjInput))

      read(12,*) nModes

    endif

    call MPI_BCAST(nModes, 1, MPI_INTEGER, root, worldComm, ierr)


    allocate(Sj(1:nModes))
    allocate(modeFreq(1:nModes))

    
    if(ionode) then

      do j = 1, nModes
        read(12,*) iDum, Sj(j), modeFreq(j) ! freq read from Sj.out is f(in Thz)*2pi
      end do

      modeFreq(:) = modeFreq(:)*Thz
        ! Convert to Hz*2pi

      close(12)

    endif

    call MPI_BCAST(modeFreq, size(modeFreq), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(Sj, size(Sj), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return 

  end subroutine readSj

!----------------------------------------------------------------------------
  subroutine readEnergy(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, EInput, dEDelta, dEPlot)

    use miscUtilities, only: ignoreNextNLinesFromFile
  
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

    character(len=300), intent(in) :: EInput
      !! Path to energy table to read

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: dEDelta(:,:)
      !! Energy for delta function
    real(kind=dp), allocatable, intent(out) :: dEPlot(:)
      !! Energy for plotting

    ! Local variables:
    integer :: iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      !! Band bounds from energy table
    integer :: iDum
      !! Dummy integer
    integer :: ibi, ibf
      !! Loop indices

    real(kind=dp) :: rDum
      !! Dummy real

    logical :: abortExecution
      !! If the program should end


    if(ionode) then
      open(12,file=trim(EInput))

      read(12,*)
      read(12,*) iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
        ! @todo Test these values against the input values

      call ignoreNextNLinesFromFile(12,6)

    endif
      
    allocate(dEDelta(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal))
    allocate(dEPlot(iBandIinit:iBandIfinal))

    if(ionode) then

      do ibf = iBandFinit, iBandFfinal
        do ibi = iBandIinit, iBandIfinal

          read(12,*) iDum, iDum, dEDelta(ibf,ibi), rDum, rDum, dEPlot(ibi) ! in Hartree

        enddo
      enddo

      dEDelta(:,:) = dEDelta(:,:)*HartreeToJ
      dEPlot(:) = dEPlot(:)*HartreeToEv

      close(12)

    endif

    call MPI_BCAST(dEDelta, size(dEDelta), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(dEPlot, size(dEPlot), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine readEnergy

!----------------------------------------------------------------------------
  subroutine readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, MifInput, matrixElement)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    character(len=300), intent(in) :: MifInput
      !! Path to zeroth-order matrix element file
      !! `allElecOverlap.isp.ik`

    ! Output variables:
    real(kind=dp), intent(out) :: matrixElement(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! Electronic matrix element

    ! Local variables:
    integer :: iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      !! Band bounds from energy table
    integer :: iDum
      !! Dummy integer
    integer :: ibi, ibf
      !! Loop indices

    real(kind=dp) :: rDum
      !! Dummy real


    if(ionode) then
      open(12,file=trim(MifInput))

      read(12,*)
      read(12,*)
      read(12,*) iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
        ! @todo Test these values against the input values
      read(12,*)

      if(order == 1) read(12,*)
        ! Ignore additional line for phonon mode 

    endif
      

    if(ionode) then

      do ibf = iBandFinit, iBandFfinal
        do ibi = iBandIinit, iBandIfinal

          read(12,*) iDum, iDum, rDum, rDum, rDum, matrixElement(ibf,ibi) ! in Hartree^2

        enddo
      enddo

      matrixElement(:,:) = matrixElement(:,:)*HartreeToJ**2

      close(12)

    endif

    return

  end subroutine readMatrixElement

!----------------------------------------------------------------------------
  function G0ExpArg(time) 

    ! Input variables:
    !integer, intent(in) :: nModes
      ! Number of phonon modes

    !real(kind=dp), intent(in) :: beta
      ! 1/kb*T
    !real(kind=dp), intent(in) :: modeFreq(nModes)
      ! Frequency for each mode
    !real(kind=dp), intent(in) :: Sj(nModes)
      ! Huang-Rhys factor for each mode
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    ! Output variables:
    complex(kind=dp) :: G0ExpArg

    ! Local variables
    integer :: j
      !! Loop index

    real(kind=dp) :: nj
      !! \(n_j\) occupation number
    real(kind=dp) :: omega
      !! Local storage of frequency for this mode


    G0ExpArg = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do j = 1, nModes

      omega = modeFreq(j)

      nj = 1.0_dp/(exp(hbar*omega*beta) - 1.0_dp)

      G0ExpArg = G0ExpArg + Sj(j)*((nj+1.0_dp)*exp(ii*omega*time) + nj*exp(-ii*omega*time) - (2.0_dp*nj + 1.0_dp))
    enddo

  end function G0ExpArg

!----------------------------------------------------------------------------
  function Aj_t(j,time)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: j
      !! Mode index
    !integer, intent(in) :: nModes
      ! Number of phonon modes

    !real(kind=dp), intent(in) :: beta
      ! 1/kb*T
    !real(kind=dp), intent(in) :: modeFreq(nModes)
      ! Frequency for each mode
    !real(kind=dp), intent(in) :: Sj(nModes)
      ! Huang-Rhys factor for each mode
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    ! Output variables:
    complex(kind=dp) :: Aj_t

    ! Local variables
    real(kind=dp) :: nj
      !! \(n_j\) occupation number
    real(kind=dp) :: omega
      !! Local storage of frequency for this mode


    omega = modeFreq(j)

    nj = 1.0_dp/(exp(hbar*omega*beta) - 1.0_dp)

    Aj_t = (2*hbar/omega)*(nj*exp(-ii*omega*time) + (nj+1)*exp(ii*omega*time) + &
            Sj(j)*(1 + nj*exp(-ii*omega*time) - (nj+1)*exp(ii*omega*time))**2)

  end function Aj_t

!----------------------------------------------------------------------------
  subroutine getAndWriteTransitionRate(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, mDim, order, dEDelta, dEPlot, gamma0, &
        matrixElement, temperature)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: mDim
      !! Size of first dimension for matrix element
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dEDelta(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! Energy for delta function
    real(kind=dp), intent(in) :: dEPlot(iBandIinit:iBandIfinal)
      !! Energy for plotting
    real(kind=dp), intent(in) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(in) :: matrixElement(mDim,iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! Electronic matrix element
    real(kind=dp), intent(in) :: temperature

    ! Local variables:
    integer :: iTime, ibi, ibf
      !! Loop indices
    integer :: updateFrequency
      !! Frequency of steps to write status update

    real(kind=dp) :: Eif
      !! Local storage of dEDelta(ibi,ibf)
    real(kind=dp) :: expPrefactor
      !! Prefactor for the exponential inside the integral.
      !! For zeroth-order, this is just the matrix element,
      !! but for first-order this is \(\sum_j M_j A_j\).
    real(kind=dp) :: t1, t2
      !! Time for each time step
    real(kind=dp) :: timer1, timer2
      !! Timers
    real(kind=dp) :: transitionRate(iBandIinit:iBandIfinal)
      !! \(Gamma_i\) transition rate

    complex(kind=dp) :: expArg_t1_base, expArg_t2_base
      !! Base exponential argument for each time step
    complex(kind=dp) :: expArg_t1, expArg_t2
      !! Exponential argument for each time step


    updateFrequency = ceiling(nStepsLocal/10.0)
    call cpu_time(timer1)

    transitionRate(:) = 0.0_dp
    do iTime = 1, nStepsLocal-2, 2
      ! Loop should not include the first step (0), 
      ! but it should include the last step (nStepsLocal-1).
      ! iTime stops at nStepsLocal-2 because nStepsLocal-1
      ! is calculated by t2 at the last step.

      if(ionode .and. (mod(iTime,updateFrequency) == 0 .or. mod(iTime+1,updateFrequency) == 0)) then

        call cpu_time(timer2)
        write(*,'(i2,"% complete with transition-rate loop. Time in loop: ",f10.2," secs")') iTime*100/nStepsLocal, timer2-timer1

      endif

      t1 = (float(iTime) + myid*float(nStepsLocal))*dt
        ! Must do this arithmetic with floats to avoid
        ! integer overflow
      expArg_t1_base = G0ExpArg(t1) - gamma0*t1

      t2 = t1 + dt
      expArg_t2_base = G0ExpArg(t2) - gamma0*t2

      do ibi = iBandIinit, iBandIfinal
        do ibf = iBandFinit, iBandFfinal

          expPrefactor = matrixElement(1,ibf,ibi)
          Eif = dEDelta(ibf,ibi)

          expArg_t1 = expArg_t1_base + ii*Eif/hbar*t1
          expArg_t2 = expArg_t2_base + ii*Eif/hbar*t2

          transitionRate(ibi) = transitionRate(ibi) + &
            Real(4.0_dp*expPrefactor*exp(expArg_t1) + 2.0_dp*expPrefactor*exp(expArg_t2))
            ! We are doing multiple sums, but they are all commutative.
            ! Here we add in the contribution to the integral at this time
            ! step from a given final state. The loop over final states 
            ! adds in the contributions from all final states. 

        enddo
      enddo
    enddo

    do ibi = iBandIinit, iBandIfinal
      do ibf = iBandFinit, iBandFfinal

        expPrefactor = matrixElement(1,ibf,ibi)
        Eif = dEDelta(ibf,ibi)

        t1 = myid*float(nStepsLocal)*dt
          ! Must do this arithmetic with floats to avoid
          ! integer overflow
        transitionRate(ibi) = transitionRate(ibi) + Real(expPrefactor*exp(G0ExpArg(t1) + ii*Eif/hbar*t1 - gamma0*t1))
          ! Add \(t_0\) that was skipped in the loop

        t2 = ((myid+1)*float(nStepsLocal) - 1.0_dp)*dt
        transitionRate(ibi) = transitionRate(ibi) - Real(expPrefactor*exp(G0ExpArg(t2) + ii*Eif/hbar*t2 - gamma0*t2)) 
          ! Subtract off last time that had a coefficient
          ! of 2 in the loop but should really have a 
          ! coefficient of 1

      enddo
    enddo

    call mpiSumDoubleV(transitionRate, worldComm)
      ! Combine results from all processes

    if(ionode) then

      open(unit=37, file=trim(outputDir)//'transitionRate.txt')

      write(37,'("# Total number of initial states, Initial states (bandI, bandF) Format : ''(3i10)''")')
      write(37,'(3i10)') iBandIfinal-iBandIinit+1, iBandIinit, iBandIfinal

      write(37,'("# Temperature (K): ", f7.1)') temperature

      write(37,'("# Initial state, dEPlot (eV), Transition rate Format : ''(i10, f10.5, ES24.15E3)''")')

      ! Multiply by prefactor for Simpson's integration method 
      ! and prefactor for time-domain integral
      if(order == 0) then
        transitionRate(:) = transitionRate(:)*(dt/3.0_dp)*(2.0_dp/(hbar*hbar))
      else if(order == 1) then
        transitionRate(:) = transitionRate(:)*(dt/3.0_dp)*(1.0_dp/(2.0_dp*hbar*hbar))
      endif        

      do ibi = iBandIinit, iBandIfinal
        write(37,'(i10, f10.5,ES24.14E3)') ibi, dEPlot(ibi), transitionRate(ibi)
      enddo

    endif

    return

  end subroutine getAndWriteTransitionRate

end module LSF0mod
