module LSFmod
  
  use constants, only: dp, HartreeToJ, HartreeToEv, eVToJ, ii, hbar, THzToHz, kB, BohrToMeter, elecMToKg
  use base, only: nKPoints, order
  use TMEmod, only: getMatrixElementFNameWPath, getMatrixElementFName, readMatrixElement
  use PhononPPMod, only: diffOmega, readSjOneFreq, readSjTwoFreq, omega, omegaPrime, Sj, SjPrime
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
  integer :: nModes
    !! Number of phonon modes
  integer :: nTransitions
    !! Total number of transitions 

  real(kind=dp) :: beta
    !! 1/kb*T
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
  real(kind=dp) :: smearingExpTolerance
    !! Tolerance for the Lorentzian-smearing
    !! exponential used to calculate max time
  real(kind=dp) :: temperature

  character(len=300) :: matrixElementDir
    !! Path to matrix element file `allElecOverlap.isp.ik`. 
    !! For first-order term, the path is just within each 
    !! subdirectory.
  character(len=300) :: MjBaseDir
    !! Path to the base directory for the first-order
    !! matrix element calculations
  character(len=300) :: outputDir
    !! Path to output transition rates
  character(len=300) :: prefix
    !! Prefix of directories for first-order matrix
    !! elements
  character(len=300) :: SjInput
    !! Path to Sj.out file
  character(len=300) :: volumeLine
    !! Volume line from overlap file to be
    !! output exactly in transition rate file


  namelist /inputParams/ energyTableDir, matrixElementDir, MjBaseDir, SjInput, temperature, hbarGamma, dt, &
                         smearingExpTolerance, outputDir, order, prefix, iSpin, diffOmega

contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iSpin, order, beta, dt, gamma0, hbarGamma, maxTime, smearingExpTolerance, temperature, &
        diffOmega, energyTableDir, matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iSpin
      !! Spin channel to use
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
    
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(out) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(out) :: SjInput
      !! Path to Sj.out file

  
    call initialize(iSpin, order, dt, hbarGamma, smearingExpTolerance, diffOmega, temperature, energyTableDir, &
          matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

    if(ionode) then

      read(5, inputParams, iostat=ierr)
        !! * Read input variables

    
      if(ierr /= 0) call exitError('LSF module', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(iSpin, order, dt, hbarGamma, smearingExpTolerance, diffOmega, temperature, energyTableDir, &
            matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

      dt = dt/THzToHz

      gamma0 = hbarGamma*1e-3*eVToJ/hbar
        ! Input expected in meV

      beta = 1.0d0/(kB*temperature)

      maxTime = -log(smearingExpTolerance)/gamma0
      write(*,'("Max time: ", ES24.15E3)') maxTime

    endif

    call MPI_BCAST(iSpin, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)
  
    call MPI_BCAST(beta, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(hbarGamma, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(gamma0, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(smearingExpTolerance, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(temperature, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(maxTime, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(diffOmega, 1, MPI_LOGICAL, root, worldComm, ierr)
  
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(matrixElementDir, len(matrixElementDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(MjBaseDir, len(MjBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(SjInput, len(SjInput), MPI_CHARACTER, root, worldComm, ierr)
    
    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(iSpin, order, dt, hbarGamma, smearingExpTolerance, diffOmega, temperature, energyTableDir, &
        matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

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
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    real(kind=dp), intent(out) :: temperature
    
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(out) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(out) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements
    character(len=300), intent(out) :: SjInput
      !! Path to Sj.out file


    iSpin = 1
    order = -1

    dt = 1d-6
    hbarGamma = 0.0_dp
    smearingExpTolerance = 0.0_dp
    temperature = 0.0_dp

    diffOmega = .false.

    energyTableDir = ''
    matrixElementDir = ''
    MjBaseDir = ''
    SjInput = ''
    outputDir = './'
    prefix = 'disp-'

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iSpin, order, dt, hbarGamma, smearingExpTolerance, diffOmega, temperature, energyTableDir, &
        matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

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
    real(kind=dp), intent(in) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    real(kind=dp), intent(in) :: temperature
    
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to read
    character(len=300), intent(in) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjBaseDir
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
    character(len=300) :: fName
      !! File name for matrix element file to get
      !! nKPoints from 

    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iSpin', iSpin, 1, 2)
    abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution 

    abortExecution = checkDoubleInitialization('dt', dt, 1.0d-10, 1.0d-4) .or. abortExecution
    abortExecution = checkDoubleInitialization('hbarGamma', hbarGamma, 0.1_dp, 20.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('smearingExpTolerance', smearingExpTolerance, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('temperature', temperature, 0.0_dp, 1500.0_dp) .or. abortExecution
      ! These limits are my best guess as to what is reasonable; they are not
      ! hard and fast, but you should think about the application of the theory
      ! to numbers outside these ranges.

    abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.'//trim(int2str(iSpin))//'.1') .or. abortExecution
    abortExecution = checkFileInitialization('SjInput', SjInput) .or. abortExecution

    write(*,'("diffOmega = ",L)') diffOmega

    if(order == 0) then 
      abortExecution = checkDirInitialization('matrixElementDir', matrixElementDir, getMatrixElementFName(1,iSpin)) .or. abortExecution

      fName = getMatrixElementFNameWPath(1,iSpin,matrixElementDir)

    else if(order == 1) then

      abortExecution = checkDirInitialization('MjBaseDir', MjBaseDir, &
            '/'//trim(prefix)//'0001/'//trim(getMatrixElementFNameWPath(1,iSpin,matrixElementDir))) .or. abortExecution
      write(*,'("prefix = ''",a,"''")') trim(prefix)
      write(*,'("matrixElementDir = ''",a,"''")') trim(matrixElementDir)

      fName = trim(MjBaseDir)//'/'//trim(prefix)//'0001/'//trim(getMatrixElementFNameWPath(1,iSpin,matrixElementDir))

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
  subroutine getAndWriteTransitionRate(nTransitions, ibi, iSpin, mDim, order, nModes, dE, gamma0, & 
          matrixElement, temperature, volumeLine)
    
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
    real(kind=dp), intent(in) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(in) :: matrixElement(mDim,nTransitions,nkPerPool)
      !! Electronic matrix element
    real(kind=dp), intent(in) :: temperature

    character(len=300), intent(in) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: iTime, iE, ikLocal, ikGlobal
      !! Loop indices
    integer :: updateFrequency
      !! Frequency of steps to write status update

    real(kind=dp) :: Eif
      !! Local storage of energy for delta function
    real(kind=dp) :: expPrefactor
      !! Prefactor for the exponential inside the integral.
      !! For zeroth-order, this is just the matrix element,
      !! but for first-order this is \(\sum_j M_j A_j\).
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

    complex(kind=dp) :: Aj(nModes)
      !! \(A_j(t)\) for all modes
    complex(kind=dp) :: expArg_base
      !! Base exponential argument for each time step
    complex(kind=dp) :: expArg
      !! Exponential argument for each time step

      
    updateFrequency = ceiling(nStepsLocal/10.0)
    call cpu_time(timer1)


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
      expArg_base = G0ExpArg(time) - gamma0*time
      if(order == 1) Aj(:) = getAj_t(time)

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
            expPrefactor = getFirstOrderPrefactor(nModes, matrixElement(:,iE,ikLocal), Aj)
          endif

          Eif = dE(1,iE,ikLocal)

          expArg = expArg_base + ii*Eif/hbar*time

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
      if(order == 0) then
        transitionRate(:,:) = transitionRate(:,:)*(dt/3.0_dp)*(2.0_dp/(hbar*hbar))
      else if(order == 1) then
        transitionRate(:,:) = transitionRate(:,:)*(dt/3.0_dp)*(1.0_dp/(2.0_dp*hbar*hbar))
      endif        

      do ikLocal = 1, nkPerPool

        ikGlobal = ikLocal+ikStart_pool-1

        open(unit=37, file=trim(outputDir)//'transitionRate.'//trim(int2str(iSpin))//"."//trim(int2str(ikGlobal)))

        write(37,'(a)') trim(volumeLine)

        write(37,'("# Total number of transitions, Initial states (bandI, bandF) Format : ''(3i10)''")')
        write(37,'(3i10)') nTransitions, ibi(1), ibi(nTransitions)

        write(37,'("# Temperature (K): ", f7.1)') temperature

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
  function G0ExpArg(time) 

    ! Input variables:
    !integer, intent(in) :: nModes
      ! Number of phonon modes

    !real(kind=dp), intent(in) :: nj(nModes)
      ! \(n_j\) occupation number
    !real(kind=dp), intent(in) :: omega(nModes)
      ! Initial-state frequency for each mode
    !real(kind=dp), intent(in) :: omegaPrime(nModes)
      ! Optional final-state frequency for each mode
    !real(kind=dp), intent(in) :: Sj(nModes)
      ! Huang-Rhys factor for each mode with omega_j
    !real(kind=dp), intent(in) :: SjPrime(nModes)
      ! Huang-Rhys factor for each mode with omega_j'
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    !logical, intent(in) :: diffOmega
      ! If initial- and final-state frequencies 
      ! should be treated as different

    ! Output variables:
    complex(kind=dp) :: G0ExpArg

    ! Local variables
    real(kind=dp) :: sinOmegaPrime(nModes)
      !! sin(omega' t/2)
    complex(kind=dp) :: Aj_t(nModes)
      !! Aj from text. See equation below.
    complex(kind=dp) :: expTimesNBarPlus1(nModes)
      !! Local storage of \(e^{i\omega_j t}(\bar{n}_j + 1)\) 
    complex(kind=dp) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed


    posExp_t(:) = exp(ii*omega(:)*time)

    expTimesNBarPlus1(:) = posExp_t(:)*(nj(:)+1)

    if(diffOmega) then
      Aj_t(:) = (expTimesNBarPlus1(:) + nj(:))/(expTimesNBarPlus1(:) - nj(:))

      sinOmegaPrime(:) = sin(omegaPrime(:)*time/2)

      G0ExpArg = -2*ii*sum(sinOmegaPrime/(ii*Aj_t*sinOmegaPrime(:)/Sj(:) - cos(omegaPrime(:)*time/2)/SjPrime(:)))
    else
      G0ExpArg = sum(Sj(:)*(expTimesNBarPlus1(:) + nj(:)/posExp_t(:) - (2.0_dp*nj(:) + 1.0_dp)))
    endif

  end function G0ExpArg

!----------------------------------------------------------------------------
  function getAj_t(time) result(Aj_t)
    
    implicit none

    ! Input variables:
    !integer, intent(in) :: nModes
      ! Number of phonon modes

    !real(kind=dp), intent(in) :: nj(nModes)
      ! \(n_j\) occupation number
    !real(kind=dp), intent(in) :: omega(nModes)
      ! Frequency for each mode
    !real(kind=dp), intent(in) :: Sj(nModes)
      ! Huang-Rhys factor for each mode
    real(kind=dp), intent(in) :: time
      !! Time at which to calculate the \(G_0(t)\) argument

    ! Output variables:
    complex(kind=dp) :: Aj_t(nModes)

    ! Local variables
    complex(kind=dp) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed


    posExp_t(:) = exp(ii*omega(:)*time)

    Aj_t(:) = (2*hbar/omega(:))*(nj(:)/posExp_t(:) + (nj(:)+1)*posExp_t(:) + &
            Sj(:)*(1 + nj(:)/posExp_t(:) - (nj(:)+1)*posExp_t(:))**2)

  end function getAj_t

!----------------------------------------------------------------------------
  function getFirstOrderPrefactor(nModes, matrixElement, Aj_t) result(prefactor)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of phonon modes

    real(kind=dp), intent(in) :: matrixElement(nModes)
      !! Matrix elements for each mode, given a
      !! choice of band combination

    complex(kind=dp), intent(in) :: Aj_t(nModes)
      !! \(A_j(t)\) for all modes

    ! Output variables:
    complex(kind=dp) :: prefactor
      !! First-order exponential prefactor


    prefactor = sum(matrixElement(:)*Aj_t(:))

  end function getFirstOrderPrefactor
  
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
