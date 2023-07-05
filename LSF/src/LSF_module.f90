module LSFmod
  
  use constants, only: dp, HartreeToJ, HartreeToEv, eVToJ, ii, hbar, THzToHz, kB, BohrToMeter, elecMToKg
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nKPoints, order
  use TMEmod, only: getMatrixElementFNameWPath, getMatrixElementFName
  use miscUtilities, only: int2strLeadZero, int2str
  use errorsAndMPI

  use energyTabulatorMod, only: energyTableDir, readEnergyTable

  implicit none 

  integer :: iTime_start, iTime_end
    !! Start and end time steps for this process
  integer :: nStepsLocal
    !! Number of time steps for each
    !! process to complete


  integer :: iSpin
    !! Spin channel to use
  integer :: nModes
    !! Number of phonon modes

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
  real(kind=dp),allocatable :: omega(:)
    !! Frequency for each mode
  real(kind=dp), allocatable :: Sj(:)
    !! Huang-Rhys factor for each mode
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


  namelist /inputParams/ iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, energyTableDir, matrixElementDir, MjBaseDir, SjInput, &
                        temperature, hbarGamma, dt, smearingExpTolerance, outputDir, order, prefix, iSpin

contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, beta, dt, gamma0, hbarGamma, maxTime, &
        smearingExpTolerance, temperature, energyTableDir, matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
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

  
    call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, dt, hbarGamma, smearingExpTolerance, temperature, energyTableDir, &
          matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

    if(ionode) then

      read(5, inputParams, iostat=ierr)
        !! * Read input variables

    
      if(ierr /= 0) call exitError('LSF module', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, dt, hbarGamma, smearingExpTolerance, temperature, energyTableDir, &
            matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

      dt = dt/THzToHz

      gamma0 = hbarGamma*1e-3*eVToJ/hbar
        ! Input expected in meV

      beta = 1.0d0/(kB*temperature)

      maxTime = -log(smearingExpTolerance)/gamma0
      write(*,'("Max time: ", ES24.15E3)') maxTime

    endif

    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)
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
  
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(matrixElementDir, len(matrixElementDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(MjBaseDir, len(MjBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(SjInput, len(SjInput), MPI_CHARACTER, root, worldComm, ierr)
    
    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, dt, hbarGamma, smearingExpTolerance, temperature, energyTableDir, &
        matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
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

    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time


    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    iSpin = 1
    order = -1

    dt = 1d-6
    hbarGamma = 0.0_dp
    smearingExpTolerance = 0.0_dp
    temperature = 0.0_dp

    energyTableDir = ''
    matrixElementDir = ''
    MjBaseDir = ''
    SjInput = ''
    outputDir = './'
    prefix = 'disp-'

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(*, '(/5X,"LSF starts on ",A9," at ",A9)') &
             cdate, ctime

      write(*, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs


    endif

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, dt, hbarGamma, smearingExpTolerance, temperature, &
        energyTableDir, matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)

    implicit none

    ! Input variables
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
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


    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9))
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 
    abortExecution = checkIntInitialization('iSpin', iSpin, 1, 2) .or. abortExecution 
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
  subroutine readSj(SjInput, nModes, omega, Sj)

    implicit none

    ! Input variables
    character(len=300), intent(in) :: SjInput
      !! Path to Sj.out file

    ! Output variables:
    integer, intent(out) :: nModes
      !! Number of phonon modes

    real(kind=dp),allocatable, intent(out) :: omega(:)
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
    allocate(omega(1:nModes))

    
    if(ionode) then

      do j = 1, nModes
        read(12,*) iDum, Sj(j), omega(j) ! freq read from Sj.out is f(in Thz)*2pi
      end do

      omega(:) = omega(:)*THzToHz
        ! Convert to Hz*2pi

      close(12)

    endif

    call MPI_BCAST(omega, size(omega), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(Sj, size(Sj), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return 

  end subroutine readSj

!----------------------------------------------------------------------------
  subroutine readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, fName, matrixElement, volumeLine)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    character(len=300), intent(in) :: fName
      !! Path to matrix element file `allElecOverlap.isp.ik`

    ! Output variables:
    real(kind=dp), intent(out) :: matrixElement(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! Electronic matrix element

    character(len=300), intent(out) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      !! Band bounds from energy table
    integer :: iDum
      !! Dummy integer
    integer :: ibi, ibf
      !! Loop indices

    real(kind=dp) :: rDum
      !! Dummy real


    if(indexInPool == 0) then
      open(12,file=trim(fName))

      read(12,*)
      read(12,*)

      read(12,'(a)') volumeLine

      read(12,*)
      read(12,*) iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
        ! @todo Test these values against the input values
      read(12,*)

      if(order == 1) read(12,*)
        ! Ignore additional line for phonon mode 

    endif
      

    if(indexInPool == 0) then

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
  subroutine getAndWriteTransitionRate(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, iSpin, mDim, order, nModes, dE, &
        gamma0, matrixElement, nj, temperature, volumeLine)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikGlobal
      !! K-point index
    integer, intent(in) :: iSpin
      !! Spin index
    integer, intent(in) :: mDim
      !! Size of first dimension for matrix element
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dE(4,iBandFinit:iBandFfinal,iBandIinit:iBandIFinal)
      !! All energy differences from energy table
    real(kind=dp), intent(in) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(in) :: matrixElement(mDim,iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! Electronic matrix element
    real(kind=dp), intent(in) :: nj(nModes)
      !! \(n_j\) occupation number
    real(kind=dp), intent(in) :: temperature

    character(len=300), intent(in) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: iTime, ibi, ibf
      !! Loop indices
    integer :: updateFrequency
      !! Frequency of steps to write status update

    real(kind=dp) :: Eif
      !! Local storage of energy for delta function
    real(kind=dp) :: expPrefactor_t1, expPrefactor_t2
      !! Prefactor for the exponential inside the integral.
      !! For zeroth-order, this is just the matrix element,
      !! but for first-order this is \(\sum_j M_j A_j\).
    real(kind=dp) :: t0
      !! Initial time for this process
    real(kind=dp) :: t1, t2
      !! Time for each time step
    real(kind=dp) :: timer1, timer2
      !! Timers
    real(kind=dp) :: transitionRate(iBandIinit:iBandIfinal)
      !! \(Gamma_i\) transition rate

    complex(kind=dp) :: Aj_t1(nModes), Aj_t2(nModes)
      !! \(A_j(t)\) for all modes
    complex(kind=dp) :: expArg_t1_base, expArg_t2_base
      !! Base exponential argument for each time step
    complex(kind=dp) :: expArg_t1, expArg_t2
      !! Exponential argument for each time step

      
    updateFrequency = ceiling(nStepsLocal/10.0)
    call cpu_time(timer1)


    t0 = myid*float(nStepsLocal)*dt

    transitionRate(:) = 0.0_dp
    do iTime = 1, nStepsLocal-2, 2
      ! Loop should not include the first step (0), 
      ! but it should include the last step (nStepsLocal-1).
      ! iTime stops at nStepsLocal-2 because nStepsLocal-1
      ! is calculated by t2 at the last step.

      if(ionode .and. (mod(iTime,updateFrequency) == 0 .or. mod(iTime+1,updateFrequency) == 0)) then

        call cpu_time(timer2)
        write(*,'("    ", i2,"% complete with transition-rate loop. Time in loop: ",f10.2," secs")') iTime*100/nStepsLocal, timer2-timer1

      endif

      t1 = t0 + float(iTime)*dt
        ! Must do this arithmetic with floats to avoid
        ! integer overflow
      expArg_t1_base = G0ExpArg(t1) - gamma0*t1
      if(order == 1) Aj_t1(:) = getAj_t(t1)

      t2 = t1 + dt
      expArg_t2_base = G0ExpArg(t2) - gamma0*t2
      if(order == 1) Aj_t2(:) = getAj_t(t2)

      do ibi = iBandIinit, iBandIfinal
        do ibf = iBandFinit, iBandFfinal

          if(order == 0) then
            expPrefactor_t1 = matrixElement(1,ibf,ibi)
            expPrefactor_t2 = matrixElement(1,ibf,ibi)
          else if(order == 1) then
            expPrefactor_t1 = getFirstOrderPrefactor(nModes, matrixElement(:,ibf,ibi), Aj_t1)
            expPrefactor_t2 = getFirstOrderPrefactor(nModes, matrixElement(:,ibf,ibi), Aj_t2)
          endif

          Eif = dE(1,ibf,ibi)

          expArg_t1 = expArg_t1_base + ii*Eif/hbar*t1
          expArg_t2 = expArg_t2_base + ii*Eif/hbar*t2

          transitionRate(ibi) = transitionRate(ibi) + &
            Real(4.0_dp*expPrefactor_t1*exp(expArg_t1) + 2.0_dp*expPrefactor_t2*exp(expArg_t2))
            ! This is the Simpson's method integration. We skip step 0 and
            ! add 2*the last point. The endpoints get adjusted below.
            !
            ! We are doing multiple sums (integral and sum over final states), 
            ! but they are all commutative. Here we add in the contribution 
            ! to the integral at this time step from a given final state. The 
            ! loop over final states adds in the contributions from all final 
            ! states. 

        enddo
      enddo
    enddo

    

    ! Need to adjust endpoints of the integration
    t1 = t0
    expArg_t1_base = G0ExpArg(t1) - gamma0*t1
    if(order == 1) Aj_t1(:) = getAj_t(t1)

    t2 = t0 + float(nStepsLocal-1)*dt
    expArg_t2_base = G0ExpArg(t2) - gamma0*t2
    if(order == 1) Aj_t2(:) = getAj_t(t2)


    do ibi = iBandIinit, iBandIfinal
      do ibf = iBandFinit, iBandFfinal

        if(order == 0) then
          expPrefactor_t1 = matrixElement(1,ibf,ibi)
          expPrefactor_t2 = matrixElement(1,ibf,ibi)
        else if(order == 1) then
          expPrefactor_t1 = getFirstOrderPrefactor(nModes, matrixElement(:,ibf,ibi), Aj_t1)
          expPrefactor_t2 = getFirstOrderPrefactor(nModes, matrixElement(:,ibf,ibi), Aj_t2)
        endif

        Eif = dE(1,ibf,ibi)

        expArg_t1 = expArg_t1_base + ii*Eif/hbar*t1
        expArg_t2 = expArg_t2_base + ii*Eif/hbar*t2

        transitionRate(ibi) = transitionRate(ibi) + Real(expPrefactor_t1*exp(expArg_t1))
          ! Add \(t_0\) that was skipped in the loop

        transitionRate(ibi) = transitionRate(ibi) - Real(expPrefactor_t2*exp(expArg_t2)) 
          ! Subtract off last time that had a coefficient
          ! of 2 in the loop but should really have a 
          ! coefficient of 1

      enddo
    enddo

    call mpiSumDoubleV(transitionRate, intraPoolComm)
      ! Combine results from all processes in pool

    if(indexInPool == 0) then

      open(unit=37, file=trim(outputDir)//'transitionRate.'//trim(int2str(iSpin))//"."//trim(int2str(ikGlobal)))

      write(37,'(a)') trim(volumeLine)

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
        write(37,'(i10, f10.5,ES24.14E3)') ibi, dE(4,iBandFinit,ibi), transitionRate(ibi)
          ! Plotting energy doesn't depend on final band, so
          ! just pick one
      enddo

      write(*, '("  Transition rate of k-point ", i4, " and spin ", i1, " written.")') ikGlobal, iSpin

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

    real(kind=dp) :: omegaj
      !! Local storage of frequency for this mode

    complex(kind=dp) :: posExp_t(nModes)
      !! Local storage of \(e^{i\omega_j t}\) for speed


    posExp_t(:) = exp(ii*omega(:)*time)

    G0ExpArg = sum(Sj(:)*((nj(:)+1.0_dp)*posExp_t(:) + nj(:)/posExp_t(:) - (2.0_dp*nj(:) + 1.0_dp)))

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
    real(kind=dp) :: omegaj
      !! Local storage of frequency for this mode

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

    ! Local variables:
    integer :: j
      !! Loop index


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


    inquire(file=trim(outputDir)//'transitionRate.'//trim(int2str(iSpin))//"."//trim(int2str(ikGlobal)), exist=fileExists)
    
  end function transitionRateFileExists

end module LSFmod
