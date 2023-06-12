module LSF0mod
  
  use constants, only: dp, HartreeToJ, HartreeToEv, eVToJ
  use errorsAndMPI

  implicit none 
  real(kind=dp),parameter :: Kb =  1.38064852d-23
  real(kind=dp),parameter :: pi= 3.14159265358979
  real(kind=dp),parameter :: tpi = 6.2831853071795864769 
  real(kind=dp),parameter :: Thz = 1.0d12
  real(kind=dp),parameter :: hbar = 1.0545718d-34
  real(kind=dp),parameter :: q_comvert = 3.920205055d50


  integer :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !! Energy band bounds for initial and final state
  integer :: nModes
    !! Number of phonon modes

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
  real(kind=dp), allocatable :: matrixElement(:,:)
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
  character(len=300) :: M0Input
    !! Path to zeroth-order matrix element file
    !! `allElecOverlap.isp.ik`
  character(len=300) :: outputDir
    !! Path to output transition rates
  character(len=300) :: SjInput
    !! Path to Sj.out file


  integer :: nstep, nw, nn

  namelist /inputParams/ iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, EInput, M0input, SjInput, &
                        temperature, nn, hbarGamma, dt, smearingExpTolerance, outputDir

contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, beta, dt, gamma0, maxTime, smearingExpTolerance, temperature, &
        EInput, M0Input, outputDir, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

    real(kind=dp), intent(out) :: beta
      !! 1/kb*T
    real(kind=dp), intent(out) :: dt
      !! Time step size
    real(kind=dp), intent(out) :: gamma0
      !! \(\gamma\) for Lorentzian smearing
    real(kind=dp), intent(out) :: maxTime
      !! Max time for integration
    real(kind=dp), intent(out) :: smearingExpTolerance
      !! Tolerance for the Lorentzian-smearing
      !! exponential used to calculate max time
    real(kind=dp), intent(out) :: temperature

    character(len=300), intent(out) :: EInput
      !! Path to energy table to read
    character(len=300), intent(out) :: M0Input
      !! Path to zeroth-order matrix element file
      !! `allElecOverlap.isp.ik`
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(out) :: SjInput
      !! Path to Sj.out file

    ! Local variables:
    real(kind=dp) :: hbarGamma
      !! \(\hbar\gamma\) for Lorentzian smearing
      !! to guarantee convergence

  
    call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, dt, hbarGamma, smearingExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

    if(ionode) then

      read(5, inputParams, iostat=ierr)
        !! * Read input variables
    
      if(ierr /= 0) call exitError('LSF0 main', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, dt, hbarGamma, smearingExpTolerance, temperature, EInput, &
            M0Input, outputDir, SjInput)

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
  
    call MPI_BCAST(beta, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(hbarGamma, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(smearingExpTolerance, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(temperature, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(maxTime, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  
    call MPI_BCAST(EInput, len(EInput), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(M0Input, len(M0Input), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(SjInput, len(SjInput), MPI_CHARACTER, root, worldComm, ierr)
    
    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, dt, hbarGamma, smearingExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

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
    character(len=300), intent(out) :: M0Input
      !! Path to zeroth-order matrix element file
      !! `allElecOverlap.isp.ik`
    character(len=300), intent(out) :: outputDir
      !! Path to store transition rates
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

    dt = 1d-6
    hbarGamma = 0.0_dp
    smearingExpTolerance = 0.0_dp
    temperature = 0.0_dp

    EInput = ''
    M0Input = ''
    SjInput = ''
    outputDir = './'

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(*, '(/5X,"LSF0 starts on ",A9," at ",A9)') &
             cdate, ctime

      write(*, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs


    endif

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, dt, hbarGamma, smearingExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

    implicit none

    ! Input variables
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

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
    character(len=300), intent(in) :: M0Input
      !! Path to zeroth-order matrix element file
      !! `allElecOverlap.isp.ik`
    character(len=300), intent(in) :: outputDir
      !! Path to store transition rates
    character(len=300), intent(in) :: SjInput
      !! Path to Sj.out file

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9))
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 

    abortExecution = checkDoubleInitialization('dt', dt, 1.0d-10, 1.0d-4) .or. abortExecution
    abortExecution = checkDoubleInitialization('hbarGamma', hbarGamma, 0.0_dp, 20.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('smearingExpTolerance', smearingExpTolerance, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('temperature', temperature, 0.0_dp, 1500.0_dp) .or. abortExecution
      ! These limits are my best guess as to what is reasonable; they are not
      ! hard and fast, but you should think about the application of the theory
      ! to numbers outside these ranges.

    abortExecution = checkFileInitialization('EInput', Einput) .or. abortExecution
    abortExecution = checkFileInitialization('M0Input', M0input) .or. abortExecution
    abortExecution = checkFileInitialization('SjInput', Sjinput) .or. abortExecution

    call system('mkdir -p '//trim(outputDir))


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif

    return 

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readSj(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, SjInput, nModes, modeFreq, Sj)

    implicit none

    ! Input variables
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

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
  subroutine readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, M0Input, matrixElement)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

    character(len=300), intent(in) :: M0Input
      !! Path to zeroth-order matrix element file
      !! `allElecOverlap.isp.ik`

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: matrixElement(:,:)
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
      open(12,file=trim(M0Input))

      read(12,*)
      read(12,*) iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
        ! @todo Test these values against the input values

    endif
      

    allocate(matrixElement(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal))

    if(ionode) then

      do ibf = iBandFinit, iBandFfinal
        do ibi = iBandIinit, iBandIfinal

          read(12,*) iDum, iDum, rDum, rDum, rDum, matrixElement(ibf,ibi) ! in Hartree^2

        enddo
      enddo

      matrixElement(:,:) = matrixElement(:,:)*HartreeToJ**2

      close(12)

    endif

    call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine readMatrixElement

function G0_t(inputt) result(G0t) 
  integer :: ifreq
  real(kind=dp) :: inputt, nj, omega, e_factor
  complex(kind=dp) G0t, tmp1, tmp2
  tmp1 = 0.0_dp
  tmp2 = 0.0_dp
  do ifreq=1, nModes
   omega=modeFreq(ifreq)
   nj=1/(exp(hbar*omega*beta)-1)
   tmp1=tmp1+(nj+1)*Sj(ifreq)*exp(cmplx(0.0,omega*inputt,dp))+nj*Sj(ifreq)*exp(cmplx(0.0,-omega*inputt,dp))-(2*nj+1)*Sj(ifreq)
   !tmp2=tmp2+(nj+1)*Sj(ifreq)/e_factor*exp(cmplx(0.0,omega*inputt,dp))+nj*Sj(ifreq)*e_factor*exp(cmplx(0.0,-omega*inputt,dp))-(2*nj+1)*Sj(ifreq)
  end do
  G0t= tmp1! - exp(tmp2)
!time-dependent line-shape-function, expand it to the first order of Sj, you will get familiar form with Huang-Rhys factor
  !!G0t=G0t*exp(-0.25*(alpha*inputt)**2)
  !narrow Gaussian to simulate delta function, what we care about is the area, we make it narrow so neighboring LSF do not interfere
end function

end module LSF0mod
