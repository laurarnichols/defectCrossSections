module capcs
  
  use constants, only: dp, HartreeToJ, HartreeToEv
  use errorsAndMPI

  implicit none 
  real(kind=dp),parameter :: Kb =  1.38064852d-23
  real(kind=dp),parameter :: pi= 3.14159265358979
  real(kind=dp),parameter :: tpi = 6.2831853071795864769 
  real(kind=dp),parameter :: Thz = 1.0d12
  real(kind=dp),parameter :: hbar = 1.0545718d-34
  real(kind=dp),parameter :: eV = 1.6021766208d-19
  real(kind=dp),parameter :: meV =  1.6021766208d-22
  real(kind=dp),parameter :: q_comvert = 3.920205055d50


  integer :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !! Energy band bounds for initial and final state

  real(kind=dp) :: beta
    !! 1/kb*T
  real(kind=dp), allocatable :: dEDelta(:,:)
    !! Energy for delta function
  real(kind=dp), allocatable :: dEPlot(:)
    !! Energy for plotting
  real(kind=dp) :: gamma0
    !! Gamma parameter to use for Lorentzian smearing
    !! to guarantee convergence
  real(kind=dp) :: gammaExpTolerance
    !! Tolerance for the gamma exponential (Lorentzian
    !! smearing) used to calculate max time
  real(kind=dp), allocatable :: matrixElement(:,:)
    !! Electronic matrix element
  real(kind=dp) :: maxTime
    !! Max time for integration
  real(kind=dp), allocatable :: Sj(:)
    !! Huang-Rhys factor for each mode
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


  real(kind=dp),allocatable :: ipfreq(:)
  real(kind=dp) :: dstep
  integer :: nstep, nw, nfreq, nn

  namelist /inputParams/ iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, EInput, M0input, SjInput, &
                        temperature, nn, gamma0, dstep, gammaExpTolerance, outputDir

contains

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, gamma0, gammaExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

    implicit none

    ! Output variables
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

    real(kind=dp), intent(out) :: gamma0
      !! Gamma parameter to use for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(out) :: gammaExpTolerance
      !! Tolerance for the gamma exponential (Lorentzian
      !! smearing) used to calculate max time
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

    gamma0 = 0.0_dp
    gammaExpTolerance = 0.0_dp
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

  end subroutine

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, gamma0, gammaExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

    implicit none

    ! Output variables
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state

    real(kind=dp), intent(in) :: gamma0
      !! Gamma parameter to use for Lorentzian smearing
      !! to guarantee convergence
    real(kind=dp), intent(in) :: gammaExpTolerance
      !! Tolerance for the gamma exponential (Lorentzian
      !! smearing) used to calculate max time
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

    abortExecution = checkDoubleInitialization('gamma0', gamma0, 0.0_dp, 20.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('gammaExpTolerance', gammaExpTolerance, 0.0_dp, 1.0_dp) .or. abortExecution
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

  end subroutine

!----------------------------------------------------------------------------
  subroutine readEnergy(dEDelta, dEPlot)
  
    implicit none

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


    open(12,file=trim(EInput))

    read(12,*)
    read(12,*) iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      ! @todo Test these values against the input values
      

    allocate(dEDelta(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal))
    allocate(dEPlot(iBandIinit:iBandIfinal))

    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal

        read(12,*) iDum, iDum, dEDelta(ibf,ibi), rDum, rDum, dEPlot(ibi) ! in Hartree

      enddo
    enddo

    dEDelta(:,:) = dEDelta(:,:)*HartreeToJ
    dePlot(:) = dEPlot(:)*HartreeToEv

    close(12)

    return

  end subroutine readEnergy

!----------------------------------------------------------------------------
  subroutine readMatrixElement(matrixElement)

    implicit none

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

    open(12,file=trim(M0Input))

    read(12,*)
    read(12,*) iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      ! @todo Test these values against the input values
      

    allocate(matrixElement(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal))

    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal

        read(12,*) iDum, iDum, rDum, rDum, rDum, matrixElement(ibf,ibi) ! in Hartree^2

      enddo
    enddo

    matrixElement(:,:) = matrixElement(:,:)*HartreeToJ**2

    close(12)

    return

  end subroutine readMatrixElement

!----------------------------------------------------------------------------
subroutine readphonon()
implicit none
integer :: ifreq, modeIndex
open(12,file=trim(SjInput))
read(12,*)nfreq
!read mode number
allocate(Sj(1:nfreq))
allocate(ipfreq(1:nfreq))

do ifreq=1,nfreq
   read(12,*)modeIndex,Sj(ifreq),ipfreq(ifreq) !freq read from Sj.out is f(in Thz)*2pi
end do
!read frequency and Sj
ipfreq=ipfreq*Thz!*tpi convert to Hz*2pi
!transform frequency to Hz
close(12)
!write(*,*) "Frequency Read"
end subroutine

end module capcs



function G0_t(inputt) result(G0t) 
use capcs
  integer :: ifreq
  real(kind=dp) :: inputt, nj, omega, e_factor
  complex(kind=dp) G0t, tmp1, tmp2
  tmp1 = 0.0_dp
  tmp2 = 0.0_dp
  do ifreq=1, nfreq
   omega=ipfreq(ifreq)
   nj=1/(exp(hbar*omega*beta)-1)
   tmp1=tmp1+(nj+1)*Sj(ifreq)*exp(cmplx(0.0,omega*inputt,dp))+nj*Sj(ifreq)*exp(cmplx(0.0,-omega*inputt,dp))-(2*nj+1)*Sj(ifreq)
   !tmp2=tmp2+(nj+1)*Sj(ifreq)/e_factor*exp(cmplx(0.0,omega*inputt,dp))+nj*Sj(ifreq)*e_factor*exp(cmplx(0.0,-omega*inputt,dp))-(2*nj+1)*Sj(ifreq)
  end do
  G0t= tmp1! - exp(tmp2)
!time-dependent line-shape-function, expand it to the first order of Sj, you will get familiar form with Huang-Rhys factor
  !!G0t=G0t*exp(-0.25*(alpha*inputt)**2)
  !narrow Gaussian to simulate delta function, what we care about is the area, we make it narrow so neighboring LSF do not interfere
end function

program captureCS
  use capcs
  use OMP_LIB

  implicit none

  ! Local variables:
  integer :: iTime, ibi, ibf
    !! Loop indices

  real(kind=dp) :: Eif
    !! Local storage of dEDelta(ibi,ibf)
  real(kind=dp) :: Mif
    !! Local storage of matrixElement(ibi,ibf)
  real(kind=dp), allocatable :: transitionRate(:)
    !! \(Gamma_i\) transition rate
  real(kind=dp) :: timerStart, timerEnd
    !! Timers

  real(kind=dp) :: dtime, inputt, t1, t2, inta, intb, transitionRateGlobal
  complex(kind=dp) :: temp,tmpa1,tmpa2,G0_t,tmpb1,tmpb2, tmpa3, tmpb3


  call cpu_time(timerStart)

  call mpiInitialization()

  call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, gamma0, gammaExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

  if(ionode) then

    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('LSF0 main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

    call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, gamma0, gammaExpTolerance, temperature, EInput, M0Input, outputDir, SjInput)

    beta = 1.0d0/Kb/temperature

    maxTime = -log(gammaExpTolerance)/gamma0
    write(*,'("Max time: ", ES24.15E3)') maxTime

    call readphonon()
    call readEnergy(dEDelta, dEPlot)
    call readMatrixElement(matrixElement)

  endif

  allocate(transitionRate(iBandIinit:iBandIfinal))
  transitionRate(:) = 0.0d0

  !!limit = 1e-4
  dtime = 1/Thz
  write(*,*)"simps"
  write(*,*)dstep
  dstep = dtime*dstep

  gamma0 = gamma0*meV/hbar
  
  nstep=200000000
  write(*,*) nstep*dstep

  transitionRateGlobal=0.0d0

!$OMP PARALLEL DO default(shared) private(inputt,t1,t2,iTime,ibi,ibf,tmpa1,tmpa2,tmpa3,tmpb1,tmpb2,tmpb3)&
& firstprivate(Eif,dstep,nstep) reduction(+:transitionRateGlobal,transitionRate) 
  do iTime = 1, nstep-1, 2
    inputt=(float(iTime))*dstep

    t1=inputt + nstep*dstep*nn
    tmpa2 = G0_t(t1)-gamma0*t1 ! Lorentz

    t2 = t1 + dstep
    tmpb2 = G0_t(t2)-gamma0*t2 ! Lorentz

    do ibi = iBandIinit, iBandIfinal
      do ibf = iBandFinit, iBandFfinal

        Mif = matrixElement(ibf,ibi)
        Eif = dEDelta(ibf,ibi)

        tmpa3 = tmpa2+cmplx(0.0,t1*Eif/hbar,dp)
        tmpb3 = tmpb2+cmplx(0.0,t2*Eif/hbar,dp)

        transitionRate(ibi) = transitionRate(ibi) + Real(4d0*Mif*exp(tmpa3) + 2d0*Mif*exp(tmpb3))
          ! We are doing multiple sums, but they are all commutative.
          ! Here we add in the contribution to the integral at this time
          ! step from a given final state. The loop over final states 
          ! adds in the contributions from all final states. 
      enddo
    enddo
  enddo

  do ibi = iBandIinit, iBandIfinal
    do ibf = iBandFinit, iBandFfinal

      Mif = matrixElement(ibf,ibi)
      Eif = dEDelta(ibf,ibi)

      t1 = nstep*dstep*nn
      t2 = (float(nstep))*dstep+nstep*dstep*nn

      transitionRate(ibi) = transitionRate(ibi) + Real(Mif*exp(G0_t(t1)+cmplx(0.0,t1*Eif/hbar,dp)-gamma0*t1)) ! Lorentz
      transitionRate(ibi) = transitionRate(ibi) - Real(Mif*exp(G0_t(t2)+cmplx(0.0,t2*Eif/hbar,dp)-gamma0*t2)) ! Lorentz

    enddo

    transitionRate(ibi) = transitionRate(ibi)*2.0d0/3.0d0*dstep
    transitionRate(ibi) = transitionRate(ibi)/hbar/hbar

    write(*,'(i10, f10.5,f7.1,i5,ES35.14E3)') ibi, dEPlot, temperature, nn, transitionRate(ibi)

  enddo

end program
