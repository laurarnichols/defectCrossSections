program LSF0main
  use LSF0mod

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

  real(kind=dp) :: inputt, t1, t2, inta, intb, transitionRateGlobal
  complex(kind=dp) :: temp,tmpa1,tmpa2,G0_t,tmpb1,tmpb2, tmpa3, tmpb3


  call cpu_time(timerStart)

  call mpiInitialization()

  call readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, beta, dt, gamma0, gammaExpTolerance, maxTime, temperature, EInput, M0Input, &
        outputDir, SjInput)

  call readSj(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, SjInput, nModes, modeFreq, Sj)

  call readEnergy(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, EInput, dEDelta, dEPlot)

  call readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, M0Input, matrixElement)


  allocate(transitionRate(iBandIinit:iBandIfinal))
  transitionRate(:) = 0.0d0


  gamma0 = gamma0*meV/hbar
  
  nstep=200000000
  write(*,*) nstep*dt

  transitionRateGlobal=0.0d0

  do iTime = 1, nstep-1, 2
    inputt=(float(iTime))*dt

    t1=inputt + nstep*dt*nn
    tmpa2 = G0_t(t1)-gamma0*t1 ! Lorentz

    t2 = t1 + dt
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

      t1 = nstep*dt*nn
      t2 = (float(nstep))*dt+nstep*dt*nn

      transitionRate(ibi) = transitionRate(ibi) + Real(Mif*exp(G0_t(t1)+cmplx(0.0,t1*Eif/hbar,dp)-gamma0*t1)) ! Lorentz
      transitionRate(ibi) = transitionRate(ibi) - Real(Mif*exp(G0_t(t2)+cmplx(0.0,t2*Eif/hbar,dp)-gamma0*t2)) ! Lorentz

    enddo

    transitionRate(ibi) = transitionRate(ibi)*2.0d0/3.0d0*dt
    transitionRate(ibi) = transitionRate(ibi)/hbar/hbar

    write(*,'(i10, f10.5,f7.1,i5,ES35.14E3)') ibi, dEPlot, temperature, nn, transitionRate(ibi)

  enddo

end program LSF0main

