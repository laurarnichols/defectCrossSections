program LSFmain

  use LSFmod
  use miscUtilities, only: int2strLeadZero

  implicit none

  ! Local variables:
  integer :: j
    !! Loop index
  integer :: mDim
    !! Size of first dimension for matrix element

  real(kind=dp) :: timerStart, timerEnd, timer1, timer2
    !! Timers

  character(len=300) :: fName
    !! File name for first-order matrix elements


  call cpu_time(timerStart)

  call mpiInitialization()


  if(ionode) write(*, '("Reading inputs: [ ] Parameters  [ ] Sj  [ ] dE  [ ] Matrix elements")')
  call cpu_time(timer1)

  call readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, beta, dt, gamma0, hbarGamma, maxTime, &
        smearingExpTolerance, temperature, EInput, MifInput, MjDir, outputDir, prefix, SjInput)


  call cpu_time(timer2)
  if(ionode) write(*, '("Reading inputs: [X] Parameters  [ ] Sj  [ ] dE  [ ] Matrix elements (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)

  call readSj(SjInput, nModes, omega, Sj)


  call cpu_time(timer2)
  if(ionode) write(*, '("Reading inputs: [X] Parameters  [X] Sj  [ ] dE  [ ] Matrix elements (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)

  call readEnergy(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, EInput, dEDelta, dEPlot)


  call cpu_time(timer2)
  if(ionode) write(*, '("Reading inputs: [X] Parameters  [X] Sj  [X] dE  [ ] Matrix elements (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)

  if(order == 0) then
    ! Read single zeroth-order matrix element

    mDim = 1
    allocate(matrixElement(mDim,iBandFinit:iBandFfinal,iBandIinit:iBandIfinal))
    
    call readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, MifInput, matrixElement(1,:,:), volumeLine)

  else if(order == 1) then
    ! Read matrix elements for all modes

    mDim = nModes
    allocate(matrixElement(mDim,iBandFinit:iBandFfinal,iBandIinit:iBandIfinal))
    
    do j = 1, nModes

      fName = trim(MjDir)//'/'//trim(prefix)//trim(int2strLeadZero(j,4))//'/'//trim(MifInput)

      call readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, trim(fName), matrixElement(j,:,:), volumeLine)
        ! The volume line will get overwritten each time, but that's
        ! okay because the volume doesn't change between the files. 

    enddo

    matrixElement(:,:,:) = matrixElement(:,:,:)/(BohrToMeter*sqrt(elecMToKg))**2
      ! Shifter program outputs dq in Bohr*sqrt(elecM), and that
      ! dq is directly used by TME to get the matrix element, so
      ! we need to convert to m*sqrt(kg). In the matrix element,
      ! dq is in the denominator and squared.

  endif

  call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
   

  call cpu_time(timer2)
  if(ionode) write(*, '("Reading inputs: [X] Parameters  [X] Sj  [X] dE  [X] Matrix elements (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)


  nStepsLocal = ceiling((maxTime/dt)/nProcs)
    ! Would normally calculate the total number of steps as
    !         nSteps = ceiling(maxTime/dt)
    ! then divide the steps across all of the processes. 
    ! However, the number of steps could be a very large
    ! integer that could cause overflow. Instead, directly
    ! calculate the number of steps that each process should 
    ! calculate. If you still get integer overflow, try 
    ! increasing the number of processes. 
    !
    ! Calculating the number of steps for each process
    ! this way will overestimate the number of steps
    ! needed, but that is okay.

  if(nStepsLocal < 0) call exitError('LSF main', 'integer overflow', 1)
    ! If there is integer overflow, the number will go to the
    ! most negative integer value available


  if(mod(nStepsLocal,2) == 0) nStepsLocal = nStepsLocal + 1
    ! Simpson's method requires the number of integration
    ! steps to be odd because a 3-point quadratic
    ! interpolation is used
   
  if(ionode) write(*,'("Each process is completing ", i15, " time steps.")') nStepsLocal


  call getAndWriteTransitionRate(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, mDim, order, nModes, dEDelta, &
        dEPlot, gamma0, matrixElement, temperature, volumeLine)

  
  deallocate(dEDelta)
  deallocate(dEPlot)
  deallocate(matrixElement)
  deallocate(nj)
  deallocate(omega)
  deallocate(Sj)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(timerEnd)
  if(ionode) write(*,'("************ LSF complete! (",f10.2," secs) ************")') timerEnd-timerStart

  call MPI_FINALIZE(ierr)

end program LSFmain

