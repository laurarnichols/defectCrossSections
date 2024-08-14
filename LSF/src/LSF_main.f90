program LSFmain

  use LSFmod

  implicit none

  ! Local variables:
  real(kind=dp) :: timerStart, timerEnd, timer1, timer2
    !! Timers


  call cpu_time(timerStart)

  call mpiInitialization('LSF')

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools


  if(ionode) write(*, '("Getting input parameters and reading necessary input files.")')
  call cpu_time(timer1)

  call readInputParams(iSpin, order, dt, gamma0, hbarGamma, maxTime, SjThresh, smearingExpTolerance, diffOmega, &
        newEnergyTable, oldFormat, rereadDq, reSortMEs, energyTableDir, matrixElementDir, MjBaseDir, njInput, outputDir, &
        PhononPPDir, prefix)


  nStepsLocal = ceiling((maxTime/dt)/nProcPerPool)
    ! Would normally calculate the total number of steps as
    !         nSteps = ceiling(maxTime/dt)
    ! then divide the steps across all of the processes in 
    ! the pool. However, the number of steps could be a very 
    ! large integer that could cause overflow. Instead, directly
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
   
  if(ionode) write(*,'("  Each process is completing ", i15, " time steps.")') nStepsLocal


  ! Distribute k-points in pools
  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)


  call readEnergyTable(iSpin, energyTableDir, nTransitions, ibi, ibf, dE)

  call readSj(diffOmega, PhononPPDir, nModes, omega, omegaPrime, Sj, SjPrime)

  allocate(nj(nModes))
  call readNj(nModes, njInput, nj)


  allocate(jReSort(nModes))

  if(ionode) then
    if(order == 1 .and. reSortMEs) call getjReSort(nModes, PhononPPDir, jReSort)
  endif

  call MPI_BCAST(jReSort, nModes, MPI_INTEGER, root, worldComm, ierr)


  call readAllMatrixElements(iSpin, nTransitions, ibi, nModes, jReSort, order, suffixLength, dE, newEnergyTable, &
          oldFormat, rereadDq, reSortMEs, matrixElementDir, MjBaseDir, PhononPPDir, prefix, mDim, matrixElement, volumeLine)

  deallocate(jReSort)


  call cpu_time(timer2)
  if(ionode) write(*, '("Reading input parameters and files complete! (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)

  if(ionode) write(*, '("Beginning transition-rate calculation")')
   

  call getAndWriteTransitionRate(nTransitions, ibi, iSpin, mDim, nModes, order, dE, dt, gamma0, & 
          matrixElement, nj, omega, omegaPrime, Sj, SjPrime, SjThresh, diffOmega, volumeLine)

  
  deallocate(dE)
  deallocate(matrixElement)
  deallocate(nj)
  deallocate(omega)
  deallocate(omegaPrime)
  deallocate(Sj)
  deallocate(SjPrime)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(timerEnd)
  if(ionode) write(*,'("************ LSF complete! (",f10.2," secs) ************")') timerEnd-timerStart

  call MPI_FINALIZE(ierr)

end program LSFmain

