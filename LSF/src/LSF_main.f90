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

  if(ionode) then
    if(.not. captured .and. nPools /= 1) call exitError('LSFmain', 'Scattering only currently set up for nPools = 1 (-nk 1)!',1)
  endif


  if(ionode) write(*, '("Getting input parameters and reading necessary input files.")')
  call cpu_time(timer1)

  call readInputParams(iSpin, order, dt, gamma0, hbarGamma, maxTime, SjThresh, smearingExpTolerance, addDeltaNj, &
        captured, diffOmega, newEnergyTable, oldFormat, rereadDq, reSortMEs, dqInput, energyTableDir, matrixElementDir, &
        MjBaseDir, njBaseInput, optimalPairsInput, outputDir, PhononPPDir, prefix)


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


  ! Distribute k-points in pools (k-point parallelization only currently
  ! used for capture.)
  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
  if(.not. captured) nkPerPool = 1
    ! Other parameters ignored for scattering. Set nkPerPool manually
    ! because it is used for array dimensions.


  call readEnergyTable(iSpin, captured, energyTableDir, nTransitions, ibi, ibf, iki, ikf, dE)

  call readSj(ibi, ibf, iki, ikf, nTransitions, captured, diffOmega, PhononPPDir, nModes, omega, &
          omegaPrime, Sj, SjPrime)

  allocate(njBase(nModes))
  call readNj(nModes, njBaseInput, njBase)


  allocate(jReSort(nModes))

  if(ionode) then
    if(order == 1 .and. reSortMEs) then
      call getjReSort(nModes, optimalPairsInput, jReSort)
    else
      jReSort = 0
    endif
  endif

  call MPI_BCAST(jReSort, nModes, MPI_INTEGER, root, worldComm, ierr)


  call readAllMatrixElements(iSpin, nTransitions, ibi, nModes, jReSort, order, suffixLength, dE, captured, newEnergyTable, &
          oldFormat, rereadDq, reSortMEs, dqInput, matrixElementDir, MjBaseDir, prefix, mDim, matrixElement, volumeLine)

  deallocate(jReSort)


  call cpu_time(timer2)
  if(ionode) write(*, '("Reading input parameters and files complete! (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)

  if(ionode) write(*, '("Beginning transition-rate calculation")')
   

  call getAndWriteTransitionRate(nTransitions, ibi, ibf, iki, ikf, iSpin, mDim, nModes, order, dE, dt, &
          gamma0, matrixElement, njBase, omega, omegaPrime, Sj, SjPrime, SjThresh, captured, diffOmega, volumeLine)

  
  deallocate(ibi,ibf,iki,ikf)
  deallocate(dE)
  deallocate(matrixElement)
  deallocate(njBase)
  deallocate(omega)
  deallocate(omegaPrime)
  deallocate(Sj)
  deallocate(SjPrime)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(timerEnd)
  if(ionode) write(*,'("************ LSF complete! (",f10.2," secs) ************")') timerEnd-timerStart

  call MPI_FINALIZE(ierr)

end program LSFmain

