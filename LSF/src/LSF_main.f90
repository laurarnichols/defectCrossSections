program LSFmain

  use LSFmod

  implicit none

  ! Local variables:
  integer :: j, ikLocal, ikGlobal
    !! Loop index
  integer :: mDim
    !! Size of first dimension for matrix element

  real(kind=dp) :: timerStart, timerEnd, timer1, timer2
    !! Timers

  character(len=300) :: fName
    !! File name for first-order matrix elements


  call cpu_time(timerStart)

  call mpiInitialization('LSF')

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools


  if(ionode) write(*, '("Pre-k-loop: [ ] Get parameters  [ ] Read Sj")')
  call cpu_time(timer1)

  call readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, order, beta, dt, gamma0, hbarGamma, maxTime, &
        smearingExpTolerance, temperature, energyTableDir, matrixElementDir, MjBaseDir, outputDir, prefix, SjInput)


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


  call cpu_time(timer2)
  if(ionode) write(*, '("Pre-k-loop: [X] Get parameters  [ ] Read Sj (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)

  call readSj(SjInput, nModes, omega, Sj)

  allocate(nj(nModes))
  nj(:) = 1.0_dp/(exp(hbar*omega(:)*beta) - 1.0_dp)


  call cpu_time(timer2)
  if(ionode) write(*, '("Pre-k-loop: [X] Get parameters  [X] Read Sj (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools


  allocate(dE(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal,4,nkPerPool))

  if(order == 0) then
    mDim = 1
  else if(order == 1) then
    mDim = nModes
  endif

  allocate(matrixElement(mDim,iBandFinit:iBandFfinal,iBandIinit:iBandIfinal,nkPerPool))

  if(indexInPool == 0) then

    dE = 0.0_dp
    matrixElement = 0.0_dp

    do ikLocal = 1, nkPerPool
    
      ikGlobal = ikLocal+ikStart_pool-1
        !! Get the global `ik` index from the local one

      if(transitionRateFileExists(ikGlobal, iSpin)) then
        write(*, '("Transition rate of k-point ", i4, " and spin ", i1, " already exists.")') ikGlobal, iSpin
      endif
      !else
        ! Skipping files that exist isn't implemented in the transition
        ! rate calculation yet, so we need to read all of the files.
        ! In the future, can just skip the ones that have already been 
        ! calculated.

        call readEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, iSpin, energyTableDir, dE(:,:,:,ikLocal))

        if(order == 0) then
          ! Read zeroth-order matrix element

          fName = getMatrixElementFNameWPath(ikGlobal, iSpin, matrixElementDir)

          call readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, fName, matrixElement(1,:,:,ikLocal), volumeLine)
            ! Includes conversion from Hartree to J

        else if(order == 1) then
          ! Read matrix elements for all modes
    
          do j = 1, nModes

            fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(j,4))//'/'//trim(getMatrixElementFNameWPath(ikGlobal,iSpin,matrixElementDir))

            call readMatrixElement(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, order, fName, matrixElement(j,:,:,ikLocal), volumeLine)
              ! Includes conversion from Hartree to J
              !
              ! The volume line will get overwritten each time, but that's
              ! okay because the volume doesn't change between the files. 

          enddo
  
        !endif
      endif
    enddo

    dE(:,:,1:3,:) = dE(:,:,1:3,:)*HartreeToJ
      ! First 3 columns are in Hartree. Last is in eV,
      ! but we want to leave it that way just for output.

    if(order == 1) then
      matrixElement(:,:,:,:) = matrixElement(:,:,:,:)/(BohrToMeter*sqrt(elecMToKg))**2
        ! Shifter program outputs dq in Bohr*sqrt(elecM), and that
        ! dq is directly used by TME to get the matrix element, so
        ! we need to convert to m*sqrt(kg). In the matrix element,
        ! dq is in the denominator and squared.
    endif

  endif


  call MPI_BCAST(dE, size(dE), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)
  call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)


  if(ionode) write(*, '("  Beginning transition-rate calculation")')
   

  call getAndWriteTransitionRate(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, iSpin, mDim, order, nModes, dE, &
        gamma0, matrixElement, nj, temperature, volumeLine)

  
  deallocate(dE)
  deallocate(matrixElement)
  deallocate(nj)
  deallocate(omega)
  deallocate(Sj)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(timerEnd)
  if(ionode) write(*,'("************ LSF complete! (",f10.2," secs) ************")') timerEnd-timerStart

  call MPI_FINALIZE(ierr)

end program LSFmain

