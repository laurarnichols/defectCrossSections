program LSFmain

  use LSFmod

  implicit none

  ! Local variables:
  integer, allocatable :: iDum1D(:)
    !! Integer to ignore input
  integer :: iDum1, iDum2
    !! Integers to ignore input
  integer :: j, ikLocal, ikGlobal
    !! Loop index
  integer :: mDim
    !! Size of first dimension for matrix element

  real(kind=dp), allocatable :: ME_tmp(:)
    !! Temporary storage of matrix element
  !real(kind=dp), allocatable :: randVal(:)
    !! Random adjustment to be made to frequencies
    !! to test sensitivity
  real(kind=dp), allocatable :: rDum2D(:,:)
    !! Dummy real to ignore input
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

  call readInputParams(iSpin, order, dt, gamma0, hbarGamma, maxTime, smearingExpTolerance, diffOmega, energyTableDir, &
        matrixElementDir, MjBaseDir, njInput, outputDir, prefix, SjInput)


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

  if(diffOmega) then
    call readSjTwoFreq(SjInput, nModes, omega, omegaPrime, Sj, SjPrime)

    ! This is the code that I used to test what difference different
    ! frequencies would have. I input the same two frequencies twice
    ! and randomly adjusted the omegaPrime. 
    !
    ! If doing a random adjustment, must do with a single process then
    ! broadcast, otherwise the random variables will be different across
    ! processes.
    !if(ionode) then
      !call random_seed()
      !allocate(randVal(nModes))
      !call random_number(randVal)
      !randVal = (randVal*2.0_dp - 1.0_dp)*0.50_dp  ! the number multiplied here is the % adjustment
      !omegaPrime(:) = omegaPrime(:)*(1.0_dp + randVal) ! this applies the random adjustment
      !omegaPrime(:) = omegaPrime(:)*(1.0_dp + 0.5_dp) ! this applies a uniform adjustment
      !deallocate(randVal)
    !endif

    ! Need to rebroadcast omegaPrime if we adjust it
    !call MPI_BCAST(omegaPrime, nModes, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  else
    allocate(omegaPrime(nModes), SjPrime(nModes))
      ! Need to allocate to avoid issues with passing variables
      ! and deallocating

    call readSjOneFreq(SjInput, nModes, omega, Sj)
  endif

  allocate(nj(nModes))
  call readNj(nModes, njInput, nj)


  call cpu_time(timer2)
  if(ionode) write(*, '("Pre-k-loop: [X] Get parameters  [X] Read Sj (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools


  if(ionode) then
    call readCaptureEnergyTable(1, iSpin, energyTableDir, ibi, ibf, nTransitions, rDum2D)
     ! Assume that band bounds and number of transitions do not depend on k-points or spin
     ! We ignore the energy here because we are only reading ibi, ibf, and nTransitions

    deallocate(rDum2D)
  endif

  call MPI_BCAST(nTransitions, 1, MPI_INTEGER, root, worldComm, ierr)
  if(.not. ionode) allocate(ibi(nTransitions))
  call MPI_BCAST(ibi, nTransitions, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(ibf, 1, MPI_INTEGER, root, worldComm, ierr)


  allocate(dE(3,nTransitions,nkPerPool))

  if(order == 0) then
    mDim = 1
  else if(order == 1) then
    mDim = nModes
  endif

  allocate(matrixElement(mDim,nTransitions,nkPerPool))
  allocate(ME_tmp(nTransitions))

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

        call readCaptureEnergyTable(ikGlobal, iSpin, energyTableDir, iDum1D, iDum1, iDum2, rDum2D)
          ! Assume that band bounds and number of transitions do not depend on k-points or spin

        dE(:,:,ikLocal) = rDum2D
        deallocate(rDum2D)

        if(order == 0) then
          ! Read zeroth-order matrix element

          fName = getMatrixElementFNameWPath(ikGlobal, iSpin, matrixElementDir)

          call readMatrixElement(nTransitions, order, fName, ME_tmp, volumeLine)
            ! Includes conversion from Hartree to J

          matrixElement(1,:,ikLocal) = ME_tmp

        else if(order == 1) then
          ! Read matrix elements for all modes
    
          do j = 1, nModes

            fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(j,4))//'/'//trim(getMatrixElementFNameWPath(ikGlobal,iSpin,matrixElementDir))

            call readMatrixElement(nTransitions, order, fName, ME_tmp, volumeLine)
              ! Includes conversion from Hartree to J
              !
              ! The volume line will get overwritten each time, but that's
              ! okay because the volume doesn't change between the files. 

            matrixElement(j,:,ikLocal) = ME_tmp

          enddo
  
        !endif
      endif
    enddo

    dE = dE*HartreeToJ

    if(order == 1) then
      matrixElement(:,:,:) = matrixElement(:,:,:)/(BohrToMeter*sqrt(elecMToKg))**2
        ! Shifter program outputs dq in Bohr*sqrt(elecM), and that
        ! dq is directly used by TME to get the matrix element, so
        ! we need to convert to m*sqrt(kg). In the matrix element,
        ! dq is in the denominator and squared.
    endif

  endif


  call MPI_BCAST(dE, size(dE), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)
  call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)


  if(ionode) write(*, '("  Beginning transition-rate calculation")')
   

  call getAndWriteTransitionRate(nTransitions, ibi, iSpin, mDim, order, nModes, dE, gamma0, & 
          matrixElement, volumeLine)

  
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

