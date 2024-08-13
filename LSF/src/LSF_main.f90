program LSFmain

  use LSFmod

  implicit none

  ! Local variables:
  integer, allocatable :: iDum1D(:)
    !! Integer to ignore input
  integer :: iDum1, iDum2
    !! Integers to ignore input
  integer :: j, ikLocal, ikGlobal, jStore
    !! Loop index
  integer :: mDim
    !! Size of first dimension for matrix element

  real(kind=dp), allocatable :: dENew(:)
    !! New energy to update matrix element; needed
    !! not to create a temporary array when passing
    !! a slice of dE
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
    !! File name to read


  call cpu_time(timerStart)

  call mpiInitialization('LSF')

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools


  if(ionode) write(*, '("Pre-k-loop: [ ] Get parameters  [ ] Read Sj")')
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


  call cpu_time(timer2)
  if(ionode) write(*, '("Pre-k-loop: [X] Get parameters  [ ] Read Sj (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)


  call readSj(diffOmega, PhononPPDir, nModes, omega, omegaPrime, Sj, SjPrime)

  allocate(nj(nModes))
  call readNj(nModes, njInput, nj)


  call cpu_time(timer2)
  if(ionode) write(*, '("Pre-k-loop: [X] Get parameters  [X] Read Sj (",f10.2," secs)")') timer2-timer1
  call cpu_time(timer1)


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools


  allocate(jReSort(nModes))

  if(ionode) then
    call readCaptureEnergyTable(1, iSpin, energyTableDir, ibi, ibf, nTransitions, rDum2D)
     ! Assume that band bounds and number of transitions do not depend on k-points or spin
     ! We ignore the energy here because we are only reading ibi, ibf, and nTransitions

    deallocate(rDum2D)

    if(order == 1 .and. reSortMEs) call getjReSort(nModes, PhononPPDir, jReSort)
  endif

  call MPI_BCAST(nTransitions, 1, MPI_INTEGER, root, worldComm, ierr)
  if(.not. ionode) allocate(ibi(nTransitions))
  call MPI_BCAST(ibi, nTransitions, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(ibf, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(jReSort, nModes, MPI_INTEGER, root, worldComm, ierr)


  allocate(dE(3,nTransitions,nkPerPool))
  allocate(dENew(nTransitions))

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

          dENew = dE(2,:,ikLocal)
          ! Call to readMatrixElement includes conversion from Hartree to J
          if(newEnergyTable) then
            call readMatrixElement(minval(ibi), maxval(ibi), nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
          else
            call readMatrixElement(-1, -1, nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
              ! dENew will be ignored here
          endif

          matrixElement(1,:,ikLocal) = ME_tmp

        else if(order == 1) then
          ! Read matrix elements for all modes
    
          do j = 1, nModes

            fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(j,suffixLength))//'/'&
                    //trim(getMatrixElementFNameWPath(ikGlobal,iSpin,matrixElementDir))

            dENew = dE(3,:,ikLocal)
            ! Call to readMatrixElement includes conversion from Hartree to J
            ! The volume line will get overwritten each time, but that's
            ! okay because the volume doesn't change between the files. 


            ! If resorting the matrix element files based on a different PhononPP output
            ! order, make sure to pass the resorted mode index if re-reading dq's
            if(reSortMEs) then
              jStore = jReSort(j)
            else 
              jStore = j
            endif


            if(newEnergyTable) then
              if(rereadDq) then
                call readMatrixElement(minval(ibi), maxval(ibi), nTransitions, order, dENew, newEnergyTable, oldFormat, &
                      fName, ME_tmp, volumeLine, jStore, PhononPPDir)
              else
                call readMatrixElement(minval(ibi), maxval(ibi), nTransitions, order, dENew, newEnergyTable, oldFormat, &
                      fName, ME_tmp, volumeLine)
              endif
            else
              if(rereadDq) then
                call readMatrixElement(-1, -1, nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine, &
                  jStore, PhononPPDir)
                  ! dENew will be ignored here
              else
                call readMatrixElement(-1, -1, nTransitions, order, dENew, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
                  ! dENew will be ignored here
              endif
            endif

            matrixElement(jStore,:,ikLocal) = ME_tmp

          enddo
  
        !endif
      endif
    enddo

  endif

  deallocate(jReSort)


  call MPI_BCAST(dE, size(dE), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)
  call MPI_BCAST(matrixElement, size(matrixElement), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)


  if(ionode) write(*, '("  Beginning transition-rate calculation")')
   

  call getAndWriteTransitionRate(nTransitions, ibi, iSpin, mDim, order, nModes, dE, gamma0, & 
          matrixElement, SjThresh, volumeLine)

  
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

