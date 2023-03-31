program shifterMain

  use shifterMod
  use generalComputations, only: direct2cart
  
  implicit none

  integer :: j 
    !! Mode index

  call cpu_time(t0)

  call mpiInitialization()

  call initialize(nAtoms, shift, phononFName, poscarFName)
    !! * Set default values for input variables and start timers

  if(ionode) then
    
    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('export main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

  endif

  call checkInitialization(nAtoms, shift, phononFName, poscarFName)

  nModes = 3*nAtoms - 3
    !! * Calculate the total number of modes

  allocate(atomPositionsDir(3,nAtoms))
  allocate(eigenvector(3,nAtoms,nModes))

  call readPOSCAR(nAtoms, poscarFName, atomPositionsDir, omega, realLattVec)

  call readPhonons(nAtoms, nModes, phononFName, eigenvector)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)

  allocate(shiftedPositions(3,nAtoms))

  do j = iModeStart, iModeEnd

    shiftedPositions = direct2cart(nAtoms, atomPositionsDir, realLattVec) + getDisplacement(j, nAtoms, nModes, eigenvector, shift)
  !  call writeShiftedPOSCAR()

  enddo

  deallocate(atomPositionsDir)
  deallocate(eigenvector)
  deallocate(shiftedPositions)

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ Shifter complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program shifterMain

