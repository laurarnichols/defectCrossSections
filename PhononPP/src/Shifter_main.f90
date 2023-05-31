program shifterMain

  use shifterMod
  use generalComputations, only: direct2cart
  use miscUtilities, only: int2strLeadZero
  
  implicit none

  integer :: j 
    !! Mode index

  call cpu_time(t0)

  call mpiInitialization()

  call initialize(shift, dqFName, phononFName, poscarFName, prefix)
    !! * Set default values for input variables and start timers

  if(ionode) then
    
    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('export main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

    call checkInitialization(shift, dqFName, phononFName, poscarFName, prefix)

  endif

  call MPI_BCAST(shift, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(dqFName, len(dqFName), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(phononFName, len(phononFName), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(poscarFName, len(poscarFName), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)

  call readPOSCAR(poscarFName, nAtoms, atomPositionsDir, omega, realLattVec)

  nModes = 3*nAtoms - 3
    !! * Calculate the total number of modes

  allocate(eigenvector(3,nAtoms,nModes))
  allocate(mass(nAtoms))

  call readPhonons(nAtoms, nModes, phononFName, eigenvector, mass)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)

  if(nModes < 10) then
    suffixLength = 1
  else if(nModes < 100) then
    suffixLength = 2
  else if(nModes < 1000) then
    suffixLength = 3
  else if(nModes < 10000) then
    suffixLength = 4
  else if(nModes < 100000) then
    suffixLength = 5
  endif

  allocate(shiftedPositions(3,nAtoms))
  allocate(displacement(3,nAtoms))
  allocate(generalizedNorm(nModes))

  generalizedNorm = 0.0_dp

  do j = iModeStart, iModeEnd

    call getDisplacement(j, nAtoms, nModes, eigenvector, mass, shift, displacement, generalizedNorm(j))

    shiftedPositions = direct2cart(nAtoms, atomPositionsDir, realLattVec) + displacement
    !shiftedPositions = getDisplacement(j, nAtoms, nModes, eigenvector, mass, shift)

    shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

    call writePOSCARNewPos(nAtoms, shiftedPositions, poscarFName, shiftedPOSCARFName, .true.)

  enddo

  deallocate(atomPositionsDir)
  deallocate(eigenvector)
  deallocate(shiftedPositions)
  deallocate(displacement)

  call writeDqs(nModes, generalizedNorm, dqFName)

  deallocate(generalizedNorm)

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ Shifter complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program shifterMain

