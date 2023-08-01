program PhononPPMain

  use PhononPPMod
  use generalComputations, only: direct2cart
  use cell, only: readPOSCAR, cartDisplacementToGeneralizedNorm, writePOSCARNewPos
  use miscUtilities, only: int2strLeadZero
  
  implicit none

  integer :: j 
    !! Mode index

  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call initialize(shift, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs)
    !! * Set default values for input variables and start timers

  if(ionode) then
    
    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('PhononPP main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

    call checkInitialization(shift, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs)

    call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
    call readPOSCAR(finalPOSCARFName, nAtomsFinal, atomPositionsDirFinal, omegaFinal, realLattVec)

    call standardizeCoordinates(nAtoms, atomPositionsDirInit)
    call standardizeCoordinates(nAtomsFinal, atomPositionsDirFinal)
  
    call checkCompatibility(nAtoms, nAtomsFinal, omega, omegaFinal, atomPositionsDirFinal, atomPositionsDirInit)

  endif

  call MPI_BCAST(shift, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(dqFName, len(dqFName), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(phononFName, len(phononFName), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(generateShiftedPOSCARs, 1, MPI_LOGICAL, root, worldComm, ierr)

  call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(atomPositionsDirFinal, size(atomPositionsDirFinal), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  nModes = 3*nAtoms - 3
    !! * Calculate the total number of modes

  allocate(eigenvector(3,nAtoms,nModes))
  allocate(mass(nAtoms))
  allocate(omegaFreq(nModes))

  call readPhonons(nAtoms, nModes, phononFName, eigenvector, mass, omegaFreq)

  deallocate(omegaFreq)

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

    displacement = getShiftDisplacement(j, nAtoms, nModes, eigenvector, mass, shift)

    generalizedNorm(j) = cartDisplacementToGeneralizedNorm(nAtoms, displacement, mass)*angToBohr*sqrt(daltonToElecM)
      !! Convert scaled displacement back to generalized
      !! coordinates and get norm
      !! @note
      !!   Input positions are in angstrom and input
      !!   masses are in amu, but the dq output is going
      !!   to our code, which uses Hartree atomic units
      !!   (Bohr and electron masses), so this value
      !!   must have a unit conversion.
      !! @endnote

    if(generateShiftedPOSCARs) then

      shiftedPositions = direct2cart(nAtoms, atomPositionsDirInit, realLattVec) + displacement

      shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

      call writePOSCARNewPos(nAtoms, shiftedPositions, initPOSCARFName, shiftedPOSCARFName, .true.)

    endif

  enddo

  deallocate(atomPositionsDirInit)
  deallocate(atomPositionsDirFinal)
  deallocate(eigenvector)
  deallocate(shiftedPositions)
  deallocate(displacement)

  call writeDqs(nModes, generalizedNorm, dqFName)

  deallocate(generalizedNorm)

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ PhononPP complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program PhononPPMain

