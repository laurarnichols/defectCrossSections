program PhononPPMain

  use PhononPPMod
  
  implicit none

  integer :: j
    !! Loop indices

  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, freqThresh, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, &
        phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs, singleDisp)

  ! Read one initial POSCAR to get the number of atoms
  if(ionode) then

    if(.not. singleDisp) initPOSCARFName = trim(CONTCARsBaseDir)//'/'//trim(int2str(iBandIinit))//'/CONTCAR'

    call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
    call standardizeCoordinates(nAtoms, atomPositionsDirInit)

  endif


  ! Broadcast atom number and cell parameters
  call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
  if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
  call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  ! Get the number of modes, read the phonons, and distribute
  ! across the proceses
  nModes = 3*nAtoms - 3

  allocate(eigenvector(3,nAtoms,nModes))
  allocate(mass(nAtoms))
  allocate(omegaFreq(nModes))

  call readPhonons(nAtoms, nModes, atomPositionsDirInit, phononFName, eigenvector, mass, omegaFreq)

  deallocate(atomPositionsDirInit)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)


  call calculateSj(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nAtoms, nModes, eigenvector, mass, omegaFreq, singleDisp, & 
        CONTCARsBaseDir, initPOSCARFName, finalPOSCARFName)


  deallocate(omegaFreq)


  call calculateShiftAndDq(nAtoms, nModes, eigenvector, mass, shift, generateShiftedPOSCARs, basePOSCARFName, dqFName, prefix)


  deallocate(mass)
  deallocate(eigenvector)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ PhononPP complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program PhononPPMain

