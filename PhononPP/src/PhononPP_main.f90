program PhononPPMain

  use PhononPPMod
  use generalComputations, only: direct2cart
  use cell, only: readPOSCAR, cartDispProjOnPhononEigsNorm, writePOSCARNewPos
  use miscUtilities, only: int2strLeadZero
  
  implicit none

  integer :: j
    !! Loop indices

  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, &
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




  if(ionode) then

  call readPOSCAR(basePOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
  call standardizeCoordinates(nAtoms, atomPositionsDirInit)

  endif

  call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
  if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
  call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  ! Get the displacement for each mode to 
  ! calculate the derivative of the wave function.
  ! Project onto the phonon eigenvectors (here,
  ! the effect is just to convert back to generalized
  ! coordinates because the displacement is already
  ! a scaled form of the eigenvectors), then (if needed)
  ! write the shifted positions and generalized displacement
  ! norms.
  allocate(shiftedPositions(3,nAtoms))
  allocate(projNorm(nModes))
  allocate(displacement(3,nAtoms))
  projNorm = 0.0_dp

  write(memoLine,'("  shift = ", ES9.2E1)') shift

  do j = iModeStart, iModeEnd

    displacement = getShiftDisplacement(nAtoms, eigenvector(:,:,j), realLattVec, mass, shift)

    projNorm(j) = cartDispProjOnPhononEigsNorm(nAtoms, displacement, eigenvector(:,:,j), mass, realLattVec)

    if(generateShiftedPOSCARs) then

      shiftedPositions = direct2cart(nAtoms, atomPositionsDirInit, realLattVec) + displacement

      shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

      call writePOSCARNewPos(nAtoms, shiftedPositions, basePOSCARFName, shiftedPOSCARFName, trim(memoLine)//' j = '//int2str(j), .true.)

    endif

  enddo

  deallocate(mass)
  deallocate(atomPositionsDirInit)
  deallocate(eigenvector)
  deallocate(shiftedPositions)
  deallocate(displacement)

  projNorm = projNorm*angToBohr*sqrt(daltonToElecM)
  call writeDqs(nModes, projNorm, dqFName)

  deallocate(projNorm)

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ PhononPP complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program PhononPPMain

