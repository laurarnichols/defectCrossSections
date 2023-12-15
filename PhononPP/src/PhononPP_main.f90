program PhononPPMain

  use PhononPPMod
  use generalComputations, only: direct2cart
  use cell, only: readPOSCAR, cartDispProjOnPhononEigsNorm, writePOSCARNewPos
  use miscUtilities, only: int2strLeadZero
  
  implicit none

  integer :: j
    !! Mode index

  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call readInputs(shift, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs)

  if(ionode) then
    call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
    call readPOSCAR(finalPOSCARFName, nAtomsFinal, atomPositionsDirFinal, omegaFinal, realLattVec)

    call standardizeCoordinates(nAtoms, atomPositionsDirInit)
    call standardizeCoordinates(nAtomsFinal, atomPositionsDirFinal)

    if(nAtoms /= nAtomsFinal) &
      call exitError('PhononPP main', 'number of atoms does not match: '//trim(int2str(nAtoms))//' '//trim(int2str(nAtomsFinal)), 1)

    if(abs(omega - omegaFinal) > 1e-8) call exitError('PhononPP main', 'volumes don''t match', 1)

  endif


  call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

  if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms), atomPositionsDirFinal(3,nAtoms))
  call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(atomPositionsDirFinal, size(atomPositionsDirFinal), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  nModes = 3*nAtoms - 3
    !! * Calculate the total number of modes

  allocate(eigenvector(3,nAtoms,nModes))
  allocate(mass(nAtoms))
  allocate(omegaFreq(nModes))

  call readPhonons(nAtoms, nModes, phononFName, eigenvector, mass, omegaFreq)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)


  allocate(projNorm(nModes))
  allocate(displacement(3,nAtoms))
  

  !> Define the displacement for the relaxation, 
  !> project onto the phonon eigenvectors, and
  !> get Sj
  call getRelaxDispAndCheckCompatibility(nAtoms, atomPositionsDirFinal, atomPositionsDirInit, displacement)
  
  projNorm = 0.0_dp
  do j = iModeStart, iModeEnd

    projNorm(j) = cartDispProjOnPhononEigsNorm(j, nAtoms, displacement, eigenvector(:,:,j), mass)

  enddo

  projNorm = projNorm*angToM*sqrt(daltonToElecM*elecMToKg)
  call calcAndWriteSj(nModes, omegaFreq, projNorm)

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

  allocate(shiftedPositions(3,nAtoms))

  !> Get the displacement for each mode to 
  !> calculate the derivative of the wave function.
  !> Project onto the phonon eigenvectors (here,
  !> the effect is just to convert back to generalized
  !> coordinates because the displacement is already
  !> a scaled form of the eigenvectors), then (if needed)
  !> write the shifted positions and generalized displacement
  !> norms.
  projNorm = 0.0_dp

  write(memoLine,'("  shift = ", ES9.2E1)') shift

  do j = iModeStart, iModeEnd

    displacement = getShiftDisplacement(nAtoms, eigenvector(:,:,j), mass, shift)

    projNorm(j) = cartDispProjOnPhononEigsNorm(j, nAtoms, displacement, eigenvector(:,:,j), mass)

    if(generateShiftedPOSCARs) then

      shiftedPositions = direct2cart(nAtoms, atomPositionsDirInit, realLattVec) + displacement

      shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

      call writePOSCARNewPos(nAtoms, shiftedPositions, initPOSCARFName, shiftedPOSCARFName, trim(memoLine)//' j = '//int2str(j), .true.)

    endif

  enddo

  deallocate(atomPositionsDirInit)
  deallocate(atomPositionsDirFinal)
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

