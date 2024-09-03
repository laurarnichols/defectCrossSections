program PhononPPMain

  use PhononPPMod
  
  implicit none

  real(kind=dp), allocatable :: rDum1d(:), rDum2d(:,:)
    !! Dummy variables to ignore input

  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call readInputs(disp2AtomInd, freqThresh, shift, SjThresh, temperature, allStatesBaseDir_relaxed, basePOSCARFName, dqFName, &
        energyTableDir, phononFName, phononPrimeFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, calcMaxDisp, &
        calcSj, diffOmega, dqEigvecsFinal, generateShiftedPOSCARs, singleDisp)


  call readPhonons(freqThresh, phononFName, nAtoms, nModes, coordFromPhon, eigenvector, mass, omega)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)

  call calcAndWriteNj(nModes, omega, temperature)


  if(diffOmega) then
    call readPhonons(freqThresh, phononPrimeFName, nAtomsPrime, nModesPrime, coordFromPhonPrime, eigenvectorPrime, rDum1d, omegaPrime)

    if(ionode) then
      if(nAtoms /= nAtomsPrime) call exitError('PhononPPMain', 'Number of atoms does not match in different phonon files.', 1)
      if(nModes /= nModesPrime) call exitError('PhononPPMain', 'Number of modes does not match in different phonon files.', 1)
    endif


    allocate(rDum2d(3,nAtoms))

    call getRelaxDispAndCheckCompatibility(nAtoms, coordFromPhon, coordFromPhonPrime, rDum2d)

    deallocate(rDum2d)

    call lineUpModes(nAtoms, nModes, eigenvector, eigenvectorPrime, omegaPrime)

  else
    allocate(omegaPrime(nModes))
      ! Need to allocate this variable either way so that 
      ! the inputs to calculateSj match expectations
  endif


  ! If we need to calculate something with the eigenvectors, determine
  ! which eigenvectors should be used.
  allocate(dqEigenvectors(3,nAtoms,nModes))
  if(calcSj .or. calcDq .or. generateShiftedPOSCARs .or. calcMaxDisp) then
    if(diffOmega .and. dqEigvecsFinal) then
      dqEigenvectors = eigenvectorPrime
    else
      dqEigenvectors = eigenvector
    endif
  endif


  deallocate(eigenvector)
  if(diffOmega) deallocate(eigenvectorPrime)


  if(calcSj) &
    call calculateSj(nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, omega, omegaPrime, SjThresh, diffOmega, singleDisp, &
            allStatesBaseDir_relaxed, energyTableDir, initPOSCARFName, finalPOSCARFName)


  deallocate(omega)
  deallocate(omegaPrime)


  if(calcDq .or. generateShiftedPOSCARs .or. calcMaxDisp) &
    call calculateShiftAndDq(disp2AtomInd, nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, shift, calcDq, calcMaxDisp, &
          generateShiftedPOSCARs, basePOSCARFName, dqFName, prefix)


  deallocate(coordFromPhon)
  deallocate(dqEigenvectors)
  deallocate(mass)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ PhononPP complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program PhononPPMain

