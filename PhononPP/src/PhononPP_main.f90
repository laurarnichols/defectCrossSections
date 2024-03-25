program PhononPPMain

  use PhononPPMod
  
  implicit none


  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call readInputs(dispInd, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, freqThresh, &
        shift, basePOSCARFName, CONTCARsBaseDir, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, & 
        calcMaxDisp, calcSj, generateShiftedPOSCARs, singleDisp)


  call readPhonons(freqThresh, phononFName, nAtoms, nModes, coordFromPhon, eigenvector, mass, omegaFreq)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)


  if(calcSj) &
    call calculateSj(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, nAtoms, nModes, &
          coordFromPhon, eigenvector, mass, omegaFreq, singleDisp, CONTCARsBaseDir, initPOSCARFName, finalPOSCARFName)


  deallocate(omegaFreq)


  if(calcDq .or. generateShiftedPOSCARs .or. calcMaxDisp) &
    call calculateShiftAndDq(dispInd, nAtoms, nModes, coordFromPhon, eigenvector, mass, shift, calcDq, calcMaxDisp, &
          generateShiftedPOSCARs, basePOSCARFName, dqFName, prefix)


  deallocate(coordFromPhon)
  deallocate(eigenvector)
  deallocate(mass)


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ PhononPP complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program PhononPPMain

