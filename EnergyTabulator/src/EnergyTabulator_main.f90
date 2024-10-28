program energyTabulatorMain

  use energyTabulatorMod
  
  implicit none

  real(kind=dp) :: t0, t2
    !! Timers


  call cpu_time(t0)

  call mpiInitialization('EnergyTabulator')

  ! Get inputs from user
  call readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, ikFinit, ikFfinal, &
        ispSelect, refBand, dENegThresh, dEZeroThresh, eCorrectTot, eCorrectEigRef, allStatesBaseDir_relaxed, &
        allStatesBaseDir_startPos, energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, &
        exportDirFinalFinal, singleStateExportDir, captured, elecCarrier, loopSpins)


  call getnSpinsAndnKPoints(exportDirEigs, nKPoints, nSpins)
    ! Assumes that all systems have the same number of spins and k-points


  if(captured) then

    ! Distribute k-points in pools
    call distributeItemsInSubgroups(myid, nKPoints, nProcs, nProcs, nProcs, ikStart, ikEnd, nkPerProc)

    call calcAndWriteCaptureEnergies(iBandIinit, iBandIfinal, iBandFinit, ispSelect, nSpins, refBand, eCorrectTot, &
        eCorrectEigRef, elecCarrier, loopSpins, energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal)
  else
    if(ionode) call calcAndWriteScatterEnergies(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, &
        ikFinit, ikFfinal, ispSelect, nSpins, refBand, dENegThresh, dEZeroThresh, eCorrectEigRef, elecCarrier, loopSpins, &
        allStatesBaseDir_relaxed, allStatesBaseDir_startPos, energyTableDir, exportDirEigs, singleStateExportDir)
  endif


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ Energy tabulator complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program energyTabulatorMain
