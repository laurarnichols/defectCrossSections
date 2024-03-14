program energyTabulatorMain

  use energyTabulatorMod
  
  implicit none

  real(kind=dp) :: t0, t2
    !! Timers


  call cpu_time(t0)

  call mpiInitialization('EnergyTabulator')

  ! Get inputs from user
  call readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef, energyTableDir, &
        exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier)


  call getnSpinsAndnKPoints(exportDirEigs, nKPoints, nSpins)
    ! Assumes that all systems have the same number of spins and k-points

  ! Distribute k-points in pools
  call distributeItemsInSubgroups(myid, nKPoints, nProcs, nProcs, nProcs, ikStart, ikEnd, nkPerProc)


  if(captured) then
    call calcAndWriteCaptureEnergies(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, nSpins, refBand, eCorrectTot, &
         eCorrectEigRef, elecCarrier, loopSpins, energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal)
  else
    call calcAndWriteScatterEnergies()
  endif


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ Energy tabulator complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program energyTabulatorMain
