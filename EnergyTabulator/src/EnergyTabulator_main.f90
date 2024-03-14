program energyTabulatorMain

  use energyTabulatorMod
  
  implicit none

  integer :: isp, ikLocal
    !! Loop indices

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


  ! Get total energies from exports of all different structures
  call getTotalEnergy(exportDirInitInit, eTotInitInit)
  call getTotalEnergy(exportDirFinalInit, eTotFinalInit)
  call getTotalEnergy(exportDirFinalFinal, eTotFinalFinal)


  do isp = 1, nSpins
    if(loopSpins .or. isp == ispSelect) then

      ! Get reference eigenvalue
      call getRefEig(isp, refBand, exportDirEigs, refEig)

      if(captured) call getRefToDefectEigDiff(iBandFinit, isp, refBand, exportDirInitInit, elecCarrier, dEEigRefDefect)

      do ikLocal = 1, nkPerProc

        call writeEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikLocal, isp, dEEigRefDefect, eCorrectTot, eCorrectEigRef, &
              eTotInitInit, eTotFinalInit, eTotFinalFinal, refEig, energyTableDir, exportDirEigs, captured, elecCarrier)

      enddo
    endif
  enddo


  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ Energy tabulator complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program energyTabulatorMain
