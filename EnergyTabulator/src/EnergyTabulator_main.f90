program energyTabulatorMain

  use energyTabulatorMod
  
  implicit none

  integer :: isp, ikLocal
    !! Loop indices

  real(kind=dp) :: t0, t2
    !! Timers


  call cpu_time(t0)

  call mpiInitialization('EnergyTabulator')

  call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef, energyTableDir, &
        exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier)
    !! * Set default values for input variables and start timers

  if(ionode) then
    
    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('energy tabulator main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

    call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef,&
          energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier, loopSpins)

  endif

  call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(refBand, 1, MPI_INTEGER, root, worldComm, ierr)

  call MPI_BCAST(ispSelect, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(loopSpins, 1, MPI_LOGICAL, root, worldComm, ierr)

  call MPI_BCAST(eCorrectTot, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(eCorrectEigRef, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

  call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(exportDirEigs, len(exportDirEigs), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(exportDirInitInit, len(exportDirInitInit), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(exportDirFinalInit, len(exportDirFinalInit), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(exportDirFinalFinal, len(exportDirFinalFinal), MPI_CHARACTER, root, worldComm, ierr)

  call MPI_BCAST(captured, 1, MPI_LOGICAL, root, worldComm, ierr)
  call MPI_BCAST(elecCarrier, 1, MPI_LOGICAL, root, worldComm, ierr)


  call getnSpinsAndnKPoints(exportDirEigs, nKPoints, nSpins)
    !! Assumes that all systems have the same number of spins and k-points

  call distributeItemsInSubgroups(myid, nKPoints, nProcs, nProcs, nProcs, ikStart, ikEnd, nkPerProc)
    !! * Distribute k-points in pools


  call getTotalEnergies(exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, eTotInitInit, eTotFinalInit, eTotFinalFinal)
    !! Get total energies from exports of all different structures


  do isp = 1, nSpins
    if(loopSpins .or. isp == ispSelect) then

      call getRefEig(isp, refBand, exportDirEigs, refEig)
        !! Get reference eigenvalue

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
