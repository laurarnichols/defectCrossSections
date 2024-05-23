program TMEmain
  use TMEmod
  
  implicit none

  
  call mpiInitialization('TME')
    !! Initialize MPI

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools
  
  if(ionode) call cpu_time(t0)


  call readInputParams(ibBra, ibKet, ispSelect, nPairs, order, phononModeJ, baselineDir, braExportDir, &
          ketExportDir, dqFName, energyTableDir, outputDir, loopSpins, overlapOnly, subtractBaseline)
    !! Read input, initialize, check that required variables were set, and
    !! distribute across processes
    

  nSys = 2

  call setUpSystemArray(nSys, braExportDir, ketExportDir, crystalSystem)


  call completePreliminarySetup(ispSelect, nSys, order, phononModeJ, loopSpins, overlapOnly, dqFName, ibi, &
        ibf, mill_local, nGVecsGlobal, nGVecsLocal, nKPoints, nSpins, nTransitions, dq_j, omega, recipLattVec, Ylm, &
        crystalSystem, pot)

  
  if(overlapOnly) then
    call getAndWriteOnlyOverlaps(nPairs, ibBra, ibKet, ispSelect, nGVecsLocal, nSpins, crystalSystem(1), crystalSystem(2), pot)
  else
    call getAndWriteCaptureMatrixElements(nTransitions, ibi, ibf, ispSelect, nGVecsLocal, nSpins, crystalSystem(1), crystalSystem(2), pot)
  endif 


  call MPI_BARRIER(worldComm, ierr)

  
  call finalizeCalculation(nSys, crystalSystem, pot)

  deallocate(crystalSystem)
  
  call MPI_FINALIZE(ierr)
  
end program TMEmain
