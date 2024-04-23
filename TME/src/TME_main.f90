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


  call readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, order, phononModeJ, baselineDir, &
          braExportDir, ketExportDir, dqFName, energyTableDir, outputDir, loopSpins, subtractBaseline)
    !! Read input, initialize, check that required variables were set, and
    !! distribute across processes
    

  nSys = 2

  call setUpSystemArray(nSys, braExportDir, ketExportDir, crystalSystem)


  call completePreliminarySetup(nSys, order, phononModeJ, dqFName, mill_local, nGVecsGlobal, nKPoints, nSpins, &
          dq_j, gCart, omega, recipLattVec, Ylm, crystalSystem, pot)


  call calcAndWrite2SysMatrixElements(ispSelect, nSpins, crystalSystem(1), crystalSystem(2), pot)


  call MPI_BARRIER(worldComm, ierr)
  if(ionode) write(*,'("Done with k loop!")')

  
  call finalizeCalculation(nSys, crystalSystem, pot)

  deallocate(crystalSystem)
  
  call MPI_FINALIZE(ierr)
  
end program TMEmain
