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


  call completePreliminarySetup(nSys, order, phononModeJ, dqFName, maxAngMom, mill_local, nGVecsGlobal, nKPoints, nSpins, &
        gCart, dq_j, omega, Ylm, crystalSystem)


  call calcAndWrite2SysMatrixElements(crystalSystem(1), crystalSystem(2))


  call MPI_BARRIER(worldComm, ierr)
  if(ionode) write(*,'("Done with k loop!")')

  
  call finalizeCalculation(nSys, crystalSystem)

  deallocate(crystalSystem)
  
  call MPI_FINALIZE(ierr)
  
end program TMEmain
