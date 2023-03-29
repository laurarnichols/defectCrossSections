program wfcExportVASPMain
  !! Export scf output from VASP to form usable for TME
  !!
  !! input:  namelist "&inputParams", with variables
  !!   VASPDir     temporary directory where VASP files reside
  !!   exportDir   output directory. 
  !! A directory "exportDir" is created and
  !! output files are put there. All the data
  !! are accessible through the ""exportDir"/input" file.
  !!
  !! <h2>Walkthrough</h2>
  !! 

  use wfcExportVASPMod
  
  implicit none

  integer :: iT
    !! Index for deallocating `pot` variables

  real(kind=dp) :: omegaPOS, realLattVecPOS(3,3)
    !! Variables from POSCAR to compare to WAVECAR

  character(len=300) :: fName
    !! File name for CONTCAR


  call cpu_time(t0)

  call mpiInitialization()

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools

  call initialize(gammaOnly, exportDir, VASPDir)
    !! * Set default values for input variables, open output file,
    !!   and start timers

  if ( ionode ) then
    
    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('export main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

    call execute_command_line('mkdir -p '//trim(exportDir))
      !! * Make the export directory
    
    mainOutputFile = trim(exportDir)//"/input"
      !! * Set the name of the main output file for export
    
    write(iostd,*) "Opening file "//trim(mainOutputFile)
    open(mainOutFileUnit, file=trim(mainOutputFile))
      !! * Open main output file

  endif

  call MPI_BCAST(exportDir, len(exportDir), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(VASPDir, len(VASPDir), MPI_CHARACTER, root, worldComm, ierr)
  call MPI_BCAST(gammaOnly, 1, MPI_LOGICAL, root, worldComm, ierr)

  if(ionode) &
    write(*, '("[ ] WAVECAR  [ ] vasprun.xml  [ ] Set up grid  [ ] POTCAR")')
  call cpu_time(t1)


  call readWAVECAR(VASPDir, realLattVec, recipLattVec, bandOccupation, omega, wfcVecCut, &
      kPosition, fftGridSize, nBands, nKPoints, nPWs1kGlobal, nRecords, nSpins, eigenE)
    !! * Read cell and wavefunction data from the WAVECAR file


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [ ] vasprun.xml  [ ] Set up grid  [ ] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call read_vasprun_xml(realLattVec, nKPoints, VASPDir, eFermi, kWeight, iType, nAtoms, nAtomsEachType, nAtomTypes)
    !! * Read the k-point weights and cell info from the `vasprun.xml` file

  allocate(atomPositionsDir(3,nAtoms))

  if(ionode) then
    fName = trim(VASPDir)//'CONTCAR'

    call readPOSCAR(nAtoms, fName, atomPositionsDir, omegaPOS, realLattVecPOS)
      !! * Get coordinates from CONTCAR

    omegaPOS = omegaPOS*angToBohr**3
    realLattVecPOS = realLattVecPOS*angToBohr

    if(abs(omegaPOS - omega) > 1e-5 .or. maxval(abs(realLattVecPOS-realLattVec)) > 1e-5) &
      call exitError('export main', 'omega and lattice vectors from POSCAR and WAVECAR not consistent', 1)

  endif

  call MPI_BCAST(atomPositionsDir, size(atomPositionsDir), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools

  call distributeItemsInSubgroups(myBgrpId, nBands, nProcPerPool, nProcPerBgrp, nBandGroups, ibStart_bgrp, ibEnd_bgrp, nbPerBgrp)
    !! * Distribute bands in groups

  call distributeItemsInSubgroups(indexInBgrp, nAtoms, nProcPerBgrp, nProcPerBgrp, nProcPerBgrp, iaStart_bgrp, iaEnd_bgrp, naPerProcInBgrp)
    !! * Distribute atoms across processes in band group


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [ ] Set up grid  [ ] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call calculateGvecs(fftGridSize, recipLattVec, gVecInCart, gIndexLocalToGlobal, gVecMillerIndicesGlobal, iMill, nGVecsGlobal, nGVecsLocal)
    !! * Calculate Miller indices and G-vectors and split
    !!   over processors


  call reconstructFFTGrid(nGVecsLocal, gIndexLocalToGlobal, nKPoints, nPWs1kGlobal, kPosition, gVecInCart, recipLattVec, wfcVecCut, gKIndexGlobal, &
      gKIndexLocalToGlobal, gKIndexOrigOrderLocal, gKSort, gToGkIndexMap, maxGIndexGlobal, maxGkVecsLocal, maxNumPWsGlobal, maxNumPWsPool, &
      nGkLessECutGlobal, nGkLessECutLocal, nGkVecsLocal)
    !! * Determine which G-vectors result in \(G+k\)
    !!   below the energy cutoff for each k-point and
    !!   sort the indices based on \(|G+k|^2\)


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [X] Set up grid  [ ] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  deallocate(gIndexLocalToGlobal)
  deallocate(gVecInCart)


  allocate(pot(nAtomTypes))

  call readPOTCAR(nAtomTypes, VASPDir, pot)
    !! * Read in pseudopotential information from POTCAR


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [X] Set up grid  [X] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call projAndWav(maxGkVecsLocal, maxNumPWsGlobal, nAtoms, nAtomTypes, nBands, nGkVecsLocal, nGVecsGlobal, nKPoints, &
      nRecords, nSpins, gKIndexOrigOrderLocal, gKSort, gVecMillerIndicesGlobal, nPWs1kGlobal, atomPositionsDir, kPosition, omega, &
      recipLattVec, exportDir, VASPDir, gammaOnly, pot)


  deallocate(nPWs1kGlobal, gKIndexOrigOrderLocal, gKSort, nGkVecsLocal, iGkStart_pool, iGkEnd_pool)


  if(ionode) &
    write(*, '("[ ] K-points  [ ] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues")')
  call cpu_time(t1)


  call writeKInfo(nBands, nKPoints, nGkLessECutGlobal, nSpins, bandOccupation, kWeight, kPosition)
    !! * Calculate ground state and write out k-point 
    !!   information to `input` file


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [ ] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  deallocate(kPosition)
  deallocate(kWeight)


  call writeGridInfo(nGVecsGlobal, nKPoints, nSpins, maxNumPWsGlobal, gKIndexGlobal, gVecMillerIndicesGlobal, nGkLessECutGlobal, maxGIndexGlobal, exportDir)
    !! * Write out grid boundaries and miller indices
    !!   for just \(G+k\) combinations below cutoff energy
    !!   in one file and all miller indices in another 
    !!   file


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)
      

  deallocate(gKIndexGlobal, gVecMillerIndicesGlobal, nGkLessECutGlobal)


  call writeCellInfo(iType, nAtoms, nBands, nAtomTypes, realLattVec, recipLattVec, atomPositionsDir)
    !! * Write out the real- and reciprocal-space lattice vectors, 
    !!   the number of atoms, the number of types of atoms, the
    !!   final atom positions, number of bands, and number of spins,
    !!   then calculate the number of atoms of each type


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [ ] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)

  
  deallocate(iType)
  deallocate(atomPositionsDir)


  call writePseudoInfo(nAtomTypes, nAtomsEachType, pot)
    !! * For each atom type, write out the element name,
    !!   number of atoms of this type, projector info,
    !!   radial grid info, and partial waves


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [X] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  deallocate(nAtomsEachType)

  if(ionode) then
    do iT = 1, nAtomTypes
      deallocate(pot(iT)%radGrid, pot(iT)%dRadGrid, pot(iT)%wae, pot(iT)%wps)
    enddo
  endif


  call writeEigenvalues(nBands, nKPoints, nSpins, eFermi, bandOccupation, eigenE)
    !! * Write Fermi energy and eigenvalues and occupations for each band


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [X] Pseudo  [X] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  deallocate(eigenE)
  deallocate(bandOccupation)

  if (ionode) close(mainOutFileUnit)

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if (ionode) write(*,'("************ VASP Export complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program wfcExportVASPMain

