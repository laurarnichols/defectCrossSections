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


  call mpiInitialization()

  call initialize(gammaOnly, exportDir, VASPDir)
    !! * Set default values for input variables, open output file,
    !!   and start timers

  if ( ionode ) then
    
    read(5, inputParams, iostat=ios)
      !! * Read input variables
    
    if (ios /= 0) call exitError('export main', 'reading inputParams namelist', abs(ios))
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
      kPosition, nBands, nKPoints, nPWs1kGlobal, nRecords, nSpins, eigenE)
    !! * Read cell and wavefunction data from the WAVECAR file


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [ ] vasprun.xml  [ ] Set up grid  [ ] POTCAR (",f10.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call distributeKpointsInPools(nKPoints)
    !! * Figure out how many k-points there should be per pool


  call read_vasprun_xml(realLattVec, nKPoints, VASPDir, atomPositionsDir, eFermi, kWeight, fftGridSize, iType, nAtoms, nAtomsEachType, nAtomTypes)
    !! * Read the k-point weights and cell info from the `vasprun.xml` file


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [ ] Set up grid  [ ] POTCAR (",f10.2," secs)")') &
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
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [X] Set up grid  [ ] POTCAR (",f10.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  deallocate(gIndexLocalToGlobal)
  deallocate(gVecInCart)


  allocate(pot(nAtomTypes))

  if (ionode) write(iostd,*) "Reading POTCAR"

  call readPOTCAR(nAtomTypes, VASPDir, pot)
    !! * Read in pseudopotential information from POTCAR


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [X] Set up grid  [X] POTCAR (",f10.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call projAndWav(fftGridSize, maxGkVecsLocal, maxNumPWsGlobal, nAtoms, nAtomTypes, nBands, nGkVecsLocal, nGVecsGlobal, nKPoints, &
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
    write(*, '("[X] K-points  [ ] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues (",f10.2," secs)")') &
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
    write(*, '("[X] K-points  [X] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues (",f10.2," secs)")') &
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
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [ ] Pseudo  [ ] Eigenvalues (",f10.2," secs)")') &
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
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [X] Pseudo  [ ] Eigenvalues (",f10.2," secs)")') &
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
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [X] Pseudo  [X] Eigenvalues (",f10.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  deallocate(eigenE)
  deallocate(bandOccupation)

  if (ionode) close(mainOutFileUnit)

  call MPI_Barrier(worldComm, ierr)
 
  if (ionode) write(iostd,*) "************ VASP Export complete! ************"

  call MPI_FINALIZE(ierr)

end program wfcExportVASPMain

