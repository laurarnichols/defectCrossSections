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

  if (ionode) write(iostd,*) "Reading WAVECAR"

  call readWAVECAR(VASPDir, realLattVec, recipLattVec, bandOccupation, omega, wfcVecCut, &
      kPosition, nBands, nKPoints, nPWs1kGlobal, nRecords, nSpins, eigenE)
    !! * Read cell and wavefunction data from the WAVECAR file

  if (ionode) write(iostd,*) "Done reading WAVECAR"


  call distributeKpointsInPools(nKPoints)
    !! * Figure out how many k-points there should be per pool


  if (ionode) write(iostd,*) "Reading vasprun.xml"

  call read_vasprun_xml(realLattVec, nKPoints, VASPDir, atomPositionsDir, eFermi, kWeight, fftGridSize, iType, nAtoms, nAtomsEachType, nAtomTypes)
    !! * Read the k-point weights and cell info from the `vasprun.xml` file

  if (ionode) write(iostd,*) "Done reading vasprun.xml"


  if (ionode) write(iostd,*) "Calculating G-vectors"

  call calculateGvecs(fftGridSize, recipLattVec, gVecInCart, gIndexLocalToGlobal, gVecMillerIndicesGlobal, iMill, nGVecsGlobal, nGVecsLocal)
    !! * Calculate Miller indices and G-vectors and split
    !!   over processors

  if (ionode) write(iostd,*) "Done calculating G-vectors"


  if (ionode) write(iostd,*) "Reconstructing FFT grid"

  call reconstructFFTGrid(nGVecsLocal, gIndexLocalToGlobal, nKPoints, nPWs1kGlobal, kPosition, gVecInCart, recipLattVec, wfcVecCut, gKIndexGlobal, &
      gKIndexLocalToGlobal, gKIndexOrigOrderLocal, gKSort, gToGkIndexMap, maxGIndexGlobal, maxGkVecsLocal, maxNumPWsGlobal, maxNumPWsPool, &
      nGkLessECutGlobal, nGkLessECutLocal, nGkVecsLocal)
    !! * Determine which G-vectors result in \(G+k\)
    !!   below the energy cutoff for each k-point and
    !!   sort the indices based on \(|G+k|^2\)

  if (ionode) write(iostd,*) "Done reconstructing FFT grid"


  deallocate(gIndexLocalToGlobal)
  deallocate(gVecInCart)


  allocate(pot(nAtomTypes))

  if (ionode) write(iostd,*) "Reading POTCAR"

  call readPOTCAR(nAtomTypes, VASPDir, pot)
    !! * Read in pseudopotential information from POTCAR

  if (ionode) write(iostd,*) "Done reading POTCAR"


  if (ionode) write(iostd,*) "Getting and writing projectors, projections, and wfc"

  call projAndWav(fftGridSize, maxGkVecsLocal, maxNumPWsGlobal, nAtoms, nAtomTypes, nBands, nGkVecsLocal, nGVecsGlobal, nKPoints, &
      nRecords, nSpins, gKIndexOrigOrderLocal, gKSort, gVecMillerIndicesGlobal, nPWs1kGlobal, atomPositionsDir, kPosition, omega, &
      recipLattVec, exportDir, VASPDir, gammaOnly, pot)

  if (ionode) write(iostd,*) "Done getting and writing projectors, projections, and wfc"

  deallocate(nPWs1kGlobal, gKIndexOrigOrderLocal, gKSort, nGkVecsLocal, iGkStart_pool, iGkEnd_pool)


  if (ionode) write(iostd,*) "Writing k-point info"

  call writeKInfo(nBands, nKPoints, nGkLessECutGlobal, nSpins, bandOccupation, kWeight, kPosition)
    !! * Calculate ground state and write out k-point 
    !!   information to `input` file

  if (ionode) write(iostd,*) "Done writing k-point info"


  deallocate(kPosition)
  deallocate(kWeight)

  if (ionode) write(iostd,*) "Writing grid info"

  call writeGridInfo(nGVecsGlobal, nKPoints, nSpins, maxNumPWsGlobal, gKIndexGlobal, gVecMillerIndicesGlobal, nGkLessECutGlobal, maxGIndexGlobal, exportDir)
    !! * Write out grid boundaries and miller indices
    !!   for just \(G+k\) combinations below cutoff energy
    !!   in one file and all miller indices in another 
    !!   file

  if (ionode) write(iostd,*) "Done writing grid info"
      

  deallocate(gKIndexGlobal, gVecMillerIndicesGlobal, nGkLessECutGlobal)

  if (ionode) write(iostd,*) "Writing cell info"

  call writeCellInfo(iType, nAtoms, nBands, nAtomTypes, realLattVec, recipLattVec, atomPositionsDir)
    !! * Write out the real- and reciprocal-space lattice vectors, 
    !!   the number of atoms, the number of types of atoms, the
    !!   final atom positions, number of bands, and number of spins,
    !!   then calculate the number of atoms of each type

  if (ionode) write(iostd,*) "Done writing cell info"

  
  deallocate(iType)
  deallocate(atomPositionsDir)

  if (ionode) write(iostd,*) "Writing pseudo info"

  call writePseudoInfo(nAtomTypes, nAtomsEachType, pot)
    !! * For each atom type, write out the element name,
    !!   number of atoms of this type, projector info,
    !!   radial grid info, and partial waves

  if (ionode) write(iostd,*) "Done writing pseudo info"

  deallocate(nAtomsEachType)

  if(ionode) then
    do iT = 1, nAtomTypes
      deallocate(pot(iT)%radGrid, pot(iT)%dRadGrid, pot(iT)%wae, pot(iT)%wps)
    enddo
  endif


  if (ionode) write(iostd,*) "Writing eigenvalues"

  call writeEigenvalues(nBands, nKPoints, nSpins, eFermi, bandOccupation, eigenE)
    !! * Write Fermi energy and eigenvalues and occupations for each band

  if (ionode) write(iostd,*) "Done writing eigenvalues"


  deallocate(eigenE)
  deallocate(bandOccupation)

  if (ionode) close(mainOutFileUnit)

  call MPI_Barrier(worldComm, ierr)
 
  if (ionode) write(iostd,*) "************ VASP Export complete! ************"

  call MPI_FINALIZE(ierr)

END PROGRAM wfcExportVASPMain

