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

  call mpiInitialization()

  call initialize(exportDir, VASPDir)
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


  if (ionode) write(iostd,*) "Reading WAVECAR"

  call readWAVECAR(VASPDir, realSpaceLatticeVectors, recipSpaceLatticeVectors, wfcECut, bandOccupation, omega, wfcVecCut, &
      kPosition, nb1max, nb2max, nb3max, nBands, maxGkNum, nKPoints, nPWs1kGlobal, nSpins, &
      eigenE)
    !! * Read cell and wavefunction data from the WAVECAR file

  if (ionode) write(iostd,*) "Done reading WAVECAR"


  call distributeKpointsInPools(nKPoints)
    !! * Figure out how many k-points there should be per pool


  if (ionode) write(iostd,*) "Calculating G-vectors"

  call calculateGvecs(nb1max, nb2max, nb3max, recipSpaceLatticeVectors, gVecInCart, gIndexLocalToGlobal, gVecMillerIndicesGlobal, nGVecsGlobal, &
      nGVecsLocal)
    !! * Calculate Miller indices and G-vectors and split
    !!   over processors

  if (ionode) write(iostd,*) "Done calculating G-vectors"


  if (ionode) write(iostd,*) "Reconstructing FFT grid"

  call reconstructFFTGrid(nGVecsLocal, gIndexLocalToGlobal, maxGkNum, nKPoints, nPWs1kGlobal, recipSpaceLatticeVectors, gVecInCart, &
      wfcVecCut, kPosition, gKIndexLocalToGlobal, gToGkIndexMap, nGkLessECutLocal, nGkLessECutGlobal, maxGIndexGlobal, maxNumPWsGlobal, maxNumPWsPool)
    !! * Determine which G-vectors result in \(G+k\)
    !!   below the energy cutoff for each k-point and
    !!   sort the indices based on \(|G+k|^2\)

  if (ionode) write(iostd,*) "Done reconstructing FFT grid"


  deallocate(gIndexLocalToGlobal)
  deallocate(gVecInCart)
  deallocate(nPWs1kGlobal)

  if (ionode) write(iostd,*) "Reading vasprun.xml"

  call read_vasprun_xml(realSpaceLatticeVectors, nKPoints, VASPDir, eFermi, kWeight, iType, nAtoms, nAtomTypes)
    !! * Read the k-point weights and cell info from the `vasprun.xml` file

  if (ionode) write(iostd,*) "Done reading vasprun.xml"


  allocate(ps(nAtomTypes))


  if (ionode) write(iostd,*) "Reading POTCAR"

  call readPOTCAR(nAtomTypes, VASPDir, ps)
    !! * Read in pseudopotential information from POTCAR

  if (ionode) write(iostd,*) "Done reading POTCAR"


  if (ionode) write(iostd,*) "Writing k-point info"

  call writeKInfo(nKPoints, maxNumPWsPool, gKIndexLocalToGlobal, nBands, nGkLessECutGlobal, nGkLessECutLocal, maxGIndexGlobal, maxNumPWsGlobal, &
      bandOccupation, kWeight, kPosition, gKIndexGlobal)
    !! * Calculate ground state and global \(G+k\) indices
    !!   and write out k-point information to `input` file

  if (ionode) write(iostd,*) "Done writing k-point info"


  deallocate(kWeight)

  if (ionode) write(iostd,*) "Writing grid info"

  call writeGridInfo(nGVecsGlobal, nKPoints, maxNumPWsGlobal, gKIndexGlobal, gVecMillerIndicesGlobal, nGkLessECutGlobal, maxGIndexGlobal, exportDir)
    !! * Write out grid boundaries and miller indices
    !!   for just \(G+k\) combinations below cutoff energy
    !!   in one file and all miller indices in another 
    !!   file

  if (ionode) write(iostd,*) "Done writing grid info"
      

  deallocate(gVecMillerIndicesGlobal)

  if (ionode) write(iostd,*) "Writing cell info"

  call writeCellInfo(iType, nAtoms, nBands, nAtomTypes, nSpins, realSpaceLatticeVectors, recipSpaceLatticeVectors, atomPositions, nAtomsEachType)
    !! * Write out the real- and reciprocal-space lattice vectors, 
    !!   the number of atoms, the number of types of atoms, the
    !!   final atom positions, number of bands, and number of spins,
    !!   then calculate the number of atoms of each type

  if (ionode) write(iostd,*) "Done writing cell info"

  
  deallocate(iType)
  deallocate(atomPositions)

  if (ionode) write(iostd,*) "Writing pseudo info"

  call writePseudoInfo(nAtomTypes, nAtomsEachType, ps)
    !! * For each atom type, write out the element name,
    !!   number of atoms of this type, projector info,
    !!   radial grid info, and partial waves

  if (ionode) write(iostd,*) "Done writing pseudo info"


  if (ionode) write(iostd,*) "Writing eigenvalues"

  call writeEigenvalues(nBands, nKPoints, nSpins, eFermi, bandOccupation, eigenE)
    !! * Write Fermi energy and eigenvalues and occupations for each band

  if (ionode) write(iostd,*) "Done writing eigenvalues"


  deallocate(eigenE)
  deallocate(bandOccupation)

  if (ionode) close(mainOutFileUnit)
 
  if (ionode) write(iostd,*) "************ VASP Export complete! ************"

END PROGRAM wfcExportVASPMain

