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
  integer :: nAtoms_
    !! Input number of atoms from CONTCAR file

  real(kind=dp) :: omegaPOS, realLattVecPOS(3,3)
    !! Variables from POSCAR to compare to WAVECAR

  character(len=300) :: fName
    !! File name for CONTCAR


  call cpu_time(t0)

  call mpiInitialization('Export')

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools

  call readInputParams(nDispkPerCoord, patternArr, energiesOnly, gammaOnly, groupForGroupVelocity, exportDir, pattern, VASPDir)
    !! * Initialize, read, check, and broadcast input parameters

  if(ionode) &
    write(*, '("[ ] WAVECAR  [ ] vasprun.xml  [ ] Set up grid  [ ] POTCAR")')
  call cpu_time(t1)


  call readWAVECAR(VASPDir, realLattVec, recipLattVec, bandOccupation, omega, wfcVecCut, &
      kPosition, fftGridSize, nBands, nKPoints, nPWs1kGlobal, nSpins, reclenWav, eigenE)
    !! * Read cell and wavefunction data from the WAVECAR file


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [ ] vasprun.xml  [ ] Set up grid  [ ] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call read_vasprun_xml(nKPoints, VASPDir, eFermi, eTot, kWeight, iType, nAtoms, nAtomsEachType, nAtomTypes)
    !! * Read the k-point weights and cell info from the `vasprun.xml` file

  if(ionode) then
    fName = trim(VASPDir)//'/CONTCAR'
    write(*,*) trim(fName)

    call readPOSCAR(fName, nAtoms_, atomPositionsDir, omegaPOS, realLattVecPOS)
      !! * Get coordinates from CONTCAR

    if(nAtoms_ /= nAtoms) call exitError('export main', 'nAtoms from vasprun.xml and '//trim(fName)//' don''t match!', 1)

    omegaPOS = omegaPOS*angToBohr**3
    realLattVecPOS = realLattVecPOS*angToBohr

    if(abs(omegaPOS - omega) > 1e-5 .or. maxval(abs(realLattVecPOS-realLattVec)) > 1e-5) &
      call exitError('export main', 'omega and lattice vectors from POSCAR and WAVECAR not consistent', 1)

  endif

  if(.not. ionode) allocate(atomPositionsDir(3,nAtoms))

  call MPI_BCAST(atomPositionsDir, size(atomPositionsDir), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools

  call distributeItemsInSubgroups(myBgrpId, nBands, nProcPerPool, nProcPerBgrp, nBandGroups, ibStart_bgrp, ibEnd_bgrp, nbPerBgrp)
    !! * Distribute bands in groups

  call distributeItemsInSubgroups(indexInBgrp, nAtoms, nProcPerBgrp, nProcPerBgrp, nProcPerBgrp, iaStart_bgrp, iaEnd_bgrp, naPerProcInBgrp)
    !! * Distribute atoms across processes in band group


  allocate(nGkLessECutGlobal(nKPoints))
  nGkLessECutGlobal(:) = 0
    ! This is needed because this array is passed to 
    ! `writeKInfo` even though it isn't used in the
    ! energies-only version. Setting this here is just
    ! for safety.

  if(.not. energiesOnly) then

    call cpu_time(t2)
    if(ionode) &
      write(*, '("[X] WAVECAR  [X] vasprun.xml  [ ] Set up grid  [ ] POTCAR (",f7.2," secs)")') &
            t2-t1
    call cpu_time(t1)

    call calculateGvecs(fftGridSize, recipLattVec, gVecInCart, gIndexLocalToGlobal, gVecMillerIndicesGlobal, iMill, nGVecsGlobal, nGVecsLocal)
      !! * Calculate Miller indices and G-vectors and split
      !!   over processors


    call reconstructFFTGrid(nGVecsLocal, gIndexLocalToGlobal, nKPoints, nPWs1kGlobal, kPosition, gVecInCart, recipLattVec, &
        wfcVecCut, gKIndexOrigOrderLocal, gKSort, maxGIndexGlobal, maxGkVecsLocal, maxNumPWsGlobal, maxNumPWsPool, &
        nGkLessECutGlobal, nGkVecsLocal)
      !! * Determine which G-vectors result in \(G+k\)
      !!   below the energy cutoff for each k-point and
      !!   sort the indices based on \(|G+k|^2\)


    deallocate(gIndexLocalToGlobal)
    deallocate(gVecInCart)
  endif


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [X] Set up grid  [ ] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  allocate(pot(nAtomTypes))

  call readPOTCAR(nAtomTypes, VASPDir, pot)
    !! * Read in pseudopotential information from POTCAR


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] WAVECAR  [X] vasprun.xml  [X] Set up grid  [X] POTCAR (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  if(.not. energiesOnly) then

    call projAndWav(maxGkVecsLocal, maxNumPWsGlobal, nAtoms, nAtomTypes, nBands, nGkVecsLocal, nGVecsGlobal, nKPoints, &
        nSpins, gKIndexOrigOrderLocal, gKSort, gVecMillerIndicesGlobal, nPWs1kGlobal, reclenWav, atomPositionsDir, kPosition, omega, &
        recipLattVec, exportDir, VASPDir, gammaOnly, pot)


    deallocate(gKIndexOrigOrderLocal, gKSort, nGkVecsLocal, iGkStart_pool, iGkEnd_pool)
  endif

  deallocate(nPWs1kGlobal)


  if(ionode) &
    write(*, '("[ ] K-points  [ ] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues")')
  call cpu_time(t1)


  call writeKInfo(nBands, nKPoints, nGkLessECutGlobal, nSpins, bandOccupation, kWeight, kPosition)
    !! * Calculate ground state and write out k-point 
    !!   information to `input` file

  deallocate(kPosition)
  deallocate(kWeight)
  deallocate(nGkLessECutGlobal)


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [ ] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call writeGridInfo(nGVecsGlobal, gVecMillerIndicesGlobal, maxGIndexGlobal, exportDir)
    !! * Write out grid boundaries and miller indices
    !!   for just \(G+k\) combinations below cutoff energy
    !!   in one file and all miller indices in another 
    !!   file

  if(.not. energiesOnly) deallocate(gVecMillerIndicesGlobal)


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [ ] Cell  [ ] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)
      

  call writeCellInfo(iType, nAtoms, nBands, nAtomTypes, realLattVec, recipLattVec, atomPositionsDir)
    !! * Write out the real- and reciprocal-space lattice vectors, 
    !!   the number of atoms, the number of types of atoms, the
    !!   final atom positions, number of bands, and number of spins,
    !!   then calculate the number of atoms of each type
  
  deallocate(iType)
  deallocate(atomPositionsDir)


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [ ] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call writePseudoInfo(nAtomTypes, nAtomsEachType, pot)
    !! * For each atom type, write out the element name,
    !!   number of atoms of this type, projector info,
    !!   radial grid info, and partial waves

  deallocate(nAtomsEachType)

  if(ionode) then
    do iT = 1, nAtomTypes
      deallocate(pot(iT)%radGrid, pot(iT)%dRadGrid, pot(iT)%wae, pot(iT)%wps)
    enddo
  endif

  deallocate(pot)


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [X] Pseudo  [ ] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  call writeEigenvalues(nBands, nKPoints, nSpins, bandOccupation, eFermi, eTot, eigenE)
    !! * Write Fermi energy and eigenvalues and occupations for each band

  if(groupForGroupVelocity) then
    call writeGroupedEigenvalues(nBands, nDispkPerCoord, nKPoints, nSpins, bandOccupation, patternArr, eigenE, pattern)
    deallocate(patternArr)
  endif

  deallocate(eigenE)
  deallocate(bandOccupation)


  call cpu_time(t2)
  if(ionode) &
    write(*, '("[X] K-points  [X] Grid  [X] Cell  [X] Pseudo  [X] Eigenvalues (",f7.2," secs)")') &
          t2-t1
  call cpu_time(t1)


  if(ionode) then
   close(mainOutFileUnit)
  endif

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if (ionode) write(*,'("************ VASP Export complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program wfcExportVASPMain

