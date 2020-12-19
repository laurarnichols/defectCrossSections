program wfcExportVASPMain
  !! Export scf output from VASP to form usable for TME
  !!
  !! input:  namelist "&inputParams", with variables
  !!   prefix      prefix of input files saved by program pwscf
  !!   QEDir       temporary directory where QE files reside
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

#ifdef __MPI
  call mpiInitialization()
#endif

  call initialize(exportDir, QEDir, VASPDir)
    !! * Set default values for input variables, open output file,
    !!   and start timers

  if ( ionode_local ) then
    
    read(5, inputParams, iostat=ios)
      !! * Read input variables
    
    if (ios /= 0) call exitError('export main', 'reading inputParams namelist', abs(ios))
      !! * Exit calculation if there's an error
    
    ios = f_mkdir_safe(trim(exportDir))
      !! * Make the export directory
    
    mainOutputFile = trim(exportDir)//"/input"
      !! * Set the name of the main output file for export
    
    write(stdout,*) "Opening file "//trim(mainOutputFile)
    open(mainout, file=trim(mainOutputFile))
      !! * Open main output file

  endif


  if (ionode_local) write(stdout,*) "Reading WAVECAR"

  call readWAVECAR(VASPDir, at_local, bg_local, ecutwfc_local, occ, omega_local, vcut_local, &
      xk_local, nb1max, nb2max, nb3max, nbnd_local, ngk_max, nkstot_local, nplane_g, nspin_local, &
      eigenE)
    !! * Read cell and wavefunction data from the WAVECAR file

  if (ionode_local) write(stdout,*) "Done reading WAVECAR"


  call distributeKpointsInPools(nkstot_local)
    !! * Figure out how many k-points there should be per pool


  if (ionode_local) write(stdout,*) "Calculating G-vectors"

  call calculateGvecs(nb1max, nb2max, nb3max, bg_local, gCart_local, ig_l2g, mill_g, ngm_g_local, &
      ngm_local)
    !! * Calculate Miller indices and G-vectors and split
    !!   over processors

  if (ionode_local) write(stdout,*) "Done calculating G-vectors"


  if (ionode_local) write(stdout,*) "Reconstructing FFT grid"

  call reconstructFFTGrid(ngm_local, ig_l2g, ngk_max, nkstot_local, nplane_g, bg_local, gCart_local, &
      vcut_local, xk_local, igk_l2g, igk_large, ngk_local, ngk_g, npw_g, npwx_g, npwx_local)
    !! * Determine which G-vectors result in \(G+k\)
    !!   below the energy cutoff for each k-point and
    !!   sort the indices based on \(|G+k|^2\)

  if (ionode_local) write(stdout,*) "Done reconstructing FFT grid"


  deallocate(ig_l2g)
  deallocate(gCart_local)
  deallocate(nplane_g)

  if (ionode_local) write(stdout,*) "Reading vasprun.xml"

  call read_vasprun_xml(at_local, nkstot_local, VASPDir, eFermi, wk_local, ityp, nat, nsp)
    !! * Read the k-point weights and cell info from the `vasprun.xml` file

  if (ionode_local) write(stdout,*) "Done reading vasprun.xml"


  allocate(ps(nsp))


  if (ionode_local) write(stdout,*) "Reading POTCAR"

  call readPOTCAR(nsp, VASPDir, ps)
    !! * Read in pseudopotential information from POTCAR

  if (ionode_local) write(stdout,*) "Done reading POTCAR"


  if (ionode_local) write(stdout,*) "Writing k-point info"

  call writeKInfo(nkstot_local, npwx_local, igk_l2g, nbnd_local, ngk_g, ngk_local, npw_g, npwx_g, &
      occ, wk_local, xk_local, igwk)
    !! * Calculate ground state and global \(G+k\) indices
    !!   and write out k-point information to `input` file

  if (ionode_local) write(stdout,*) "Done writing k-point info"


  deallocate(wk_local)

  if (ionode_local) write(stdout,*) "Writing grid info"

  call writeGridInfo(ngm_g_local, nkstot_local, npwx_g, igwk, mill_g, ngk_g, npw_g, exportDir)
    !! * Write out grid boundaries and miller indices
    !!   for just \(G+k\) combinations below cutoff energy
    !!   in one file and all miller indices in another 
    !!   file

  if (ionode_local) write(stdout,*) "Done writing grid info"
      

  deallocate(mill_g)

  if (ionode_local) write(stdout,*) "Writing cell info"

  call writeCellInfo(ityp, nat, nbnd_local, nsp, nspin_local, at_local, bg_local, tau, nnTyp)
    !! * Write out the real- and reciprocal-space lattice vectors, 
    !!   the number of atoms, the number of types of atoms, the
    !!   final atom positions, number of bands, and number of spins,
    !!   then calculate the number of atoms of each type

  if (ionode_local) write(stdout,*) "Done writing cell info"

  
  deallocate(ityp)
  deallocate(tau)

  if (ionode_local) write(stdout,*) "Writing pseudo info"

  call writePseudoInfo(nsp, nnTyp, ps)
    !! * For each atom type, write out the element name,
    !!   number of atoms of this type, projector info,
    !!   radial grid info, and partial waves

  if (ionode_local) write(stdout,*) "Done writing pseudo info"


  if (ionode_local) write(stdout,*) "Writing eigenvalues"

  call writeEigenvalues(nbnd_local, nkstot_local, nspin_local, eFermi, occ, eigenE)
    !! * Write Fermi energy and eigenvalues and occupations for each band

  if (ionode_local) write(stdout,*) "Done writing eigenvalues"


  deallocate(eigenE)
  deallocate(occ)

  if (ionode_local) close(mainout)
 
  if (ionode_local) write(stdout,*) "************ VASP Export complete! ************"

END PROGRAM wfcExportVASPMain

