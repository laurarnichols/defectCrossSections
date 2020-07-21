program wfcExportVASPMain
  !! Export scf output from VASP to form usable for TME
  !!
  !! input:  namelist "&inputParams", with variables
  !!   prefix       prefix of input files saved by program pwscf
  !!   outdir       temporary directory where files resides
  !!   exportDir    output directory. 
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

  call initialize()
    !! * Set default values for input variables, open output file,
    !!   and start timers

  if ( ionode_local ) then
    
    call input_from_file()
      !! @todo Remove this once extracted from QE #end @endtodo
    
    read(5, inputParams, iostat=ios)
      !! * Read input variables
    
    if (ios /= 0) call exitError('export main', 'reading inputParams namelist', abs(ios))
      !! * Exit calculation if there's an error

    outdir = QEDir
      !! * Convert the QEDir to outdir to match what QE uses
      !! @todo Remove this once extracted from QE #end @endtodo
    
    ios = f_mkdir_safe(trim(exportDir))
      !! * Make the export directory
    
    mainOutputFile = trim(exportDir)//"/input"
      !! * Set the name of the main output file for export
    
    write(stdout,*) "Opening file "//trim(mainOutputFile)
    open(mainout, file=trim(mainOutputFile))
      !! * Open main output file

  endif

  tmp_dir = trim(outdir)
  CALL mp_bcast( outdir, root, world_comm_local )
  CALL mp_bcast( tmp_dir, root, world_comm_local )
  CALL mp_bcast( prefix, root, world_comm_local )
    !! @todo Remove this once extracted from QE #end @endtodo

  CALL read_file
  CALL openfil_pp
    !! * Read QE input files
    !! @todo Remove this once extracted from QE #end @endtodo
  
  call readWAVECAR(VASPDir, at_local, bg_local, ecutwfc_local, occ, omega_local, vcut_local, &
      xk_local, nb1max, nb2max, nb3max, nbnd_local, nkstot_local, nplane, npmax, nspin_local)
    !! * Read cell and wavefunction data from the WAVECAR file

  stop

  call distributeKpointsInPools(nkstot_local)
    !! * Figure out how many k-points there should be per pool

  call calculateGvecs(nb1max, nb2max, nb3max, npmax, bg_local, gCart_local, ig_l2g, itmp_g, &
      ngm_g_local, ngm_local)
    !! * Calculate Miller indices and G-vectors and split
    !!   over processors

  call reconstructFFTGrid(ngm_local, ig_l2g, nkstot_local, gCart_local, vcut_local, xk_local, &
      igk_l2g, igk_large, ngk_local, ngk_g, npw_g, npwx_g, npwx_local, nplane, npmax)
    !! * Determine which G-vectors result in \(G+k\)
    !!   below the energy cutoff for each k-point and
    !!   sort the indices based on \(|G+k|^2\)

  deallocate(ig_l2g)
  deallocate(gCart_local)
  deallocate(nplane)

  call writeKInfo(nkstot_local, npwx_local, igk_l2g, nbnd_local, ngk_g, ngk_local, npw_g, npwx_g, &
      occ, xk_local, igwk)

  deallocate(occ)

  CALL write_export (mainOutputFile, exportDir)

  CALL stop_pp
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
 
END PROGRAM wfcExportVASPMain

