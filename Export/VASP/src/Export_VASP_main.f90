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
      !! @todo Remove this once extracted from QE @endtodo
    
    read(5, inputParams, iostat=ios)
      !! * Read input variables
    
    if (ios /= 0) call exitError ('export main', 'reading inputParams namelist', abs(ios) )
      !! * Exit calculation if there's an error

    outdir = QEDir
      !! * Convert the QEDir to outdir to match what QE uses;
      !!   eventually will get rid of both directories, but 
      !!   needed for now
    
    ios = f_mkdir_safe( trim(exportDir) )
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

  CALL read_file
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
  CALL openfil_pp
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
  
  call readWAVECAR(VASPDir, at_local, bg_local, ecutwfc_local, omega_local, vcut_local, xk_local, &
      ikEnd, ikStart, nb1max, nb2max, nb3max, nbnd_local, nk_Pool, nkstot_local, nplane, npmax, &
      nspin_local)
    !! * Read data from the WAVECAR file

  call calculateGvecs(ikEnd, ikStart, nb1max, nb2max, nb3max, nk_Pool, nkstot_local, nplane, npmax, &
      bg_local, vcut_local, xk_local, igk_l2g, igk_large, itmp_g, ngk_local, ngk_g, ngm_local, &
      ngm_g_local, npw_g, npwx_g, npwx_local)

  deallocate(nplane)

  CALL write_export (mainOutputFile, exportDir)

  CALL stop_pp
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
 
END PROGRAM wfcExportVASPMain

