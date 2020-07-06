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
  use wvfct, only : ecutwfc, npwx
  
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
  
  call readWAVECAR(VASPDir, nspin_local, ecutwfc_local, vcut_local, at_local, &
        nkstot_local, nbnd_local, omega_local, bg_local, xk_local, ngm_g_local, &
        ngm_local, mill_local, gCart_local, nk_Pool)
    !! * Read data from the WAVECAR file

  stop

  call reconstructMainGrid(nk_Pool, nkstot_local, xk_local, igk_l2g, itmp_g)

  CALL write_export (mainOutputFile, exportDir)

  CALL stop_pp
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
 
END PROGRAM wfcExportVASPMain

