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

  if ( ionode ) then
    
    read(5, inputParams, iostat=ios)
      !! * Read input variables
    
    if (ios /= 0) call exitError ('export main', 'reading inputParams namelist', abs(ios) )
      !! * Exit calculation if there's an error
    
    ios = f_mkdir_safe( trim(exportDir) )
      !! * Make the export directory
    
    mainOutputFile = trim(exportDir)//"/input"
      !! * Set the name of the main output file for export
    
    call readWAVECAR()
      !! * Read data from the WAVECAR file

    write(stdout,*) "Opening file "//trim(mainOutputFile)
    open(mainout, file=trim(mainOutputFile))
      !! * Open main output file

  endif

  !CALL mp_bcast( outdir, root, world_comm )
  !CALL mp_bcast( exportDir, root, world_comm )
  !CALL mp_bcast( prefix, root, world_comm )

  CALL read_file
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
  CALL openfil_pp
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
  

  call distributeKpointsInPools()

  call reconstructMainGrid()

  CALL write_export (mainOutputFile, exportDir)

  CALL stop_pp
    !! @todo Figure out what this subroutine does and what can be moved here @endtodo
 
END PROGRAM wfcExportVASPMain

