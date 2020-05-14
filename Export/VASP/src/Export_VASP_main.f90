!-----------------------------------------------------------------------
program wfcExportVASPMain
  !-----------------------------------------------------------------------
  !
  ! writes PWSCF data for postprocessing purposes in XML format using IOTK lib
  ! Wave-functions are collected and written using IO_BASE module.
  !
  ! input:  namelist "&inputpp", with variables
  !   prefix       prefix of input files saved by program pwscf
  !   outdir       temporary directory where files resides
  !   exportDir    output directory. A directory
  !                "exportDir" is created and
  !                output files are put there. All the data
  !                are accessible through the ""exportDir"/input" file.
  !


  USE wrappers,      ONLY : f_mkdir_safe
  USE pwcom

  USE io_global, ONLY : ionode, ionode_id
  USE io_files,  ONLY : prefix, tmp_dir, outdir
  USE ions_base, ONLY : ntype => nsp
  USE iotk_module
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : kunit
  USE mp_world,  ONLY: world_comm
  USE mp,        ONLY: mp_bcast
  USE environment,   ONLY : environment_start
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ik, i, kunittmp, ios
  !
  real(kind = dp), parameter :: ryToHartree = 0.5_dp
  !
  CHARACTER(len=256) :: pp_file, exportDir
  LOGICAL :: writeWFC
  !
  NAMELIST /inputpp/ prefix, outdir, exportDir, writeWFC
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW_EXPORT' )
  !
  !   set default values for variables in namelist
  !
  prefix = ''
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  exportDir = './Export'
  !
  writeWFC  = .true.        ! gdb : by default the wavefunctions are needed, 
                            !       this gives the user the ability not to write the wavefunctions
  !
  !    Reading input file
  !
  IF ( ionode ) THEN
    !
    CALL input_from_file ( )
    !
    READ(5, inputpp, IOSTAT=ios)
    !
    IF (ios /= 0) CALL errore ('pw_export', 'reading inputpp namelist', abs(ios) )
    !
    ios = f_mkdir_safe( trim(exportDir) )
    !
    pp_file = trim(exportDir)//"/input"
    !
    !
  ENDIF
  !
  ! ... Broadcasting variables
  !
  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( outdir, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  !

#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  CALL write_export (pp_file, exportDir, kunittmp )

  CALL stop_pp
 
END PROGRAM pw_export_for_TME

