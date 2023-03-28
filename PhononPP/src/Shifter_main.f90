program shifterMain

  use shifterMod
  
  implicit none


  call cpu_time(t0)

  call mpiInitialization()

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if (ionode) write(iostd,'("************ VASP Export complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program wfcExportVASPMain

