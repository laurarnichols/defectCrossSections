program shifterMain

  use shifterMod
  
  implicit none

  integer :: j 
    !! Mode index

  call cpu_time(t0)

  call mpiInitialization()

  call initialize(shift, phononFName, poscarFName)
    !! * Set default values for input variables and start timers

  if(ionode) then
    
    read(5, inputParams, iostat=ierr)
      !! * Read input variables
    
    if(ierr /= 0) call exitError('export main', 'reading inputParams namelist', abs(ierr))
      !! * Exit calculation if there's an error

  endif

  call checkInitialization(shift, phononFName, poscarFName)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeStart, nModesLocal)
    !! @todo Get modes from somewhere @endtodo

  do j = iModeStart, iModeEnd

    call shiftMode()
    call writeShiftedPOSCAR()

  enddo

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if (ionode) write(iostd,'("************ Shifter complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program wfcExportVASPMain

