program PhononPPMain

  use PhononPPMod
  use generalComputations, only: direct2cart
  use cell, only: readPOSCAR, cartDispProjOnPhononEigsNorm, writePOSCARNewPos
  use miscUtilities, only: int2strLeadZero
  
  implicit none

  integer :: j, ibi, ibf
    !! Loop indices

  call cpu_time(t0)

  call mpiInitialization('PhononPP')

  call readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, &
        phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs, singleDisp)

  ! Read one initial POSCAR to get the number of atoms
  if(ionode) then

    if(.not. singleDisp) initPOSCARFName = trim(CONTCARsBaseDir)//'/'//trim(int2str(iBandIinit))//'/CONTCAR'

    call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
    call standardizeCoordinates(nAtoms, atomPositionsDirInit)

  endif


  ! Broadcast atom number and cell parameters
  call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
  call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

  if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
  call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


  ! Get the number of modes, read the phonons, and distribute
  ! across the proceses
  nModes = 3*nAtoms - 3

  allocate(eigenvector(3,nAtoms,nModes))
  allocate(mass(nAtoms))
  allocate(omegaFreq(nModes))

  call readPhonons(nAtoms, nModes, phononFName, eigenvector, mass, omegaFreq)

  call distributeItemsInSubgroups(myid, nModes, nProcs, nProcs, nProcs, iModeStart, iModeEnd, nModesLocal)


  if(singleDisp) then

    SjFName = 'Sj.out'
    call getSingleDisp(nAtoms, nModes, atomPositionsDirInit, eigenvector, mass, omega, omegaFreq, finalPOSCARFName, SjFName)

    deallocate(atomPositionsDirInit)

  else

    ! I do a loop over all of the bands. Depending on the band selection,
    ! this could result in duplicate pairs (e.g., .1.2 and .2.1), but I 
    ! can't think of a general way to exlude these pairs without a significant
    ! amount of work to track the pairs that have already been calcualted.
    ! Any solution I can think of would not hold for both hole and electron
    ! capture and/or different ranges of bands selected by the user. For
    ! now, I just have the code output all of the duplicate pairs and the
    ! user/LSF code can use whatever they need from that.
    do ibi = iBandIinit, iBandIfinal

      ! Get the initial positions for this band (unless we did
      ! it earlier)
      if(ibi /= iBandIinit) then
        if(ionode) then
          initPOSCARFName = trim(CONTCARsBaseDir)//'/'//trim(int2str(ibi))//'/CONTCAR'

          call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
          call standardizeCoordinates(nAtoms, atomPositionsDirInit)
        endif

        if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
        call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      endif

      do ibf = iBandFinit, iBandFfinal

        if(ibf /= ibi) then
          ! Set the final positions and output file names
          finalPOSCARFName = trim(CONTCARsBaseDir)//'/'//trim(int2str(ibf))//'/CONTCAR'
          SjFName = 'Sj.'//trim(int2str(ibi))//'.'//trim(int2str(ibf))//'.out'

          ! Get displacement and output ranked Sj for a single pair of states
          call getSingleDisp(nAtoms, nModes, atomPositionsDirInit, eigenvector, mass, omega, omegaFreq, finalPOSCARFName, SjFName)

        endif
      enddo

      ! Positions are allocated in readPOSCAR, so we need to 
      ! deallocate them after every loop
      deallocate(atomPositionsDirInit)
    enddo

  endif


  deallocate(omegaFreq)




  if(nModes < 10) then
    suffixLength = 1
  else if(nModes < 100) then
    suffixLength = 2
  else if(nModes < 1000) then
    suffixLength = 3
  else if(nModes < 10000) then
    suffixLength = 4
  else if(nModes < 100000) then
    suffixLength = 5
  endif




  if(generateShiftedPOSCARs) then
    ! Read base POSCAR to shift from
    if(ionode) then

    call readPOSCAR(basePOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
    call standardizeCoordinates(nAtoms, atomPositionsDirInit)

    endif

    if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
    call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

  endif


  ! Get the displacement for each mode to 
  ! calculate the derivative of the wave function.
  ! Project onto the phonon eigenvectors (here,
  ! the effect is just to convert back to generalized
  ! coordinates because the displacement is already
  ! a scaled form of the eigenvectors), then (if needed)
  ! write the shifted positions and generalized displacement
  ! norms.
  allocate(shiftedPositions(3,nAtoms))
  allocate(projNorm(nModes))
  allocate(displacement(3,nAtoms))
  projNorm = 0.0_dp

  write(memoLine,'("  shift = ", ES9.2E1)') shift

  do j = iModeStart, iModeEnd

    displacement = getShiftDisplacement(nAtoms, eigenvector(:,:,j), mass, shift)

    projNorm(j) = cartDispProjOnPhononEigsNorm(j, nAtoms, displacement, eigenvector(:,:,j), mass)

    if(generateShiftedPOSCARs) then

      shiftedPositions = direct2cart(nAtoms, atomPositionsDirInit, realLattVec) + displacement

      shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

      call writePOSCARNewPos(nAtoms, shiftedPositions, basePOSCARFName, shiftedPOSCARFName, trim(memoLine)//' j = '//int2str(j), .true.)

    endif

  enddo

  deallocate(mass)
  deallocate(atomPositionsDirInit)
  deallocate(eigenvector)
  deallocate(shiftedPositions)
  deallocate(displacement)

  projNorm = projNorm*angToBohr*sqrt(daltonToElecM)
  call writeDqs(nModes, projNorm, dqFName)

  deallocate(projNorm)

  call MPI_Barrier(worldComm, ierr)
 
  call cpu_time(t2)
  if(ionode) write(*,'("************ PhononPP complete! (",f10.2," secs) ************")') t2-t0

  call MPI_FINALIZE(ierr)

end program PhononPPMain

