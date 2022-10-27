program transitionMatrixElements
  !
  use mpi
  use declarations
  !
  implicit none
  !
  real(kind = dp) :: t1, t2
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  !
  allocate ( nPWsI(0:numprocs-1), nPWsF(0:numprocs-1) )
  !
  if ( myid == root ) then
    !
    call cpu_time(t0)
    !
    ! Reading input, initializing and checking all variables of the calculation.
    !
    call readInput()
    !
    call readPWsSet()
    !
    allocate ( Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nKptsPC) )
    allocate ( paw_SDKKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_PsiPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_SDPhi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_fi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( eigvI (iBandIinit:iBandIfinal), eigvF (iBandFinit:iBandFfinal) )
    !
    !
    Ufi(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    allocate ( counts(0:numprocs-1), displmnt(0:numprocs-1) )
    !
    call distributePWsToProcs(numOfGvecs, numprocs)
    !
    nPWsI(:) = 0
    nPWsF(:) = 0
    !
    do i = 0, numprocs - 1
      nPWsI(i) = 1 + sum(counts(:i-1))
      nPWsF(i) = sum(counts(:i))
    enddo
    !  
  endif
  !
  call MPI_BCAST(iBandIinit,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iBandFinit,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(nKptsPC,     1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(nProjsPC,    1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nProjsSD,    1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nBands,      1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nSpins,      1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(numOfPWs,    1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(numOfGvecs,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(nPWsI, numprocs, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nPWsF, numprocs, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate ( gvecs(3, numOfGvecs) )
  call MPI_BCAST(gvecs, size(gvecs), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(numOfTypesPC, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate ( atomsPC(numOfTypesPC) )
  !
  call MPI_BCAST(JMAX, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  do i = 1, numOfTypesPC
    !
    call MPI_BCAST(atomsPC(i)%numOfAtoms, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atomsPC(i)%lMax,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atomsPC(i)%lmMax,      1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atomsPC(i)%nMax,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atomsPC(i)%iRc,        1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    !
    if ( myid /= root ) then 
      allocate( atomsPC(i)%lps(atomsPC(i)%lMax) )
      allocate( atomsPC(i)%r  (atomsPC(i)%nMax) )
      allocate( atomsPC(i)%rab(atomsPC(i)%nMax) )
      allocate( atomsPC(i)%F(atomsPC(i)%iRc, atomsPC(i)%lMax ) )
      allocate( atomsPC(i)%F1(atomsPC(i)%iRc, atomsPC(i)%lMax, atomsPC(i)%lMax ) )
      allocate( atomsPC(i)%bes_J_qr ( 0:JMAX, atomsPC(i)%iRc) )
    endif
    !
    call MPI_BCAST(atomsPC(i)%lps, size(atomsPC(i)%lps), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atomsPC(i)%r,   size(atomsPC(i)%r),   MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atomsPC(i)%rab, size(atomsPC(i)%rab), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atomsPC(i)%F,   size(atomsPC(i)%F),   MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atomsPC(i)%F1,  size(atomsPC(i)%F1),  MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atomsPC(i)%bes_J_qr, size(atomsPC(i)%bes_J_qr), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  enddo
  !
  call MPI_BCAST(nIonsPC, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate( posIonPC(3, nIonsPC), TYPNIPC(nIonsPC) )
  call MPI_BCAST(TYPNIPC,  size(TYPNIPC),  MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(posIonPC, size(posIonPC), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(numOfTypes, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate ( atoms(numOfTypes) )
  !
  do i = 1, numOfTypes
    !
    call MPI_BCAST(atoms(i)%numOfAtoms, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atoms(i)%lMax,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atoms(i)%lmMax,      1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atoms(i)%nMax,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(atoms(i)%iRc,        1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    !
    if ( myid /= root ) then 
      allocate( atoms(i)%lps(atoms(i)%lMax) )
      allocate( atoms(i)%r(atoms(i)%nMax) )
      allocate( atoms(i)%rab(atoms(i)%nMax) )
      allocate( atoms(i)%F(atoms(i)%iRc, atoms(i)%lMax ) )
      allocate( atoms(i)%F1(atoms(i)%iRc, atoms(i)%lMax, atoms(i)%lMax ) )
      allocate( atoms(i)%bes_J_qr( 0:JMAX, atoms(i)%iRc) )
    endif
    !
    call MPI_BCAST(atoms(i)%lps, size(atoms(i)%lps), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atoms(i)%r, size(atoms(i)%r), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atoms(i)%rab, size(atoms(i)%rab), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atoms(i)%F, size(atoms(i)%F), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atoms(i)%F1, size(atoms(i)%F1), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(atoms(i)%bes_J_qr, size(atoms(i)%bes_J_qr), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  enddo
  !
  call MPI_BCAST(nIonsSD, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate( posIonSD(3, nIonsSD), TYPNISD(nIonsSD) )
  call MPI_BCAST(TYPNISD,  size(TYPNISD),  MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(posIonSD, size(posIonSD), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  allocate ( paw_id(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
  !
  do ik = 1, nKptsPC
    !
    if ( myid == root ) then
      !
      tmes_file_exists = .false.
      call checkIfCalculated(ik,tmes_file_exists)
      !
    endif
    !
    call MPI_BCAST(tmes_file_exists, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    !
    if ( .not.tmes_file_exists ) then
      !
      allocate ( cProjPC(nProjsPC, nBands, nSpins) )
      allocate ( cProjSD(nProjsSD, nBands, nSpins) )
      !
      if ( myid == root ) then
        !
        write(iostd, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4)') ik, nKptsPC
        flush(iostd)
        !
        write(iostd, *)
        write(iostd, '("    Plane waves part begun.")')
        write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> begun.")')
        call cpu_time(t1)
        allocate( wfcPC (numOfPWs, iBandIinit:iBandIfinal), wfcSD (numOfPWs, iBandFinit:iBandFfinal ) )
        !
        call calculatePWsOverlap(ik)
        !
        call cpu_time(t2)
        write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
        !
        write(iostd, '("    Plane waves part done in", f10.2, " secs.")') t2-t1
        write(iostd, *)
        write(iostd, '("    PAW part begun.")')
        !
        call cpu_time(t1)
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> begun.")')
        flush(iostd)
        !
        call readProjectionsPC(ik)
        !
        allocate ( cProjBetaPCPsiSD(nProjsPC, nBands, nSpins) )
        call projectBetaPCwfcSD(ik)
        !
        deallocate ( wfcSD )
        !
        call pawCorrectionPsiPC()
        !
        deallocate ( cProjBetaPCPsiSD )
        !
        call cpu_time(t2)
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> done in", f10.2, " secs.")') t2-t1
        call cpu_time(t1)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> begun.")')
        flush(iostd)
        !
        call readProjectionsSD(ik)
        !
        allocate ( cProjBetaSDPhiPC(nProjsSD, nBands, nSpins) )
        call projectBetaSDwfcPC(ik)
        !
        deallocate ( wfcPC )
        !
        call pawCorrectionSDPhi()
        deallocate ( cProjBetaSDPhiPC )
        !
        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
        !
        call cpu_time(t1)
        write(iostd, '("      <\\vec{k}|PAW_PC> begun.")')
        flush(iostd)
        !
      endif
      !
      call MPI_BCAST(cProjPC, size(cProjPC), MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(cProjSD, size(cProjSD), MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
      !
      allocate ( pawKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nPWsI(myid):nPWsF(myid)) )
      !
      call pawCorrectionKPC()
      !
      if ( myid == root ) then
        call cpu_time(t2)
        write(iostd, '("      <\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
        !
        call cpu_time(t1)
        write(iostd, '("      <PAW_SD|\\vec{k}> begun.")')
        flush(iostd)
        !
      endif
      !
      allocate ( pawSDK(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nPWsI(myid):nPWsF(myid) ) )
      !
      call pawCorrectionSDK()
      !
      if ( myid == root) then
        !
        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\vec{k}> done in", f10.2, " secs.")') t2-t1
        !
        call cpu_time(t1)
        write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> begun.")')
        flush(iostd)
        !
      endif
      !
      paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      !
      do ibi = iBandIinit, iBandIfinal
        !
        do ibf = iBandFinit, iBandFfinal   
          paw_id(ibf,ibi) = sum(pawSDK(ibf,ibi,:)*pawKPC(ibf,ibi,:))
        enddo
        !
      enddo
      !
      if ( myid == root ) paw_SDKKPC(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      !
      CALL MPI_REDUCE(paw_id, paw_SDKKPC, size(paw_id), MPI_DOUBLE_COMPLEX, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      !
      if ( myid == root ) then
        !
        call cpu_time(t2)
        write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
        !
        Ufi(:,:,ik) = Ufi(:,:,ik) + paw_SDPhi(:,:) + paw_PsiPC(:,:) + paw_SDKKPC(:,:)*16.0_dp*pi*pi/omega
        !
        call writeResults(ik)
        !
      endif
      !
      deallocate ( cProjPC, pawKPC )
      deallocate ( cProjSD, pawSDK )
      !
    else
      !
      if ( myid == root ) call readUfis(ik)
      ! 
    endif
    !
  enddo
  !
  if ( allocated ( paw_id ) ) deallocate ( paw_id )
  if ( myid == root ) then
    if ( allocated ( paw_PsiPC ) ) deallocate ( paw_PsiPC )
    if ( allocated ( paw_SDPhi ) ) deallocate ( paw_SDPhi )
  endif
  !
  ! Calculating Vfi
  !
  if ( myid == root ) then
    !
    if (calculateVfis ) call calculateVfiElements()
    !
    ! Finalize Calculation
    ! 
    call finalizeCalculation()
    !
  endif
  !
  call MPI_FINALIZE(ierr)
  !
end program transitionMatrixElements
