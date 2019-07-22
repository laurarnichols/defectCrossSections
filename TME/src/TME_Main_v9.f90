program transitionMatrixElements
  !! @todo Finish documentation for main program @endtodo
  !!
  !! <h2>Walkthrough</h2>
  !!
  use mpi
  use TMEModule
    !! * Use pre-built mpi library and declarations module that 
    !! is defined in TME_Module_v28.f90
  !
  implicit none
  !
  real(kind = dp) :: t1, t2
    !! * Declare start and end times
  !
  call MPI_INIT(ierr)
    !! * Initialize MPI environment
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    !! * Determine the rank or ID of the calling process
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    !! * Determine the size of the MPI pool (i.e., the number of processes)
  !
  allocate ( nPWsI(0:numprocs-1), nPWsF(0:numprocs-1) )
  !
  if ( myid == root ) then
    !! * If this is the root process
    !!    * Call [[TMEModule(module):initializeCalculation(subroutine)]] 
    !!      to start timer and set default values
    !!    * Call [[TMEModule(module):readInput(subroutine)]] to read program input and
    !!      QE Export output
    !!    * Call [[TMEModule(module):readPWsSet(subroutine)]] to read g vectors from `mgrid` file
    !!    * Initialize all values in `Ufi` matrix to complex double zero
    !!    * Figure out how many g vectors/plane waves to give each process
    !!    * Initialize the number of initial and final plane waves to zero for each process
    !!    * For each process, calculate the initial (before this process) and final (after this process)
    !!      number of plane waves
    !
    call initializeCalculation(solidDefect, perfectCrystal, elementsPath, VFisOutput, ki, kf, eBin, &
                               iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, calculateVFis, t0)
    ! 
    call readInput(perfectCrystal, solidDefect, elementsPath, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, &
                       ki, kf, calculateVfis, VfisOutput)
    !
    call readPWsSet()
    !
    !> @todo Figure out if need to allocate space for arrays so soon @endtodo
    allocate ( counts(0:numprocs-1) )!, displmnt(0:numprocs-1) )
    allocate ( Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, perfectCrystal%nKpts) )
    allocate ( paw_SDKKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_PsiPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_SDPhi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_fi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( eigvI (iBandIinit:iBandIfinal), eigvF (iBandFinit:iBandFfinal) )
    !
    Ufi(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    call distributePWsToProcs(solidDefect%numOfGvecs, numprocs)
      !! @todo Figure out if SD and PC `numOfGvecs` should be the same @endtodo
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
  !--------------------------------------------------------------------------------------------------------
  !> Broadcast variables from root process to all other processes, allocating space as needed
  !
  call MPI_BCAST(iBandIinit,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iBandFinit,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(perfectCrystal%nKpts,     1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(perfectCrystal%nProjs,    1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%nProjs,    1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%nBands,      1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%nSpins,      1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%numOfPWs,    1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%numOfGvecs,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(nPWsI, numprocs, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nPWsF, numprocs, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  call MPI_BCAST(perfectCrystal%numOfTypes, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(JMAX, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate ( gvecs(3, solidDefect%numOfGvecs) )
  call MPI_BCAST(gvecs, size(gvecs), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  if ( myid /= root ) allocate ( perfectCrystal%atoms(perfectCrystal%numOfTypes) )
  !
  do i = 1, perfectCrystal%numOfTypes
    !
    call MPI_BCAST(perfectCrystal%atoms(i)%numOfAtoms, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%numProjs,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%lmMax,      1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%nMax,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%iRc,        1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    !
    if ( myid /= root ) then 
      allocate( perfectCrystal%atoms(i)%projAngMom(perfectCrystal%atoms(i)%numProjs) )
      allocate( perfectCrystal%atoms(i)%r  (perfectCrystal%atoms(i)%nMax) )
      allocate( perfectCrystal%atoms(i)%rab(perfectCrystal%atoms(i)%nMax) )
      allocate( perfectCrystal%atoms(i)%F(perfectCrystal%atoms(i)%iRc, perfectCrystal%atoms(i)%numProjs ) )
      allocate(perfectCrystal%atoms(i)%F1(perfectCrystal%atoms(i)%iRc, perfectCrystal%atoms(i)%numProjs, &
               perfectCrystal%atoms(i)%numProjs))
      allocate( perfectCrystal%atoms(i)%bes_J_qr ( 0:JMAX, perfectCrystal%atoms(i)%iRc) )
    endif
    !
    call MPI_BCAST(perfectCrystal%atoms(i)%projAngMom, size(perfectCrystal%atoms(i)%projAngMom), &
                    MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%r,   size(perfectCrystal%atoms(i)%r),   MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%rab, size(perfectCrystal%atoms(i)%rab), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%F,   size(perfectCrystal%atoms(i)%F),   MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%F1,  size(perfectCrystal%atoms(i)%F1),  MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(perfectCrystal%atoms(i)%bes_J_qr, size(perfectCrystal%atoms(i)%bes_J_qr), &
                   MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  enddo
  !
  call MPI_BCAST(perfectCrystal%nIons, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate( perfectCrystal%posIon(3, perfectCrystal%nIons), perfectCrystal%atomTypeIndex(perfectCrystal%nIons) )
  call MPI_BCAST(perfectCrystal%atomTypeIndex,  size(perfectCrystal%atomTypeIndex),  MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(perfectCrystal%posIon, size(perfectCrystal%posIon), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  call MPI_BCAST(solidDefect%numOfTypes, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate ( solidDefect%atoms(solidDefect%numOfTypes) )
  !
  do i = 1, solidDefect%numOfTypes
    !
    call MPI_BCAST(solidDefect%atoms(i)%numOfAtoms, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(solidDefect%atoms(i)%numProjs,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(solidDefect%atoms(i)%lmMax,      1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(solidDefect%atoms(i)%nMax,       1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(solidDefect%atoms(i)%iRc,        1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    !
    if ( myid /= root ) then 
      allocate( solidDefect%atoms(i)%projAngMom(solidDefect%atoms(i)%numProjs) )
      allocate( solidDefect%atoms(i)%r(solidDefect%atoms(i)%nMax) )
      allocate( solidDefect%atoms(i)%rab(solidDefect%atoms(i)%nMax) )
      allocate( solidDefect%atoms(i)%F(solidDefect%atoms(i)%iRc, solidDefect%atoms(i)%numProjs ) )
      allocate( solidDefect%atoms(i)%F1(solidDefect%atoms(i)%iRc, solidDefect%atoms(i)%numProjs, solidDefect%atoms(i)%numProjs ) )
      allocate( solidDefect%atoms(i)%bes_J_qr( 0:JMAX, solidDefect%atoms(i)%iRc) )
    endif
    !
    call MPI_BCAST(solidDefect%atoms(i)%projAngMom, size(solidDefect%atoms(i)%projAngMom), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(solidDefect%atoms(i)%r, size(solidDefect%atoms(i)%r), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(solidDefect%atoms(i)%rab, size(solidDefect%atoms(i)%rab), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(solidDefect%atoms(i)%F, size(solidDefect%atoms(i)%F), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(solidDefect%atoms(i)%F1, size(solidDefect%atoms(i)%F1), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(solidDefect%atoms(i)%bes_J_qr, size(solidDefect%atoms(i)%bes_J_qr), MPI_DOUBLE_PRECISION, &
                   root,MPI_COMM_WORLD,ierr)
  enddo
  !
  call MPI_BCAST(solidDefect%nIons, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  !
  if ( myid /= root ) allocate( solidDefect%posIon(3, solidDefect%nIons), solidDefect%atomTypeIndex(solidDefect%nIons) )
  call MPI_BCAST(solidDefect%atomTypeIndex,  size(solidDefect%atomTypeIndex),  MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(solidDefect%posIon, size(solidDefect%posIon), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  allocate ( paw_id(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
  !
  !--------------------------------------------------------------------------------------------------------
  !
  do ik = 1, perfectCrystal%nKpts
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
      allocate ( perfectCrystal%cProj(perfectCrystal%nProjs, solidDefect%nBands, solidDefect%nSpins) )
      allocate ( solidDefect%cProj(solidDefect%nProjs, solidDefect%nBands, solidDefect%nSpins) )
      !
      if ( myid == root ) then
        !
        write(iostd, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4)') ik, perfectCrystal%nKpts
        flush(iostd)
        !
        write(iostd, *)
        write(iostd, '("    Plane waves part begun.")')
        write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> begun.")')
        call cpu_time(t1)
        allocate( perfectCrystal%wfc (solidDefect%numOfPWs, iBandIinit:iBandIfinal), &
                  solidDefect%wfc (solidDefect%numOfPWs, iBandFinit:iBandFfinal ) )
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
        call readProjections(ik, perfectCrystal)
        !
        allocate ( cProjBetaPCPsiSD(perfectCrystal%nProjs, solidDefect%nBands, solidDefect%nSpins) )
        call projectBeta(ik, perfectCrystal, solidDefect)
        !
        deallocate ( solidDefect%wfc )
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
        call readProjections(ik, solidDefect)
        !
        allocate ( cProjBetaSDPhiPC(solidDefect%nProjs, solidDefect%nBands, solidDefect%nSpins) )
        call projectBeta(ik, solidDefect, perfectCrystal)
        !
        deallocate ( perfectCrystal%wfc )
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
        !do ibi = iBandIinit, iBandIfinal
        !  !
        !  do ibf = ibi, ibi
        !    paw = paw_SDPhi(ibf,ibi) + paw_PsiPC(ibf,ibi)
        !    write(iostd,'(" paw ", 2i4, 6ES14.5E3)') ibi, ibf, paw_SDPhi(ibf,ibi), paw_PsiPC(ibf,ibi), paw
        !  enddo
        !  !
        !  flush(iostd)
        !  !
        !enddo
        !
        !call pawCorrection()
        !do ibi = iBandIinit, iBandIfinal
        !  !
        !  do ibf = ibi, ibi
        !    paw = paw_SDPhi(ibf,ibi) + paw_PsiPC(ibf,ibi) + paw_fi(ibf,ibi)
        !    write(iostd,'(" paw ", 2i4, 6f15.10)') ibi, ibf, Ufi(ibf,ibi,ik), paw, Ufi(ibf,ibi,ik) + paw
        !  enddo
        !  !
        !  flush(iostd)
        !  !
        !enddo
        !
      endif
      !
      call MPI_BCAST(perfectCrystal%cProj, size(perfectCrystal%cProj), MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(solidDefect%cProj, size(solidDefect%cProj), MPI_DOUBLE_COMPLEX, root, MPI_COMM_WORLD, ierr)
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
        Ufi(:,:,ik) = Ufi(:,:,ik) + paw_SDPhi(:,:) + paw_PsiPC(:,:) + paw_SDKKPC(:,:)*16.0_dp*pi*pi/solidDefect%omega
          !! @todo Figure out if should be solid defect volume or pristine @endtodo
          !! @todo Are pristine and solid defect volume the same? @endtodo
        !
        call writeResults(ik)
        !
        !write(iostd,*)'--------------------------------------------------------------------------------------------'
        !
        !do ibi = iBandIinit, iBandIfinal
        !  !
        !  do ibf = iBandFinit, iBandFfinal
        !    !paw = paw_SDPhi(ibf,ibi) + paw_PsiPC(ibf,ibi) + paw_SDKKPC(ibf,ibi)*16.0_dp*pi*pi/solidDefect%omega
        !    !write(iostd,'(" paw ", 2i4, 6f15.10)') ibi, ibf, Ufi(ibf, ibi, ik), paw, Ufi(ibf, ibi, ik) + paw
        !    !Ufi(ibf, ibi, ik) = Ufi(ibf, ibi, ik) + paw
        !    write(iostd,'(" Ufi ", 2i4, 2ES24.15E3)') ibi, ibf, Ufi(ibf, ibi, ik)
        !  enddo
        !  !
        !  flush(iostd)
        !  !
        !enddo
        !
        !write(iostd,*)'--------------------------------------------------------------------------------------------'
        !flush(iostd)
        !
      endif
      !
      deallocate ( perfectCrystal%cProj, pawKPC )
      deallocate ( solidDefect%cProj, pawSDK )
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
