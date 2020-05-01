program transitionMatrixElements
  !! @todo Add detailed math derivation and summary for main program @endtodo
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
  real(kind = dp) :: t1
    !! Start time used in various timers
  real(kind = dp) :: t2
    !! End time used in various timers
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
                               calculateVFis, t0)
    ! 
    call readInput(perfectCrystal, solidDefect, elementsPath, ki, kf, calculateVfis, VfisOutput)
    !
    call readPWsSet()
    !
    !> @todo Figure out if need to allocate space for arrays so soon @endtodo
    allocate ( counts(0:numprocs-1) )!, displmnt(0:numprocs-1) )
    allocate ( Ufi(solidDefect%iBandL:solidDefect%iBandH, &
                   perfectCrystal%iBandL:perfectCrystal%iBandH, perfectCrystal%nKpts) )
    allocate ( paw_SDKKPC(solidDefect%iBandL:solidDefect%iBandH, &
                          perfectCrystal%iBandL:perfectCrystal%iBandH) )
    allocate ( perfectCrystal%paw_Wfc(solidDefect%iBandL:solidDefect%iBandH, &
                                      perfectCrystal%iBandL:perfectCrystal%iBandH) )
    allocate ( solidDefect%paw_Wfc(solidDefect%iBandL:solidDefect%iBandH, &
                                   perfectCrystal%iBandL:perfectCrystal%iBandH) )
    allocate ( paw_fi(solidDefect%iBandL:solidDefect%iBandH, perfectCrystal%iBandL:perfectCrystal%iBandH) )
    allocate ( eigvI (perfectCrystal%iBandL:perfectCrystal%iBandH) )
    allocate ( eigvF (solidDefect%iBandL:solidDefect%iBandH) )
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
  !> * Broadcast variables from root process to all other processes, allocating space as needed
  !
  call MPI_BCAST(perfectCrystal%iBandL,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(perfectCrystal%iBandH, 1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%iBandL,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(solidDefect%iBandH, 1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
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
    call MPI_BCAST(perfectCrystal%atoms(i)%iRAugMax,        1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    !
    if ( myid /= root ) then 
      allocate( perfectCrystal%atoms(i)%projAngMom(perfectCrystal%atoms(i)%numProjs) )
      allocate( perfectCrystal%atoms(i)%r  (perfectCrystal%atoms(i)%nMax) )
      allocate( perfectCrystal%atoms(i)%rab(perfectCrystal%atoms(i)%nMax) )
      allocate( perfectCrystal%atoms(i)%F(perfectCrystal%atoms(i)%iRAugMax, perfectCrystal%atoms(i)%numProjs ) )
      allocate(perfectCrystal%atoms(i)%F1(perfectCrystal%atoms(i)%iRAugMax, perfectCrystal%atoms(i)%numProjs, &
               perfectCrystal%atoms(i)%numProjs))
      allocate( perfectCrystal%atoms(i)%bes_J_qr ( 0:JMAX, perfectCrystal%atoms(i)%iRAugMax) )
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
    call MPI_BCAST(solidDefect%atoms(i)%iRAugMax,        1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    !
    if ( myid /= root ) then 
      allocate( solidDefect%atoms(i)%projAngMom(solidDefect%atoms(i)%numProjs) )
      allocate( solidDefect%atoms(i)%r(solidDefect%atoms(i)%nMax) )
      allocate( solidDefect%atoms(i)%rab(solidDefect%atoms(i)%nMax) )
      allocate( solidDefect%atoms(i)%F(solidDefect%atoms(i)%iRAugMax, solidDefect%atoms(i)%numProjs ) )
      allocate( solidDefect%atoms(i)%F1(solidDefect%atoms(i)%iRAugMax, &
                solidDefect%atoms(i)%numProjs, solidDefect%atoms(i)%numProjs ) )
      allocate( solidDefect%atoms(i)%bes_J_qr( 0:JMAX, solidDefect%atoms(i)%iRAugMax) )
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
  allocate ( paw_id(solidDefect%iBandL:solidDefect%iBandH, &
                    perfectCrystal%iBandL:perfectCrystal%iBandH) )
  !
  !--------------------------------------------------------------------------------------------------------
  !
  do ik = 1, perfectCrystal%nKpts
    !! * For each k point
    !!    * If I'm root, check if the matrix elements have already 
    !!      been calculated and send the result to all other processes
    !!    * If the matrix elements haven't already been 
    !!      calculated
    !!       * Allocate space for the projections `cProj`
    !!       * If I'm root
    !!          * Start a timer and output that starting to
    !!            calculate overlap
    !!          * Allocate space for the wavefunctions
    !!          * Read wavefunctions and calculate 
    !!            `U_fi`\(= \langle\Phi_f|\Psi_i\rangle\)
    !!          * Stop timer and output how long calculating overlap 
    !!            took
    !!          * Output that starting projector augmented wave 
    !!            portion
    !!          * Start a timer and report that doing PAW for 
    !!            perfect crystal
    !!          * Read \(\langle\beta_{\Psi}|\Psi\rangle\)
    !!            from [[pw_export_for_tme(program)]]
    !!          * Allocate space for cross projection
    !!            \(\langle\beta_{\Psi}|\Phi\rangle\)
    !!          * Calculate cross projection of \(|\Phi\rangle\)
    !!          * Deallocate space for solid defect wavefunction
    !!          * Calculate the 2nd and 3rd terms in C3 from paper
    !!          * Deallocate space for the cross projection of 
    !!            \(|\Phi\rangle\)
    !!          * Stop timer and output that finished PAW for perfect
    !!            crystal
    !!          * Start timer and output that started PAW for solid 
    !!            defect
    !!          * Read \(\langle\beta_{\Phi}|\Phi\rangle\) from 
    !!            [[pw_export_for_tme(program)]]
    !!          * Allocate space for cross projection 
    !!            \(\langle\beta_{\Phi}|\Psi\rangle\)
    !!          * Calculate cross projection of \(|\Psi\rangle\)
    !!          * Deallocate space for the perfect crystal wavefunction
    !!          * Calculate the 4th and 5th terms in C3 from paper
    !!          * Deallocate space for the cross projection of 
    !!            \(|\Psi\rangle\)
    !!          * Stop timer and output that finished PAW for solid 
    !!            defect
    !!          * Start a timer and output that started k projections 
    !!            for perfect crystal
    !!       * Send projectors to all other processes
    !!       * Allocate space for k projections for perfect crystal
    !!       * Calculate k projections for perfect crystal
    !!       * If I'm root
    !!          * Stop timer and output that done with k projection for 
    !!            perfect crystal
    !!          * Start timer and output that started k projection for 
    !!            solid defect
    !!       * Allocate space for k projections for solid defect
    !!       * Calculate k projections for solid defect
    !!       * If I'm root
    !!          * Stop timer and output that done with k projection for 
    !!            solid defect
    !!          * Start timer and output that combining k projections
    !!       * For each initial and final band, sum the product of the 
    !!         k projections of the solid defect and perfect crystal 
    !!         to get the last term in equation C3 in the paper
    !!       * Allocate space for `paw_SDKKPC` if root and sum 
    !!         individual sums from processes into a single total
    !!       * If I'm root
    !!          * Stop timer and output that done with summing k
    !!            projections
    !!          * Add PAW corrections to initially calculated overlap
    !!          * Output transition matrix elements for a given k point
    !!       * Deallocate space for the projections
    !!    * If the transition matrix elements for this k point have 
    !!      already been calculated and I'm root, read in the existing
    !!      output file
    !
    if ( myid == root ) then
      ! If I'm root, check if the matrix elements have
      ! already been calculated
      !
      tmes_file_exists = .false.
      call checkIfCalculated(ik,tmes_file_exists)
      !
    endif
    !
    call MPI_BCAST(tmes_file_exists, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
      ! Send the result to all other processes
    !
    if ( .not.tmes_file_exists ) then
      ! If the matrix elements haven't already been calculated
      !
      allocate ( perfectCrystal%cProj(perfectCrystal%nProjs, solidDefect%nBands, solidDefect%nSpins) )
      allocate ( solidDefect%cProj(solidDefect%nProjs, solidDefect%nBands, solidDefect%nSpins) )
        ! Allocate space for the projections `cProj`
      !
      if ( myid == root ) then
        ! If I'm root
        !
        write(iostd, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4)') ik, perfectCrystal%nKpts
        flush(iostd)
        !
        write(iostd, *)
        write(iostd, '("    Plane waves part begun.")')
        write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> begun.")')
        call cpu_time(t1)
          ! Start a timer and output that starting to calculate overlap
        !
        allocate( perfectCrystal%wfc (solidDefect%numOfPWs, perfectCrystal%iBandL:perfectCrystal%iBandH), &
                  solidDefect%wfc (solidDefect%numOfPWs, solidDefect%iBandL:solidDefect%iBandH ) )
          ! Allocate space for the wavefunctions
        !
        call calculatePWsOverlap(ik)
          ! Read wavefunctions and calculate `U_fi`\(= \langle\Phi_f|\Psi_i\rangle\)
        !
        call cpu_time(t2)
        write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> done in", f10.2, " secs.")') t2-t1
        write(iostd, '("    Plane waves part done in", f10.2, " secs.")') t2-t1
          ! Stop timer and output how long calculating overlap took
        !
        flush(iostd)
        write(iostd, *)
        write(iostd, '("    PAW part begun.")')
          ! Output that starting projector augmented wave portion
        !
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> begun.")')
        flush(iostd)
        call cpu_time(t1)
          ! Start a timer and report that doing PAW for perfect crystal
        !
        call readProjections(ik, perfectCrystal)
          ! Read \(\langle\beta_{\Psi}|\Psi\rangle\) from [[pw_export_for_tme(program)]]
        !
        allocate ( perfectCrystal%cCrossProj(perfectCrystal%nProjs, solidDefect%nBands, solidDefect%nSpins) )
          ! Allocate space for cross projection \(\langle\beta_{\Psi}|\Phi\rangle\)
        !
        call projectBeta(ik, perfectCrystal, solidDefect)
          ! Calculate cross projection of \(|\Phi\rangle\)
        !
        deallocate ( solidDefect%wfc )
          ! Deallocate space for solid defect wavefunction
        !
        call pawCorrectionWfc(perfectCrystal)
          ! Calculate the 2nd and 3rd terms in C3 from paper
        !
        deallocate ( perfectCrystal%cCrossProj )
          ! Deallocate space for the cross projection of \(|\Phi\rangle\)
        !
        call cpu_time(t2)
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> done in", f10.2, " secs.")') t2-t1
          ! Stop timer and output that finished PAW for perfect crystal
        !
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> begun.")')
        flush(iostd)
        call cpu_time(t1)
          ! Start timer and output that started PAW for solid defect
        !
        call readProjections(ik, solidDefect)
          ! Read \(\langle\beta_{\Phi}|\Phi\rangle\) from [[pw_export_for_tme(program)]]
        !
        allocate ( solidDefect%cCrossProj(solidDefect%nProjs, solidDefect%nBands, solidDefect%nSpins) )
          ! Allocate space for cross projection \(\langle\beta_{\Phi}|\Psi\rangle\)
        !
        call projectBeta(ik, solidDefect, perfectCrystal)
          ! Calculate cross projection of \(|\Psi\rangle\)
        !
        deallocate ( perfectCrystal%wfc )
          ! Deallocate space for the perfect crystal wavefunction
        !
        call pawCorrectionWfc(solidDefect)
          ! Calculate the 4th and 5th terms in C3 from paper
        !
        deallocate ( solidDefect%cCrossProj )
          ! Deallocate space for the cross projection of \(|\Psi\rangle\)
        !
        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
          ! Stop timer and output that finished PAW for solid defect
        !
        write(iostd, '("      <\\vec{k}|PAW_PC> begun.")')
        flush(iostd)
        call cpu_time(t1)
          ! Start a timer and output that started k projections for perfect crystal
        !
        !do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
        !  !
        !  do ibf = ibi, ibi
        !    paw = solidDefect%paw_Wfc(ibf,ibi) + perfectCrystal%paw_Wfc(ibf,ibi)
        !    write(iostd,'(" paw ", 2i4, 6ES14.5E3)') ibi, ibf, solidDefect%paw_Wfc(ibf,ibi), perfectCrystal%paw_Wfc(ibf,ibi), paw
        !  enddo
        !  !
        !  flush(iostd)
        !  !
        !enddo
        !
        !call pawCorrection()
        !do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
        !  !
        !  do ibf = ibi, ibi
        !    paw = solidDefect%paw_Wfc(ibf,ibi) + perfectCrystal%paw_Wfc(ibf,ibi) + paw_fi(ibf,ibi)
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
        ! Send projectors to all other processes
      !
      allocate ( perfectCrystal%pawK(solidDefect%iBandL:solidDefect%iBandH, &
                                     perfectCrystal%iBandL:perfectCrystal%iBandH, nPWsI(myid):nPWsF(myid)) )
        ! Allocate space for k projections for perfect crystal
      !
      call pawCorrectionK(perfectCrystal)
        ! Calculate k projections for perfect crystal
      !
      if ( myid == root ) then
        ! If I'm root
        !
        call cpu_time(t2)
        write(iostd, '("      <\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
          ! Stop timer and output that done with k projection for perfect crystal
        !
        write(iostd, '("      <PAW_SD|\\vec{k}> begun.")')
        flush(iostd)
        call cpu_time(t1)
          ! Start timer and output that started k projection for solid defect  
        !
      endif
      !
      allocate ( solidDefect%pawK(solidDefect%iBandL:solidDefect%iBandH, &
                                  perfectCrystal%iBandL:perfectCrystal%iBandH, nPWsI(myid):nPWsF(myid) ) )
        ! Allocate space for k projections for solid defect
      !
      call pawCorrectionK(solidDefect)
        ! Calculate k projections for solid defect
      !
      if ( myid == root) then
        ! If I'm root
        !
        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\vec{k}> done in", f10.2, " secs.")') t2-t1
          ! Stop timer and output that done with k projection for solid defect
        !
        write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> begun.")')
        flush(iostd)
        call cpu_time(t1)
          ! Start timer and output that combining k projections
        !
      endif
      !
      paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      !
      do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
        !
        do ibf = solidDefect%iBandL, solidDefect%iBandH   
          ! For each initial and final band, sum the product of the k projections
          ! of the solid defect and perfect crystal to get the last term in 
          ! equation C3 in the paper
          !
          paw_id(ibf,ibi) = sum(solidDefect%pawK(ibf,ibi,:)*perfectCrystal%pawK(ibf,ibi,:))
          !
        enddo
        !
      enddo
      !
      if ( myid == root ) paw_SDKKPC(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      CALL MPI_REDUCE(paw_id, paw_SDKKPC, size(paw_id), MPI_DOUBLE_COMPLEX, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        ! Allocate space for `paw_SDKKPC` if root and sum individual sums from 
        ! processes into a single total
      !
      if ( myid == root ) then
        ! If I'm root
        !
        call cpu_time(t2)
        write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
          ! Stop timer and output that done with summing k projections
        !
        Ufi(:,:,ik) = Ufi(:,:,ik) + solidDefect%paw_Wfc(:,:) + perfectCrystal%paw_Wfc(:,:) + &
                      paw_SDKKPC(:,:)*16.0_dp*pi*pi/solidDefect%omega
          ! Add PAW corrections to initially calculated overlap
        !
        call writeResults(ik)
          ! Output transition matrix elements for a given k point
        !
        !write(iostd,*)'--------------------------------------------------------------------------------------------'
        !
        !do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
        !  !
        !  do ibf = solidDefect%iBandL, solidDefect%iBandH
        !    !paw = solidDefect%paw_Wfc(ibf,ibi) + perfectCrystal%paw_Wfc(ibf,ibi) + paw_SDKKPC(ibf,ibi)*16.0_dp*pi*pi/solidDefect%omega
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
      deallocate ( perfectCrystal%cProj, perfectCrystal%pawK )
      deallocate ( solidDefect%cProj, solidDefect%pawK )
        ! Deallocate space for the projections
      !
    else
      ! If the transition matrix elements for this k point have already been 
      ! calculated and I'm root, read in the existing output file
      !
      if ( myid == root ) call readUfis(ik)
      ! 
    endif
    !
  enddo
  !
  if ( allocated ( paw_id ) ) deallocate ( paw_id )
  if ( myid == root ) then
    if ( allocated ( perfectCrystal%paw_Wfc ) ) deallocate ( perfectCrystal%paw_Wfc )
    if ( allocated ( solidDefect%paw_Wfc ) ) deallocate ( solidDefect%paw_Wfc )
  endif
    !! * Deallocate space for the PAW corrections
  !
  !
  if ( myid == root ) then
  !! * If my ID is root, calculate `VfiElements` if needed and finalize
  !!   calculation
    !
    if (calculateVfis ) call calculateVfiElements()
    !
    call finalizeCalculation()
    !
  endif
  !
  call MPI_FINALIZE(ierr)
  !
end program transitionMatrixElements
