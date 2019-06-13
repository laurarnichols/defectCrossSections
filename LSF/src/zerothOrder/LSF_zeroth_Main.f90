program lineShapeFunction
  !
  ! Pull in modules
  use mpi
  use lsf
  !
  implicit none
  !
  ! Define an integer for ????
  integer :: lll, iPhonon
  !
  character(len=2) :: charI
  !
  ! Initialize mpi and set up processes
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  !
  ! If root process
  if ( myid == root ) then
    !
    ! Start a timer
    call cpu_time(ti)
    !
    ! Read input, check all variables needed and initialize the calculation.
    call readInputs()
    !
    call computeGeneralizedDisplacements()
    !
    call computeVariables()
    !
    call initializeLSF()
    !
  endif
  !
  ! Broadcast calculation parameters to all processes
  call MPI_BCAST(nModes   ,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(maximumNumberOfPhonons,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(minimumNumberOfPhonons,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nEnergies,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(deltaE,     1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  ! For all processes that aren't the root
  if ( myid /= root ) then
    ! Allocate space for arrays
    allocate( phonF(nModes), x(nModes), Sj(nModes), coth(nModes), wby2kT(nModes) )
    allocate( besOrderNofModeM(0:maximumNumberOfPhonons + 1, nModes) )
    !allocate( Vfis(-nEnergies:nEnergies) )
  endif
  !
  ! Broadcast arrays to all processes
  call MPI_BCAST( phonF, size(phonF), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( x, size(x), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( Sj, size(Sj), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( coth, size(coth), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( wby2kT, size(wby2kT), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( besOrderNofModeM, size(besOrderNofModeM), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !call MPI_BCAST( Vfis, size(Vfis), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  ! Allocate space for arrays
  allocate( lsfVsEbyBands(-nEnergies:nEnergies) )
  allocate( iEbinsByBands(-nEnergies:nEnergies) )
  !
  allocate( pj(nModes) )
  !
  if ( 3 <= maximumNumberOfPhonons ) then
    !
    allocate( iModeIs(0:numprocs-1) )
    allocate( iModeFs(0:numprocs-1) )
    !
    allocate ( iEbinsByPhonons(-nEnergies:nEnergies), lsfVsEbyPhonons(-nEnergies:nEnergies) )
    !
  end if
  !
  do iPhonon = minimumNumberOfPhonons, MIN0(maximumNumberOfPhonons,4)
    if ( ( ( iPhonon == 1 .or. iPhonon == 2 ) .and. myid == root ) .or. iPhonon > 2 ) then
      if ( iPhonon > 2 ) then
        !
        iModeIs(:) =  0
        iModeFs(:) = -1
        !
      endif
      !
      lsfVsEbyBands(:) = 0.0_dp
      iEbinsByBands(:) = 0
      ! 
      if ( myid == root ) then
        call cpu_time(t1)
        !
        call lsfMbyOneBand(iPhonon)
        !
        if ( iPhonon > 1 ) then
          !
          call lsfMbyTwoBands(iPhonon)
          !
        else if ( iPhonon > 2 ) then
          !
          call parallelIsFsBy3()
          !
        endif
        !
      endif
      !
      if ( iPhonon > 2 ) then
        !
        call MPI_BCAST(iModeIs, size(iModeIs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iModeFs, size(iModeFs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
        !
        call lsfMbyThreeBands(iPhonon)
        !
      endif
      !
      if ( iPhonon > 3 ) then
        !
        iModeIs(:) =  0
        iModeFs(:) = -1
        !
        if ( myid == root ) call parallelIsFsBy4()
        !
        call MPI_BCAST(iModeIs, size(iModeIs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(iModeFs, size(iModeFs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
        !
        call lsfDeterministicFourPhononsByFourBands()
        !
      endif
      !
      if ( iPhonon > 2 ) then
        !
        iEbinsByPhonons = 0
        lsfVsEbyPhonons = 0.0_dp
        !
        CALL MPI_REDUCE(iEbinsByBands, iEbinsByPhonons, size(iEbinsByBands), MPI_INTEGER, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        CALL MPI_REDUCE(lsfVsEbyBands, lsfVsEbyPhonons, size(lsfVsEbyPhonons), &
                                                                    MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        !
      endif
      if ( myid == root ) then
        !
        call cpu_time(t2)
        !
        write(iostd,'(i2, " modes, time needed :," , f10.2, " secs.")') iPhonon, t2-t1
        flush(iostd)
        !
        if ( iPhonon > 2 ) then
          ! calculate the DOS and update the total lsfVsE
          !
          call calculateDE(iPhonon, iEbinsByPhonons, de)
          lsfVsE(:) = lsfVsE(:) + lsfVsEbyPhonons(:)/de
          !
        else
          !
          ! calculate the DOS and update the total lsfVsE
          !
          call calculateDE(iPhonon, iEbinsByBands, de)
          lsfVsE(:) = lsfVsE(:) + lsfVsEbyBands(:)/de
          !
        endif
        !
        write(iostd,*) 'DE', iPhonon,  de
        flush(iostd)
        !
        charI = ''
        write(charI, "(i2.2)") iPhonon
        !
        open(1, file='lsfVsEwithUpTo'//trim(charI)//'phonons', status='unknown')
        !
        write(1,'("#", i10, " energies", i5, " phonons")') nEnergies, iPhonon
        !
        do iE = -nEnergies, nEnergies
          E = real(iE, dp)*deltaE
          !
          if ( iPhonon < 3 ) then
            !
            write(1,'(F16.8,2E18.6e3)') E*HartreeToEv, lsfVsE(iE), lsfVsEbyBands(iE)/de
            !
          else
            !
            write(1,'(F16.8,2E18.6e3)') E*HartreeToEv, lsfVsE(iE), lsfVsEbyPhonons(iE)/de
            !
          endif
          !
        enddo
        !
        close(1)
        !
      endif
      !
    endif 
    !
  enddo
  !
  if ( maximumNumberOfPhonons >= 5 ) then
    !
    open(unit=un, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
    !
    if ( myid == root ) then
      !
      if (istat /= 0) then
        !
        write(iostd, *) 'File "/dev/urandom" not found! A pseudo random generator will be used!'
        !
      else
        !
        write(iostd, *) 'File "/dev/urandom" will be used to generate real random numbers!'
        !
      endif
      !
      flush(iostd)
      !
    endif
    !
    if (istat /= 0) close(un)
    !
    ! allocate ( lsfbyPhononsPerProc(-nEnergies:nEnergies) )
    !
    if ( minimumNumberOfPhonons < 6 ) minimumNumberOfPhonons = 5
    do m = minimumNumberOfPhonons, maximumNumberOfPhonons
      !
      lsfVsEbyBands(:) = 0.0_dp
      iEbinsByBands(:) = 0
      !
      iModeIs(:) =  0
      iModeFs(:) = -1
      !
      if ( myid == root ) then
        !
        call lsfMbyOneBand(m)
        call lsfMbyTwoBands(m)
        !
        call parallelIsFsBy3()
        !
      endif
      !
      call MPI_BCAST(iModeIs, size(iModeIs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iModeFs, size(iModeFs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      !
      call lsfMbyThreeBands(m)
      !
      !iEbinsByPhonons = 0
      lsfVsEbyPhonons = 0.0_dp
      !
      !CALL MPI_REDUCE(iEbinsByBands, iEbinsByPhonons, size(iEbinsByBands), MPI_INTEGER, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(lsfVsEbyBands, lsfVsEbyPhonons, size(lsfVsEbyPhonons), &
                                                                  MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      !
      if (istat /= 0) CALL init_random_seed() 
      !
      iModeIs(:) =  0
      iModeFs(:) = -1
      !
      if ( myid == root ) then
        !
        iMint = int(nMC/numprocs)
        iMmod = mod(nMC,numprocs)
        !
        iModeIs(0) = 1
        iModeFs(numprocs-1) = nMC
        !
        do i = numprocs - 1, 1, -1
          !
          iModeIs(i) = i*iMint + 1
          !
          if ( iMmod > 0 ) then
            !
            iModeIs(i) = iModeIs(i) + iMmod
            iMmod = iMmod - 1
            !
          endif
          !
          iModeFs(i-1) = iModeIs(i) - 1
          !
        enddo
        !
      endif
      !
      call MPI_BCAST(iModeIs, size(iModeIs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iModeFs, size(iModeFs), MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
      !
      !write(6,*) myid, iModeIs(myid), iModeFs(myid)
      !
      do l = 4, m
        !
        !write(iostd,*) "---------------------------------"
        !write(iostd,*) m, " by ", l
        !flush(iostd)
        !
        times = 1.0_dp
        mi = l-1
        !
        do ni = m - 1, m - l + 1, -1
          !
          times = times*dble(ni)/dble(mi)
          mi = mi - 1
          !
        enddo
        !
        allocate( pj0s(int(times + 1.e-3_dp), l) )
        !
        pj0s(:,:) = 0
        !
        !write(6,*) 'distrubutePhononsInBands', m, l, times, int(times + 1.e-3_dp)
        call distrubutePhononsInBands(m, l) 
        !
        allocate( pms( 0:2**l-1, 0:l-1 ) )
        !
        pms(:,:) = 0
        !
        call calculatePlusMinusStates(l)
        !
        lsfVsEbyBands(:) = 0.0_dp
        !
        call lsfWithMphonons(m, l, int(times + 1.e-3_dp))
        !
        lsfbyPhononsPerProc(:) = 0.0_dp
        CALL MPI_REDUCE(lsfVsEbyBands, lsfbyPhononsPerProc, size(lsfbyPhononsPerProc), &
                                                            MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)

        if ( myid == root ) then
          !
          weight = nModes
          !
          do iMode = 2, l
            !
            weight = weight*(nModes - iMode + 1)/iMode
            !
          enddo
          !
          write(iostd, 101) m, l, times*weight
          write(iostd, 102) m, l, real(nMC, dp)
          write(iostd, 103) m, l, times*real(nMC, dp)
          write(iostd, 104) weight/real(nMC, dp)
          flush(iostd)
          !
          lsfVsEbyPhonons(:) = lsfVsEbyPhonons(:) + lsfbyPhononsPerProc(:)*weight/real(nMC, dp)
          !
          !do iE = -nEnergies, nEnergies
          !  write(1000 + 10*m + l, *) real(iE, dp)*deltaE*HartreeToEv, lsfVsEbyBands(iE)*(weight/real(nMC, dp))
          !enddo
          !close(1000 + 10*m + l)
          !
        endif
        !
        deallocate( pj0s, pms )
        !
      enddo
      !
      iEbinsByPhonons = 0
      CALL MPI_REDUCE(iEbinsByBands, iEbinsByPhonons, size(iEbinsByBands), MPI_INTEGER, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      !
      if ( myid == root ) then
        ! 
        call calculateDE(m, iEbinsByPhonons, de)
        lsfVsE(:) = lsfVsE(:) + lsfVsEbyPhonons(:)/de
        !
        write(iostd,*) 'DE', m,  de
        flush(iostd)
        !
        if ( m < 10 ) then
          !
          write(fn,'("lsfVsEwithUpTo", i1, "phonons")') m
          !
        elseif ( m < 100 ) then
          !
          write(fn,'("lsfVsEwithUpTo", i2, "phonons")') m
          !
        elseif ( m < 1000 ) then
          !
          write(fn,'("lsfVsEwithUpTo", i3, "phonons")') m
          !
        else
          !
          write(fn,'("lsfVsEwithUpTo", i4, "phonons")') m
          !
        endif
        !
        open(unit=5000, file=trim(fn), status='unknown')
        !
        !write(5000,'("# ", i5, " phonons")') m
        write(5000,'("#", i10, " energies", i5, " phonons")') nEnergies, m
        !
        do iE = -nEnergies, nEnergies
          !
          E = real(iE, dp)*deltaE
          !
          !vg = 1.0_dp
          !if (E > 0.0_dp) vg = sqrt(2.0_dp*E)
          write(5000,'(F16.8,2E18.6e3)') E*HartreeToEv, lsfVsE(iE), lsfVsEbyPhonons(iE)/de
          !twoPi*abCM**2*volume*Vfis(iE)*lsfVsE(iE)/vg
          !
          !write(5000, *) real(iE, dp)*deltaE*HartreeToEv, lsfVsE(iE), lsfVsEbyPhonons(iE)/de
          !
        enddo
        !
        close(5000)
        !
      endif
      !
    enddo
    !
    if (istat == 0) close(un)
    !
  endif
  !
  if ( myid == root ) then
    !
    call writeLSFandCrossSection()
    !
    call cpu_time(tf)
    !
    write(iostd,'(" Time needed: ", f10.2, " secs.")') tf-ti
    !
  endif
  !
  101 format("   Total number of configurations of ", i4, " phonons by ", i4, " bands : ", E20.10E3)
  102 format("   Total number of configurations of ", i4, " phonons by ", i4, " bands sampled : ", E20.10E3)
  103 format("   Total number of configurations of ", i4, " phonons by ", i4, " bands calculated : ", E20.10E3)
  104 format("   Each sampled configuration will be weighted by : ", E20.10E3)
  !
  deallocate( lsfVsEbyBands, iEbinsByBands, pj )
  !
  if ( 3 <= maximumNumberOfPhonons ) then
    !
    deallocate( iModeIs, iModeFs )
    !
    deallocate ( iEbinsByPhonons, lsfVsEbyPhonons )
    !
  end if
  !
  call MPI_FINALIZE(ierr)
  !
end program lineShapeFunction
