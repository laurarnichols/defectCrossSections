program lineShapeFunction
  !
  use mpi
  use lsf
  use generalComputations
    !! Include the `generalComputations` module
    !! for call to `computeGeneralizedDisplacements`
    !! and `computeVariables`
  !
  implicit none
  !
  integer :: lll, iPhonon
  !
  character(len=2) :: charI
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    !! * Initialize mpi and set up processes
  !
  if ( myid == root ) then
    !! * If root process
    !!    * Start a timer
    !!    * Read input, check all variables needed and initialize the calculation
    !!    * Calculate \(\delta q_j\)
    !!    * Compute main parts of equations 42 and 43 in paper
    !!      to make whole formula more manageable
    !!    * Initialize or read LSF from file
    !
    call cpu_time(ti)
      ! Start a timer
    !
    call readInputs()
      ! Read input, check all variables needed and initialize the calculation
    !
    allocate( genCoord(nModes) )
    !
    call computeGeneralizedDisplacements(nOfqPoints, nModes, genCoord, nAtoms, atomM, phonD, atomD)
      ! Calculate \(\delta q_j\)
    !
    deallocate( atomM, phonD, atomD )
    !
    allocate( x(nModes), Sj(nModes), coth(nModes), wby2kT(nModes), s2L(nModes) )
    allocate( besOrderNofModeM(0:maximumNumberOfPhonons + 1, nModes) )
    !
    call computeVariables(x, Sj, coth, wby2kT, phonF, genCoord, kT, s2L, nModes, maximumNumberOfPhonons, &
                          besOrderNofModeM)
      ! Compute main parts of equations 42 and 43 in paper
      ! to make whole formula more manageable
    !
    deallocate( genCoord )
    !
    call initializeLSF()
      ! Initialize or read LSF from file
    !
  endif
  !
  !-----------------------------------------------------------------------------------------------------------
  !> * Allocate space for variables and broadcast to all processes so can move forward
  !>   with calculation
  call MPI_BCAST(nModes   ,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(maximumNumberOfPhonons,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(minimumNumberOfPhonons,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nEnergies,  1, MPI_INTEGER,root,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(deltaE,     1, MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
  if ( myid /= root ) then
    !
    allocate( phonF(nModes), x(nModes), Sj(nModes), coth(nModes), wby2kT(nModes) )
    allocate( besOrderNofModeM(0:maximumNumberOfPhonons + 1, nModes) )
    !allocate( Vfis(-nEnergies:nEnergies) )
    !
  endif
  !
  call MPI_BCAST( phonF, size(phonF), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( x, size(x), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( Sj, size(Sj), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( coth, size(coth), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( wby2kT, size(wby2kT), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  call MPI_BCAST( besOrderNofModeM, size(besOrderNofModeM), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !call MPI_BCAST( Vfis, size(Vfis), MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
  !
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
  !-----------------------------------------------------------------------------------------------------------
  !
  do iPhonon = minimumNumberOfPhonons, MIN0(maximumNumberOfPhonons,4)
    !! * For each possible phonon number less than 5
    !!    * Calculate the line shape function explicitly by 
    !!      summing contributions from using different numbers of bands 
    !!    * For more than 2 phonons, split up the possible configurations 
    !!      among the processes for speed then sum
    !!    * Calculate minimum energy bin size such that no bin is empty
    !!    * Output resulting line shape function
    !! @todo Redo the loop for less than 5 phonons to be more clear and streamlined @endtodo 
    !
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
          !! @todo Change this to have `size(iEbinsByBands)` @endtodo
        CALL MPI_REDUCE(lsfVsEbyBands, lsfVsEbyPhonons, size(lsfVsEbyPhonons), &
                                                                    MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
        !
      endif
        !! @todo Add `else iEbinsByPhonons = iEbinsByBands` to remove if below @endtodo
        !! @todo Maybe change variable names to be clearer @endtodo
      !
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
            !! @todo Figure out how getting `de` is "calculating DOS" and if not where DOS is @endtodo
            !! @todo Figure out why DOS isn't in sum as in formula @endtodo
          !
          lsfVsE(:) = lsfVsE(:) + lsfVsEbyPhonons(:)/de
            ! Add so that can continue from file if needed
          !
        else
          !
          ! calculate the DOS and update the total lsfVsE
          !
          call calculateDE(iPhonon, iEbinsByBands, de)
          !
          lsfVsE(:) = lsfVsE(:) + lsfVsEbyBands(:)/de
            ! Add so that can continue from file if needed
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
    allocate ( lsfbyPhononsPerProc(-nEnergies:nEnergies) )
    !
    if ( minimumNumberOfPhonons < 6 ) minimumNumberOfPhonons = 5
    do m = minimumNumberOfPhonons, maximumNumberOfPhonons
      !! * For each possible phonon number greater than or equal to 5
      !!    * Calculate the contribution to the line shape function of splitting
      !!      up the phonons in 1-3 bands explicitly, splitting up the possible 
      !!      configurations among the processes for speed then summing all 
      !!      contributions
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
        ! If I'm root, figure out how the Monte Carlo
        ! steps should be split up among the processes
        !! @todo Move the behavior of splitting up Monte Carlo steps to a subroutine @endtodo
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
      do l = 4, m
        !! @todo Move this to a subroutine @endtodo
        !
        !write(iostd,*) "---------------------------------"
        !write(iostd,*) m, " by ", l
        !flush(iostd)
        !
        !> @todo Replace this with `binomialCoefficient(kPhonons-1, kPhonons-nBands)` @endtodo
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
        !
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
    !! * Write out the configuration with the most phonons
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
