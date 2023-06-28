
module lsf
  !
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
  integer, parameter :: int32 = selected_int_kind(5)
  integer, parameter :: int64 = selected_int_kind(15)
  integer, parameter :: iostd = 16, un = 3
  integer, parameter :: root = 0
  !
  real(kind = dp), parameter ::           pi = 3.1415926535897932_dp
  real(kind = dp), parameter ::        twopi = 2.0_dp*pi
  real(kind = dp), parameter ::         abCM = 0.529177219217e-8_dp
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
  real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
  !
  integer(kind = int32) :: myid, numprocs, ios, istat, ierr
  integer :: iMode, l, m, nMC, nProcMax
  integer :: iMint, iMmod, i, printsteps, iE, ni, mi
  integer :: nAtoms, nOfqPoints, nModes, minimumNumberOfPhonons, maximumNumberOfPhonons, nEnergies
  !
  real(kind = dp) :: ti, tf, t1, t2
  real(kind = dp) :: weight, times, de, E ! , vg
  real(kind = dp) :: temperature, maxEnergy, deltaE, kT ! , volume
  !
  integer, allocatable :: iModeIs(:), iModeFs(:)
  integer, allocatable :: pj(:), pj0s(:,:), pms(:,:), s2L(:)
  integer, allocatable :: iEbinsByBands(:), iEbinsByPhonons(:)
  !
  real(kind = dp), allocatable :: atomD(:,:), atomM(:), phonQ(:,:), phonF(:), genCoord(:), Mjs(:,:)
  real(kind = dp), allocatable :: wby2kT(:), phonD(:,:,:,:), x(:), Sj(:), coth(:), besOrderNofModeM(:,:)
  real(kind = dp), allocatable :: lsfVsEbyBands(:), lsfVsE(:), lsfVsEbyPhonons(:), lsfbyPhononsPerProc(:)
  !
  integer :: modes
  !
  character(len = 6), parameter :: output = 'status'
!  character(len = 256) :: MjsInput, PhononsInput, crossSectionOutput, fn, continueLSFfromFile
  character(len = 256) :: MjsInput, phononsInput, fn, continueLSFfromFile, equilibriumAtomicPositions
  !
  logical :: file_exists
  !
!  namelist /elphscat/ MjsInput, PhononsInput, crossSectionOutput, temperature, maxEnergy, continueLSFfromFile, volume, &
  namelist /lsfInput/ MjsInput, equilibriumAtomicPositions, phononsInput, &
                      continueLSFfromFile, maximumNumberOfPhonons, nMC, &
                      temperature, maxEnergy, &
                      modes
  !
  !
contains
  !
  !
  subroutine readInputs()
    !
    implicit none
    !
    ! Check if file output exists. If it does, delete it.
    !
    inquire(file = output, exist = file_exists)
    if ( file_exists ) then
      open (unit = 11, file = output, status = "old")
      close(unit = 11, status = "delete")
    endif
    !
    ! Open new output file.
    !
    open (iostd, file = output, status='new')
    !
    call initialize()
    !
    READ (5, lsfInput, iostat = ios)
    !
    call checkAndUpdateInput()
    !
    return
    !
  end subroutine readInputs
  !
  !
  subroutine initializeLSF()
    !
    implicit none
    !
    real(kind = dp) :: dummyD
    integer :: dummyI
    character (len = 1) :: dummyC1
    character (len = 8) :: dummyC8
    character (len = 9) :: dummyC9
    !
    allocate ( lsfVsE(-nEnergies:nEnergies) )
    !
    minimumNumberOfPhonons = 1
    lsfVsE(:) = 0.0_dp
    !
    if (continueLSFfromFile /= '') then
      !
      inquire(file = trim(continueLSFfromFile), exist = file_exists)
      if ( file_exists ) then
        !
        open(unit = 11, file = trim(continueLSFfromFile), status = "old")
        !
        read(11,'(a1, i10, a9, i5, a8)') dummyC1, dummyI, dummyC9, minimumNumberOfPhonons, dummyC8
        ! 
        minimumNumberOfPhonons = minimumNumberOfPhonons + 1
        write(iostd, '(" Minimum number of phonons : ", i5)') minimumNumberOfPhonons
        !
        do iE = -nEnergies, nEnergies
          read(11,*) dummyD, lsfVsE(iE)
        enddo
        close(11)
        !
      endif
    endif
    !
    return
    !
  end subroutine initializeLSF
  !
  !
  subroutine initialize()
    !
    implicit none
    !
    MjsInput = ''
    phononsInput = ''
    !crossSectionOutput = ''
    temperature = -1.0_dp
    deltaE = -1.0_dp
    minimumNumberOfPhonons =  1
    maximumNumberOfPhonons = -1
    nMC = -1
    modes = -1
    !
    return
    !
  end subroutine initialize
  !
  !
  subroutine checkAndUpdateInput()
    !
    implicit none
    !
    logical :: abortExecution = .false.
    !
    if ( MjsInput == '' ) then
      write(iostd, '(" MjsInput is not defined!")')
      write(iostd, '(" A default value will be used. MjsInput = VfisVsE")')
      MjsInput = 'VfisVsE'
    else
      write(iostd, '(" Mjs elements input : ", a)') trim(MjsInput)
    endif
    !
    if ( phononsInput == '' ) then
      write(iostd, '(" PhononsInput is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Phonons input : ", a)') trim(PhononsInput)
    endif
    !
    if ( temperature < 0.0_dp ) then
      write(iostd, '(" Variable temperature has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Tempetature : ", f10.2, " Kelvin.")') temperature
      kT = temperature*8.6173324d-5*eVToHartree
    endif
    !
!    if ( deltaE < 0 ) then
!      write(iostd, '(" Variable deltaE has not been set.")')
!    else
!      write(iostd, '(" DeltaE : ", f10.5, " eV.")') deltaE
!      deltaE = deltaE*eVToHartree
!    endif
    !
    if ( maximumNumberOfPhonons < 0 ) then
      write(iostd, '(" Variable maximumNumberOfPhonons has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Maximum number of phonons : ", i5)') maximumNumberOfPhonons
    endif
    if ( nMC  < 0 ) then
      write(iostd, '(" Variable nMC has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Number of Monte Carlo steps : ", i15)') nMC
    endif
    !
    if ( modes < 0 ) then
      write(iostd, '(" Variable modes has not been set.")')
      write(iostd, '(" All modes will be used.")')
    endif
    !
    if ( abortExecution ) then
      write(iostd, '(" *************************** ")')
      write(iostd, '(" * Program stops!          * ")')
      write(iostd, '(" *************************** ")')
      stop
    endif
    !
    maxEnergy = 10.0_dp*eVToHartree
    !
    nEnergies = 5040 ! 2520 ! 10080
    deltaE = maxEnergy/real(nEnergies, dp)
    !
    write(iostd,*) 'nEnergies', nEnergies
    write(iostd,*) 'maxEnergy', maxEnergy, 'deltaE', deltaE
    !
    flush(iostd)
    !
    return
    !
  end subroutine checkAndUpdateInput
  !
  !
  subroutine readPhonons()
    !
    implicit none
    !
    integer :: iAtom, iMode, iq
    real(kind = dp) :: dummyD, freqInTHz
    !
    CHARACTER :: dummyC
    !
    !write(6,*) trim(PhononsInput)
    open(1, file=trim(PhononsInput), status="old")
    !
    read(1,*) nOfqPoints, nAtoms
    !
    nModes = 3*nAtoms - 3
    !
    write(iostd, '(" Number of atoms    : ", i5)') nAtoms
    write(iostd, '(" Number of q-Points : ", i5)') nOfqPoints
    write(iostd, '(" Number of modes    : ", i5)') nModes
    if ( modes < 0 ) then
      modes = nModes
      write(iostd, '(" Number of modes to be used : ", i5)') modes
    endif
    flush(iostd)
    !
    read (1,*)
    !
    allocate( atomD(3,nAtoms), atomM(nAtoms) )
    !
    atomD = 0.0_dp
    atomM = 0.0_dp
    !
    do iAtom = 1, nAtoms
      read(1,*) atomD(1,iAtom), atomD(2,iAtom), atomD(3,iAtom), atomM(iAtom)
    enddo
    !
    read(1,*)
    !
    allocate( phonQ(3,nOfqPoints), phonF(nModes), phonD(3,nAtoms,nModes,nOfqPoints) )
    !
    phonQ = 0.0_dp
    phonF = 0.0_dp
    phonD = 0.0_dp
    !
    do iq = 1, nOfqPoints
      !
      read (1,*) dummyC, dummyC, dummyC, phonQ(1,iq), phonQ(2,iq), phonQ(3,iq), dummyC
      !
      do iMode = 1, nModes
        !
        read(1,*)
        !
        read(1,*) freqInTHz, dummyC, dummyD, dummyC, dummyD, dummyC, dummyD, dummyC
        phonF(iMode) = dble(freqInTHz)*THzToHartree 
        !
        read(1,*) dummyC, dummyC, dummyC, dummyC, dummyC, dummyC
        !
        do iAtom = 1, nAtoms
          !
          read (1,*) dummyD, dummyD, dummyD, phonD(1,iAtom,iMode,iq), phonD(2,iAtom,iMode,iq), phonD(3,iAtom,iMode,iq)
          !
        enddo
        !
      enddo
      !
    enddo
    !
    close(1)
    !
    flush(iostd)
    !
    return
    !
  end subroutine readPhonons
  !
  !
  subroutine computeGeneralizedDisplacements()
    !
    implicit none
    !
    integer :: iq, iMode, iAtom
    !
    allocate( genCoord(nModes) )
    !
    do iq = 1, nOfqPoints
      !
      do iMode = 1, nModes
        !
        genCoord(iMode) = 0.0_dp
        !
        do iAtom = 1, nAtoms
          !
          genCoord(iMode) = genCoord(iMode) + sqrt(1822.88833218_dp*atomM(iAtom))*sum(phonD(:,iAtom,iMode,iq)*atomD(:,iAtom))
          !
        enddo
        !
      enddo
      !
    enddo
    !
!    open(11, file='generalizedDisplacements', status='unknown')
!    !
!    write(11, '("#Mode, frequency (eV),        genCoord(Mode),     genCoord(Mode)^2")')
!    !
!    do iMode = 1, nModes
!     write(11, '(i4,1x,3E20.10E3)') iMode, phonF(iMode)*1.0e3_dp*HartreeToEv, genCoord(iMode), genCoord(iMode)*genCoord(iMode)
!    enddo
!    !
!    close(11)
    !
    deallocate( atomM, phonD, atomD )
    !
    return
    !
  end subroutine computeGeneralizedDisplacements
  !
  !
  subroutine computeVariables()
    !
    implicit none
    !
    integer :: i, j, nm, nb
    !
    real(kind = dp), allocatable :: bi(:), di(:), bk(:), dk(:)
    !
    allocate( x(nModes), Sj(nModes), coth(nModes), wby2kT(nModes) )
    !
    x = 0.0_dp
    Sj = 0.0_dp
    coth = 0.0_dp
    wby2kT = 0.0_dp
    !
    Sj(:) = 0.5_dp*phonF(:)*genCoord(:)*genCoord(:)
    wby2kT(:) = phonF(:)/(2.0_dp*kT)
    coth(:) = cosh(wby2kT(:))/sinh(wby2kT(:))
    x(:) = Sj(:)/sinh(wby2kT(:))
    !
    allocate( s2L(nModes) )
    s2L(:) = 0
    !
    do iMode = 1, nModes
      s2L(iMode) = iMode
    enddo
    !
    call arrangeLargerToSmaller()
    !
    open(11, file='modes', status='unknown')
    !
    write(11, '("#Mode, frequency (eV),        genCoord(Mode),     genCoord(Mode)^2,  Sj/sinh(wby2kT)")')
    !
    do iMode = 1, nModes
      !
    write(11,'(i4,1x,4E20.10E3)') s2L(iMode), phonF(s2L(iMode))*1.0e3_dp*HartreeToEv, & 
                                     genCoord(s2L(iMode)), genCoord(s2L(iMode))**2, x(s2L(iMode))
      !
    enddo
    !
    close(11)
    !
    deallocate( genCoord )
    !
    nb = maximumNumberOfPhonons
    allocate( besOrderNofModeM(0:nb + 1, nModes) )
    allocate( bi(0:nb + 1), di(0:nb + 1) )
    allocate( bk(0:nb + 1), dk(0:nb + 1) )
    !
    do j = 1, nModes
      !
      bi(:) = 0.0_dp
      !
      nm = nb + 1
      call iknb(nb + 1, x(j), nm, bi) ! , di, bk, dk)
      !
      do i = 0, nb + 1
        !
        besOrderNofModeM(i,j) = bi(i)
        !
      enddo
      !
      !write(6,*) j, x(j) !, (besOrderNofModeM(i,j), i = 0, 5) ! nb + 1) !, phonF(j)
      !
    enddo
    !
    return
    !
  end subroutine computeVariables
  !
  !
  subroutine arrangeLargerToSmaller()
    !
    implicit none
    !
    integer :: i, iMode
    !
    real(kind = dp), allocatable :: temp(:)
    real(kind = dp) :: tmpr
    integer :: tmpi
    !
    allocate( temp(nModes) )
    !
    temp(:) = 0.0_dp
    temp(:) = x(:)
    !
    do iMode = 1, nModes
      !
      do i = 1, nModes-1
        !
        if ( temp(i) < temp(i+1) ) then ! exp(wby2kT(i))*bessi(1,x(i)) < exp(wby2kT(i+1))*bessi(1,x(i+1)) ) then
                                        ! if ( exp(wby2kT(i))*bessi(1,x(i)) < exp(wby2kT(i+1))*bessi(1,x(i+1)) ) then
          !
!          tmpi = s2L(i)
!          s2L(i) = s2L(i+1)
!          s2L(i+1) = tmpi
          !
          tmpr = temp(i)
          temp(i) = temp(i+1)
          temp(i+1) = tmpr
          !

          tmpr = Sj(i)
          Sj(i) = Sj(i+1)
          Sj(i+1) = tmpr
          !
          tmpr = x(i)
          x(i) = x(i+1)
          x(i+1) = tmpr
          !
          tmpr = coth(i)
          coth(i) = coth(i+1)
          coth(i+1) = tmpr
          !
          tmpr = wby2kT(i)
          wby2kT(i) = wby2kT(i+1)
          wby2kT(i+1) = tmpr
          !
          tmpr = phonF(i)
          phonF(i) = phonF(i+1)
          phonF(i+1) = tmpr
          !
        endif
        !
      enddo
      !
    enddo
    !
    deallocate ( temp )
    !
    return
    !
  end subroutine arrangeLargerToSmaller
  !
  !  
  subroutine readMjs()
    !
    implicit none
    !
    integer :: i, iE0, iE, numOfMjs
    real(kind = dp) :: dummyD1, dummyD2, Ee, MjOfE, MjOfE0, eBin, DHifMin, eifMin, volume
    character(len =  1) :: dummyC1
    character(len = 32) :: dummyC32
    character(len = 35) :: dummyC35
    character(len = 256) :: modeFolder
    !
    allocate ( Mjs(modes,-nEnergies:nEnergies) )
    !
    Mjs(:,:) = 0.0_dp
    !
    do iMode = 1, modes
      !
      if ( s2L(iMode) < 10 ) then
        write(modeFolder, '("mode_", i1, "/", a)') s2L(iMode), trim(MjsInput)
      else if ( s2L(iMode) < 100 ) then
        write(modeFolder, '("mode_", i2, "/", a)') s2L(iMode), trim(MjsInput)
      else if ( s2L(iMode) < 1000 ) then
        write(modeFolder, '("mode_", i3, "/", a)') s2L(iMode), trim(MjsInput)
      else if ( s2L(iMode) < 10000 ) then
        write(modeFolder, '("mode_", i4, "/", a)') s2L(iMode), trim(MjsInput)
      endif
      !

      !
      inquire(file= trim(modeFolder), exist = file_exists)
      if ( .not.file_exists ) then
        !
        write(iostd, '(" File : ", a, " does not exist!")') trim(modeFolder)
        !
      else
        !
        write(iostd, '(" Reading file : ", a)') trim(modeFolder)
        !
        open(1, file=trim(modeFolder), status="old")
        !
        read(1, *)
        read(1, '(a32, ES24.15E3, a35)') dummyC32, volume, dummyC35
        read(1, '(a32, ES24.15E3, a35)') dummyC32, eifMin, dummyC35
        read(1, '(a32, ES24.15E3, a35)') dummyC32, DHifMin, dummyC35
        read(1, '(a32, ES24.15E3, a35)') dummyC32, eBin, dummyC35
        read(1, *)
        !
        read(1, '(i10)') numOfMjs
        !
        read(1, '(3ES24.15E3)' ) Ee, MjOfE0, dummyD1
        !
        Mjs(iMode,1) = MjOfE0
        !energy(1) = Ee
        !
        iE = int(Ee/de) + 1
        !
        do i = 2, numOfMjs
          !
          iE0 = iE ! int(energy(i-1)/deltaE) + 1 !  iE
          read(1, '(3ES24.15E3)') Ee, MjOfE, dummyD2
!          energy(i) = Ee
          iE = int(Ee/de) + 1
          !Vfis(iE0:iE) = VfiOfE0
          Mjs(iMode, i) = MjOfE
          !VfiOfE0 = VfiOfE
          !
          write(6,*) iMode, Mjs(iMode,i)
        enddo
        !
        close(1)
        !
      endif
      !
    enddo
    !
    !do iE = 0, numOfVfis ! -nEnergies, nEnergies
    !  write(44,*) energy(iE)*HartreeToEv, Vfis(iE), lsf(iE)
    !enddo
    !
    !close(44)
    !
    return
    !
  end subroutine readMjs
  !
  !
  subroutine lsfDeterministicFourPhononsByFourBands()
    !
    implicit none
    !
    integer :: ic
    integer :: iMode1, iMode2, iMode3, iMode4
    integer :: pm1, pm2, pm3, pm4
    !
    real(kind = dp) :: t1, t2
    !
    if ( myid == root ) then
      write(iostd,*) 'Four modes'
      flush(iostd)
    endif
    !
    ! Four modes
    ! 
    call cpu_time(t1)
    !
    ic = 0
    do iMode1 = iModeIs(myid), iModeFs(myid)
      do iMode2 = iMode1+1, modes - 2
        do iMode3 = iMode2+1, modes - 1
          do iMode4 = iMode3+1, modes
            !
            do pm1 = -1, 1, 2
              do pm2 = -1, 1, 2
                do pm3 = -1, 1, 2
                  do pm4 = -1, 1, 2
                    !
                    pj(:) = 0
                    pj(s2L(iMode1)) = pm1
                    pj(s2L(iMode2)) = pm2
                    pj(s2L(iMode3)) = pm3
                    pj(s2L(iMode4)) = pm4
                    !
                    call lsfOfConfigurationPj()
                    !
                    ic = ic + 1
                    !
                  enddo
                enddo
              enddo
            enddo
            !
          enddo
        enddo
      enddo
    enddo
    !
    call cpu_time(t2)
    !
    return
    !
  end subroutine lsfDeterministicFourPhononsByFourBands
  !
  !
  subroutine lsfOfConfigurationPj()
    !
    implicit none
    !
    integer :: iE, j
    !
    real(kind = dp) :: E, Fj, prodFj, sumOverj, besPj, besRatio
    !
    prodFj = 1.0_dp
    sumOverj = 0.0_dp
    do j = 1, modes
      !
      Fj = 1.0_dp
      besPj = besOrderNofModeM(abs(pj(s2L(j))), s2L(j))
      if ( pj(s2L(j)) > 0 ) then
        if ( besPj > 1.0e-15_dp ) then 
          Fj = exp(pj(s2L(j))*wby2kT(s2L(j)) - Sj(s2L(j))*coth(s2L(j)))*besPj
        else 
          Fj = 0.0_dp
        endif
      else
        Fj = exp(pj(s2L(j))*wby2kT(s2L(j)) - Sj(s2L(j))*coth(s2L(j)))*besPj
      endif
      prodFj = prodFj * Fj
      !
      besRatio = 0.5_dp*x(s2L(j))/(abs(pj(s2L(j)))+1)
      if  ( besPj > 1.0e-15_dp ) besRatio = besOrderNofModeM(abs(pj(s2L(j)))+1, s2L(j))/besPj
      sumOverj = sumOverj + (abs(pj(s2L(j))) + x(s2L(j))*besRatio)
      !
    enddo
    !
    E = sum(pj(:)*phonF(:))
    iE = 0
    if ( abs(E) > 1.0e-6_dp*evToHartree ) then
      iE = int(abs(E)/deltaE) + 1
      if ( E < 0.0_dp ) iE = -iE
    endif
    !
    iEbinsByBands(iE) = iEbinsByBands(iE) + 1
    lsfVsEbyBands(iE) = lsfVsEbyBands(iE) + prodFj*sumOverj
    !
    return
    !
  end subroutine lsfOfConfigurationPj
  !
  !
  subroutine calculatePlusMinusStates(l)
    !
    implicit none
    !
    integer, intent(in) :: l
    !
    integer :: iDes, other(0:l-1)
    !
    do iDes = 0, 2**l - 1
      !
      other(:) = 0
      !
      call decimalToOther(iDes, l, 2, other)
      !
      pms(iDes,:) = other(:)
      !
    enddo
    !
    return
    !
  end subroutine calculatePlusMinusStates
  !
  !
  subroutine distrubutePhononsInBands(m, l)
    !
    implicit none
    !
    integer, intent(in) :: m, l
    !
    integer :: i, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12
    !
    if ( l == 1 ) then
      !
      pj0s(1,1) = m
      !
    else if ( l == m - 1 ) then
      !
      do i = 1, l
        pj0s(i,:) = 1
        pj0s(i,i) = m-(l-1)
      enddo
      !
    else if ( l == m ) then
      !
      pj0s(1,:) = 1
      !
    else if ( l == 2 ) then
      !
      do i = 1, m - 1
        !
        pj0s(i, 1) = i
        pj0s(i, 2) = m-i
        !
      enddo
      !
    else if ( l == 3 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            !
            if ( i1 + i2 + i3 == m ) then
              !
              pj0s(i, 1) = i1
              pj0s(i, 2) = i2
              pj0s(i, 3) = i3
              !
              i = i + 1
              !
            endif
            !
          enddo
        enddo
      enddo
      !
      !write(6,*) 'l = 3, i = ', i - 1
      !
    else if ( l == 4 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              !
              if ( i1 + i2 + i3 + i4 == m ) then
                !
                pj0s(i, 1) = i1
                pj0s(i, 2) = i2
                pj0s(i, 3) = i3
                pj0s(i, 4) = i4
                !
                i = i + 1
                !
              endif
              !
            enddo
          enddo
        enddo
      enddo
      !write(6,*) 'l = 4, i = ', i - 1
      !
    else if ( l == 5 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                !
                if ( i1 + i2 + i3 + i4 + i5 == m ) then
                  !
                  pj0s(i, 1) = i1
                  pj0s(i, 2) = i2
                  pj0s(i, 3) = i3
                  pj0s(i, 4) = i4
                  pj0s(i, 5) = i5
                  !
                  i = i + 1
                  !
                endif
                !
              enddo
            enddo
          enddo
        enddo
      enddo
      !
      !write(6,*) 'l = 5, i = ', i - 1
      !
    else if ( l == 6 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !
      !write(6,*) 'l = 6, i = ', i - 1
      !
    else if ( l == 7 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                do i7 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 + i7 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    pj0s(i, 7) = i7
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !
      !write(6,*) 'l = 7, i = ', i
      !
    else if ( l == 8 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                do i7 = 1, m - (l - 1)
                do i8 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    pj0s(i, 7) = i7
                    pj0s(i, 8) = i8
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !
      !write(6,*) 'l = 8, i = ', i
      !
    else if ( l == 9 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                do i7 = 1, m - (l - 1)
                do i8 = 1, m - (l - 1)
                do i9 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    pj0s(i, 7) = i7
                    pj0s(i, 8) = i8
                    pj0s(i, 9) = i9
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
                enddo
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !
    else if ( l == 10 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                do i7 = 1, m - (l - 1)
                do i8 = 1, m - (l - 1)
                do i9 = 1, m - (l - 1)
                do i10 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    pj0s(i, 7) = i7
                    pj0s(i, 8) = i8
                    pj0s(i, 9) = i9
                    pj0s(i, 10) = i10
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
                enddo
                enddo
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !
    else if ( l == 11 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                do i7 = 1, m - (l - 1)
                do i8 = 1, m - (l - 1)
                do i9 = 1, m - (l - 1)
                do i10 = 1, m - (l - 1)
                do i11 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 + i11 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    pj0s(i, 7) = i7
                    pj0s(i, 8) = i8
                    pj0s(i, 9) = i9
                    pj0s(i, 10) = i10
                    pj0s(i, 11) = i11
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
                enddo
                enddo
                enddo
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      !
    else if ( l == 12 ) then
      !
      i = 1
      do i1 = 1, m - (l - 1)
        do i2 = 1, m - (l - 1)
          do i3 = 1, m - (l - 1)
            do i4 = 1, m - (l - 1)
              do i5 = 1, m - (l - 1)
                do i6 = 1, m - (l - 1)
                do i7 = 1, m - (l - 1)
                do i8 = 1, m - (l - 1)
                do i9 = 1, m - (l - 1)
                do i10 = 1, m - (l - 1)
                do i11 = 1, m - (l - 1)
                do i12 = 1, m - (l - 1)
                  !
                  if ( i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 + i11 + i12 == m ) then
                    !
                    pj0s(i, 1) = i1
                    pj0s(i, 2) = i2
                    pj0s(i, 3) = i3
                    pj0s(i, 4) = i4
                    pj0s(i, 5) = i5
                    pj0s(i, 6) = i6
                    pj0s(i, 7) = i7
                    pj0s(i, 8) = i8
                    pj0s(i, 9) = i9
                    pj0s(i, 10) = i10
                    pj0s(i, 11) = i11
                    pj0s(i, 12) = i12
                    !
                    i = i + 1
                    !
                  endif
                  !
                enddo
                enddo
                enddo
                enddo
                enddo
                enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    endif
    !
    return
    !
  end subroutine distrubutePhononsInBands
  !
  !
  subroutine lsfWithMphonons(m, l, tTimes)
    !
    implicit none
    !
    integer, intent(in) :: m, l, tTimes
    !
    integer :: ii
    !
    integer :: iMC, iM, i, pick, j, picks(l), iE, iDes, iRand, steps
    !
    real(kind = dp) :: E, Fj, prodFj, sumOverj, besPj, besRatio, randy
    !
    logical :: picked
    !
    if (myid == root) then
      write(iostd,'(i4," phonons by", i3, " bands started.")') m, l
      flush(iostd)
    endif
    !
    printSteps = int( (iModeFs(myid) - iModeIs(myid) + 1.01_dp)/10 )
    !
    do iMC = iModeIs(myid), iModeFs(myid)
      !
      if ( mod(iMC-iModeIs(myid)+1, printSteps) == 0 ) then
        if (myid == root) then
          steps = iModeFs(myid) - iModeIs(myid) + 1
          write(iostd,'(i4," phonons by", i3," bands.", i12," over ",i12," MC iters per processor done.")') m, l, iMC, steps
          flush(iostd)
        endif
      endif
      !
      picks(:) = 0
      !
      if (istat == 0) then
        !
        do iM = 1, l
          picked = .false.

 10       read(un) iRand
          iRand = mod(abs(iRand), modes) + 1
          do i = 1, iM-1
            if (picks(i) == iRand) goto 10
          enddo
          picks(iM) = iRand
          !
        enddo
        !
      else
        !
        do iM = 1, l
          picked = .false.
          !
 11       CALL RANDOM_NUMBER(randy)
          !
          pick = int( modes*randy ) + 1
          do i = 1, l
            if ( pick .eq. picks(i) ) picked = .true.
          enddo
          if ( picked ) goto 11
          picks(iM) = pick
        enddo
        !
      endif
      !
      do ii = 1, tTimes
        !
        do iDes = 0, 2**l - 1
          !
          pj(:) = 0
          !
          do iM = 1, l
            pj(s2L(picks(iM))) = pj0s(ii,iM)*(-1)**(pms(iDes,iM-1))
          enddo
          !
          if ( abs(sum(abs(pj(picks(:)))) - m) > 0 ) then
            if (myid == root) then 
              write(iostd,*) 'ERROR', m, sum(abs(pj(s2L(picks(:))))), pj(s2L(picks(:)))
              do iM = 1, l
                if (abs(pj(picks(iM))) < 1 ) then
                  write(iostd, *) 'ERROR 1', picks(iM)
                  write(iostd, *) 'ERROR 2', pj(s2L(picks(iM)))
                  flush(iostd)
                endif
              enddo
            endif
          endif
          !
          prodFj = 1.0_dp
          sumOverj = 0.0_dp
          !
          do j = 1, modes
            !
            Fj = 1.0_dp
            besPj = besOrderNofModeM(abs(pj(s2L(j))), s2L(j))
            if ( pj(s2L(j)) > 0 ) then
              if ( besPj > 1.0e-15_dp ) then
                Fj = exp(pj(s2L(j))*wby2kT(s2L(j)) - Sj(s2L(j))*coth(s2L(j)))*besPj
              else
                Fj = 0.0_dp
              endif
            else
              Fj = exp(pj(s2L(j))*wby2kT(s2L(j)) - Sj(s2L(j))*coth(s2L(j)))*besPj
            endif
            !
            prodFj = prodFj * Fj
            !
            besRatio = 0.5_dp*x(s2L(j))/(abs(pj(s2L(j)))+1)
            if  ( besPj > 1.0e-15_dp ) besRatio = besOrderNofModeM(abs(pj(s2L(j)))+1, s2L(j))/besPj
            !
            sumOverj = sumOverj + (abs(pj(s2L(j))) + x(s2L(j))*besRatio)
            !
          enddo
          !
          E = sum(pj(:)*phonF(:))
          !
          iE = 0
          if ( abs(E) > 1.0e-6_dp*evToHartree ) then
            iE = int(abs(E)/deltaE) + 1
            if ( E < 0.0_dp ) iE = -iE
          endif
          !
          iEbinsByBands(iE) = iEbinsByBands(iE) + 1
          lsfVsEbyBands(iE) = lsfVsEbyBands(iE) + prodFj*sumOverj
          !
        enddo
        !
      enddo
      !
    enddo
    !
    if (myid == root) then
      write(iostd,'("---------------------------------------------")')
      write(iostd,'(i4," phonons by", i3, " bands done.")') m, l
      flush(iostd)
    endif
    !
    return
    !
  end subroutine lsfWithMphonons
  !
  !
  subroutine decimalToOther(iDec, n, iBase, other)
    !
    implicit none
    !
    integer, intent(in) :: n, iBase
    integer :: iDec, m
    integer :: other(0:n-1), j
    !
    m = iDec
    do j = n-1, 1, -1
      other(j) = int(m/(iBase**j))
      m = mod(iDec, iBase**j)
    enddo
    other(0) = mod(m,iBase)
    !
    return
    !
  end subroutine decimalToOther
  !
  !
  subroutine calculateDE(maxM, iEbins, de)
    !
    implicit none
    !
    integer, intent(in) :: maxM, iEbins(-nEnergies: nEnergies)
    real(dp), intent(out) :: de
    !
    integer :: iE, j, ic, ib, iEmMax, nSteps, jMax, iEstep
    !
    integer, allocatable :: tmpB(:), iEsteps(:)
    !
    logical :: empty
    !
    allocate ( tmpB(nEnergies) )
    !
    ic = 1
    do j = 1, nEnergies
      if (mod(nEnergies,j) == 0) then
        tmpB(ic) = int((dble(nEnergies) + 1.e-8_dp)/j)
        ic = ic + 1
      endif
    enddo
    !
    nSteps = ic - 1
    allocate( iEsteps(nSteps) )
    iEsteps(:) = tmpB(nSteps:1:-1)
    deallocate(tmpB)
    !
    iEmMax = int(maxM*maxval(phonF(:))/deltaE) + 1
    !
    j = 1
    do while ( ( iEmMax > iEsteps(j) ) .and. (j < nSteps) ) 
      j = j + 1
    enddo
    !
    jMax = j - 1
    if (jMax > nSteps) jMax = nSteps
    !
    empty = .true.
    j = jMax
    do while ( ( empty .eqv. .true. ) .and. (j > 1) )
      !
      empty = .true.
      iEstep = iEsteps(j)
      do iE = 1, iEmMax-1, iEstep
        ib = sum( iEbins(iE:iE+iEstep-1) )
        if ( ib < 1 ) then
          empty = .false.
        endif
      enddo
      j = j - 1
      !
    enddo
    !
    j = j + 2
    !
    iEstep = iEsteps(j)
    de = deltaE*real(iEstep, dp)
    !
    deallocate( iEsteps )
    !
    return
    !
  end subroutine calculateDE
  ! 
  !
  subroutine lsfMbyOneBand(m)
    !
    implicit none
    !
    integer, intent(in) :: m
    !
    integer :: iMode1, pm1
    !
    real(dp) :: t1, t2
    !
    call cpu_time(t1)
    !
    do iMode1 = 1, modes
      !
      do pm1 = -m, m, 2*m
        !
        pj(:) = 0
        pj(s2L(iMode1)) = pm1
        !
        call lsfOfConfigurationPj()
        !
      enddo
      !
    enddo
    !
    call cpu_time(t2)
    !
    write(iostd,'(" LSF of: ", i4, " phonons using one band done in: ", f10.2, " secs.")') m, t2-t1
    flush(iostd)
    !
    return
    !
  end subroutine lsfMbyOneBand
  !
  !
  subroutine lsfMbyTwoBands(m)
    !
    implicit none
    !
    integer, intent(in) :: m
    !
    integer :: iMode1, iMode2, pm1, pm2, l
    !
    real(dp) :: t1, t2
    !
    call cpu_time(t1)
    !
    do l = 1, m - 1
      !
      do iMode1 = 1, modes - 1
        do iMode2 = iMode1+1, modes
          !
          do pm1 = -l, l, 2*l
            do pm2 = -(m-l), (m-l), 2*(m-l)
              !
              pj(:) = 0
              pj(s2L(iMode1)) = pm1
              pj(s2L(iMode2)) = pm2
              !
              call lsfOfConfigurationPj()
              !
            enddo
          enddo
          !
        enddo
      enddo
      !
    enddo
    !
    call cpu_time(t2)
    !
    if ( myid == root ) then
      write(iostd,'(" LSF of: ", i4, " phonons using two bands done in: ", f10.2, " secs.")') m, t2-t1
      flush(iostd)
    endif
    !
    return
    !
  end subroutine lsfMbyTwoBands
  !
  !
  subroutine lsfMbyThreeBands(m)
    !
    implicit none
    !
    integer, intent(in) :: m
    !
    real(dp) :: t1, t2, times3
    integer :: iMode1, iMode2, iMode3, ni, mi, iDes, ii
    !
    call cpu_time(t1)
    !
    times3 = 1.0_dp
    mi = 2
    do ni = m - 1, m - 3 + 1, -1
      times3 = times3*dble(ni)/dble(mi)
      mi = mi - 1
    enddo
    !
    allocate( pj0s(int(times3 + 1.e-3_dp), 3) )
    pj0s(:,:) = 0
    !
    call distrubutePhononsInBands(m, 3)
    !
    allocate( pms( 0:2**3-1, 0:3-1 ) )
    pms(:,:) = 0
    !
    call calculatePlusMinusStates( 3 )
    !
    do iMode1 = iModeIs(myid), iModeFs(myid)
      do iMode2 = iMode1+1, modes - 1
        do iMode3 = iMode2+1, modes
          !
          do ii = 1, int(times3 + 1.e-3_dp)
            !
            do iDes = 0, 2**3 - 1
              !
              pj(:) = 0
              !
              pj(s2L(iMode1)) = pj0s(ii,1)*(-1)**(pms(iDes,1-1))
              pj(s2L(iMode2)) = pj0s(ii,2)*(-1)**(pms(iDes,2-1))
              pj(s2L(iMode3)) = pj0s(ii,3)*(-1)**(pms(iDes,3-1))
              !
              call lsfOfConfigurationPj()
              !
            enddo
            !
          enddo
          !
        enddo
      enddo
    enddo
    !
    deallocate ( pj0s, pms )
    !    
    call cpu_time(t2)
    !
    if ( myid == root ) then
      write(iostd,'(" LSF of: ", i4, " phonons using three bands done in: ", f10.2, " secs.")') m, t2-t1
      flush(iostd) 
    endif
    !
    return
    !
  end subroutine lsfMbyThreeBands
  !
  !
  subroutine writeLSFandCrossSection()
    !
    implicit none
    ! 
    integer :: iE
    real(kind = dp) :: E!, vg
    !
    open(1, file='lsfVsE', status='unknown')
    !
    do iE = -nEnergies, nEnergies
      !
      E = real(iE, dp)*deltaE 
      !
      !vg = 1.0_dp
      !if (E > 0.0_dp) vg = sqrt(2.0_dp*E)
      !
      !write(1,'(F16.8,2E18.6e3)') E*HartreeToEv, lsfVsE(iE), twoPi*abCM**2*volume*Vfis(iE)*lsfVsE(iE)/vg
      write(1,'(F16.8,E18.6e3)') E*HartreeToEv, lsfVsE(iE)!, twoPi*abCM**2*volume*Vfis(iE)*lsfVsE(iE)/vg
      !
    enddo
    !
    close(1)
    !
    return
    !
  end subroutine writeLSFandCrossSection
  !
  !
  subroutine init_random_seed()
    !               
    implicit none
    !
    integer(kind = int32), allocatable :: seed(:)
    integer(kind = int32) :: n !, i, n, dt(8), pid
    integer :: t
    !
    call random_seed(size = n)
    !
    allocate(seed(n))
    !
    ! Fallback to XOR:ing the current time and pid. The PID is
    ! useful in case one launches multiple instances of the same
    ! program in parallel.
    !
    call system_clock(t)
    !
    seed = 5347
    !
!    if (t == 0) then
!      call date_and_time(values=dt)
!      t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
!          + dt(2) * 31 * 24 * 60 * 60 * 1000 &
!          + dt(3) * 24 * 60 * 60 * 1000 &
!          + dt(5) * 60 * 60 * 1000 &
!          + dt(6) * 60 * 1000 + dt(7) * 1000 &
!          + dt(8)
!    end if
!    pid = getpid()
!    t = ieor(t, int(pid, kind(t)))
!    do i = 1, n
!      seed(i) = lcg(t)
!    end do
!    !
    call random_seed(put=seed)
    !
  end subroutine init_random_seed
  !
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  !
  integer function lcg(s)
    !
    integer :: s
    !
    if (s == 0) then
      s = 104729
    else
      !s = mod(s, 4294967296)
      s = mod(s, 4294967)
    end if
    !
    !s = mod(s * 279470273, 4294967291)
    s = mod(s * 279470273, 4294967)
    !
    lcg = int(mod(s, huge(0)), kind(0))
    !
  end function lcg
  !
  !

      subroutine iknb ( n, x, nm, bi) !, di, bk, dk )
!
! Modified : when x < 10^(-15) return the limiting value for small argument [ I_n(x) ~ (x/2)^n Gamma(n+1) ]
! 
!c*********************************************************************72
!c
!cc IKNB compute Bessel function In(x) and Kn(x).
!c
!c  Discussion:
!c
!c    Compute modified Bessel functions In(x) and Kn(x),
!c    and their derivatives.
!c
!c  Licensing:
!c
!c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!c    they give permission to incorporate this routine into a user program 
!c    provided that the copyright is acknowledged.
!c
!c  Modified:
!c
!c    17 July 2012
!c
!c  Author:
!c
!c    Shanjie Zhang, Jianming Jin
!c
!c  Reference:
!c
!c    Shanjie Zhang, Jianming Jin,
!c    Computation of Special Functions,
!c    Wiley, 1996,
!c    ISBN: 0-471-11963-6,
!c    LC: QA351.C45.
!c
!c  Parameters:
!c
!c    Input, integer N, the order of In(x) and Kn(x).
!c
!c    Input, double precision X, the argument.
!c
!c    Output, integer NM, the highest order computed.
!c
!c    Output, double precision BI(0:N), DI(0:N), BK(0:N), DK(0:N),
!c    the values of In(x), In'(x), Kn(x), Kn'(x).
!c
      implicit none

      integer, intent(in) :: n

!      double precision :: a0
      double precision :: bi(0:n)
!      double precision :: bkl
      double precision :: bs
      double precision :: el
      double precision :: f
      double precision :: f0
      double precision :: f1
!      double precision :: g
!      double precision :: g0
!      double precision :: g1
      integer :: k
!      integer :: k0
!      integer :: l
      integer :: m, ik
!      integer :: msta1
!      integer :: msta2
      integer :: nm
      double precision :: pi
!      double precision :: r
      double precision :: s0
      double precision :: sk0
!      double precision :: vt
      double precision :: x, ifact

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      nm = n

      if ( x .le. 1.0D-15 ) then
        do k = 0, n
          ifact = 1.0_dp
          do ik = 2, k
            ifact = ifact*ik
          enddo
          bi(k) = (0.5_dp*x)**k/ifact
        end do
        return
      end if

      if ( n .eq. 0 ) then
        nm = 1
      end if

      m = msta1 ( x, 200 )
      if ( m .lt. nm ) then
        nm = m
      else
        m = msta2 ( x, nm, 15 )
      end if

      bs = 0.0D+00
      sk0 = 0.0D+00
      f0 = 0.0D+00
      f1 = 1.0D-100
      do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 + f0
        if ( k .le. nm ) then
          bi(k) = f
        end if
        if ( k .ne. 0 .and. k .eq. 2 * int ( k / 2 ) ) then
          sk0 = sk0 + 4.0D+00 * f / k
        end if
        bs = bs + 2.0D+00 * f
        f0 = f1
        f1 = f
      end do

      s0 = exp ( x ) / ( bs - f )
      do k = 0, nm
        bi(k) = s0 * bi(k)
      end do

      return
      end SUBROUTINE iknb
  !
  !

  SUBROUTINE iknb2(n,x,nm,bi,di,bk,dk) 
    !
    !    ============================================================ 
    !    Purpose: Compute modified Bessel functions In(x) and Kn(x), 
    !             and their derivatives 
    !    Input:   x --- Argument of In(x) and Kn(x) ( 0 ó x ó 700 ) 
    !             n --- Order of In(x) and Kn(x) 
    !    Output:  BI(n) --- In(x) 
    !             DI(n) --- In'(x) 
    !             BK(n) --- Kn(x) 
    !             DK(n) --- Kn'(x) 
    !             NM --- Highest order computed 
    !    Routines called: 
    !             MSTA1 and MSTA2 for computing the starting point 
    !             for backward recurrence 
    !    =========================================================== 
    !   
    INTEGER, INTENT(IN)     :: n 
    REAL (dp), INTENT(IN)   :: x 
    INTEGER, INTENT(OUT)    :: nm 
    REAL (dp), INTENT(OUT)  :: bi(0:n) 
    REAL (dp), INTENT(OUT)  :: di(0:n) 
    REAL (dp), INTENT(OUT)  :: bk(0:n) 
    REAL (dp), INTENT(OUT)  :: dk(0:n) 
    ! 
    REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.5772156649015329_dp 
    REAL (dp)  :: a0, bkl, bs, f, f0, f1, g, g0, g1, r, s0, sk0, vt 
    INTEGER    :: k, k0, l, m 
    ! 
    nm = n 
    IF (x <= 1.0D-50) THEN 
      DO  k = 0, n 
        bi(k) = 0.0D0 
        di(k) = 0.0D0 
        bk(k) = 1.0D+300 
        dk(k) = -1.0D+300 
      END DO 
      bi(0) = 1.0D0 
      di(1) = 0.5D0 
      RETURN 
    END IF 
    IF (n == 0) nm = 1 
    m = msta1(x, 200)
    IF (m < nm) THEN 
      nm = m 
    ELSE 
      m = msta2(x, nm, 15)
    END IF 
    !write(6,*)'mmmmmmmmm', m
    bs = 0.0D0 
    sk0 = 0.0D0 
    f0 = 0.0D0 
    f1 = 1.0D-100 
    DO  k = m, 0, -1 
      f = 2*(k+1)/x*f1 + f0 
      IF (k <= nm) bi(k) = f 
      IF (k /= 0 .AND. k == 2*INT(k/2)) sk0 = sk0 + 4.0D0 * f / k 
      bs = bs + 2.0D0 * f 
      f0 = f1 
      f1 = f 
    END DO 
    !s0 = EXP(x) / (bs-f) 
    !write(6,*) f, f1
    s0 = EXP(x) / (bs-f1) 
    bi(0:nm) = s0 * bi(0:nm) 
    IF (x <= 8.0D0) THEN 
      bk(0) = -(LOG(0.5D0*x)+el) * bi(0) + s0 * sk0 
      bk(1) = (1.0D0/x-bi(1)*bk(0)) / bi(0) 
    ELSE 
      a0 = SQRT(pi/(2.0D0*x)) * EXP(-x) 
      k0 = 16 
      IF (x >= 25.0) k0 = 10 
      IF (x >= 80.0) k0 = 8 
      IF (x >= 200.0) k0 = 6 
      DO  l = 0, 1 
        bkl = 1.0D0 
        vt = 4 * l 
        r = 1.0D0 
        DO  k = 1, k0 
          r = 0.125D0 * r * (vt - (2*k-1)**2) / (k*x) 
          bkl = bkl + r 
        END DO 
        bk(l) = a0 * bkl 
      END DO 
    END IF 
    g0 = bk(0) 
    g1 = bk(1) 
    DO  k = 2, nm 
      g = 2*(k-1)/x*g1 + g0 
      bk(k) = g 
      g0 = g1 
      g1 = g 
    END DO 
    di(0) = bi(1) 
    dk(0) = -bk(1) 
    DO  k = 1, nm 
      di(k) = bi(k-1) - k / x * bi(k) 
      dk(k) = -bk(k-1) - k / x * bk(k) 
    END DO 
    RETURN 
    !
  END SUBROUTINE iknb2 
  !
  !
  FUNCTION msta1(x, mp) RESULT(fn_val) 
    !
    !       =================================================== 
    !       Purpose: Determine the starting point for backward 
    !                recurrence such that the magnitude of 
    !                Jn(x) at that point is about 10^(-MP) 
    !       Input :  x     --- Argument of Jn(x) 
    !                MP    --- Value of magnitude 
    !       Output:  MSTA1 --- Starting point 
    !       =================================================== 
    ! 
    REAL (dp), INTENT(IN)  :: x 
    INTEGER, INTENT(IN)    :: mp 
    INTEGER                :: fn_val 
    ! 
    REAL (dp)  :: a0, f, f0, f1 
    INTEGER    :: it, n0, n1, nn 
    ! 
    a0 = ABS(x) 
    n0 = INT(1.1*a0) + 1 
    f0 = envj(n0,a0) - mp 
    n1 = n0 + 5 
    f1 = envj(n1,a0) - mp 
    DO  it = 1, 20 
      nn = n1 - int((n1-n0)/(1.0_dp - f0/f1))
      f = envj(nn,a0) - mp 
      IF (ABS(nn-n1) < 1) EXIT 
      n0 = n1 
      f0 = f1 
      n1 = nn 
      f1 = f 
    END DO 
    ! 
    fn_val = nn 
    !
    RETURN
    !
  END FUNCTION msta1 
  !
  ! 
  FUNCTION msta2(x, n, mp) RESULT(fn_val) 
    !
    !       =================================================== 
    !       Purpose: Determine the starting point for backward 
    !                recurrence such that all Jn(x) has MP 
    !                significant digits 
    !       Input :  x  --- Argument of Jn(x) 
    !                n  --- Order of Jn(x) 
    !                MP --- Significant digit 
    !       Output:  MSTA2 --- Starting point 
    !       =================================================== 
    !
    REAL (dp), INTENT(IN)  :: x 
    INTEGER, INTENT(IN)    :: n 
    INTEGER, INTENT(IN)    :: mp 
    INTEGER                :: fn_val 
    ! 
    REAL (dp)  :: a0, ejn, f, f0, f1, hmp, obj 
    INTEGER    :: it, n0, n1, nn 
    ! 
    a0 = ABS(x) 
    hmp = 0.5_dp * mp 
    ejn = envj(n, a0) 
    IF (ejn <= hmp) THEN 
      obj = mp 
      n0 = INT(1.1*a0) 
    ELSE 
      obj = hmp + ejn 
      n0 = n 
    END IF 
    !!!!!!!!
    if ( n0 < 1 ) n0 = 1
    !!!!!!!!
    f0 = envj(n0,a0) - obj 
    n1 = n0 + 5 
    f1 = envj(n1,a0) - obj 
    !
    DO  it = 1, 20 
      nn = n1 - int((n1-n0)/(1.0_dp - f0/f1))
      f = envj(nn, a0) - obj 
      IF (ABS(nn-n1) < 1) EXIT 
      n0 = n1 
      f0 = f1 
      n1 = nn 
      f1 = f 
    END DO 
    !
    fn_val = nn + 10 
    !
    RETURN
    ! 
  END FUNCTION msta2 
  !
  ! 
  FUNCTION envj(n, x) RESULT(fn_val) 
    ! 
    INTEGER, INTENT(IN)    :: n 
    REAL (dp), INTENT(IN)  :: x 
    REAL (dp)              :: fn_val 
    !
    fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n) 
    !
    RETURN 
    !
  END FUNCTION envj 
  !
  !
  subroutine parallelIsFsBy3()
    !
    implicit none
    !
    real(dp) :: parTotal, parTotal2, totalStates, states, averageStatesPerProc
    integer :: i, iState, iproc
    !
    totalStates = modes/6.0_dp
    totalStates = totalStates*(modes-1)*(modes-2)
    !
    if ( modes > numprocs ) then
      !
      nProcMax = numprocs
      !
      parTotal = 0
      iproc = 0
      iModeIs(iproc) = 1
      iState = modes
      !
      do while ( (iState > 1) .and. (nProcMax - iproc > 1) )
        !
        states = real(iState-1, dp)*real(iState-2, dp)/2.0_dp
        parTotal = parTotal + states
        parTotal2 = parTotal + real(iState-2, dp)*real(iState-3, dp)/2.0_dp 
        !
        averageStatesPerProc = totalStates/(nProcMax - iproc) - 0.01_dp
        if ( ( parTotal > averageStatesPerProc ) .or. ( parTotal2 > averageStatesPerProc ) ) then
          iModeFs(iproc) = modes - iState + 1
          iproc = iproc + 1
          iModeIs(iproc) = modes - iState + 2
          totalStates = totalStates - parTotal
          parTotal = 0
        endif
        !
        iState = iState - 1
        !
      enddo
      !
      iModeFs(iproc) = modes - 2
      !
    else
      !
      nProcMax = modes - 2
      !
      do i = 0, nProcMax - 1
        iModeIs(i) = i+1
        iModeFs(i) = i+1
      enddo
    endif
    !
    write(6,*) iModeIs(:), iModeFs(:)
    !
    return 
    !
  end subroutine parallelIsFsBy3
  !
  !
  subroutine parallelIsFsBy4()
    !
    implicit none
    !
    real(dp) :: parTotal, parTotal2, totalStates, states, averageStatesPerProc
    integer :: i, iState, iproc
    !
    totalStates = modes/24.0_dp
    totalStates = totalStates*(modes-1)*(modes-2)*(modes-3)
    !
    if ( modes > numprocs ) then
      !
      nProcMax = numprocs
      !
      parTotal = 0
      iproc = 0
      iModeIs(iproc) = 1
      iState = modes
      !
      do while ( (iState > 1) .and. (nProcMax - iproc > 1) )
        !
        states = real(iState-1, dp)*real(iState-2, dp)*real(iState-3, dp)/6.0_dp
        parTotal = parTotal + states
        parTotal2 = parTotal + real(iState-2, dp)*real(iState-3, dp)*real(iState-4, dp)/6.0_dp 
        !
        averageStatesPerProc = totalStates/(nProcMax - iproc) - 0.01_dp
        if ( ( parTotal > averageStatesPerProc ) .or. ( parTotal2 > averageStatesPerProc ) ) then
          iModeFs(iproc) = modes - iState + 1
          iproc = iproc + 1
          iModeIs(iproc) = modes - iState + 2
          totalStates = totalStates - parTotal
          parTotal = 0
        endif
        !
        iState = iState - 1
        !
      enddo
      !
      iModeFs(iproc) = modes - 2
      !
    else
      !
      nProcMax = modes - 2
      !
      do i = 0, nProcMax - 1
        iModeIs(i) = i+1
        iModeFs(i) = i+1
      enddo
    endif
    !
    return
    !
  end subroutine parallelIsFsBy4
  !
  !
end module lsf
