
module lsf
  !
  use constants
  !
  implicit none
  !
  integer, parameter :: int32 = selected_int_kind(5)
  integer, parameter :: int64 = selected_int_kind(15)
  integer, parameter :: un = 3
  integer, parameter :: root = 0
  !
  character(len = 6), parameter :: output = 'status'
  !
  !
  integer(kind = int32) :: ierr
  integer(kind = int32) :: ios
  integer(kind = int32) :: istat
  integer(kind = int32) :: myid
  integer(kind = int32) :: numprocs
  integer :: i
  integer :: iE
  integer :: iMint
  integer :: iMmod
  integer :: iMode
  integer :: l
  integer :: m
  integer :: maximumNumberOfPhonons
  integer :: mi
  integer :: minimumNumberOfPhonons
  integer :: nAtoms
  integer :: nEnergies
  integer :: ni
  integer :: nMC
  integer :: nModes
  integer :: nOfqPoints
  integer :: nProcMax
  integer :: printsteps
  !
  real(kind = dp) :: de
  real(kind = dp) :: deltaE
  real(kind = dp) :: E 
  real(kind = dp) :: kT
  real(kind = dp) :: maxEnergy
  real(kind = dp) :: t1
  real(kind = dp) :: t2
  real(kind = dp) :: tf
  real(kind = dp) :: ti
  real(kind = dp) :: temperature
  real(kind = dp) :: times
  !real(kind = dp) :: vg
  !real(kind = dp) :: volume
  real(kind = dp) :: weight
  !
  character(len = 256) :: continueLSFfromFile
  !character(len = 256) :: crossSectionOutput
  character(len = 256) :: fn
  character(len = 256) :: phononsInputFormat
  character(len = 256) :: phononsInput
  !character(len = 256) :: VfisInput
  !
  logical :: file_exists
  !
  !
  integer, allocatable :: iEbinsByBands(:)
  integer, allocatable :: iEbinsByPhonons(:)
  integer, allocatable :: iModeFs(:)
  integer, allocatable :: iModeIs(:)
  integer, allocatable :: pj(:)
  integer, allocatable :: pj0s(:,:)
  integer, allocatable :: pms(:,:)
  integer, allocatable :: s2L(:)
  !
  real(kind = dp), allocatable :: atomD(:,:)
  real(kind = dp), allocatable :: atomM(:)
  real(kind = dp), allocatable :: besOrderNofModeM(:,:)
  real(kind = dp), allocatable :: coth(:)
  real(kind = dp), allocatable :: genCoord(:)
  real(kind = dp), allocatable :: lsfbyPhononsPerProc(:)
  real(kind = dp), allocatable :: lsfVsE(:)
  real(kind = dp), allocatable :: lsfVsEbyBands(:)
  real(kind = dp), allocatable :: lsfVsEbyPhonons(:)
  real(kind = dp), allocatable :: phonD(:,:,:,:)
  real(kind = dp), allocatable :: phonF(:)
  real(kind = dp), allocatable :: phonQ(:,:)
  real(kind = dp), allocatable :: Sj(:)
  !real(kind = dp), allocatable :: Vfis(:)
  real(kind = dp), allocatable :: wby2kT(:)
  real(kind = dp), allocatable :: x(:)
  !
!  namelist /elphscat/ VfisInput, PhononsInput, temperature, maxEnergy, continueLSFfromFile, volume, &
  namelist /lsfInput/ phononsInput, phononsInputFormat, temperature, &
                      continueLSFfromFile, maximumNumberOfPhonons, nMC
  !
  !
contains
  !
  !
  subroutine readInputs()
    !! Read input parameters and read phonon output
    !!
    !! <h2>Walkthrough</h2>
    !!
    use readInputFiles
      !! * Include `readInputFiles` module for reading 
      !! phonons input
    !
    implicit none
    !
    !> * Check if output file exists; if it does delete it
    inquire(file = output, exist = file_exists)
    if ( file_exists ) then
      open (unit = 11, file = output, status = "old")
      close(unit = 11, status = "delete")
    endif
    !
    open (iostd, file = output, status='new')
      !! * Open new output file
    !
    call initialize()
      !! * Set default values of input parameters
    !
    READ (5, lsfInput, iostat = ios)
      !! * Read input parameters
    !
    call checkAndUpdateInput()
      !! * Check if input parameters were updated and do some basic checks
    !
    !> * Read the phonons output from QE or VASP
    if ( trim(phononsInputFormat) == 'VASP' ) then
      !
      call readPhonons(phononsInput, nOfqPoints, nAtoms, nModes, atomD, atomM, phonQ, phonF, phonD)
      !
    else if ( trim(phononsInputFormat) == 'QE' ) then
      !
      call readPhononsQE()
      !
    else 
      !
      write(iostd, '(" Unknown phonons input format : ", (a) )') trim(phononsInputFormat)
      write(iostd, '(" Phonons input format implemened are : ''VASP'' and ''QE''")')
      write(iostd, '(" Program stops!")')
      stop
      !
    endif
    !
!    call readVfis()
    !
    return
    !
  end subroutine readInputs
  !
  !
  subroutine initializeLSF()
    !! Allocate and initialize `lsfVsE` and 
    !! `minimumNumberOfPhonons`. If a file was
    !! given to continue from, read in both 
    !! variables from the file.
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
    if (trim(continueLSFfromFile) /= '') then
      !
      inquire(file = trim(continueLSFfromFile), exist = file_exists)
      if ( file_exists ) then
        !
        open(unit = 11, file = trim(continueLSFfromFile), status = "old")
        !
        read(11,'(a1, i10, a9, i5, a8)') dummyC1, dummyI, dummyC9, minimumNumberOfPhonons, dummyC8
        ! 
        minimumNumberOfPhonons = minimumNumberOfPhonons + 1
          !! @todo Figure out why increase `minimumNumberOfPhonons` by 1
        !
        write(iostd, '(" Minimum number of phonons : ", i5)') minimumNumberOfPhonons
        !
        do iE = -nEnergies, nEnergies
          !
          read(11,*) dummyD, lsfVsE(iE)
          !
        enddo
        !
        close(11)
        !
      endif
      !
    endif
    !
    return
    !
  end subroutine initializeLSF
  !
  !
  subroutine initialize()
    !! Set default values for input parameters
    !
    implicit none
    !
    !VfisInput = ''
    phononsInput = ''
    phononsInputFormat = ''
    temperature = -1.0_dp
    minimumNumberOfPhonons =  1
    maximumNumberOfPhonons = -1
    nMC = -1
    !
    return
    !
  end subroutine initialize
  !
  !
  subroutine checkAndUpdateInput()
    !! Check that the input variables don't still have their default
    !! values. The program will abort here if:
    !!   * `phononsInput` is undefined
    !!   * `phononsInputFormat` is undefined
    !!   * `temperature` is undefined
    !!   * `maximumNumberOfPhonons` is undefined
    !!   * number of Monte Carlo steps (`nMc`) is not
    !!     set and `maximumNumberOfPhonons` is greater
    !!     than 4
    !!
    !! This subroutine also sets `maxEnergy` and `nEnergies`
    !! and calculates `kT` and `deltaE`.
    !
    implicit none
    !
    logical :: abortExecution = .false.
    !
    if ( trim(phononsInput) == '' ) then
      write(iostd, '(" PhononsInput is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Phonons input : ", a)') trim(phononsInput)
    endif
    !
    if ( trim(phononsInputFormat) == '' ) then
      write(iostd, '(" PhononsInputFormat is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Phonons input format : ", a)') trim(phononsInputFormat)
    endif
    !
    if ( temperature < 0.0_dp ) then
      write(iostd, '(" Variable temperature has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Tempetature : ", f10.2, " Kelvin.")') temperature
      kT = temperature*8.6173324e-5_dp*eVToHartree
    endif
    !
    if ( maximumNumberOfPhonons < 0 ) then
      write(iostd, '(" Variable maximumNumberOfPhonons has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Maximum number of phonons : ", i5)') maximumNumberOfPhonons
    endif
    !
    if ( nMC < 0 ) then 
      if ( maximumNumberOfPhonons > 4 ) then
        write(iostd, '(" For calculations with configurations with more than 4 phonon modes ")')
        write(iostd, '(" the number of Monte Carlo steps ''nMC'' must be set.")')
        abortExecution = .true.
      endif
    else
      if ( maximumNumberOfPhonons > 4 ) then
        write(iostd, '(" Number of Monte Carlo steps : ", i15)') nMC
      else 
        write(iostd, '(" The number of Monte Carlo steps ''nMC'' is set to : ", i15, " but")') nMC
        write(iostd, '(" will not be used since the Monte Carlo sheme is used for calculations")')
        write(iostd, '(" with configurations with more that 4 phonon modes.")')
      endif
    endif
    !
    if ( abortExecution ) then
      write(iostd, '(" *************************** ")')
      write(iostd, '(" * Program stops!          * ")')
      write(iostd, '(" * Please check the input. * ")')
      write(iostd, '(" *************************** ")')
      stop
    endif
    !
    maxEnergy = 10.0_dp*eVToHartree
    !
    nEnergies = 5040 ! 2520 ! 10080
    deltaE = maxEnergy/real(nEnergies, dp)
    !
    write(iostd, '(" The resolution in energy is :", f10.2, " meV.")') deltaE*1000.0_dp*HartreeToEv
    !
    flush(iostd)
    !
    return
    !
  end subroutine checkAndUpdateInput
  !
  !
  subroutine readPhononsQE()
    !
    implicit none
    !
    integer :: iAtom, iMode, iq
    real(kind = dp) :: dummyD, freqInTHz
    !
    CHARACTER :: dummyC
    !
    !write(6,*) trim(phononsInput)
    open(1, file=trim(phononsInput), status="old")
    !
    read(1,*) nOfqPoints, nAtoms, nModes
    !
    write(iostd, '(" Number of atoms : ", i5)') nAtoms
    write(iostd, '(" Number of q-Points : ", i5)') nOfqPoints
    write(iostd, '(" Number of modes : ", i5)') nModes
    flush(iostd)
    !
    read (1,*)
    !
    allocate( atomD(3,nAtoms), atomM(nAtoms) )
    !
    atomD = 0.0_dp
    atomM = 0.0_dp
    !The unit is in Bohr
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
      read(1,*)
      !
      do iMode = 1, nModes
        !
        read(1,*) dummyC, dummyC, dummyC, dummyC, freqInTHz
        phonF(iMode) = dble(freqInTHz)*THzToHartree
        !
        do iAtom = 1, nAtoms
          !
          read(1,*) dummyC, phonD(1,iAtom,iMode,iq), dummyD, phonD(2,iAtom,iMode,iq), dummyD, phonD(3,iAtom,iMode,iq), dummyC
          write(6,*) phonD(1,iAtom,iMode,iq),phonD(2,iAtom,iMode,iq),phonD(3,iAtom,iMode,iq)
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
  end subroutine readPhononsQE
  !
  !  
!  subroutine readVfis()
!    !
!    implicit none
!    !
!    integer :: i, iE0, iE, dummyI, nEVfi
!    real(kind = dp) :: dummyD, E, VfiOfE, VfiOfE0
!    character :: dummyC
!    !
!    open(1, file=trim(VfisInput), status="old")
!    !
!    read(1, *) dummyC, nEVfi
!    !
!    allocate ( Vfis(-nEnergies:nEnergies) )
!    !
!    Vfis = 0.0_dp
!    !
!    read(1, '(d22.14,i5,4d22.14)' ) E, dummyI, dummyD, VfiOfE0, dummyD
!    !
!    E = E*eVToHartree 
!    iE = int(E/deltaE) + 1
!    !
!    do i = 1, nEVfi - 1
!      !
!      iE0 = iE
!      read(1, '(d22.14,i5,4d22.14)' ) E, dummyI, dummyD, VfiOfE, dummyD
!      E = E*eVToHartree
!      iE = int(E/deltaE) + 1
!      Vfis(iE0:iE) = VfiOfE0
!      VfiOfE0 = VfiOfE
!      !
!    enddo
!    !
!    close(1)
!    !
!    !do iE = -nEnergies, nEnergies
!    !  write(44,*) real(iE, dp)*deltaE*HartreeToEv, Vfis(iE)
!    !enddo
!    !
!    return
!    !
!  end subroutine readVfis
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
      do iMode2 = iMode1+1, nModes - 2
        !write(iostd, '("myID ", i15, " iMode1 = ", i5, " / ", i5, " in each proc. iMode2 =", i5, " / ", i5)') myid, & 
        !           iMode1-iModeIs(myid)+1, iModeFs(myid)-iModeIs(myid)+1,iMode2 - iMode1+1 + 1, nModes - 2 - iMode1+1 + 1
        do iMode3 = iMode2+1, nModes - 1
          do iMode4 = iMode3+1, nModes
            !
            do pm1 = -1, 1, 2
              do pm2 = -1, 1, 2
                do pm3 = -1, 1, 2
                  do pm4 = -1, 1, 2
                    !
                    pj(:) = 0
                    pj(iMode1) = pm1
                    pj(iMode2) = pm2
                    pj(iMode3) = pm3
                    pj(iMode4) = pm4
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
    !! Calculate a portion of the line shape function \(F\)
    !! for a given set \(\{p_j\}\)
    !!
    !! <h2>Background</h2>
    !! 
    !! The line shape function is defined in the paper as 
    !! \[F = \dfrac{1}{\Omega_k}\sum_{p_j}\left[\left(\prod_{j=1}^{M}F_j\right)
    !!       \sum_{j=1}^M\left(p_j + \dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}
    !!       \dfrac{I_{p_j+1}\left[\dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}\right]}{I_{p_j}\left[\dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}\right]}\right)
    !!       D(\omega_j)\right]\]
    !! This function calculates 
    !! \[\sum_{p_j}\left[\left(\prod_{j=1}^{M}F_j\right) 
    !!   \sum_{j=1}^M\left(p_j + \dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}
    !!   \dfrac{I_{p_j+1}\left[\dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}\right]}{I_{p_j}\left[\dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}\right]}\right)\right]\]
    !! The \(F_j\) terms are calculated based on (42) from the paper
    !! \[F_j = \exp\left[\dfrac{p_j\hbar\omega_j}{2kT} - 
    !!   S_j\coth\left(\dfrac{\hbar\omega_j}{2kT}\right)\right]
    !!   I_{p_j}\left[\dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}\right]\]
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: iE, j
    !
    real(kind = dp) :: E, Fj, prodFj, sumOverj, besPj, besRatio
    !
    prodFj = 1.0_dp
    sumOverj = 0.0_dp
    do j = 1, nModes
      !! * For each phonon mode
      !!    * Get the modified Bessel function \(I_{p_j}\)
      !!    * If the number of additional phonons \(p_j > 0\)
      !!      and \(I_{p_j} > 10^{-15}\), \(F_j = 0\)
      !!    * Otherwise, calculate \(F_j\) from (42)
      !!    * Multiply \(F_j\) on running product to get \(\prod_{j=1}^{M}F_j\)
      !!    * Calculate \(I_{p_j+1}(x)I_{p_j}(x)\)
      !!    * Add results to running total to get innermost sum in (43) in paper
      !!      (see Background for details)
      !
      Fj = 1.0_dp
      !
      besPj = besOrderNofModeM(abs(pj(j)), j)
      !
      !> @todo Change this to merge if statements @endtodo
      !if ( pj(j) > 0 .and. besPj > 1.0e-15_dp ) then
      !  !
      !  Fj = 0.0_dp
      !  !
      if ( pj(j) > 0 ) then
        !
        if ( besPj > 1.0e-15_dp ) then 
          !
          Fj = exp(pj(j)*wby2kT(j) - Sj(j)*coth(j))*besPj
          !
        else 
          !
          Fj = 0.0_dp
            !! @todo Figure out why don't just exit here because will be multiplying by 0 @endtodo
          !
        endif
        !
      else
        !
        Fj = exp(pj(j)*wby2kT(j) - Sj(j)*coth(j))*besPj
        !
      endif
      !
      prodFj = prodFj * Fj
      !
      besRatio = 0.5_dp*x(j)/(abs(pj(j))+1)
        ! Small argument approximation:
        ! \[I_{p_j}(x) \approx (x/2)^{p_j}/\Gamma(p_j+1)\]
        ! which means 
        ! \[I_{p_j+1}(x)/I_{p_j}(x) \approx \dfrac{(x/2)^{p_j+1}}{(x/2)^{p_j}}\dfrac{\Gamma(p_j+1)}{\Gamma(p_j+2)}
        !   = x/2\dfrac{p_j!}{(p_j+1)!} = x/2(n+1)\]
      !
      if  ( besPj > 1.0e-15_dp ) besRatio = besOrderNofModeM(abs(pj(j))+1, j)/besPj
        !! @todo Redo `besRatio` if statement to be more clear that it is if/else @endtodo
      !
      sumOverj = sumOverj + (abs(pj(j)) + x(j)*besRatio)
      !
    enddo
    !
    E = sum(pj(:)*phonF(:))
      !! * Calculate the energy gained by the extra phonons
    !
    !> * Calculate the energy index
    iE = 0
    if ( abs(E) > 1.0e-6_dp*evToHartree ) then
      iE = int(abs(E)/deltaE) + 1
      if ( E < 0.0_dp ) iE = -iE
    endif
    !
    iEbinsByBands(iE) = iEbinsByBands(iE) + 1
      !! * Increment the number of energies in the calculated bin
    !
    lsfVsEbyBands(iE) = lsfVsEbyBands(iE) + prodFj*sumOverj
      !! * Combine terms to get the middle portion of (43) in the paper (see Background for details)
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
          iRand = mod(abs(iRand), nModes) + 1
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
          pick = int( nModes*randy ) + 1
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
            pj(picks(iM)) = pj0s(ii,iM)*(-1)**(pms(iDes,iM-1))
          enddo
          !
          if ( abs(sum(abs(pj(picks(:)))) - m) > 0 ) then
            if (myid == root) then 
              write(iostd,*) 'ERROR', m, sum(abs(pj(picks(:)))), pj(picks(:))
              do iM = 1, l
                if (abs(pj(picks(iM))) < 1 ) then
                  write(iostd, *) 'ERROR 1', picks(iM)
                  write(iostd, *) 'ERROR 2', pj(picks(iM))
                  flush(iostd)
                endif
              enddo
            endif
          endif
          !
          prodFj = 1.0_dp
          sumOverj = 0.0_dp
          !
          do j = 1, nModes
            !
            Fj = 1.0_dp
            besPj = besOrderNofModeM(abs(pj(j)), j)
            if ( pj(j) > 0 ) then
              if ( besPj > 1.0e-15_dp ) then
                Fj = exp(pj(j)*wby2kT(j) - Sj(j)*coth(j))*besPj
              else
                Fj = 0.0_dp
              endif
            else
              Fj = exp(pj(j)*wby2kT(j) - Sj(j)*coth(j))*besPj
            endif
            !
            prodFj = prodFj * Fj
            !
            besRatio = 0.5_dp*x(j)/(abs(pj(j))+1)
            if  ( besPj > 1.0e-15_dp ) besRatio = besOrderNofModeM(abs(pj(j))+1, j)/besPj
            !
            sumOverj = sumOverj + (abs(pj(j)) + x(j)*besRatio)
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
  subroutine lsfMbyOneBand(M)
    !! Calculate the line shape function for 
    !! various sets \(\{p_j\}\) where all \(p_j\)
    !! are zero except one that has values
    !! \(-M, 0, M, 2M\) 
    !
    implicit none
    !
    integer, intent(in) :: M
      !! Number of phonons
    !
    integer :: iMode1, pm1
    !
    real(dp) :: t1, t2
    !
    call cpu_time(t1)
    !
    do iMode1 = 1, nModes
      !
      do pm1 = -M, M, 2*M
        !
        pj(:) = 0
        pj(iMode1) = pm1
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
    !! Calculate the line shape function for 
    !! various sets \(\{p_j\}\) where all \(p_j\)
    !! are zero except two that have values
    !! ??
    !!
    !! @todo Figure out how came up with limits @endtodo
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
      do iMode1 = 1, nModes - 1
        do iMode2 = iMode1+1, nModes
          !
          do pm1 = -l, l, 2*l
            do pm2 = -(m-l), (m-l), 2*(m-l)
              !
              pj(:) = 0
              pj(iMode1) = pm1
              pj(iMode2) = pm2
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
      do iMode2 = iMode1+1, nModes - 1
        do iMode3 = iMode2+1, nModes
          !
          do ii = 1, int(times3 + 1.e-3_dp)
            !
            do iDes = 0, 2**3 - 1
              !
              pj(:) = 0
              !
              pj(iMode1) = pj0s(ii,1)*(-1)**(pms(iDes,1-1))
              pj(iMode2) = pj0s(ii,2)*(-1)**(pms(iDes,2-1))
              pj(iMode3) = pj0s(ii,3)*(-1)**(pms(iDes,3-1))
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
  subroutine parallelIsFsBy3()
    !
    implicit none
    !
    real(dp) :: parTotal, parTotal2, totalStates, states, averageStatesPerProc
    integer :: i, iState, iproc
    !
    totalStates = nModes/6.0_dp
    totalStates = totalStates*(nModes-1)*(nModes-2)
    !
    !write(iostd, *) 'totalStates', totalStates
    if ( nModes > numprocs ) then
      !
      nProcMax = numprocs
      !
      parTotal = 0
      iproc = 0
      iModeIs(iproc) = 1
      iState = nModes
      !
      do while ( (iState > 1) .and. (nProcMax - iproc > 1) )
        !
        states = real(iState-1, dp)*real(iState-2, dp)/2.0_dp
        parTotal = parTotal + states
        parTotal2 = parTotal + real(iState-2, dp)*real(iState-3, dp)/2.0_dp 
        !
        averageStatesPerProc = totalStates/(nProcMax - iproc) - 0.01_dp
        if ( ( parTotal > averageStatesPerProc ) .or. ( parTotal2 > averageStatesPerProc ) ) then
          iModeFs(iproc) = nModes - iState + 1
          iproc = iproc + 1
          iModeIs(iproc) = nModes - iState + 2
          totalStates = totalStates - parTotal
          parTotal = 0
        endif
        !
        iState = iState - 1
        !
      enddo
      !
      iModeFs(iproc) = nModes - 2
      !
    else
      !
      nProcMax = nModes - 2
      !
      do i = 0, nProcMax - 1
        iModeIs(i) = i+1
        iModeFs(i) = i+1
      enddo
    endif
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
    totalStates = nModes/24.0_dp
    totalStates = totalStates*(nModes-1)*(nModes-2)*(nModes-3)
    !
    !write(iostd, *) 'totalStates', totalStates
    if ( nModes > numprocs ) then
      !
      nProcMax = numprocs
      !
      parTotal = 0
      iproc = 0
      iModeIs(iproc) = 1
      iState = nModes
      !
      do while ( (iState > 1) .and. (nProcMax - iproc > 1) )
        !
        states = real(iState-1, dp)*real(iState-2, dp)*real(iState-3, dp)/6.0_dp
        parTotal = parTotal + states
        parTotal2 = parTotal + real(iState-2, dp)*real(iState-3, dp)*real(iState-4, dp)/6.0_dp 
        !
        averageStatesPerProc = totalStates/(nProcMax - iproc) - 0.01_dp
        if ( ( parTotal > averageStatesPerProc ) .or. ( parTotal2 > averageStatesPerProc ) ) then
          iModeFs(iproc) = nModes - iState + 1
          iproc = iproc + 1
          iModeIs(iproc) = nModes - iState + 2
          totalStates = totalStates - parTotal
          parTotal = 0
        endif
        !
        iState = iState - 1
        !
      enddo
      !
      iModeFs(iproc) = nModes - 2
      !
    else
      ! 
      nProcMax = nModes - 2
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
