
module MjModule
  !
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
  integer, parameter :: int32 = selected_int_kind(5)
  integer, parameter :: int64 = selected_int_kind(15)
  integer, parameter :: iostd = 16, un = 3
  !
  real(kind = dp), parameter ::         abCM = 0.529177219217e-8_dp
  real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
  real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
  real(kind = dp), parameter ::           pi = 3.1415926535897932_dp
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  real(kind = dp), parameter ::        twopi = 2.0_dp*pi
  !
  character(len = 6), parameter :: output = 'status'
  !
  !
  integer :: ios
  integer :: modeF
  integer :: modeI
  integer :: nAtoms
  integer :: nModes
  integer :: nOfqPoints
  integer :: qPoint
  !
  real(kind = dp) :: kT
  real(kind = dp) :: maxDisplacement
  real(kind = dp) :: t1
  real(kind = dp) :: t2
  real(kind = dp) :: temperature
  real(kind = dp) :: tf
  real(kind = dp) :: ti
  !
  character(len = 256) :: equilibriumAtomicPositions
  character(len = 256) :: newAtomicPositions
  character(len = 256) :: phononsInput
  character(len = 256) :: QEInput
  !
  logical :: file_exists
  logical :: readQEInput
  !
  !
  integer, allocatable :: s2L(:)
  !
  real(kind = dp), allocatable :: atomD(:,:)
  real(kind = dp), allocatable :: atomM(:)
  real(kind = dp), allocatable :: atomPosition(:,:)
  real(kind = dp), allocatable :: besOrderNofModeM(:,:)
  real(kind = dp), allocatable :: coth(:)
  real(kind = dp), allocatable :: genCoord(:)
  real(kind = dp), allocatable :: newAtomicPosition(:,:)
  real(kind = dp), allocatable :: phonD(:,:,:,:)
  real(kind = dp), allocatable :: phonF(:)
  real(kind = dp), allocatable :: phonQ(:,:)
  real(kind = dp), allocatable :: Sj(:)
  real(kind = dp), allocatable :: wby2kT(:)
  real(kind = dp), allocatable :: x(:)
  !
  character(len = 2), allocatable :: elements(:)
  !
  namelist /MjInput/ QEInput, phononsInput, temperature, equilibriumAtomicPositions, modeI, modeF, qPoint, maxDisplacement
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
    implicit none
    !
    !> * Check if file output exists; if it does, delete it
    inquire(file = output, exist = file_exists)
    !
    if ( file_exists ) then
      !
      open (unit = 11, file = output, status = "old")
      !
      close(unit = 11, status = "delete")
      !
    endif
    !
    open (iostd, file = output, status='new')
      !! * Open new output file
    !
    call initialize()
      !! * Set default values of input parameters
    !
    READ (5, MjInput, iostat = ios)
      !! * Read input parameters
    !
    call checkAndUpdateInput()
      !! * Check if input parameters were updated 
      !!   and do some basic checks
    !
    call readPhonons()
      !! * Read the phonons output from QE or VASP
    !
    call readAtomicPositions()
      !! * Read the equilibrium atomic positions
    !
    return
    !
  end subroutine readInputs
  !
  !
  subroutine initialize()
    !! Set default values for input parameters
    !
    implicit none
    !
    QEInput = ''
    phononsInput = ''
    equilibriumAtomicPositions = ''
    temperature = -1.0_dp
    maxDisplacement = -1.0_dp
    modeI = -1
    modeF = -1
    !
    return
    !
  end subroutine initialize
  !
  !
  subroutine checkAndUpdateInput()
    !! Check that the input variables don't still have their default
    !! values. The program will abort here if:
    !! * `equilibriumAtomicPositions` is not defined
    !! * `phononsInput` is not defined
    !! * `temperature` is not defined
    !! * `modeI` or `modeF` is not defined
    !! * `modeF < modeI`
    !! * `maxDisplacement` is not defined
    !
    implicit none
    !
    logical :: abortExecution = .false.
    !
    write(iostd, *)
    !
    if ( equilibriumAtomicPositions == '' ) then
      write(iostd, '(" equilibriumAtomicPositions is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Equilibrium Atomic Positions input : ", a)') trim(equilibriumAtomicPositions)
    endif
    !
    if ( phononsInput == '' ) then
      write(iostd, '(" PhononsInput is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Phonons input : ", a)') trim(PhononsInput)
    endif
    !
    if ( QEInput == '' ) then
      write(iostd, '(" QEInput is not defined!")')
      readQEInput = .false.
    else
      readQEInput = .true.
      write(iostd, '(" Quantum Espresso input : ", a)') trim(QEInput)
    endif
    !
    if ( temperature < 0.0_dp ) then
      write(iostd, '(" Variable temperature has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Temperature : ", f10.2, " Kelvin.")') temperature
      kT = temperature*8.6173324d-5*eVToHartree
    endif
    !
    if ( modeI < 0 ) then
      write(iostd, '(" Variable modeI has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Initial mode : ", i5)') modeI 
    endif
    !   
    if ( modeF < 0 ) then
      write(iostd, '(" Variable modeF has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Final mode : ", i5)') modeF
    endif
    !
    if ( modeF < modeI ) then
      write(iostd, '(" Final mode is set smaller than initial one!")')
      abortExecution = .true.
    endif
    !
    if ( maxDisplacement < 0 ) then
      write(iostd, '(" Variable maxDisplacement has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Maximum atomic displacement in each direction : ", f15.10)') maxDisplacement
    endif
    !
    if ( abortExecution ) then
      write(iostd, '(" *************************** ")')
      write(iostd, '(" * Program stops!          * ")')
      write(iostd, '(" *************************** ")')
      stop
    endif
    !
    return
    !
  end subroutine checkAndUpdateInput
  !
  !
  subroutine readPhonons()
    !! Read the number of atoms and q points and
    !! get the phonon information like frequency 
    !! and displacements
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: iAtom
      !! Loop index over atoms
    integer :: iMode
      !! Loop index over phonon modes
    integer :: iq
      !! Loop index over q points
    !
    real(kind = dp) :: dummyD
      !! Dummy variable to ignore input
    real(kind = dp) :: freqInTHz 
      !! Input frequency in THz
    !
    character :: dummyC
      !! Dummy variable to ignore input
    !
    open(1, file=trim(phononsInput), status="old")
      !! * Open `phononsInput` file
    !
    read(1,*) nOfqPoints, nAtoms
      !! * Read in the number of q points and number of atoms
    !
    write(iostd, '(" Number of atoms : ", i5)') nAtoms
    write(iostd, '(" Number of q-Points : ", i5)') nOfqPoints
    flush(iostd)
      !! * Write the number of atoms and q points to the output file
    !
    nModes = 3*nAtoms - 3
      !! * Calculate the number of phonon modes
    !
    read (1,*)
      !! * Ignore the next line as it is blank
    !
    allocate( atomD(3,nAtoms), atomM(nAtoms) )
    !
    atomD = 0.0_dp
    atomM = 0.0_dp
    !
    do iAtom = 1, nAtoms
      !! * For each atom, read in the displacement (either pristine-defect 
      !!   or defect-pristine) and the atom mass
      !
      read(1,*) atomD(1,iAtom), atomD(2,iAtom), atomD(3,iAtom), atomM(iAtom)
      !
    enddo
    !
    read(1,*)
      !! * Ignore the next line as it is blank
    !
    allocate( phonQ(3,nOfqPoints), phonF(nModes), phonD(3,nAtoms,nModes,nOfqPoints) )
    !
    phonQ = 0.0_dp
    phonF = 0.0_dp
    phonD = 0.0_dp
    !
    do iq = 1, nOfqPoints
      !! * For each q point
      !!    * Read in the coordinates in \(q\)-space
      !!    * For each phonon mode
      !!      * Read in the phonon frequency in THz
      !!      * Convert the frequency to Hartree
      !!      * Read in the atom displacements
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
      !! * Close `phononsInput` file
    !
    flush(iostd)
    !
    return
    !
  end subroutine readPhonons
  !
  !
  subroutine readAtomicPositions()
    !! Read in the element and equilibrium position for 
    !! each atom
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer :: iAtom
      !! Loop index over atoms
    !
    open(1, file=trim(equilibriumAtomicPositions), status="old")
      !! * Open the `equilibriumAtomicPositions` file
    !
    allocate( elements(nAtoms), atomPosition(3,nAtoms) )
    !
    atomPosition(:,:) = 0.0_dp
    !
    do iAtom = 1, nAtoms
      !! * For each atom, read in the element and equilibrium position
      !
      read(1,*) elements(iAtom), atomPosition(1,iAtom), atomPosition(2,iAtom), atomPosition(3,iAtom)
      !
    enddo
    !
    close(1)
      !! * Close the `equilibriumAtomicPositions` file
    !
    return
    !
  end subroutine readAtomicPositions
  !
  !
  subroutine computeGeneralizedDisplacements()
    !! Calculate the generalized displacements
    !! by dotting the phonon displacements with 
    !! the atom displacements
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: iAtom
      !! Loop index over atoms
    integer :: iMode
      !! Loop index over phonon modes
    integer :: iq
      !! Loop index over q points
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
          !! * For each q point, mode, and atom combination, calculate
          !!   the generalized displacement as 
          !!   \[\sum_{\text{mode}} \sqrt{1823m}\mathbf{\Delta r}_{\text{phonon}}\cdot\mathbf{\Delta r}_{\text{atom}}\]
          !
          genCoord(iMode) = genCoord(iMode) + sqrt(1822.88833218_dp*atomM(iAtom))*sum(phonD(:,iAtom,iMode,iq)*atomD(:,iAtom))
          !
        enddo
        !
      enddo
      !
    enddo
    !
    deallocate( atomM, atomD )
    !
    return
    !
  end subroutine computeGeneralizedDisplacements
  !
  !
  subroutine computeVariables()
    !! Calculate biggest portions of equations 42 and 43 to make
    !! entire equation more manageable
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: i, j, iMode, nm, 
    integer :: nb
      !! Final phonon mode
    integer :: iMode
      !! Loop index over phonon modes
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
      !! * Calculate 
      !!   \[S_j = \dfrac{\omega_j}{2\hbar}N\delta q_j^2\] 
      !!   where \(\omega_i\rightarrow\)`phonF(j)`, 
      !!   \(\delta q_j\rightarrow\)`genCoord(j)`, and \(N\)
      !!   is the number of atoms per supercell
      !! @note This is equation 44 in the paper @endnote
      !! @todo Figure out why there is no \(N\) in this equation in the code @endtodo
    !
    wby2kT(:) = phonF(:)/(2.0_dp*kT)
      !! * Calculate \(\hbar\omega_j/2kT\) that is the argument 
      !!   to hyperbolic trig functions
      !! @note This is mainly from equations 42 and 43 in the paper @endnote
    !
    coth(:) = cosh(wby2kT(:))/sinh(wby2kT(:))
      !! * Calculate \(\coth(\hbar\omega_j/2kT)\) 
      !! @note This is in equation 42 in the paper @endnote
      !! @todo Figure out if this needs to be another variable @endtodo
    !
    x(:) = Sj(:)/sinh(wby2kT(:))
      !! * Calculate the argument to the modified Bessel functions
      !!   \[\dfrac{S_j}{\sinh(\hbar\omega_j/2kT)}\]
    !
    allocate( s2L(nModes) )
    s2L(:) = 0
    !
    do iMode = 1, nModes
      !! * Create an array of the indices that will be rearranged
      !
      s2L(iMode) = iMode
      !
    enddo
    !
    call arrangeLargerToSmaller()
      !! * Rearrange the indices based on ordering the arguments (`x`)
      !!   largest to smallest
    !
    open(11, file='modes', status='unknown')
    !
    write(11, '("#Mode, frequency (eV),        genCoord(Mode),     genCoord(Mode)^2,  Sj/sinh(wby2kT)")')
      !! @todo Frequency should actually be in eV/\(\hbar\). Check that that's the case. @endtodo
    !
    do iMode = 1, nModes
      !! * Write out the index, \(\omega_j\) in eV, \(\delta q_j\), \(\delta q_j^2\), 
      !!   and \(S_j/\sinh(\hbar\omega_j/2kT)\) for each phonon mode 
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
    nb = modeF
      !! @todo Figure out why this is assigned to another variable @endtodo
    !
    allocate( besOrderNofModeM(0:nb + 1, nModes) )
    allocate( bi(0:nb + 1), di(0:nb + 1) )
    allocate( bk(0:nb + 1), dk(0:nb + 1) )
      !! @todo Figure out if still need `di`, `bk`, and `dk` @endtodo
    !
    do j = 1, nModes
      !
      bi(:) = 0.0_dp
      !
      nm = nb + 1
        !! @todo Figure out if need to have this in loop. Why change `nm` in `iknb`? @endtodo
      !
      call iknb(nb + 1, x(j), nm, bi) ! , di, bk, dk)
        !! @todo Figure out if should send `nm` as it is immediately modified and not used here @endtodo
      !
      do i = 0, nb + 1
        !! * Store the modified Bessel function for each mode and order
        !! @todo Possibly change `besOrderNofModeM` to `modBesOrderNofModeM` @endtodo
        !! @todo Figure out why this loop is here. Can not just pass `besOrderNofModeM(:,j)` @endtodo
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
    !! Sort `s2L` based on descending order
    !! of `x`
    !!
    !! @todo Change this to a more efficient algorithm @endtodo
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
        if ( temp(i) < temp(i+1) ) then 
          !
          tmpi = s2L(i)
          s2L(i) = s2L(i+1)
          s2L(i+1) = tmpi
          !
          tmpr = temp(i)
          temp(i) = temp(i+1)
          temp(i+1) = tmpr
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
  subroutine displaceAtoms()
    !
    implicit none
    !
    integer :: istat, iAtom, iMode, iRand
    real(kind = dp) :: ran, norm
    !
    allocate ( newAtomicPosition(3, nAtoms) )
    !
    open(unit=12, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
    !
    write(iostd, *)
    !
    do iMode = modeI, modeF
      !
      write(iostd, '(" Calculating new atomic positions for mode :", i10)') s2L(iMode)
      !
      do iAtom = 1, nAtoms
        !
        norm = sqrt(sum(phonD(:,iAtom,s2L(iMode),qPoint)**2)) 
        !
        read(12) iRand
        ran = mod(abs(iRand),10000000)/1.0e7_dp
        !
        newAtomicPosition(:,iAtom) = atomPosition(:,iAtom) + maxDisplacement*ran*phonD(:,iAtom,s2L(iMode),qPoint)/norm
        !
      enddo
      !
    enddo
    !
    close(12)
    !
  end subroutine displaceAtoms
  !
  !
  subroutine writeNewAtomicPositions()
    !
    implicit none
    !
    integer :: iAtom, iMode
    !
    do iMode = modeI, modeF
      !
      write(iostd, '(" Writing new atomic positions for mode :", i10)') s2L(iMode)
      !
      if ( s2L(iMode) < 10 ) then
        write(newAtomicPositions, '("newPositionsForMode", i1)') s2L(iMode)
      else if ( s2L(iMode) < 100 ) then
        write(newAtomicPositions, '("newPositionsForMode", i2)') s2L(iMode)
      else if ( s2L(iMode) < 1000 ) then
        write(newAtomicPositions, '("newPositionsForMode", i3)') s2L(iMode)
      else if ( s2L(iMode) < 10000 ) then
        write(newAtomicPositions, '("newPositionsForMode", i4)') s2L(iMode)
      else
        newAtomicPositions = 'newPositions'
      endif
      !
      open(21, file=trim(newAtomicPositions), status='unknown')
      !
      do iAtom = 1, nAtoms
        !
        !write(6, '(i4, f15.12, 3f15.5, " | ", 3f15.5)') iAtom, maxDisplacement*ran, &
        !          atomPosition(:,iAtom), atomPosition(:,iAtom) + maxDisplacement*ran*phonD(:,iAtom,s2L(iMode),qPoint)/norm
        write(21,*) elements(iAtom), newAtomicPosition(:,iAtom)
      enddo
      !
      close(21)
      !
    enddo
    !
  end subroutine writeNewAtomicPositions
  !
  !
  subroutine exportQEInput()
    !
    integer :: iAtom
      !! Loop index over atoms
    integer :: iMode
      !! Loop index over phonon modes
    character(len = 256) :: line, fn, modeFolder, mkDir, s2LStr
    !
    do iMode = modeI, modeF
      !
      call int2str(s2L(iMode), s2LStr)
      !
      write(modeFolder, '("mode_", a)') trim(s2LStr)
      !
      !> If the folder for this mode doesn't exist, create it
      inquire(file= trim(modeFolder), exist = file_exists)
      if ( .not.file_exists ) then
        !
        write(mkDir, '("mkdir -p ", a)') trim(modeFolder)
        call system(mkDir)
        !
      endif
      !
      fn = trim(QEInput)
      fn = fn(INDEX(QEInput, '/', BACK=.TRUE.):INDEX(QEInput, '.in')-1)
      !
      write(iostd, '(" Writing new QE input file for mode :", i10)') s2L(iMode)
      !
      write(fn, '(a, "_mode", a, ".in")') trim(fn), trim(s2LStr)
      !
      fn = trim(modeFolder)//"/"//trim(fn)
      !
      open(2, file=trim(fn), status="unknown")
      !
      open(1, file=trim(QEInput), status="old")
      !
      do 
        read(1,'(a)', END = 100) line
        write(2,'(a)') trim(line)
        if ( INDEX(line, 'ATOMIC_POSITIONS') /= 0 ) then
          do iAtom = 1, nAtoms
            read(1,'(a)') line
            write(2,*) elements(iAtom), newAtomicPosition(:,iAtom)
          enddo
        endif
      enddo
100   continue
      !
      close(1)
      close(2)
      !
    enddo
    !
  end subroutine exportQEInput



!  subroutine readMjs()
!    !
!    implicit none
!    !
!    integer :: i, iE0, iE, dummyI, nEMjs
!    real(kind = dp) :: dummyD, E, MjsOfE, MjOfE0
!    character :: dummyC
!    !
!    open(1, file=trim(MjsInput), status="old")
!    !
!    read(1, *) dummyC, nEMjs
!    !
!    allocate ( Mjs(-nEnergies:nEnergies) )
!    !
!    Mjs = 1.0_dp
!    !
!!    read(1, '(d22.14,i5,4d22.14)' ) E, dummyI, dummyD, MjsOfE0, dummyD
!!    !
!!    E = E*eVToHartree 
!!    iE = int(E/deltaE) + 1
!!    !
!!    do i = 1, nEMjs - 1
!!      !
!!      iE0 = iE
!!      read(1, '(d22.14,i5,4d22.14)' ) E, dummyI, dummyD, MjsOfE, dummyD
!!      E = E*eVToHartree
!!      iE = int(E/deltaE) + 1
!!      Mjs(iE0:iE) = MjsOfE0
!!      MjsOfE0 = MjsOfE
!!      !
!!    enddo
!!    !
!!    close(1)
!!    !
!!    !do iE = -nEnergies, nEnergies
!!    !  write(44,*) real(iE, dp)*deltaE*HartreeToEv, Mjs(iE)
!!    !enddo
!!    !
!    return
!    !
!  end subroutine readMjs
  !
  !
  subroutine iknb ( n, x, nm, bi) !, di, bk, dk )
    !! author: Shanjie Zhang, Jianming Jin
    !! date: 17 July 2012
    !
    ! Modified : when x < 10^(-15) return the limiting value for small argument [ I_n(x) ~ (x/2)^n Gamma(n+1) ]
    ! 
    !*********************************************************************72
    !
    !! IKNB compute Bessel function In(x) and Kn(x).
    !!
    !!  <h2>Discussion</h2>
    !!
    !!    Compute modified Bessel functions In(x) and Kn(x),
    !!    and their derivatives.
    !!
    !!  <h2>Licensing</h2>
    !!
    !!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
    !!    they give permission to incorporate this routine into a user program 
    !!    provided that the copyright is acknowledged.
    !!
    !!
    !!  <h2>Reference</h2>
    !!
    !!    Shanjie Zhang, Jianming Jin,
    !!    Computation of Special Functions,
    !!    Wiley, 1996,
    !!    ISBN: 0-471-11963-6,
    !!    LC: QA351.C45.
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: n
      !! Order of \(I_n(x)\) and \(K_n(x)\)
    integer, intent(inout) :: nm
      !! The highest order computed
    integer :: ik
    integer :: k
!    integer :: k0
!    integer :: l
    integer :: m
!    integer :: msta1
!    integer :: msta2
    !
!    double precision :: a0
    double precision :: bi(0:n)
      !! \(I_n(x)\)
!    double precision :: bkl
    double precision :: bs
    double precision :: el
    double precision :: f
    double precision :: fact
    double precision :: f0
    double precision :: f1
!    double precision :: g
!    double precision :: g0
!    double precision :: g1
    double precision :: pi
!    double precision :: r
    double precision :: s0
    double precision :: sk0
!    double precision :: vt
    double precision :: x
      !! The argument
    !
    pi = 3.141592653589793D+00
    el = 0.5772156649015329D+00
    nm = n
    !
    if ( x .le. 1.0D-15 ) then
      !! * If \(x < 10^{-15}\), use the limiting value for  a small argument 
      !!   \[I_n(x) \approx \left(\dfrac{x}{2}\right)^2\Gamma(n+1)\]
      !!   where \(\Gamma(n+1) = n!\) to calculate multiple orders of the 
      !!   modified Bessel function (up to \(n\)) for a single value of \(x\)
      !
      do k = 0, n
        ! For each order
        !
        fact = 1.0_dp
        !
        do ik = 2, k
          ! Calculate the factorial
          !
          fact = fact*ik
          !
        enddo
        !
        bi(k) = (0.5_dp*x)**k/fact
          ! Calculate \(I_n(x)\)
        !
      enddo
      !
      return
      !
    endif
    !
    if ( n .eq. 0 ) then
      !
      nm = 1
      !
    end if
    !
    m = msta1 ( x, 200 )
      !! * Otherwise, get the starting order \(m\) such
      !!   that \(J_n(x)\approx10^{-200}\)
    !
    if ( m .lt. nm ) then
      !! * If \(m\) is less than the max order to be 
      !!   calculated, set the max to \(m\)
      !
      nm = m
      !
    else
      !! * Otherwise, set \(m\) to the starting order such
      !!   that all \(J_{nm}(x)\) have 15 significant digits
      !
      m = msta2 ( x, nm, 15 )
      !
    end if
    !
    bs = 0.0D+00
    sk0 = 0.0D+00
    f0 = 0.0D+00
    f1 = 1.0D-100
    !
    do k = m, 0, -1
    !
      f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 + f0
      !
      if ( k .le. nm ) then
        !
        bi(k) = f
        !
      end if
      !
      if ( k .ne. 0 .and. k .eq. 2 * int ( k / 2 ) ) then
        !
        sk0 = sk0 + 4.0D+00 * f / k
        !
      end if
      !
      bs = bs + 2.0D+00 * f
      !
      f0 = f1
      !
      f1 = f
      !
    enddo
    !
    s0 = exp ( x ) / ( bs - f )
    !
    do k = 0, nm
      !
      bi(k) = s0 * bi(k)
      !
    end do
    !
    return
  end subroutine iknb
  !
  !
  FUNCTION msta1(x, mp) RESULT(fn_val) 
    !! Determine the starting point for backward 
    !! recurrence such that the magnitude of 
    !! \(J_n(x)\) at that point is about 
    !! \(10^{-\text{mp}}\) 
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: mp 
      !! Magnitude
    integer :: fn_val 
      !! Starting point
    integer :: it
      !! Loop index
    integer    :: n0, n1, nn
      !! Order of Bessel function
    ! 
    real(kind = dp), intent(in) :: x 
      !! Argument of \(J_n(x)\)
    real(kind = dp) :: a0
      !! \(|x|\)
    real(kind = dp) :: f, f0, f1
      !! \(-\log(J_n(x)) - \text{mp}\)
    ! 
    a0 = abs(x) 
    !
    n0 = int(1.1*a0) + 1 
    !
    f0 = envj(n0,a0) - mp 
      !! * Get initial guess for \(f = -log(J_{n_0}(a_0)) - \text{mp}\)
    !
    n1 = n0 + 5 
    !
    f1 = envj(n1,a0) - mp 
      !! * Get next guess for \(f\)
    !
    do  it = 1, 20 
      !! * Recursively minimize \(f\)
      !
      nn = n1 - int((n1-n0)/(1.0_dp - f0/f1))
      !
      f = envj(nn,a0) - mp 
      !
      if (abs(nn-n1) < 1) exit 
      !
      n0 = n1 
      !
      f0 = f1 
      !
      n1 = nn 
      !
      f1 = f 
      !
    enddo 
    ! 
    fn_val = nn 
    !
    return
    !
  end function msta1 
  !
  ! 
  function msta2(x, n, mp) result(fn_val) 
    !! Determine the starting point for backward 
    !! recurrence such that all \(J_n(x)\) has mp
    !! significant digits 
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: n 
      !! Order of \(J_n(x)\)
    integer, intent(in) :: mp 
      !! Significant digit
    integer :: fn_val 
      !! Starting point
    integer :: it, n0, n1, nn 
    ! 
    real(kind = dp), intent(in) :: x 
      !! Argument of \(J_n(x)\)
    real(kind = dp) :: a0, ejn, f, f0, f1, hmp, obj 
    ! 
    a0 = ABS(x) 
    ! 
    hmp = 0.5_dp * mp 
    ! 
    ejn = envj(n, a0) 
    ! 
    if (ejn <= hmp) then
      ! 
      obj = mp 
      ! 
      n0 = int(1.1*a0) 
      ! 
    else
      ! 
      obj = hmp + ejn 
      ! 
      n0 = n 
      ! 
    endif 
    !
    if ( n0 < 1 ) n0 = 1
    !
    f0 = envj(n0,a0) - obj 
    !
    n1 = n0 + 5 
    !
    f1 = envj(n1,a0) - obj 
    !
    do  it = 1, 20 
      !
      nn = n1 - int((n1-n0)/(1.0_dp - f0/f1))
      !
      f = envj(nn, a0) - obj 
      !
      if (abs(nn-n1) < 1) exit
      !
      n0 = n1 
      !
      f0 = f1 
      !
      n1 = nn 
      !
      f1 = f 
      !
    enddo
    !
    fn_val = nn + 10 
    !
    return
    ! 
  end function msta2 
  !
  ! 
  function envj(n, x) result(fn_val) 
    !! Estimates \(-\log(J_n(x))\) from the estimate
    !! \[J_n(x) \approx \dfrac{1}{\sqrt{2\pi n}}\left(\dfrac{ex}{2n}\right)^n\]
    ! 
    implicit none
    !
    integer, intent(in) :: n 
      !! Order of Bessel function
    real(kind = dp), intent(in) :: x 
      !! Argument
    real(kind = dp) :: fn_val 
      !! \(-\log(J_n(x))\)
    !
    fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n) 
      ! \(6.28 = 2\pi\) and \(1.36 = e/2\)
    !
    RETURN 
    !
  END FUNCTION envj 
  !
  !
!=====================================================================================================
! Utility functions that simplify the code and may be used multiple times
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  function wasRead(inputVal, variableName, usage, abortExecution) 
    !! Determine if an input variable still has the default value.
    !! If it does, output an error message and possibly set the program
    !! to abort. Not all variables would cause the program to abort,
    !! so this program assumes that if you pass in the logical `abortExecution`
    !! then the variable is required and causes the program to abort 
    !! if missing.
    !!
    !! I could not find a clean way to allow this function to receive
    !! different types of variables (integer, real, character, etc.), so
    !! I made the argument be an integer so that each type could be sent
    !! in a different way. Each case is set up so that the value is tested to
    !! see if it is less than zero to determine if the variable still has
    !! its default value
    !!
    !! * For strings, the default value is `''`, so pass in 
    !! `LEN(trim(variable))-1` as this should be less than zero if
    !! the string still has the default value and greater than or equal 
    !! to zero otherwise
    !! * For integers the default values are less than zero, so just pass as is 
    !! * Real variables also have a negative default value, so just pass the
    !! value cast from real to integer
    !!
    implicit none
    !
    integer, intent(in) :: inputVal
      !! Value to compare with 0 to see if a variable has been read;
    !
    character(len=*), intent(in) :: variableName
      !! Name of the variable used in output message
    character(len=*), intent(in) :: usage
      !! Example of how the variable can be used
    !
    logical, optional, intent(inout) :: abortExecution
      !! Optional logical for if the program should be aborted 
    logical :: wasRead
      !! Whether or not the input variable was read from the input file;
      !! this is the return value
    !
    !! <h2>Walkthrough</h2>
    !!
    wasRead = .true.
      !! * Default return value is true
    !
    if ( inputVal < 0) then
      !! * If the input variable still has the default value
      !!    * output an error message
      !!    * set the program to abort if that variable was sent in
      !!    * set the return value to false to indicate that the 
      !!      variable wasn't read
      !
      write(iostd, *)
      write(iostd, '(" Variable : """, a, """ is not defined!")') variableName
      write(iostd, '(" usage : ", a)') usage
      if(present(abortExecution)) then
        !
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
        !
      endif 
      !
      wasRead = .false.
      !
    endif
    !
    return
    !
  end function wasRead
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine int2str(integ, string)
    !! Write a given integer to a string, using only as many digits as needed
    !
    implicit none
    integer :: integ
    character(len = 300) :: string
    !
    if ( integ < 10 ) then
      write(string, '(i1)') integ
    else if ( integ < 100 ) then
      write(string, '(i2)') integ
    else if ( integ < 1000 ) then
      write(string, '(i3)') integ
    else if ( integ < 10000 ) then
      write(string, '(i4)') integ
    endif
    !
    string = trim(string)
    !
    return
    !
  end subroutine int2str
  !
end module MjModule
