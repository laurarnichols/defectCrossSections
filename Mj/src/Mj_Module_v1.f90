
module MjModule
  !
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
  integer, parameter :: int32 = selected_int_kind(5)
  integer, parameter :: int64 = selected_int_kind(15)
  integer, parameter :: iostd = 16, un = 3
  !
  real(kind = dp), parameter ::           pi = 3.1415926535897932_dp
  real(kind = dp), parameter ::        twopi = 2.0_dp*pi
  real(kind = dp), parameter ::         abCM = 0.529177219217e-8_dp
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
  real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
  !
  integer :: nAtoms, nOfqPoints, nModes
  !
  integer :: ios
  !
  real(kind = dp) :: ti, tf, t1, t2
  real(kind = dp) :: temperature, kT
  !
  integer, allocatable :: s2L(:)
  !
  real(kind = dp), allocatable :: atomD(:,:), atomM(:), phonQ(:,:), phonF(:), genCoord(:)
  real(kind = dp), allocatable :: atomPosition(:,:), newAtomicPosition(:,:)
  real(kind = dp), allocatable :: wby2kT(:), phonD(:,:,:,:), x(:), Sj(:), coth(:), besOrderNofModeM(:,:)
  !
  real(kind = dp) ::  maxDisplacement
  !
  integer :: modeI, modeF, qPoint
  !
  character(len = 2), allocatable :: elements(:)
  character(len = 6), parameter :: output = 'status'
  character(len = 256) :: phononsInput, equilibriumAtomicPositions, newAtomicPositions, QEInput
  !
  logical :: file_exists, readQEInput
  !
  namelist /MjInput/ QEInput, phononsInput, temperature, equilibriumAtomicPositions, modeI, modeF, qPoint, maxDisplacement
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
    READ (5, MjInput, iostat = ios)
    !
    call checkAndUpdateInput()
    !
    call readPhonons()
    !
    call readAtomicPositions()
    !
    return
    !
  end subroutine readInputs
  !
  !
  subroutine initialize()
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
    !
    implicit none
    !
    integer :: iAtom, iMode, iq
    real(kind = dp) :: dummyD, freqInTHz
    !
    CHARACTER :: dummyC
    !
    open(1, file=trim(phononsInput), status="old")
    !
    read(1,*) nOfqPoints, nAtoms
    !
    write(iostd, '(" Number of atoms : ", i5)') nAtoms
    write(iostd, '(" Number of q-Points : ", i5)') nOfqPoints
    flush(iostd)
    !
    nModes = 3*nAtoms - 3
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
  subroutine readAtomicPositions()
    !
    implicit none
    !
    integer :: iAtom
    !
    open(1, file=trim(equilibriumAtomicPositions), status="old")
    !
    allocate( elements(nAtoms), atomPosition(3,nAtoms) )
    !
    atomPosition(:,:) = 0.0_dp
    !
    do iAtom = 1, nAtoms
      read(1,*) elements(iAtom), atomPosition(1,iAtom), atomPosition(2,iAtom), atomPosition(3,iAtom)
    enddo
    !
    close(1)
    !
    return
    !
  end subroutine readAtomicPositions
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
    deallocate( atomM, atomD )
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
    integer :: i, j, iMode, nm, nb
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
    nb = modeF
    !
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
    temp(:) = x(:) ! exp(wby2kT(:) - Sj(:)*coth(:))*besOrderNofModeM(1,:)
      ! x(:)
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
    integer :: iAtom, iMode
    character(len = 256) :: line, fn, modeFolder, mkDir
    !
    do iMode = modeI, modeF
      !
      if ( s2L(iMode) < 10 ) then
        write(modeFolder, '("mode_", i1)') s2L(iMode)
      else if ( s2L(iMode) < 100 ) then
        write(modeFolder, '("mode_", i2)') s2L(iMode)
      else if ( s2L(iMode) < 1000 ) then
        write(modeFolder, '("mode_", i3)') s2L(iMode)
      else if ( s2L(iMode) < 10000 ) then
        write(modeFolder, '("mode_", i4)') s2L(iMode)
      endif
      !
      inquire(file= trim(modeFolder), exist = file_exists)
      if ( .not.file_exists ) then
        write(mkDir, '("mkdir -p ", a)') trim(modeFolder)
        call system(mkDir)
      endif
      !
      fn = trim(QEInput)
      fn = fn(INDEX(QEInput, '/', BACK=.TRUE.):INDEX(QEInput, '.in')-1)
      !
      write(iostd, '(" Writing new QE input file for mode :", i10)') s2L(iMode)
      !
      if ( s2L(iMode) < 10 ) then
        write(fn, '(a, "_mode", i1, ".in")') trim(fn), s2L(iMode)
      else if ( s2L(iMode) < 100 ) then
        write(fn, '(a, "_mode", i2, ".in")') trim(fn), s2L(iMode)
      else if ( s2L(iMode) < 1000 ) then
        write(fn, '(a, "_mode", i3, ".in")') trim(fn), s2L(iMode)
      else if ( s2L(iMode) < 10000 ) then
        write(fn, '(a, "_mode", i4, ".in")') trim(fn), s2L(iMode)
      endif
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
end module MjModule
