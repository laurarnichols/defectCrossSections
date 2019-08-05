module generalComputations
  !
  use constants
  !
  implicit none
  !
  contains
  !
  subroutine computeGeneralizedDisplacements(nOfqPoints, nModes, genCoord, nAtoms, atomM, phonD, atomD)
    !! Calculate the generalized displacements
    !! by dotting the phonon displacements with
    !! the atom displacements
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: nAtoms
      !! Number of atoms in system
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: nOfqPoints
      !! Number of q points
    integer :: iAtom
      !! Loop index over atoms
    integer :: iMode
      !! Loop index over phonon modes
    integer :: iq
      !! Loop index over q points
    !
    real(kind = dp), intent(in) :: atomD(3,nAtoms)
      !! Equilibrium displacements in defect versus pristine
    real(kind = dp), intent(in) :: atomM(nAtoms)
      !! Atom masses
    real(kind = dp), intent(out) :: genCoord(nModes)
      !! Generalized coordinates \(\delta q_j\)
    real(kind = dp), intent(in) :: phonD(3,nAtoms,nModes,nOfqPoints)
      !! Phonon displacements
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
    return
    !
  end subroutine computeGeneralizedDisplacements                                                                                                                              
  !  
  subroutine computeVariables(x, Sj, coth, wby2kT, phonF, genCoord, kT, s2L, nModes, maximumNumberOfPhonons, besOrderNofModeM)
    !! Calculate biggest portions of equations 42 and 43 to make
    !! entire equation more manageable
    !!
    !! <h2>Walkthrough</h2>
    !!
    use miscUtilities
      !! Include `miscUtilities` module for call to
      !! `arrangeLargerToSmaller`
    !
    implicit none
    !
    integer, intent(in) :: maximumNumberOfPhonons
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(out) :: s2L(nModes)
    integer :: i, j, nm
    integer :: nb
      !! Final phonon mode
    integer :: iMode
      !! Loop index over phonon modes
    !
    real(kind = dp), intent(out) :: besOrderNofModeM(0:maximumNumberOfPhonons+1,nModes)
    real(kind = dp), intent(out) :: coth(nModes)
      !! Hyperbolic cotangent
    real(kind = dp), intent(in) :: genCoord(nModes)
      !! Generalized coordinates \(\delta q_j\)
    real(kind = dp), intent(in) :: kT
    real(kind = dp), intent(in) :: phonF(nModes)
      !! Phonon frequencies in Hartree
    real(kind = dp), intent(out) :: Sj(nModes)
      !! \(S_j\) in equation 44 in paper
    real(kind = dp), intent(out) :: wby2kT(nModes)
      !! \(\omega/2kT\)
    real(kind = dp), intent(out) :: x(nModes)
      !! Argument to modified Bessel function
    real(kind = dp), allocatable :: bi(:)
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
    s2L(:) = 0
    !
    do iMode = 1, nModes
      !! * Create an array of the indices that will be rearranged
      !
      s2L(iMode) = iMode
      !
    enddo
    !
    call arrangeLargerToSmaller(nModes, x, s2L)
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
    nb = maximumNumberOfPhonons
      !! @todo Figure out why this is assigned to another variable @endtodo
    !
    allocate( bi(0:nb + 1) )
      !! @todo Figure out if still need `di`, `bk`, and `dk` @endtodo
    !
    do j = 1, nModes
      !
      bi(:) = 0.0_dp
      !
      nm = nb + 1
        !! @todo Figure out if need to have this in loop. Why change `nm` in `iknb`? @endtodo
      !
      call iknb(nb + 1, x(j), nm, bi)
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
    deallocate( bi )
    !
    return
    !
  end subroutine computeVariables
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
  function msta1(x, mp) result(fn_val) 
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
    fn_val = 0.5_dp * log10(6.28_dp*n) - n * log10(1.36_dp*x/n) 
      ! \(6.28 = 2\pi\) and \(1.36 = e/2\)
    !
    return
    !
  end function envj 
  !
end module generalComputations
