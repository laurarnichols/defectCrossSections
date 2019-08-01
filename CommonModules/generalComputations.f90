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
end module generalComputations
