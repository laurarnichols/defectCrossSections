module readInputFiles
  !
  use constants
  !
  implicit none
  !
  contains
  !
  subroutine readPhonons(phononsInput, nOfqPoints, nAtoms, nModes, phonQ, phonF, phonD)
    !! Read the number of atoms and q points and
    !! get the phonon information like frequency 
    !! and displacements
    !!
    !! @note This subroutine is only meant for QE 5.3.0 phonon output @endnote
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: nOfqPoints
      !! Number of q points
    integer, intent(in) :: nAtoms
      !! Number of atoms in system
    integer, intent(out) :: nModes
      !! Number of phonon modes
    integer :: iAtom
      !! Loop index over atoms
    integer :: iMode
      !! Loop index over phonon modes
    integer :: iq
      !! Loop index over q points
    !
    real(kind = dp), allocatable, intent(out) :: phonD(:,:,:,:)
      !! Phonon displacements
    real(kind = dp), allocatable, intent(out) :: phonF(:)
      !! Phonon frequencies
    real(kind = dp), allocatable, intent(out) :: phonQ(:,:)
      !! Phonon q point coordinates
    real(kind = dp) :: dummyD
      !! Dummy variable to ignore input
    real(kind = dp) :: freqInTHz 
      !! Input frequency in THz
    !
    character(len = 256), intent(in) :: phononsInput
      !! QE phonon output file name
    character :: dummyC
      !! Dummy variable to ignore input
    !
    open(1, file=trim(phononsInput), status="old")
      !! * Open `phononsInput` file
    !
    !
    nModes = 3*nAtoms - 3
      !! * Calculate the number of phonon modes
    write(iostd, '(" Number of atoms : ", i5)') nAtoms
    write(iostd, '(" Number of q-Points : ", i5)') nOfqPoints
    write(iostd, '(" Number of modes : ", i5)') nModes
    flush(iostd)
      !! * Write the number of atoms, q points, and modes to the output file
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
      read(1,*) 
      read(1,*) 
        !! Ignore the first 2 lines
      !
      read (1,*) dummyC, dummyC, phonQ(1,iq), phonQ(2,iq), phonQ(3,iq)
      !
      read(1,*)
        !! Skip the asterisk separator
      !
      do iMode = 1, nModes
        !
        read(1,*) dummyC, dummyC, dummyC, dummyC, freqInTHz, dummyC, dummyC, dummyD, dummyC
        phonF(iMode) = dble(freqInTHz)*THzToHartree 
        !
        do iAtom = 1, nAtoms
          !
          read (1,*) dummyC, dummyD, dummyD, dummyD, phonD(1,iAtom,iMode,iq), phonD(2,iAtom,iMode,iq), &
                     phonD(3,iAtom,iMode,iq), dummyC
          !
        enddo
        !
      enddo
      !
      read(1,*)
        !! Skip the asterisk separator
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
end module readInputFiles
