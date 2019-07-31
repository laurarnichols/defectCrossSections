module readInputFiles
  !
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
    !! Used to make reals double precision
  integer, parameter :: iostd = 16
    !! Unit number for output file
  !
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  !
  contains
  !
  subroutine readPhonons(phononsInput, nOfqPoints, nAtoms, nModes, atomD, atomM, phonQ, phonF, phonD)
    !! Read the number of atoms and q points and
    !! get the phonon information like frequency 
    !! and displacements
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(out) :: nOfqPoints
      !! Number of q points
    integer, intent(out) :: nAtoms
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
    real(kind = dp), allocatable, intent(out) :: atomD(:,:)
      !! Atom displacements when comparing defective and perfect crystals
    real(kind = dp), allocatable, intent(out) :: atomM(:)
      !! Atom masses
    real(kind = dp), allocatable, intent(out) :: phonD(:,:,:,:)
    real(kind = dp), allocatable, intent(out) :: phonF(:)
    real(kind = dp), allocatable, intent(out) :: phonQ(:,:)
    real(kind = dp) :: dummyD
      !! Dummy variable to ignore input
    real(kind = dp) :: freqInTHz 
      !! Input frequency in THz
    !
    character(len = 256), intent(in) :: phononsInput
      !! QE/VASP phonon output file name
    character :: dummyC
      !! Dummy variable to ignore input
    !
    open(1, file=trim(phononsInput), status="old")
      !! * Open `phononsInput` file
    !
    read(1,*) nOfqPoints, nAtoms
      !! * Read in the number of q points and number of atoms
    !
    nModes = 3*nAtoms - 3
      !! * Calculate the number of phonon modes
    write(iostd, '(" Number of atoms : ", i5)') nAtoms
    write(iostd, '(" Number of q-Points : ", i5)') nOfqPoints
    write(iostd, '(" Number of modes : ", i5)') nModes
    flush(iostd)
      !! * Write the number of atoms, q points, and modes to the output file
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
end module readInputFiles
