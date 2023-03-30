module cell

  use constants, only: dp, pi
  use errorsAndMPI
  
  implicit none

  integer :: nAtoms
    !! Number of atoms
  integer :: nAtomTypes
    !! Number of types of atoms
  integer, allocatable :: nAtomsEachType(:)
    !! Number of atoms of each type

  real(kind=dp), allocatable :: atomPositionsDir(:,:)
    !! Atom positions in direct coordinates
  real(kind=dp) :: omega
    !! Volume of unit cell
  real(kind=dp) :: realLattVec(3,3)
    !! Real space lattice vectors
  real(kind=dp) :: recipLattVec(3,3)
    !! Reciprocal lattice vectors

  contains
      
!----------------------------------------------------------------------------
  subroutine readPOSCAR(nAtoms, poscarFName, atomPositionsDir, omega, realLattVec)

    use generalComputations, only: cart2direct

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    character(len=300), intent(in) :: poscarFName
      !! File name for POSCAR

    ! Output variables:
    real(kind=dp), intent(out) :: atomPositionsDir(3,nAtoms)
      !! Atom positions in direct coordinates
    real(kind=dp), intent(out) :: omega
      !! Volume of unit cell
    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors

    ! Local variables:
    integer :: ix, jx, ia
      !! Loop indices

    real(kind=dp) :: invLattVec(3,3)
      !! Inverse of real space lattice vectors
    real(kind=dp) :: pos(3,nAtoms)
      !! Positions as read from file
    real(kind=dp) :: scaleParam
      !! Scale parameter for lattice vectors

    logical :: fileExists
      !! If the POSCAR file exists

    character(len=1) :: charSwitch
      !! Single character for testing input
    character(len=300) :: line
      !! Line to read from file


    !! Make sure the POSCAR file exists
    inquire(file=trim(poscarFName), exist=fileExists)

    if(.not. fileExists) call exitError('readPOSCAR', 'POSCAR file does not exist', 1)

    open(unit=15, file=trim(poscarFName))

    read(15,*) 
      !! Ignore first line

    read(15,*) scaleParam
      !! Get scaling parameter. Technically, this can be 
      !! 1 or 3 items, but I am only going to write this
      !! for 1 item right now because the NITEMS subroutine
      !! in VASP is larger than I want to deal with right
      !! now. 

    !> Read real-space lattice vectors
    do ix = 1,3
      read(15,*) (realLattVec(jx,ix), jx=1,3)
    enddo

    call calculateOmega(realLattVec, omega)

    if(scaleParam < 0._dp) then
      !! User can also give desired volume as a negative number. 
      !! Scale is then set the achieve the desired volume.

      scaleParam = (abs(scaleParam)/abs(omega))**(1._dp/3._dp)

    endif

    realLattVec = realLattVec*scaleParam
    call calculateOmega(realLattVec, omega)
      !! Recalculate volume based on scaled lattice vectors

    if(omega < 0) call exitError('readPOSCARHead', 'volume is less than zero', 1)

    
    read(15,'(A)') line
      !! 6th line contains either the number of ions or their type

    read(line,'(A1)') charSwitch
    if(.not. (charSwitch >= '0' .and. charSwitch <= '9')) read(15,*) 
      !! Determine if the line contains numbers or characters
      !! (i.e., if the optional atom-name line exists) and read
      !! the next line if it does. Ignore both the atom-name 
      !! line and the atom-number line as we get this info from
      !! the vasprun.xml file because it is more structured and
      !! easier to read.

    read(15,'(A1)') charSwitch
    if(charSwitch == 'S' .or. charSwitch == 's') read(15,'(A1)') charSwitch
      !! Ignore selective dynamics line, if it exists. Read new character
      !! switch that sets if the coordinates are direct or Cartesian,
      !! but do the processing of the coordinates after reading.

    !> Read in atom positions
    do ia = 1, nAtoms

      read(15,*) (pos(ix,ia), ix=1,3)

    enddo


    if(charSwitch == 'K' .or. charSwitch == 'k' .or. &
       charSwitch == 'C' .or. charSwitch == 'c') then

      write(*,*) "Positions in Cartesian coordinates"

      pos = pos*scaleParam

      atomPositionsDir = cart2direct(nAtoms, pos, realLattVec)
        !! Convert Cartesian coordinates to direct

    else

      write(*,*) "Positions in direct coordinates"

      atomPositionsDir = pos

    endif


    close(15)
  
    return

  end subroutine readPOSCAR

!----------------------------------------------------------------------------
  subroutine calculateOmega(realLattVec, omega)
    !! Calculate the cell volume as \(a_1\cdot a_2\times a_3\)

    use generalComputations, only: vcross

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: realLattVec(3,3)
      !! Real space lattice vectors


    ! Output variables:
    real(kind=dp), intent(out) :: omega
      !! Volume of unit cell


    ! Local variables:
    real(kind=dp) :: vtmp(3)
      !! \(a_2\times a_3\)


    call vcross(realLattVec(:,2), realLattVec(:,3), vtmp)

    omega = sum(realLattVec(:,1)*vtmp(:))

    return
  end subroutine calculateOmega

!----------------------------------------------------------------------------
  subroutine getReciprocalVectors(realLattVec, omega, recipLattVec)
    !! Calculate the reciprocal lattice vectors from the real-space
    !! lattice vectors and the cell volume

    use generalComputations, only: vcross

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: omega
      !! Volume of unit cell


    ! Output variables:
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors


    ! Local variables:
    integer :: i
      !! Loop index
    

    call vcross(2.0d0*pi*realLattVec(:,2)/omega, realLattVec(:,3), recipLattVec(:,1))
      ! \(b_1 = 2\pi/\Omega a_2\times a_3\)
    call vcross(2.0d0*pi*realLattVec(:,3)/omega, realLattVec(:,1), recipLattVec(:,2))
      ! \(b_2 = 2\pi/\Omega a_3\times a_1\)
    call vcross(2.0d0*pi*realLattVec(:,1)/omega, realLattVec(:,2), recipLattVec(:,3))
      ! \(b_3 = 2\pi/\Omega a_1\times a_2\)


    return
  end subroutine getReciprocalVectors

end module cell
