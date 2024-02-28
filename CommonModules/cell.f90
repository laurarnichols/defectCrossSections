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
  function cartDispProjOnPhononEigsNorm(nAtoms, displacement, eigenvector, mass, realLattVec) result(projNorm)
    !! Project the displacement onto the phonon eigenvectors
    !! and return the norm of that projection across all atoms

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: displacement(3,nAtoms)
      !! Displacements for each atom for this mode
    real(kind=dp), intent(in) :: eigenvector(3,nAtoms)
      !! Eigenvectors for each atom for a single mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Masses of atoms
    real(kind=dp), intent(in) :: realLattVec(3,3)
      !! Real space lattice vectors

    ! Output variables:
    real(kind=dp) :: projNorm
      !! Norm of displacement projected onto phonon
      !! eigenvectors (in generalized coordinates)


    ! Local variables:
    integer :: ia
      !! Loop index

    real(kind=dp) :: eig(3)
      !! Displacement vector for single atom
    real(kind=dp) :: eigDotEig
      !! Dot product across all atoms for calculating
      !! projection
    real(kind=dp) :: genDisp(3)
      !! Displacement in generalized coordinates
    real(kind=dp) :: genDispDotEig
      !! Dot product across all atoms for calculating
      !! projection
    real(kind=dp) :: proj(3,nAtoms)
      !! Projection of generalized displacement 
      !! onto phonon eigenvector


    genDispDotEig = 0.0_dp
    eigDotEig = 0.0_dp
    do ia = 1, nAtoms

      genDisp = matmul(realLattVec,displacement(:,ia)*sqrt(mass(ia)))

      eig = eigenvector(:,ia)

      genDispDotEig = genDispDotEig + dot_product(genDisp, eig)

      eigDotEig = eigDotEig + dot_product(eig,eig)

    enddo


    proj(:,:) = genDispDotEig/eigDotEig*eigenvector(:,:)

    projNorm = 0.0_dp
    do ia = 1, nAtoms

      projNorm = projNorm + dot_product(proj(:,ia), proj(:,ia))

    enddo

    projNorm = sqrt(projNorm)

  end function cartDispProjOnPhononEigsNorm

!----------------------------------------------------------------------------
  subroutine writePOSCARNewPos(nAtoms, atomPositions, initFName, newFName, memoLine_, cart_)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: atomPositions(3,nAtoms)
      !! Atom positions

    character(len=300), intent(in) :: initFName
      !! File name for initial POSCAR
    character(len=300), intent(in) :: newFName
      !! File name for POSCAR with new positions
    character(len=300), optional :: memoLine_
      !! Input memo line to add to first note line

    logical, optional :: cart_
      !! Input for whether or not to write coordinates in Cartesian

    ! Local variables:
    integer :: ix, ia, iLine
      !! Loop indices

    character(len=300) :: memoLine
      !! Memo line to add to first note line

    logical :: cart
      !! Whether or not to write coordinates in Cartesian
    logical :: fileExists
      !! If the POSCAR file exists

    character(len=1) :: charSwitch
      !! Single character for testing input
    character(len=300) :: line
      !! Line to read from file



    memoLine = ''
    if(present(memoLine_)) memoLine = memoLine_
      !! Set optional argument memoLine with a default empty string

    cart = .false.
    if(present(cart_)) cart = cart_
      !! Set optional argument cart with default value false


    !> Make sure the POSCAR file exists
    inquire(file=trim(initFName), exist=fileExists)

    if(.not. fileExists) call exitError('writePOSCARNewPos', 'POSCAR file does not exist', 1)

    open(unit=15, file=trim(initFName))
    open(unit=20, file=trim(newFName))


    read(15,'(A)') line 
    write(20,'(A)') trim(line)//trim(memoLine)
      !! Add optional memoLine to first comment line
  
    !> Write next 5 lines as-is
    do iLine = 2, 6
      read(15,'(A)') line 
      write(20,'(A)') trim(line) 
    enddo


    read(line,'(A1)') charSwitch
    if(.not. (charSwitch >= '0' .and. charSwitch <= '9')) then
      !! Determine if the 6th line contains numbers or characters
      !! (i.e., if the optional atom-name line exists) and read
      !! the next line if it does. Ignore both the atom-name 
      !! line and the atom-number line as we get this info from
      !! the vasprun.xml file because it is more structured and
      !! easier to read.

      read(15,'(A)') line
      write(20,'(A)') trim(line)
    
    endif

    
    if(cart) then
      write(20,'("Cartesian")')
    else
      write(20,'("Direct")')
    endif

    do ia = 1, nAtoms

      write(20,'(3f20.16)') (atomPositions(ix,ia), ix=1,3) 

    enddo

    close(15)
    close(20)
  
    return

  end subroutine writePOSCARNewPOS
      
!----------------------------------------------------------------------------
  subroutine readPOSCAR(poscarFName, nAtoms, atomPositionsDir, omega, realLattVec)

    use generalComputations, only: cart2direct

    implicit none

    ! Input variables:
    character(len=300), intent(in) :: poscarFName
      !! File name for POSCAR

    ! Output variables:
    integer, intent(out) :: nAtoms
      !! Number of atoms

    real(kind=dp), allocatable, intent(out) :: atomPositionsDir(:,:)
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
    real(kind=dp) :: pos(3,7000)
      !! Positions as read from file
    real(kind=dp) :: scaleParam
      !! Scale parameter for lattice vectors

    logical :: fileExists
      !! If the POSCAR file exists

    character(len=1) :: charSwitch
      !! Single character for testing input
    character(len=300) :: line
      !! Line to read from file


    !> Make sure the POSCAR file exists
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

    if(omega < 0) call exitError('readPOSCAR', 'volume is less than zero', 1)

    
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
    ia = 0
    do while(.true.)

      read(15,'(a)') line
      if(trim(line) == '') goto 200
        ! Read line to see if it's blank, which would
        ! indicate the end of the positions section

      ia = ia + 1
      read(line,*,end=200) (pos(ix,ia), ix=1,3)
        ! If the line isn't blank, increment the atom
        ! number and read the position for this atom 

    enddo

200 continue
    nAtoms = ia

    allocate(atomPositionsDir(3,nAtoms))

    if(charSwitch == 'K' .or. charSwitch == 'k' .or. &
       charSwitch == 'C' .or. charSwitch == 'c') then

      pos = pos*scaleParam

      atomPositionsDir = cart2direct(nAtoms, pos(:,1:nAtoms), realLattVec)
        !! Convert Cartesian coordinates to direct

    else

      atomPositionsDir = pos(:,1:nAtoms)

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
