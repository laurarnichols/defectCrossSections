module PhononPPMod
  
  use constants, only: dp, angToBohr, angToM, daltonToElecM, elecMToKg, THzToHz, pi, hbar
  use cell, only: nAtoms, omega, realLattVec
  use errorsAndMPI
  use mpi

  implicit none

  ! Variables that should not be passed as arguments:
  integer :: iModeStart, iModeEnd
    !! Start and end mode for this process

  real(kind=dp) :: t0, t1, t2
    !! Timers

  ! Variables that should be passed as arguments:
  integer :: nAtomsFinal
    !! Number of atoms from final-state POSCAR
  integer :: nModes
    !! Number of phonon modes
  integer :: nModesLocal
    !! Local number of phonon modes
  integer :: suffixLength
    !! Length of shifted POSCAR file suffix

  real(kind=dp), allocatable :: atomPositionsDirInit(:,:), atomPositionsDirFinal(:,:)
    !! Atom positions in initial and final relaxed positions
  real(kind=dp), allocatable :: displacement(:,:)
    !! Atom displacements in angstrom
  real(kind=dp), allocatable :: eigenvector(:,:,:)
    !! Eigenvectors for each atom for each mode
  real(kind=dp), allocatable :: projNorm(:)
    !! Generalized norms after displacement
  real(kind=dp), allocatable :: mass(:)
    !! Masses of atoms
  real(kind=dp) :: omegaFinal
    !! Volume of final-state supercell
  real(kind=dp), allocatable :: omegaFreq(:)
    !! Frequency for each mode
  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors
  real(kind=dp), allocatable :: shiftedPositions(:,:)
    !! Positions after shift along eigenvector

  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: phononFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: initPOSCARFName, finalPOSCARFName
    !! File name for POSCAR for relaxed initial and final charge states
  character(len=300) :: prefix
    !! Prefix for shifted POSCARs
  character(len=300) :: shiftedPOSCARFName
    !! File name for shifted POSCAR

  logical :: generateShiftedPOSCARs
    !! If shifted POSCARs should be generated


  namelist /inputParams/ initPOSCARFName, finalPOSCARFName, phononFName, prefix, nAtoms, shift, dqFName, generateShiftedPOSCARs


  contains

!----------------------------------------------------------------------------
  subroutine initialize(shift, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs)
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    ! Output variables:
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states
    character(len=300), intent(out) :: prefix
      !! Prefix for shifted POSCARs

    logical, intent(out) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated


    dqFName = 'dq.txt'
    initPOSCARFName = 'POSCAR_init'
    finalPOSCARFName = 'POSCAR_final'
    phononFName = 'mesh.yaml'
    prefix = 'ph_POSCAR'

    shift = 0.01_dp

    generateShiftedPOSCARs = .true.

    return

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(shift, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, generateShiftedPOSCARs)

    implicit none

    ! Input variables:
    real(kind=dp), intent(inout) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(in) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states
    character(len=300), intent(in) :: prefix
      !! Prefix for shifted POSCARs

    logical, intent(in) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkFileInitialization('initPOSCARFName', initPOSCARFName)
    abortExecution = checkFileInitialization('finalPOSCARFName', finalPOSCARFName) .or. abortExecution
    abortExecution = checkFileInitialization('phononFName', phononFName) .or. abortExecution
    abortExecution = checkStringInitialization('prefix', prefix) .or. abortExecution
    abortExecution = checkDoubleInitialization('shift', shift, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkStringInitialization('dqFName', dqFName) .or. abortExecution


    write(*,'("generateShiftedPOSCARs = ''",L1,"''")') generateShiftedPOSCARs


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif
    
    return

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine standardizeCoordinates(nAtoms, atomPositionsDir)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms in system

    ! Output variables:
    real(kind=dp), intent(inout) :: atomPositionsDir(3,nAtoms)
      !! Atom positions in direct coordinates

    ! Local variables:
    integer :: ia, ix
      !! Loop indices


    where(atomPositionsDir(:,:) < 0)
      atomPositionsDir(:,:) = 1.0_dp + atomPositionsDir(:,:)
    end where

    return

  end subroutine standardizeCoordinates

!----------------------------------------------------------------------------
  subroutine checkCompatibility(nAtoms1, nAtoms2, omega1, omega2, atomPositionsDir1, atomPositionsDir2)

    use miscUtilities, only: int2str

    implicit none

    ! Input variables
    integer, intent(in) :: nAtoms1, nAtoms2
      !! Number of atoms

    real(kind=dp), intent(in) :: omega1, omega2
      !! Cell volume

    ! Output variables:
    real(kind=dp), intent(inout) :: atomPositionsDir1(3,nAtoms1), atomPositionsDir2(3,nAtoms2)
      !! Atom positions

    ! Local variables:
    integer :: ia, ix
      !! Loop indices

    real(kind=dp) :: dispDir, dispDir1, dispDir2
      !! Displacement in direct coordinates

    logical :: abortExecution
      !! Whether or not to abort the execution


    if(nAtoms1 /= nAtoms2) &
      call exitError('checkCompatibility', 'number of atoms does not match: '//trim(int2str(nAtoms1))//' '//trim(int2str(nAtoms2)), 1)

    if(abs(omega1 - omega2) > 1e-8) call exitError('checkCompatibility', 'volumes don''t match', 1)


    abortExecution = .false.
    
    do ia = 1, nAtoms
      do ix = 1, 3

        dispDir1 = atomPositionsDirFinal(ix,ia) - atomPositionsDirInit(ix,ia)
        dispDir2 = 1.0_dp - atomPositionsDirFinal(ix,ia) - atomPositionsDirInit(ix,ia)

        if(abs(dispDir1) <= abs(dispDir2)) then
          dispDir = dispDir1
        else
          dispDir = dispDir2
          atomPositionsDirFinal(ix,ia) = 1.0_dp - atomPositionsDirFinal(ix,ia)
        endif

        if(abs(dispDir) > 0.2) then
          abortExecution = .true.

          write(*,'("Possible atoms out of order: ", i4, 3f10.5)') ia, dispDir
        endif

      enddo
    enddo

    if(abortExecution) call exitError('checkCompatibility', 'atoms don''t seem to be in the same order', 1)

    return

  end subroutine checkCompatibility

!----------------------------------------------------------------------------
  subroutine readPhonons(nAtoms, nModes, phononFName, eigenvector, mass, omegaFreq)

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: nModes
      !! Number of modes

    character(len=300), intent(in) :: phononFName
      !! File name for mesh.yaml phonon file

    ! Output variables:
    real(kind=dp), intent(out) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), intent(out) :: mass(nAtoms)
      !! Masses of atoms
    real(kind=dp), intent(out) :: omegaFreq(nModes)
      !! Frequency for each mode

    ! Local variables:
    integer :: j, ia, ix
      !! Loop index
    integer :: nAtomsFromPhon
      !! Number of atoms from phonon file

    real(kind=dp), allocatable :: coordFromPhon(:,:)
      !! Corodinates from phonon file
    real(kind=dp) :: qPos(3)
      !! Phonon q position

    character(len=300) :: line
      !! Line read from file

    logical :: fileExists
      !! If phonon file exists


    inquire(file=phononFName, exist=fileExists)

    if(.not. fileExists) call exitError('readPhonons', 'Phonon file '//trim(phononFName)//' does not exist', 1)

    open(57, file=phononFName)

    line = getFirstLineWithKeyword(57,'natom')
    read(line(7:len(trim(line))),*) nAtomsFromPhon

    line = getFirstLineWithKeyword(57,'points')
      !! Ignore everything next until you get to points line

    allocate(coordFromPhon(3,nAtomsFromPhon))

    do ia = 1, nAtomsFromPhon

      read(57,'(A)') ! Ignore symbol 
      read(57,'(A)') line
      read(line(17:len(trim(line))),*) coordFromPhon(:,ia)
      read(57,'(a7,f)') line, mass(ia)
        !! Read mass

    enddo

    call checkCompatibility(nAtoms, nAtomsFromPhon, omega, omega, atomPositionsDirInit, coordFromPhon)

    deallocate(coordFromPhon)

    line = getFirstLineWithKeyword(57,'q-position')
      !! Ignore everything next until you get to q-position line

    read(line(16:len(trim(line))),*) qPos
      !! Read in the q position

    if(qPos(1) > 1e-8 .or. qPos(2) > 1e-8 .or. qPos(3) > 1e-8) &
      call exitError('readPhonons', 'Code assumes phonons have no momentum (at q=0)!',1)

    read(57,'(A)') line
    read(57,'(A)') line
    read(57,'(A)') line
      !! Ignore next 3 lines

    do j = 1, nModes+3

      read(57,'(A)') ! Ignore mode number

      read(57,'(A)') line
      if(j > 3) read(line(15:len(trim(line))),*) omegaFreq(j-3)

      read(57,'(A)') ! Ignore eigenvector section header

      do ia = 1, nAtoms

        read(57,'(A)') line ! Ignore atom number

        do ix = 1, 3

          read(57,'(A)') line
          if(j > 3) read(line(10:len(trim(line))),*) eigenvector(ix,ia,j-3)
            !! Only store the eigenvectors after the first 3 modes 
            !! because those are the acoustic modes and we only want
            !! the optical modes.

        enddo

      enddo
      
    enddo

    omegaFreq(:) = omegaFreq(:)*2*pi

    close(57)

    return

  end subroutine readPhonons
  
!----------------------------------------------------------------------------
  subroutine calcAndWriteSj(nModes, omegaFreq, projNorm)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: omegaFreq(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(inout) :: projNorm(nModes)
      !! Generalized norms after displacement

    ! Local variables:
    integer :: j
      !! Loop index


    call mpiSumDoubleV(projNorm, worldComm)
      !! * Get the generalized-displacement norms
      !!   from all processes

    if(ionode) then

      open(60, file="./Sj.out")

      write(60,'(1i7)') nModes

      do j = 1, nModes

        write(60,'(1i7, 2ES24.15E3)') j, projNorm(j)**2*omegaFreq(j)*THzToHz/(2*hbar), omegaFreq(j)

      enddo

      close(60)

    endif


    return

  end subroutine calcAndWriteSj

!----------------------------------------------------------------------------
  function getShiftDisplacement(nAtoms, eigenvector, mass, shift) result(displacement)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: eigenvector(3,nAtoms)
      !! Eigenvectors for each atom for this mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Masses of atoms
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors

    ! Output variables:
    real(kind=dp) :: displacement(3,nAtoms)
      !! Displacements for each atom for this mode

    ! Local variables:
    integer :: ia
      !! Loop indices

    real(kind=dp) :: eig(3)
      !! Eigenvector for single mode and atom
    real(kind=dp) :: cartNorm
      !! Norm of eigenvectors in Cartesian coordinates


    !> Convert from generalized to Cartesian coordinates
    !> and scale displacement based on norm of entire
    !> displacement vector
    cartNorm = 0.0_dp
    do ia = 1, nAtoms

      eig = eigenvector(:,ia)/sqrt(mass(ia))

      displacement(:,ia) = eig

      cartNorm = cartNorm + dot_product(eig,eig)

    enddo

    cartNorm = sqrt(cartNorm)

    displacement = displacement*shift/cartNorm
      !! Scale the displacement so that the 
      !! norm of the entire vector is `shift`

    return

  end function getShiftDisplacement

!----------------------------------------------------------------------------
  subroutine writeDqs(nModes, projNorm, dqFName)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(inout) :: projNorm(nModes)
      !! Generalized norms after displacement

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms

    ! Local variables:
    integer :: j
      !! Loop index


    call mpiSumDoubleV(projNorm, worldComm)
      !! * Get the generalized-displacement norms
      !!   from all processes

    if(ionode) then

      open(60, file=dqFName)

      write(60,'("# Norm of generalized displacement vectors after scaling Cartesian displacement Format: ''(1i7, 1ES24.15E3)''")')

      do j = 1, nModes

        write(60,'(1i7, 1ES24.15E3)') j, projNorm(j)
          !! Output norm to file

      enddo

      close(60)

    endif
      
    return

  end subroutine writeDqs

end module PhononPPMod
