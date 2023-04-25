module shifterMod
  
  use constants, only: dp, angToBohr, daltonToElecM
  use cell
  use errorsAndMPI
  use mpi

  implicit none

  ! Variables that should not be passed as arguments:
  integer :: iModeStart, iModeEnd
    !! Start and end mode for this process

  real(kind=dp) :: t0, t1, t2
    !! Timers

  ! Variables that should be passed as arguments:
  integer :: nModes
    !! Number of phonon modes
  integer :: nModesLocal
    !! Local number of phonon modes
  integer :: suffixLength
    !! Length of shifted POSCAR file suffix

  real(kind=dp), allocatable :: displacement(:,:)
    !! Atom displacements in angstrom
  real(kind=dp), allocatable :: eigenvector(:,:,:)
    !! Eigenvectors for each atom for each mode
  real(kind=dp), allocatable :: generalizedNorm(:)
    !! Generalized norms after displacement
  real(kind=dp), allocatable :: mass(:)
    !! Masses of atoms
  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors
  real(kind=dp), allocatable :: shiftedPositions(:,:)
    !! Positions after shift along eigenvector

  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: phononFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: poscarFName
    !! File name for POSCAR
  character(len=300) :: prefix
    !! Prefix for shifted POSCARs
  character(len=300) :: shiftedPOSCARFName
    !! File name for shifted POSCAR


  namelist /inputParams/ poscarFName, phononFName, prefix, nAtoms, shift, dqFName


  contains

!----------------------------------------------------------------------------
  subroutine initialize(nAtoms, shift, dqFName, phononFName, poscarFName, prefix)
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    ! Input variables:
    !integer, intent(in) :: nProcs
      ! Number of processes


    ! Output variables:
    integer, intent(out) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: poscarFName
      !! File name for POSCAR
    character(len=300), intent(out) :: prefix
      !! Prefix for shifted POSCARs


    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time

    nAtoms = -1
    dqFName = 'dq.txt'
    poscarFName = 'POSCAR'
    phononFName = 'mesh.yaml'
    prefix = 'ph_POSCAR'
    shift = 0.01_dp

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(*, '(/5X,"Phonon post-processing: POSCAR shifter starts on ",A9," at ",A9)') &
             cdate, ctime

      write(*, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs


    endif

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(nAtoms, shift, dqFName, phononFName, poscarFName, prefix)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(inout) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(in) :: poscarFName
      !! File name for POSCAR
    character(len=300), intent(in) :: prefix
      !! Prefix for shifted POSCARs

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkFileInitialization('poscarFName', poscarFName)
    abortExecution = checkFileInitialization('phononFName', phononFName) .or. abortExecution
    abortExecution = checkStringInitialization('prefix', prefix) .or. abortExecution
    abortExecution = checkIntInitialization('nAtoms', nAtoms, 0, 5000) .or. abortExecution
    abortExecution = checkDoubleInitialization('shift', shift, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkStringInitialization('dqFName', dqFName) .or. abortExecution


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif
    
    return

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readPhonons(nAtoms, nModes, phononFName, eigenvector, mass)

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

    ! Local variables:
    integer :: j, ia, ix
      !! Loop index

    real(kind=dp) :: qPos(3)
      !! Phonon q position

    character(len=300) :: line
      !! Line read from file

    logical :: fileExists
      !! If phonon file exists


    inquire(file=phononFName, exist=fileExists)

    if(.not. fileExists) call exitError('readPhonons', 'Phonon file '//trim(phononFName)//' does not exist', 1)

    open(57, file=phononFName)

    line = getFirstLineWithKeyword(57,'points')
      !! Ignore everything next until you get to points line

    do ia = 1, nAtoms

      read(57,'(A)') ! Ignore symbol 
      read(57,'(A)') ! Ignore coordinates
      read(57,'(a7,f)') line, mass(ia)
        !! Read mass

    enddo
    line = getFirstLineWithKeyword(57,'q-position')
      !! Ignore everything next until you get to q-position line

    read(line(16:len(trim(line))-1),*) qPos
      !! Read in the q position

    if(qPos(1) > 1e-8 .or. qPos(2) > 1e-8 .or. qPos(3) > 1e-8) &
      call exitError('readPhonons', 'Code assumes phonons have no momentum (at q=0)!',1)

    read(57,'(A)') line
    read(57,'(A)') line
    read(57,'(A)') line
      !! Ignore next 3 lines

    do j = 1, nModes+3

      read(57,'(A)') ! Ignore mode number
      read(57,'(A)') ! Ignore frequency
      read(57,'(A)') ! Ignore eigenvector section header

      do ia = 1, nAtoms

        read(57,'(A)') line ! Ignore atom number

        do ix = 1, 3

          read(57,'(A)') line
          if(j > 3) read(line(10:len(trim(line))-1),*) eigenvector(ix,ia,j-3)
            !! Only store the eigenvectors after the first 3 modes 
            !! because those are the acoustic modes and we only want
            !! the optical modes.

        enddo

      enddo
      
    enddo

    close(57)

    return

  end subroutine readPhonons

!----------------------------------------------------------------------------
  subroutine getDisplacement(j, nAtoms, nModes, eigenvector, mass, shift, displacement, generalizedNorm_j)

    implicit none

    ! Input variables:
    integer, intent(in) :: j
      !! Phonon mode index
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Masses of atoms
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors

    ! Output variables:
    real(kind=dp), intent(out) :: displacement(3,nAtoms)
      !! Displacements for each atom for this mode
    real(kind=dp), intent(out) :: generalizedNorm_j
      !! Norm of eigenvectors in generalized coordinates
      !! after being scaled by `shift` in Cartesian space
      !! for a given mode j

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

      eig = eigenvector(:,ia,j)/sqrt(mass(ia))

      displacement(:,ia) = eig*shift

      cartNorm = cartNorm + dot_product(eig,eig)

    enddo

    cartNorm = sqrt(cartNorm)

    displacement = displacement/cartNorm
      !! @note
      !!  Input positions and shift magnitude are in
      !!  angstrom, and the output positions are to
      !!  POSCARs for VASP calculations, so the units
      !!  should still be angstrom.
      !! @endnote

    generalizedNorm_j = cartDisplacementToGeneralizedNorm(nAtoms, displacement, mass)*angToBohr*sqrt(daltonToElecM)
      !! Convert scaled displacement back to generalized
      !! coordinates and get norm
      !! @note
      !!   Input positions are in angstrom and input
      !!   masses are in amu, but the dq output is going
      !!   to our code, which uses Hartree atomic units
      !!   (Bohr and electron masses), so this value
      !!   must have a unit conversion.
      !! @endnote

    return

  end subroutine getDisplacement

!----------------------------------------------------------------------------
  subroutine writeDqs(nModes, generalizedNorm, dqFName)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(inout) :: generalizedNorm(nModes)
      !! Generalized norms after displacement

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms

    ! Local variables:
    integer :: j
      !! Loop index


    call mpiSumDoubleV(generalizedNorm, worldComm)
      !! * Get the generalized-displacement norms
      !!   from all processes

    if(ionode) then

      open(60, file=dqFName)

      write(60,'("# Norm of generalized displacement vectors after scaling Cartesian displacement Format: ''(1i7, 1ES24.15E3)''")')

      do j = 1, nModes

        write(60,'(1i7, 1ES24.15E3)') j, generalizedNorm(j)
          !! Output norm to file

      enddo

      close(60)

    endif
      
    return

  end subroutine writeDqs

end module shifterMod
