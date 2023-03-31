module shifterMod
  
  use constants, only: dp
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

  real(kind=dp), allocatable :: eigenvector(:,:,:)
    !! Eigenvectors for each atom for each mode
  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors
  real(kind=dp), allocatable :: shiftedPositions(:,:)
    !! Positions after shift along eigenvector

  character(len=300) :: phononFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: poscarFName
    !! File name for POSCAR


  namelist /inputParams/ poscarFName, phononFName, nAtoms, shift


  contains

!----------------------------------------------------------------------------
  subroutine initialize(nAtoms, shift, phononFName, poscarFName)
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

    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: poscarFName
      !! File name for POSCAR


    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time

    nAtoms = -1
    poscarFName = 'POSCAR'
    phononFName = 'mesh.yaml'
    shift = 0.01_dp

    call date_and_time(cdate, ctime)

    if(ionode) then

      write(*, '(/5X,"Phonon post-processing: POSCAR shifter starts on ",A9," at ",A9)') &
             cdate, ctime

      write(*, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nProcs


    endif

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(nAtoms, shift, phononFName, poscarFName)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(inout) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(in) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(in) :: poscarFName
      !! File name for POSCAR

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution
    logical :: fileExists
      !! If a file exists


    if(trim(poscarFName) == '' ) then

      write(*,*)
      write(*,'(" Variable : ""poscarFName"" is not defined!")')
      write(*,'(" usage : poscarFName = ''./POSCAR''")')
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    else
      
      inquire(file=trim(poscarFName), exist=fileExists)
      
      if(fileExists .eqv. .false.) then

        write(*,'(" File : ", a, " , does not exist!")') trim(poscarFName)
        write(*,'(" This variable is mandatory and thus the program will not be executed!")')

        abortExecution = .true.

      endif
    endif
    

    write(*, '("poscarFName = ''", a, "''")') trim(poscarFName)


    if(trim(phononFName) == '' ) then

      write(*,*)
      write(*,'(" Variable : ""phononFName"" is not defined!")')
      write(*,'(" usage : phononFName = ''./mesh.yaml''")')
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    else
      
      inquire(file=trim(phononFName), exist=fileExists)
      
      if(fileExists .eqv. .false.) then

        write(*,'(" File : ", a, " , does not exist!")') trim(phononFName)
        write(*,'(" This variable is mandatory and thus the program will not be executed!")')

        abortExecution = .true.

      endif
    endif


    write(*, '("phononFName = ''", a, "''")') trim(phononFName)

    if(nAtoms < 0) then

      write(*,*)
      write(*,'(" Variable : ""nAtoms"" is not defined!")')
      write(*,'(" usage : nAtoms = 100")')
      write(*,'(" This variable is mandatory and thus the program will not be executed!")')

      abortExecution = .true.

    endif


    if(shift < 0.0_dp ) then

      shift = 0.01_dp ! Angstrom?

      write(*,'(" Variable : ""shift"" is less than zero!")')
      write(*,'(" usage : shift = 0.01")')
      write(*,'(" A default value of 0.01 A will be used !")')

    else if(shift > 0.2_dp) then
      ! This limit is arbitrary. Seems like a good number for right now.

      write(*,'(" Variable : ""shift"" is too large!")')
      write(*,'(" usage : shift = 0.01")')
      write(*,'(" Re-run with a value less than 0.2 A!")')

      abortExecution = .true.

    endif
    
    write(*, '("shift = ", f8.4, " (A)")') shift
    

    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif
    
    return

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readPhonons(nAtoms, nModes, phononFName, eigenvector)

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

    line = getFirstLineWithKeyword(57,'q-position')
      !! Ignore everything until you get to q-position line

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
  function getDisplacement(j, nAtoms, nModes, eigenvector, shift) result(displacement)

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
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors

    ! Output variables:
    real(kind=dp) :: displacement(3,nAtoms)
      !! Displacements for each atom for this mode
    real(kind=dp) :: eig(3)
      !! Eigenvector for single mode and atom

    ! Local variables:
    integer :: ia
      !! Loop indices


    do ia = 1, nAtoms

      eig = eigenvector(:,ia,j)

      displacement(:,ia) = eig*shift/sqrt(dot_product(eig, eig))

    enddo

  end function getDisplacement

end module shifterMod
