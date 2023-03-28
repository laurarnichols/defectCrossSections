module shifterMod
  
  use constants, only: dp
  use errorsAndMPI
  use mpi

  implicit none

  ! Parameters:
  integer, parameter :: root = 0
    !! ID of the root node

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

  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors

  character(len=300) :: phononFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: poscarFName
    !! File name for POSCAR


  namelist /inputParams/ poscarFName, phononFName, shift


  contains

!----------------------------------------------------------------------------
  subroutine initialize(shift, phononFName, poscarFName)
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
  subroutine checkInitialization(shift, phononFName, poscarFName)

    implicit none

    ! Input variables:
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: poscarFName
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
      write(iostd, '(" Program stops!")')
      stop
    endif
    
    return

  end subroutine checkInitialization

end module shifterMod
