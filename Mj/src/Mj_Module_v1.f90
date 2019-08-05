
module MjModule
  !
  use constants
  !
  implicit none
  !
  integer, parameter :: int32 = selected_int_kind(5)
  integer, parameter :: int64 = selected_int_kind(15)
  integer, parameter :: un = 3
  !
  character(len = 6), parameter :: output = 'status'
  !
  !
  integer :: ios
  integer :: modeF
  integer :: modeI
  integer :: nAtoms
    !! Number of atoms
  integer :: nModes
    !! Number of phonon modes
  integer :: nOfqPoints
  integer :: qPoint
    !! Read in from input file, but no default value set
    !! @todo Make sure default value is set for `qPoint` @endtodo
  !
  real(kind = dp) :: kT
  real(kind = dp) :: maxDisplacement
  real(kind = dp) :: t1
  real(kind = dp) :: t2
  real(kind = dp) :: temperature
  real(kind = dp) :: tf
  real(kind = dp) :: ti
  !
  character(len = 256) :: equilibriumAtomicPositions
  character(len = 256) :: newAtomicPositions
  character(len = 256) :: phononsInput
  character(len = 256) :: QEInput
  !
  logical :: file_exists
  logical :: readQEInput
  !
  !
  integer, allocatable :: s2L(:)
  !
  real(kind = dp), allocatable :: atomD(:,:)
  real(kind = dp), allocatable :: atomM(:)
  real(kind = dp), allocatable :: atomPosition(:,:)
  real(kind = dp), allocatable :: besOrderNofModeM(:,:)
  real(kind = dp), allocatable :: coth(:)
  real(kind = dp), allocatable :: genCoord(:)
  real(kind = dp), allocatable :: newAtomicPosition(:,:)
  real(kind = dp), allocatable :: phonD(:,:,:,:)
  real(kind = dp), allocatable :: phonF(:)
  real(kind = dp), allocatable :: phonQ(:,:)
  real(kind = dp), allocatable :: Sj(:)
  real(kind = dp), allocatable :: wby2kT(:)
  real(kind = dp), allocatable :: x(:)
  !
  character(len = 2), allocatable :: elements(:)
  !
  namelist /MjInput/ QEInput, phononsInput, temperature, equilibriumAtomicPositions, modeI, modeF, qPoint, maxDisplacement
  !
  !
contains
  !
  !
  subroutine readInputs()
    !! Read input parameters and read phonon output
    !!
    !! <h2>Walkthrough</h2>
    !!
    use readInputFiles
    !
    implicit none
    !
    !> * Check if file output exists; if it does, delete it
    inquire(file = output, exist = file_exists)
    !
    if ( file_exists ) then
      !
      open (unit = 11, file = output, status = "old")
      !
      close(unit = 11, status = "delete")
      !
    endif
    !
    open (iostd, file = output, status='new')
      !! * Open new output file
    !
    call initialize()
      !! * Set default values of input parameters
    !
    READ (5, MjInput, iostat = ios)
      !! * Read input parameters
    !
    call checkAndUpdateInput()
      !! * Check if input parameters were updated 
      !!   and do some basic checks
    !
    call readPhonons(phononsInput, nOfqPoints, nAtoms, nModes, atomD, atomM, phonQ, phonF, phonD)
      !! * Read the phonons output from QE or VASP
    !
    call readAtomicPositions()
      !! * Read the equilibrium atomic positions
    !
    return
    !
  end subroutine readInputs
  !
  !
  subroutine initialize()
    !! Set default values for input parameters
    !
    implicit none
    !
    QEInput = ''
    phononsInput = ''
    equilibriumAtomicPositions = ''
    temperature = -1.0_dp
    maxDisplacement = -1.0_dp
    modeI = -1
    modeF = -1
    !
    return
    !
  end subroutine initialize
  !
  !
  subroutine checkAndUpdateInput()
    !! Check that the input variables don't still have their default
    !! values. The program will abort here if:
    !! * `equilibriumAtomicPositions` is not defined
    !! * `phononsInput` is not defined
    !! * `temperature` is not defined
    !! * `modeI` or `modeF` is not defined
    !! * `modeF < modeI`
    !! * `maxDisplacement` is not defined
    !
    implicit none
    !
    logical :: abortExecution = .false.
    !
    write(iostd, *)
    !
    if ( equilibriumAtomicPositions == '' ) then
      write(iostd, '(" equilibriumAtomicPositions is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Equilibrium Atomic Positions input : ", a)') trim(equilibriumAtomicPositions)
    endif
    !
    if ( phononsInput == '' ) then
      write(iostd, '(" PhononsInput is not defined!")')
      abortExecution = .true.
    else
      write(iostd, '(" Phonons input : ", a)') trim(PhononsInput)
    endif
    !
    if ( QEInput == '' ) then
      write(iostd, '(" QEInput is not defined!")')
      readQEInput = .false.
    else
      readQEInput = .true.
      write(iostd, '(" Quantum Espresso input : ", a)') trim(QEInput)
    endif
    !
    if ( temperature < 0.0_dp ) then
      write(iostd, '(" Variable temperature has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Temperature : ", f10.2, " Kelvin.")') temperature
      kT = temperature*8.6173324d-5*eVToHartree
    endif
    !
    if ( modeI < 0 ) then
      write(iostd, '(" Variable modeI has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Initial mode : ", i5)') modeI 
    endif
    !   
    if ( modeF < 0 ) then
      write(iostd, '(" Variable modeF has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Final mode : ", i5)') modeF
    endif
    !
    if ( modeF < modeI ) then
      write(iostd, '(" Final mode is set smaller than initial one!")')
      abortExecution = .true.
    endif
    !
    if ( maxDisplacement < 0 ) then
      write(iostd, '(" Variable maxDisplacement has not been set.")')
      abortExecution = .true.
    else
      write(iostd, '(" Maximum atomic displacement in each direction : ", f15.10)') maxDisplacement
    endif
    !
    if ( abortExecution ) then
      write(iostd, '(" *************************** ")')
      write(iostd, '(" * Program stops!          * ")')
      write(iostd, '(" *************************** ")')
      stop
    endif
    !
    return
    !
  end subroutine checkAndUpdateInput
  !
  subroutine readAtomicPositions()
    !! Read in the element and equilibrium position for 
    !! each atom
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer :: iAtom
      !! Loop index over atoms
    !
    open(1, file=trim(equilibriumAtomicPositions), status="old")
      !! * Open the `equilibriumAtomicPositions` file
    !
    allocate( elements(nAtoms), atomPosition(3,nAtoms) )
    !
    atomPosition(:,:) = 0.0_dp
    !
    do iAtom = 1, nAtoms
      !! * For each atom, read in the element and equilibrium position
      !
      read(1,*) elements(iAtom), atomPosition(1,iAtom), atomPosition(2,iAtom), atomPosition(3,iAtom)
      !
    enddo
    !
    close(1)
      !! * Close the `equilibriumAtomicPositions` file
    !
    return
    !
  end subroutine readAtomicPositions
  !
  !
  subroutine displaceAtoms()
    !! For each mode, generate random displacements for the atoms
    !! based on the parameters `maxDisplacement` and `phonD`
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: istat, iAtom, iMode, iRand
    real(kind = dp) :: ran, norm
    !
    allocate ( newAtomicPosition(3, nAtoms) )
    !
    open(unit=12, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
    !
    write(iostd, *)
    !
    do iMode = modeI, modeF
      !! * For each mode, generate random displacements for the atoms 
      !!   based on the parameters `maxDisplacement` and `phonD`
      !
      write(iostd, '(" Calculating new atomic positions for mode :", i10)') s2L(iMode)
        !! @todo Figure out if expect `modeI` and `modeF` to represent index of magnitude of argument `x` @endtodo
      !
      do iAtom = 1, nAtoms
        !
        norm = sqrt(sum(phonD(:,iAtom,s2L(iMode),qPoint)**2)) 
        !
        read(12) iRand
        ran = mod(abs(iRand),10000000)/1.0e7_dp
        !
        newAtomicPosition(:,iAtom) = atomPosition(:,iAtom) + maxDisplacement*ran*phonD(:,iAtom,s2L(iMode),qPoint)/norm
        !
      enddo
      !
    enddo
    !
    close(12)
    !
  end subroutine displaceAtoms
  !
  !
  subroutine writeNewAtomicPositions()
    !! Write `newAtomicPosition`s for each mode
    !
    use miscUtilities
    !
    implicit none
    !
    integer :: iAtom, iMode
    !
    character(len = 300) :: s2LStr
    !
    do iMode = modeI, modeF
      !
      write(iostd, '(" Writing new atomic positions for mode :", i10)') s2L(iMode)
      !
      if ( s2L(iMode) < 10000 ) then
        call int2str(s2L(iMode), s2LStr)
        !
        write(newAtomicPositions, '("newPositionForMode", a)') trim(s2LStr)
        !
      else
        !
        newAtomicPositions = 'newPositions'
        !
      endif
      !
      open(21, file=trim(newAtomicPositions), status='unknown')
      !
      do iAtom = 1, nAtoms
        !
        !write(6, '(i4, f15.12, 3f15.5, " | ", 3f15.5)') iAtom, maxDisplacement*ran, &
        !          atomPosition(:,iAtom), atomPosition(:,iAtom) + maxDisplacement*ran*phonD(:,iAtom,s2L(iMode),qPoint)/norm
        write(21,*) elements(iAtom), newAtomicPosition(:,iAtom)
      enddo
      !
      close(21)
      !
    enddo
    !
  end subroutine writeNewAtomicPositions
  !
  !
  subroutine exportQEInput()
    !! Create QE input files for all different modes
    !! by copying all of the `QEInput` except the 
    !! `newAtomicPosition`s for each mode
    !!
    !! <h2>Walkthrough</h2>
    !!
    use miscUtilities
    !
    implicit none
    !
    integer :: iAtom
      !! Loop index over atoms
    integer :: iMode
      !! Loop index over phonon modes
    character(len = 300) :: line, fn, modeFolder, mkDir, s2LStr
    !
    do iMode = modeI, modeF
      !! * For each mode between `modeI` and `modeF`
      !!   * If the folder for this mode doesn't already 
      !!     exist, make it
      !!   * Copy the `QEInput` file into the new folder, 
      !!     changing the positions to be the `newAtomicPosition`
      !
      call int2str(s2L(iMode), s2LStr)
      !
      write(modeFolder, '("mode_", a)') trim(s2LStr)
      !
      inquire(file= trim(modeFolder), exist = file_exists)
      if ( .not.file_exists ) then
        !
        write(mkDir, '("mkdir -p ", a)') trim(modeFolder)
        call system(mkDir)
        !
      endif
      !
      fn = trim(QEInput)
      fn = fn(INDEX(QEInput, '/', BACK=.TRUE.):INDEX(QEInput, '.in')-1)
      !
      write(iostd, '(" Writing new QE input file for mode :", i10)') s2L(iMode)
      !
      write(fn, '(a, "_mode", a, ".in")') trim(fn), trim(s2LStr)
      !
      fn = trim(modeFolder)//"/"//trim(fn)
      !
      open(2, file=trim(fn), status="unknown")
      !
      open(1, file=trim(QEInput), status="old")
      !
      do 
        !! @todo Make this loop more clear @endtodo
        !
        read(1,'(a)', END = 100) line
        !
        write(2,'(a)') trim(line)
        !
        if ( INDEX(line, 'ATOMIC_POSITIONS') /= 0 ) then
          !
          do iAtom = 1, nAtoms
            !
            read(1,'(a)') line
            !
            write(2,*) elements(iAtom), newAtomicPosition(:,iAtom)
            !
          enddo
          !
        endif
        !
      enddo
100   continue
      !
      close(1)
      close(2)
      !
    enddo
    !
  end subroutine exportQEInput



!  subroutine readMjs()
!    !
!    implicit none
!    !
!    integer :: i, iE0, iE, dummyI, nEMjs
!    real(kind = dp) :: dummyD, E, MjsOfE, MjOfE0
!    character :: dummyC
!    !
!    open(1, file=trim(MjsInput), status="old")
!    !
!    read(1, *) dummyC, nEMjs
!    !
!    allocate ( Mjs(-nEnergies:nEnergies) )
!    !
!    Mjs = 1.0_dp
!    !
!!    read(1, '(d22.14,i5,4d22.14)' ) E, dummyI, dummyD, MjsOfE0, dummyD
!!    !
!!    E = E*eVToHartree 
!!    iE = int(E/deltaE) + 1
!!    !
!!    do i = 1, nEMjs - 1
!!      !
!!      iE0 = iE
!!      read(1, '(d22.14,i5,4d22.14)' ) E, dummyI, dummyD, MjsOfE, dummyD
!!      E = E*eVToHartree
!!      iE = int(E/deltaE) + 1
!!      Mjs(iE0:iE) = MjsOfE0
!!      MjsOfE0 = MjsOfE
!!      !
!!    enddo
!!    !
!!    close(1)
!!    !
!!    !do iE = -nEnergies, nEnergies
!!    !  write(44,*) real(iE, dp)*deltaE*HartreeToEv, Mjs(iE)
!!    !enddo
!!    !
!    return
!    !
!  end subroutine readMjs
  !
  !
!=====================================================================================================
! Utility functions that simplify the code and may be used multiple times
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  function wasRead(inputVal, variableName, usage, abortExecution) 
    !! Determine if an input variable still has the default value.
    !! If it does, output an error message and possibly set the program
    !! to abort. Not all variables would cause the program to abort,
    !! so this program assumes that if you pass in the logical `abortExecution`
    !! then the variable is required and causes the program to abort 
    !! if missing.
    !!
    !! I could not find a clean way to allow this function to receive
    !! different types of variables (integer, real, character, etc.), so
    !! I made the argument be an integer so that each type could be sent
    !! in a different way. Each case is set up so that the value is tested to
    !! see if it is less than zero to determine if the variable still has
    !! its default value
    !!
    !! * For strings, the default value is `''`, so pass in 
    !! `LEN(trim(variable))-1` as this should be less than zero if
    !! the string still has the default value and greater than or equal 
    !! to zero otherwise
    !! * For integers the default values are less than zero, so just pass as is 
    !! * Real variables also have a negative default value, so just pass the
    !! value cast from real to integer
    !!
    implicit none
    !
    integer, intent(in) :: inputVal
      !! Value to compare with 0 to see if a variable has been read;
    !
    character(len=*), intent(in) :: variableName
      !! Name of the variable used in output message
    character(len=*), intent(in) :: usage
      !! Example of how the variable can be used
    !
    logical, optional, intent(inout) :: abortExecution
      !! Optional logical for if the program should be aborted 
    logical :: wasRead
      !! Whether or not the input variable was read from the input file;
      !! this is the return value
    !
    !! <h2>Walkthrough</h2>
    !!
    wasRead = .true.
      !! * Default return value is true
    !
    if ( inputVal < 0) then
      !! * If the input variable still has the default value
      !!    * output an error message
      !!    * set the program to abort if that variable was sent in
      !!    * set the return value to false to indicate that the 
      !!      variable wasn't read
      !
      write(iostd, *)
      write(iostd, '(" Variable : """, a, """ is not defined!")') variableName
      write(iostd, '(" usage : ", a)') usage
      if(present(abortExecution)) then
        !
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
        !
      endif 
      !
      wasRead = .false.
      !
    endif
    !
    return
    !
  end function wasRead
  !
end module MjModule
