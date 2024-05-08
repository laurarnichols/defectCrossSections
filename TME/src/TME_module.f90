! NOTE: This file has different programming styles because
! I haven't been able to fully clean this file up yet.
! Please excuse the mess.
module TMEmod
  use constants, only: dp, pi, eVToHartree, ii
  use miscUtilities, only: int2str, int2strLeadZero
  use energyTabulatorMod, only: energyTableDir, readCaptureEnergyTable
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nKPoints, nSpins, order, ispSelect, loopSpins

  use errorsAndMPI
  use mpi
  
  implicit none


  integer, allocatable :: mill_local(:,:)
    !! Local Miller indices
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors
  integer :: nGkVecsLocal
    !! Local number of G+k vectors on this processor
  integer :: nSys
    !! Number of systems
  integer :: phononModeJ
    !! Index of phonon mode for the calculation
    !! of \(M_j\) (only for order=1)

  real(kind=dp) :: dq_j
    !! \(\delta q_j) for displaced wave functions
    !! (only order = 1)
  real(kind=dp) :: omega
    !! Supercell volume
  real(kind=dp) :: recipLattVec(3,3)
    !! Reciprocal-space lattice vectors
  real(kind = dp) :: t1, t2
    !! For timing different processes

  complex(kind=dp), allocatable :: Ufi(:,:)
    !! All-electron overlap
  complex(kind=dp), allocatable :: Ylm(:,:)
    !! Spherical harmonics
  
  character(len=300) :: baselineDir
    !! Directory for baseline overlap to optionally
    !! be subtracted for the first-order term
  character(len=300) :: braExportDir, ketExportDir
    !! Path to export dirs
  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: outputDir
    !! Path to where matrix elements should be output

  logical :: subtractBaseline
    !! If baseline should be subtracted from first-order
    !! overlap for increased numerical accuracy in the 
    !! derivative

  type atomInfo
    integer, allocatable :: angMom(:)
      !! Angular momentum of projectors
    integer :: iRAugMax
      !! Max index of augmentation sphere
    integer :: lmMax
      !! Total number of nlm channels
    integer :: nChannels
      !! Number of l channels;
      !! also number of projectors
    integer :: nMax
      !! Number of radial grid points

    real(kind=dp), allocatable :: bes_J_qr(:,:)
      !! Needed for PAW
    real(kind=dp), allocatable :: dRadGrid(:)
      !! Derivative of radial grid
    real(kind=dp), allocatable :: F(:,:)
      !! Needed for PAW
    real(kind=dp), allocatable :: F1bra(:,:,:)
      !! Needed for PAW (bra version)
    real(kind=dp), allocatable :: F1ket(:,:,:)
      !! Needed for PAW (ket version)
    real(kind=dp), allocatable :: F2(:,:,:)
      !! Needed for PAW
    real(kind=dp), allocatable :: FI(:,:)
      !! Radial integration with Bessel function (?)
    real(kind=dp), allocatable :: radGrid(:)
      !! Radial grid points
    real(kind=dp), allocatable :: wae(:,:)
      !! AE wavefunction
    real(kind=dp), allocatable :: wps(:,:)
      !! PS wavefunction

    character(len=2) :: element
      !! Name of the element
  end type atomInfo

  ! Define a type to match the export code for the variables
  ! that come straight from the POTCAR file
  type potcar
    integer :: maxAngMom
      !! Maximum angular momentum of the projectors
    integer :: maxNAtomTypes
      !! Maximum number of atom types across all systems

    type(atomInfo), allocatable :: atom(:)
      !! Track the pseudopotential information specific to 
      !! each of the atoms
  end type potcar

  type(potcar) :: pot
    !! Global pseudopotential variable

  ! Define a type for each of the crystal inputs
  type :: crystal
    integer, allocatable :: iType(:)
      !! Atom type index
    integer :: nAtoms
      !! Number of atoms
    integer, allocatable :: nAtomsEachType(:)
      !! Number of atoms of each type
    integer :: nAtomTypes
      !! Number of types of atoms
    integer :: nKPoints
      !! Number of k-points
    integer :: nGVecsGlobal
      !! Global number of G-vectors
    integer :: nProj
      !! Number of projectors across all atom types
    integer, allocatable :: nPWs1kGlobal(:)
      !! Global number of PWs at each k-point
    integer :: nSpins
      !! Number of spins

    real(kind=dp), allocatable :: atomPositionsCart(:,:)
      !! Position of atoms in cartesian coordinates
    real(kind=dp) :: omega
      !! Supercell volume
    real(kind=dp) :: recipLattVec(3,3)
      !! Reciprocal-space lattice vectors

    complex(kind=dp), allocatable :: beta(:,:)
      !! Projectors
    complex(kind = dp), allocatable :: crossProjection(:)
      !! Cross projections with projectors from this system
      !! and the wave function from the other
    complex(kind=dp), allocatable :: exp_iGDotR(:,:)
      !! e^(iG*R)
    complex(kind = dp), allocatable :: pawK(:)
      !! PAW k correction
    complex(kind = dp) :: pawWfc
      !! PAW wave function correction
    complex(kind=dp), allocatable :: projection(:)
      !! Projections
    complex*8, allocatable :: wfc(:)
      !! Wave function coefficients

    character(len=2), allocatable :: element(:)
      !! Element names for atoms in the system
    character(len=300) :: exportDir
      !! Path to Export directory
    character(len=300) :: ID
      !! How this system should be identified in output
    character(len=3) :: sysType
      !! Type of system for overlap:
      !! bra or ket
  end type crystal

  type(crystal) :: braSys, ketSys
    !! Define variables for the system used as the
    !! <bra| and the |ket> in the overlap matrix 
    !! element
  type(crystal), allocatable :: crystalSystem(:)
    !! Array containing all crystal systems
  
  real(kind = dp) t0, tf
  
  
contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, order, phononModeJ, baselineDir, &
          braExportDir, ketExportDir, dqFName, energyTableDir, outputDir, loopSpins, subtractBaseline)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(out) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(out) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term
    character(len=300), intent(out) :: braExportDir, ketExportDir
      !! Path to export dirs
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(out) :: outputDir
      !! Path to where matrix elements should be output

    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel
    logical, intent(out) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative


    namelist /TME_Input/ ketExportDir, braExportDir, outputDir, energyTableDir, &
                         iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, &
                         order, dqFName, phononModeJ, subtractBaseline, baselineDir, &
                         ispSelect
    

    if(ionode) then
    
      call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, order, phononModeJ, baselineDir, &
              braExportDir, ketExportDir, dqFName, energyTableDir, outputDir, subtractBaseline)
    
      read(5, TME_Input, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputParams', 'reading TME_Input namelist', abs(ierr))
    
      call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, order, phononModeJ, &
              baselineDir, braExportDir, ketExportDir, dqFName, energyTableDir, outputDir, subtractBaseline, loopSpins)

    endif

    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(ispSelect, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(loopSpins, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(phononModeJ, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(subtractBaseline, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(braExportDir, len(braExportDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(ketExportDir, len(ketExportDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(baselineDir, len(baselineDir), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)
    
    return
    
  end subroutine readInputParams
  
!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, order, phononModeJ, baselineDir, &
          braExportDir, ketExportDir, dqFName, energyTableDir, outputDir, subtractBaseline)
    
    implicit none

    ! Output variables:
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(out) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(out) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term
    character(len=300), intent(out) :: braExportDir, ketExportDir
      !! Path to export dirs
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(out) :: outputDir
      !! Path to where matrix elements should be output

    logical, intent(out) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative
    

    braExportDir = ''
    ketExportDir = ''
    energyTableDir = ''
    dqFName = ''
    outputDir = './TMEs'
    baselineDir = ''

    order = -1

    ispSelect = -1

    phononModeJ = -1
    
    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    
    subtractBaseline = .false.
    
    return
    
  end subroutine initialize
  
!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, order, phononModeJ, &
          baselineDir, braExportDir, ketExportDir, dqFName, energyTableDir, outputDir, subtractBaseline, loopSpins)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(in) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(in) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term
    character(len=300), intent(in) :: braExportDir, ketExportDir
      !! Path to export dirs
    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(in) :: outputDir
      !! Path to where matrix elements should be output

    logical, intent(in) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative

    ! Output variables:
    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel
    
    ! Local variables
    logical :: abortExecution
      !! If program should stop
    
    
    write(*,'("Inputs: ")')
    
    abortExecution = checkDirInitialization('braExportDir', braExportDir, 'input')
    abortExecution = checkDirInitialization('ketExportDir', ketExportDir, 'input') .or. abortExecution
    abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.1.1') .or. abortExecution
    abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution

    if(ispSelect < 1 .or. ispSelect > 2) then
      write(*,*) "No valid choice for spin channel selection given. Looping over spin."
      loopSpins = .true.
    else
      write(*,'("Only exporting spin channel ", i2)') ispSelect
      loopSpins = .false.
    endif

    if(order == 1) then
      abortExecution = checkFileInitialization('dqFName', dqFName) .or. abortExecution
      abortExecution = checkIntInitialization('phononModeJ', phononModeJ, 1, int(1e9)) .or. abortExecution

      if(subtractBaseline) then
        if(ispSelect == 2) then
          abortExecution = checkDirInitialization('baselineDir', baselineDir, 'allElecOverlap.2.1') .or. abortExecution
        else
          abortExecution = checkDirInitialization('baselineDir', baselineDir, 'allElecOverlap.1.1') .or. abortExecution
        endif
      endif
    endif

    abortExecution = checkStringInitialization('outputDir', outputDir) .or. abortExecution
    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 


    call system('mkdir -p '//trim(outputDir))

    
    if(abortExecution) then
      write(*,'(" Program stops!")')
      stop
    endif
    
    return
    
  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine setUpSystemArray(nSys, braExportDir, ketExportDir, crystalSystem)

    implicit none

    ! Input variables:
    integer, intent(in) :: nSys
      !! Number of systems

    character(len=300), intent(in) :: braExportDir, ketExportDir
      !! Path to export dirs

    ! Output variables
    type(crystal), allocatable :: crystalSystem(:)
      !! Array containing all crystal systems


    allocate(crystalSystem(nSys))

    crystalSystem(1)%sysType = 'bra'
    crystalSystem(1)%ID = 'bra'
    crystalSystem(1)%exportDir = braExportDir

    crystalSystem(2)%sysType = 'ket'
    crystalSystem(2)%ID = 'ket'
    crystalSystem(2)%exportDir = ketExportDir

    return

   end subroutine setUpSystemArray

!----------------------------------------------------------------------------
  subroutine completePreliminarySetup(nSys, order, phononModeJ, dqFName, mill_local, nGVecsGlobal, nKPoints, nSpins, &
        dq_j, omega, recipLattVec, Ylm, crystalSystem, pot)

    implicit none

    ! Input variables:
    integer, intent(in) :: nSys
      !! Number of systems to read files for
    integer, intent(in) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(in) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms

    ! Output variables:
    integer, allocatable, intent(out) :: mill_local(:,:)
      !! Local Miller indices
    integer, intent(out) :: nGVecsGlobal
      !! Number of global G-vectors (tested to be consistent
      !! across all systems)
    integer, intent(out) :: nKPoints
      !! Number of k-points (tested to be consistent
      !! across all systems)
    integer, intent(out) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp), intent(out) :: dq_j
      !! \(\delta q_j) for displaced wave functions
      !! (only order = 1)
    real(kind=dp), intent(out) :: omega
      !! Cell volume (tested to be consistent
      !! across all systems)
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal-space lattice vectors

    complex(kind=dp), allocatable, intent(out) :: Ylm(:,:)
      !! Spherical harmonics

    type(crystal) :: crystalSystem(nSys)
      !! Array containing all crystal systems

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables:
    integer :: isys
      !! Loop variable


    if(ionode) write(*, '("Pre-k-loop: [ ] Read inputs  [ ] Read full PW grid  [ ] Set up tables ")')
    call cpu_time(t1)

    ! Initialize global variables to be tracked to make sure they
    ! are consistent across all of the systems
    nKPoints = -1
    nGVecsGlobal = -1
    nSpins = 0
    omega = -1.0

    do isys = 1, nSys

      call readInputFileSkipPseudo(nGVecsGlobal, nKPoints, nSpins, omega, crystalSystem(isys))

    enddo

    recipLattVec = crystalSystem(1)%recipLattVec
      ! Could, in principle, check that all of the systems are consistent,
      ! but I don't want to do that right now. Right now, the code just 
      ! assumes that they are consistent.


    call getGlobalPseudo(nSys, crystalSystem, pot)
    
    pot%maxAngMom = 2*pot%maxAngMom + 1


    if(order == 1) call readDqFile(phononModeJ, dqFName, dq_j)


    call cpu_time(t2)
    if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [ ] Read full PW grid  [ ] Set up tables (",f10.2," secs)")') t2-t1
    call cpu_time(t1)


    call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
      !! * Distribute k-points in pools

    call distributeItemsInSubgroups(indexInPool, nGVecsGlobal, nProcPerPool, nProcPerPool, nProcPerPool, iGStart_pool, &
            iGEnd_pool, nGVecsLocal)
      !! * Distribute G-vectors across processes in pool


    allocate(mill_local(3,nGVecsLocal))

    call getFullPWGrid(nGVecsLocal, nGVecsGlobal, crystalSystem(1), mill_local)
      ! We have already verified at this point that the number of 
      ! G-vectors is the same in all of the systems, so it doesn't
      ! matter which system we read the PW grid from.

  
    call cpu_time(t2)
    if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [X] Read full PW grid  [ ] Set up tables (",f10.2," secs)")') t2-t1
    call cpu_time(t1)


    allocate(Ylm((pot%maxAngMom+1)**2,nGVecsLocal))

    call setUpTables(nGVecsLocal, mill_local, nSys, recipLattVec, pot, Ylm, crystalSystem)
    ! Also allocate this table for the braSys so that there isn't an error when deallocating

    deallocate(mill_local)

  
    call cpu_time(t2)
    if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [X] Read full PW grid  [X] Set up tables (",f10.2," secs)")') t2-t1
    call cpu_time(t1)

    return

  end subroutine completePreliminarySetup
  
!----------------------------------------------------------------------------
  subroutine readInputFileSkipPseudo(nGVecsGlobal, nKPoints, nSpins, omega, sys)

    use miscUtilities, only: getFirstLineWithKeyword, ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
    integer, intent(inout) :: nGVecsGlobal
      !! Number of global G-vectors (tested to be consistent
      !! across all systems)
    integer, intent(inout) :: nKPoints
      !! Number of k-points (tested to be consistent
      !! across all systems)
    integer, intent(inout) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp) :: omega
      !! Cell volume (tested to be consistent
      !! across all systems)

    type(crystal) :: sys
      !! The crystal system

    ! Local variables:
    integer :: iDum
      !! Dummy integer
    integer :: ik, iA, ix, iT
      !! Loop indices
    integer :: nBands
      !! Number of bands
    integer :: lmMax
      !! Needed to calculate nProj

    real(kind=dp) :: t1, t2 
      !! Timers
    real(kind=dp) :: rDum
      !! Dummy real variable

    character(len=300) :: inputFName
      !! File name for the input file 
    character(len=300) :: textDum
      !! Dummy text
    
    
    if(ionode) then
      write(*,'(" Reading ",a," input file")') trim(sys%ID)
      call cpu_time(t1)
    
      inputFName = trim(trim(sys%exportDir)//'/input')
    
      open(50, file=trim(inputFName), status = 'old')
    
      read(50,*)
      read(50,*) sys%omega


      if(omega < 0) then
        omega = sys%omega
      else if(abs(omega - sys%omega) > 1e-8) then
        call exitError('readInput', 'volumes don''t match', 1)
      endif

    endif

    call MPI_BCAST(sys%omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)


    if(ionode) then
      read(50,*)
      read(50,*) sys%nGVecsGlobal


      if(nGVecsGlobal < 0) then
        nGVecsGlobal = sys%nGVecsGlobal
      else if(sys%nGVecsGlobal /= nGVecsGlobal) then
        call exitError('readInput', 'number of G vecs in system '//trim(sys%ID)//' does not match', 1)
      end if

    endif

    call MPI_BCAST(sys%nGVecsGlobal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nGVecsGlobal, 1, MPI_INTEGER, root, worldComm, ierr)


    if(ionode) then

      read(50,*)
      read(50,*) ! fftGridSize(1:3)
    
      read(50,*)
      read(50,'(i10)') sys%nSpins

      if(sys%nSpins > nSpins) nSpins = sys%nSpins

    end if

    call MPI_BCAST(sys%nSpins, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nSpins, 1, MPI_INTEGER, root, worldComm, ierr)


    if(ionode) then
    
      read(50,*) 
      read(50,'(i10)') sys%nKPoints


      if(nKPoints < 0) then
        nKPoints = sys%nKPoints
      else if(sys%nKPoints /= nKPoints) then
        call exitError('readInput', 'number of k-points in system '//trim(sys%ID)//' does not match', 1)
      end if
    endif

    call MPI_BCAST(sys%nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    

    allocate(sys%nPWs1kGlobal(sys%nKPoints))
    
    if(ionode) then

      read(50,*) 
    
      do ik = 1, sys%nKPoints
      
        read(50,'(2i10,4ES24.15E3)') iDum, sys%nPWs1kGlobal(ik), rDum, rDum, rDum, rDum
      
      enddo

    endif
    
    call MPI_BCAST(sys%nPWs1kGlobal, sys%nKPoints, MPI_INTEGER, root, worldComm, ierr)

    if(ionode) then

      call ignoreNextNLinesFromFile(50, 5)

      read(50,'(a5, 3ES24.15E3)') textDum, sys%recipLattVec(1:3,1)
      read(50,'(a5, 3ES24.15E3)') textDum, sys%recipLattVec(1:3,2)
      read(50,'(a5, 3ES24.15E3)') textDum, sys%recipLattVec(1:3,3)
    
      read(50,*)
      read(50,'(i10)') sys%nAtoms
    
      read(50,*)
      read(50,'(i10)') sys%nAtomTypes
    
    endif

    call MPI_BCAST(sys%recipLattVec, size(sys%recipLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(sys%nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(sys%nAtomTypes, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(sys%atomPositionsCart(3,sys%nAtoms), sys%iType(sys%nAtoms))


    if(ionode) then
    
      read(50,*) 
      do iA = 1, sys%nAtoms
        read(50,'(i10, 3ES24.15E3)') sys%iType(iA), (sys%atomPositionsCart(ix,iA), ix=1,3)
      enddo
    
      read(50,*)
      read(50,'(i10)') nBands

      if(iBandIfinal > nBands .or. iBandFfinal > nBands) &
        call exitError('readInputFile', 'band limits outside the number of bands in the system '//trim(int2str(nBands)), 1)
        ! Only need to test these bands because we tested in
        ! the `checkInitialization` subroutine to make sure
        ! that the `initial` bands are lower than the `final`
        ! bands

    endif

    call MPI_BCAST(sys%iType,  size(sys%iType),  MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(sys%atomPositionsCart, size(sys%atomPositionsCart), MPI_DOUBLE_PRECISION,root,worldComm,ierr)
    
    allocate(sys%nAtomsEachType(sys%nAtomTypes), sys%element(sys%nAtomTypes))

    if(ionode) then

      sys%nProj = 0
      do iT = 1, sys%nAtomTypes

        textDum = getFirstLineWithKeyword(50,'Element')
        read(50,*) sys%element(iT)
      
        read(50,*)
        read(50,'(i10)') sys%nAtomsEachType(iT)


        textDum = getFirstLineWithKeyword(50,'Number of channels')
        read(50,'(i10)') lmMax
      
        sys%nProj = sys%nProj + sys%nAtomsEachType(iT)*lmMax

      enddo

      close(50)

    endif

    call MPI_BCAST(sys%nAtomsEachType, sys%nAtomTypes, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(sys%element, size(sys%element), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(sys%nProj, 1, MPI_INTEGER, root, worldComm, ierr)

   if(ionode) then
    
      call cpu_time(t2)
      write(*,'(" Reading ",a," input file done in:                ", f10.2, " secs.")') trim(sys%ID), t2-t1
      write(*,*)

    endif
    
    return
    
  end subroutine readInputFileSkipPseudo

!----------------------------------------------------------------------------
  subroutine getGlobalPseudo(nSys, crystalSystem, pot)

    implicit none

    ! Input variables:
    integer, intent(in) :: nSys
      !! Number of systems to read files for

    type(crystal) :: crystalSystem(nSys)
      !! Array containing all crystal systems

    ! Output variables:
    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information
      
    ! Local variables:
    integer :: isys, iT
      !! Loop indices
    integer :: isys_maxNAtomTypes
      !! Index of system with the maximum number of atom
      !! types


    if(ionode) then

      ! Get the maximum number of atom types across all systems
      ! This will be the one to get the pseudopotential info from.
      pot%maxNAtomTypes = 0
      do isys = 1, nSys
        if(crystalSystem(isys)%nAtomTypes > pot%maxNAtomTypes) then
          pot%maxNAtomTypes = crystalSystem(isys)%nAtomTypes
          isys_maxNAtomTypes = isys
        endif
      enddo


      ! Make sure that the atom types are consistent across all systems
      do iT = 1, pot%maxNAtomTypes
        do isys = 1, nSys
          if(.not. (isys == isys_maxNAtomTypes .or. iT > crystalSystem(isys)%nAtomTypes)) then
            if(trim(crystalSystem(isys)%element(iT)) /= trim(crystalSystem(isys_maxNAtomTypes)%element(iT))) then
              write(*,'("Element mismatch detected in atom ", i5," in systems ",a," and ",a,": ",a," ",a)') &
                  iT, trim(crystalSystem(isys)%ID), trim(crystalSystem(isys_maxNAtomTypes)%ID), &
                  trim(crystalSystem(isys)%element(iT)), trim(crystalSystem(isys_maxNAtomTypes)%element(iT))
              call exitError('completePreliminarySetup', 'Elements must be consistent across all systems!', 1)
            endif
          endif
        enddo
      enddo

    endif

    call MPI_BCAST(pot%maxNAtomTypes, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(isys_maxNAtomTypes, 1, MPI_INTEGER, root, worldComm, ierr)

    call readPseudoFromInputFile(crystalSystem(isys_maxNAtomTypes), pot)


    return

  end subroutine getGlobalPseudo

!----------------------------------------------------------------------------
  subroutine readPseudoFromInputFile(sys, pot)

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    type(crystal) :: sys
      !! The crystal system

    ! Output variables:
    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables:
    integer :: iDum
      !! Dummy integer
    integer :: irm
      !! Store iRAugMax for single loop
    integer :: iT, ip, ir, ip1, ip2
      !! Loop indices

    real(kind=dp), allocatable :: aepsDiff1(:), aepsDiff2(:)
      !! Difference between wae and wps for different channels

    character(len=300) :: inputFName
      !! File name for the input file 
    character(len=300) :: textDum
      !! Dummy text


    if(ionode) then
      write(*,'(" Reading pseudo information from ",a," input file")') trim(sys%ID)
      call cpu_time(t1)
    
      inputFName = trim(trim(sys%exportDir)//'/input')
    
      open(50, file=trim(inputFName), status = 'old')

      textDum = getFirstLineWithKeyword(50,'Number of Bands')
      read(50,*) 
    
    endif


    allocate(pot%atom(pot%maxNAtomTypes))
    

    pot%maxAngMom = 0
    do iT = 1, sys%nAtomTypes
      ! Could also use pot%maxNAtomTypes as they should be the same
      
      if(ionode) then

        read(50,*) 
        read(50,*) pot%atom(iT)%element
      
        read(50,*)
        read(50,*)

        read(50,*)
        read(50,'(i10)') pot%atom(iT)%nChannels              ! number of projectors

      endif

      call MPI_BCAST(pot%atom(iT)%element, len(pot%atom(iT)%element), MPI_CHARACTER, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%nChannels, 1, MPI_INTEGER, root, worldComm, ierr)
      

      allocate(pot%atom(iT)%angMom(pot%atom(iT)%nChannels))

      if(ionode) then

        read(50,*)

        do ip = 1, pot%atom(iT)%nChannels

          read(50,'(2i10)') pot%atom(iT)%angMom(ip), iDum
          if(pot%atom(iT)%angMom(ip) > pot%maxAngMom) pot%maxAngMom = pot%atom(iT)%angMom(ip)

        enddo

      endif

      call MPI_BCAST(pot%atom(iT)%angMom, pot%atom(iT)%nChannels, MPI_INTEGER, root, worldComm, ierr)


      if(ionode) then
      
        read(50,*)
        read(50,'(i10)') pot%atom(iT)%lmMax
      
        read(50,*)
        read(50,'(2i10)') pot%atom(iT)%nMax, pot%atom(iT)%iRAugMax

      endif
    
      call MPI_BCAST(pot%atom(iT)%lmMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%nMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%iRAugMax, 1, MPI_INTEGER, root, worldComm, ierr)


      allocate(pot%atom(iT)%radGrid(pot%atom(iT)%nMax))
      
      if(ionode) then
      
        allocate(pot%atom(iT)%dRadGrid(pot%atom(iT)%nMax))

        read(50,*)

        do ir = 1, pot%atom(iT)%nMax
          read(50,'(2ES24.15E3)') pot%atom(iT)%radGrid(ir), pot%atom(iT)%dRadGrid(ir)
        enddo
       

        allocate(pot%atom(iT)%wae(pot%atom(iT)%nMax, pot%atom(iT)%nChannels))
        allocate(pot%atom(iT)%wps(pot%atom(iT)%nMax, pot%atom(iT)%nChannels))
      
        read(50,*)
        do ip = 1, pot%atom(iT)%nChannels
          do ir = 1, pot%atom(iT)%nMax
            read(50,'(2ES24.15E3)') pot%atom(iT)%wae(ir,ip), pot%atom(iT)%wps(ir,ip) 
          enddo
        enddo
        
      endif


      allocate(pot%atom(iT)%F(pot%atom(iT)%iRAugMax, pot%atom(iT)%nChannels))
      allocate(pot%atom(iT)%F1bra(pot%atom(iT)%iRAugMax, pot%atom(iT)%nChannels, pot%atom(iT)%nChannels))
      allocate(pot%atom(iT)%F1ket(pot%atom(iT)%iRAugMax, pot%atom(iT)%nChannels, pot%atom(iT)%nChannels))
      allocate(pot%atom(iT)%F2(pot%atom(iT)%iRAugMax, pot%atom(iT)%nChannels, pot%atom(iT)%nChannels))

      
      if(ionode) then
        
        irm = pot%atom(iT)%iRAugMax
        allocate(aepsDiff1(irm), aepsDiff2(irm))

        pot%atom(iT)%F = 0.0_dp
        pot%atom(iT)%F1bra = 0.0_dp
        pot%atom(iT)%F1ket = 0.0_dp
        pot%atom(iT)%F2 = 0.0_dp
      
        do ip1 = 1, pot%atom(iT)%nChannels

          aepsDiff1 = pot%atom(iT)%wae(1:irm,ip1) - pot%atom(iT)%wps(1:irm,ip1)

          pot%atom(iT)%F(1:irm,ip1)= aepsDiff1(:)*pot%atom(iT)%radGrid(1:irm)*pot%atom(iT)%dRadGrid(1:irm)
        
          do ip2 = 1, pot%atom(iT)%nChannels

            aepsDiff2 = pot%atom(iT)%wae(1:irm,ip2) - pot%atom(iT)%wps(1:irm,ip2)

            pot%atom(iT)%F1bra(1:irm,ip2,ip1) = pot%atom(iT)%wps(1:irm,ip2)*aepsDiff1(:)*pot%atom(iT)%dRadGrid(1:irm)

            pot%atom(iT)%F1ket(1:irm,ip2,ip1) = pot%atom(iT)%wps(1:irm,ip1)*aepsDiff2(:)*pot%atom(iT)%dRadGrid(1:irm)
          
            pot%atom(iT)%F2(1:irm,ip2,ip1) = aepsDiff1(:)*aepsDiff2(:)*pot%atom(iT)%dRadGrid(1:irm)

          enddo
        enddo

        deallocate(pot%atom(iT)%wae)
        deallocate(pot%atom(iT)%wps)
        deallocate(pot%atom(iT)%dRadGrid)
        deallocate(aepsDiff1, aepsDiff2)

      endif

      call MPI_BCAST(pot%atom(iT)%radGrid, size(pot%atom(iT)%radGrid), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%F, size(pot%atom(iT)%F), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%F1bra, size(pot%atom(iT)%F1bra), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%F1ket, size(pot%atom(iT)%F1ket), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(pot%atom(iT)%F2, size(pot%atom(iT)%F2), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      
    enddo

    call MPI_BCAST(pot%maxAngMom, 1, MPI_INTEGER, root, worldComm, ierr)

    if(ionode) then
    
      close(50)
    
      call cpu_time(t2)
      write(*,'(" Reading pseudo information from ",a," input file done in:                ", f10.2, " secs.")') trim(sys%ID), t2-t1
      write(*,*)

    endif

    return

  end subroutine readPseudoFromInputFile

!----------------------------------------------------------------------------
  subroutine readDqFile(phononModeJ, dqFName, dq_j)

    use miscUtilities, only: ignoreNextNLinesFromFile

    implicit none

    ! Input variables:
    integer, intent(in) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms

    ! Output variables:
    real(kind=dp), intent(out) :: dq_j
      !! \(\delta q_j) for displaced wave functions
      !! (only order = 1)

    ! Local variables:
    integer :: iDum
      !! Dummy integer

    
    if(ionode) then

      open(30,file=trim(dqFName))
      call ignoreNextNLinesFromFile(30, 1+phononModeJ-1)
        ! Ignore header and all modes before phononModeJ
      read(30,*) iDum, dq_j
      close(30)

    endif

    call MPI_BCAST(dq_j, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine readDqFile

!----------------------------------------------------------------------------
  subroutine getFullPWGrid(nGVecsLocal, nGVecsGlobal, sys, mill_local)
    !! Read full PW grid from mgrid file
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors

    type(crystal) :: sys
       !! Arbitrary system to read PW grid from

    ! Output variables:
    integer, intent(out) :: mill_local(3,nGVecsLocal)
      !! Integer coefficients for G-vectors
    
    ! Local variables:
    integer :: displacement(nProcPerPool)
      !! Offset from beginning of array for
      !! scattering/gathering G-vectors 
    integer, allocatable :: gVecMillerIndicesGlobal(:,:)
      !! Integer coefficients for G-vectors on all processors
    integer :: gVecsLocalX(nGVecsLocal), gVecsLocalY(nGVecsLocal), gVecsLocalZ(nGVecsLocal)
      !! Arrays to hold x, y, and z components of local
      !! G-vectors to make scattering simpler
    integer :: iDum
      !! Ignore dummy integer
    integer :: sendCount(nProcPerPool)
      !! Number of items to send/recieve to/from each process
    integer :: ig, ix
      !! Loop index
    

    if(ionode) then
      
      allocate(gVecMillerIndicesGlobal(nGVecsGlobal,3))

      open(72, file=trim(sys%exportDir)//"/mgrid")
        !! Read full G-vector grid from defect folder.
        !! This assumes that the grids are the same.
    
      read(72, * )
      read(72, * )
    
      do ig = 1, nGVecsGlobal

        read(72, '(4i10)') iDum, (gVecMillerIndicesGlobal(ig,ix),ix=1,3)

      enddo
    
      close(72)

    else
      allocate(gVecMillerIndicesGlobal(1,3))
        ! Must include 3 for first dimension so that scatter
        ! statement is not out of bounds
    endif


    if(myPoolId == 0) then

      sendCount = 0
      sendCount(indexInPool+1) = nGVecsLocal
      call mpiSumIntV(sendCount, intraPoolComm)
        !! * Put the number of G+k vectors on each process
        !!   in a single array per band group

      displacement = 0
      displacement(indexInPool+1) = iGStart_pool-1
      call mpiSumIntV(displacement, intraPoolComm)
        !! * Put the displacement from the beginning of the array
        !!   for each process in a single array per band group

      call MPI_SCATTERV(gVecMillerIndicesGlobal(:,1), sendCount, displacement, MPI_INTEGER, gVecsLocalX(1:nGVecsLocal), nGVecsLocal, &
          MPI_INTEGER, 0, intraPoolComm, ierr)
      call MPI_SCATTERV(gVecMillerIndicesGlobal(:,2), sendCount, displacement, MPI_INTEGER, gVecsLocalY(1:nGVecsLocal), nGVecsLocal, &
          MPI_INTEGER, 0, intraPoolComm, ierr)
      call MPI_SCATTERV(gVecMillerIndicesGlobal(:,3), sendCount, displacement, MPI_INTEGER, gVecsLocalZ(1:nGVecsLocal), nGVecsLocal, &
          MPI_INTEGER, 0, intraPoolComm, ierr)


      mill_local(1,:) = gVecsLocalX(:)
      mill_local(2,:) = gVecsLocalY(:)
      mill_local(3,:) = gVecsLocalZ(:)

    endif

    deallocate(gVecMillerIndicesGlobal)


    call MPI_BCAST(mill_local, size(mill_local), MPI_INTEGER, root, interPoolComm, ierr)
    
    return
    
  end subroutine getFullPWGrid

!----------------------------------------------------------------------------
  subroutine setUpTables(nGVecsLocal, mill_local, nSys, recipLattVec, pot, Ylm, crystalSystem)

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors
    integer, intent(in) :: mill_local(3,nGVecsLocal)
      !! Miller indices for local G-vectors
    integer, intent(in) :: nSys
      !! Number of systems to read files for

    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal-space lattice vectors

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Output variables:
    complex(kind=dp), intent(out) :: Ylm((pot%maxAngMom+1)**2,nGVecsLocal)
      !! Spherical harmonics

    type(crystal) :: crystalSystem(nSys)
      !! Array containing all crystal systems

    ! Local variables:
    integer :: ig, iT, iR, isys, iL
      !! Loop indices
    integer :: L
      !! Angular momentum for this channel

    real(kind=dp) :: gCart(3)
      !! G-vectors in Cartesian coordinates
    real(kind=dp) :: gUnit(3)
      !! Unit G-vector
    real(kind=dp) :: JL(0:pot%maxAngMom)
      !! Bessel_j temporary variable
    real(kind=dp) :: q
      !! Magnitude of G-vector


    Ylm = cmplx(0.0_dp, 0.0_dp, kind = dp)

    do iT = 1, pot%maxNAtomTypes

      allocate(pot%atom(iT)%bes_J_qr(0:pot%maxAngMom, pot%atom(iT)%iRAugMax))
      pot%atom(iT)%bes_J_qr = 0.0_dp

      allocate(pot%atom(iT)%FI(pot%atom(iT)%nChannels, nGVecsLocal))
      pot%atom(iT)%FI = 0.0_dp

    enddo
    

    do isys = 1, nSys
      allocate(crystalSystem(isys)%exp_iGDotR(crystalSystem(isys)%nAtoms, nGVecsLocal))
      crystalSystem(isys)%exp_iGDotR = 0.0_dp
    enddo


    do ig = 1, nGVecsLocal

      gCart(:) = matmul(recipLattVec, mill_local(:,ig))

      do isys = 1, nSys
        call getExpiGDotR(ig, gCart, crystalSystem(isys))
      enddo

      q = sqrt(dot_product(gCart(:),gCart(:)))

      gUnit(:) = gCart(:)
      if(abs(q) > 1.0e-6_dp) gUnit = gUnit/q
        !! Get unit vector for Ylm calculation

      call getYlm(gUnit, pot%maxAngMom, Ylm(:,ig))
        !! Calculate all the needed spherical harmonics

      do iT = 1, pot%maxNAtomTypes

        do iR = 1, pot%atom(iT)%iRAugMax

          JL = 0.0_dp
          call bessel_j(q*pot%atom(iT)%radGrid(iR), pot%maxAngMom, JL) ! returns the spherical bessel at qr point

          pot%atom(iT)%bes_J_qr(:,iR) = JL(:)

        enddo

        do iL = 1, pot%atom(iT)%nChannels
          L = pot%atom(iT)%angMom(iL)

          pot%atom(iT)%FI(iL,ig) = dot_product(pot%atom(iT)%bes_J_qr(L,:),pot%atom(iT)%F(:,iL))
            ! radial part integration F contains dRadGrid
        enddo

      enddo

    enddo

    do iT = 1, pot%maxNAtomTypes
      deallocate(pot%atom(iT)%bes_J_qr)
    enddo

    return

  end subroutine setUpTables

!----------------------------------------------------------------------------
  subroutine getExpiGDotR(ig, gCart, sys)

    implicit none

    ! Input variables:
    integer, intent(in) :: ig
      !! G-vector index

    real(kind=dp), intent(in) :: gCart(3)
      !! Single G vector in Cartesian coordinates

    ! Output variables:
    type(crystal) :: sys
      !! System to get G+R vectors for

    ! Local variables:
    integer :: iA
      !! Loop index
    integer :: sign_i
      !! Sign in front of i

    real(kind=dp) :: GDotR
      !! G * R


    if(trim(sys%sysType) == 'bra') then
      sign_i = 1
    else if(trim(sys%sysType) == 'ket') then
      sign_i = -1
    endif
        

    do iA = 1, sys%nAtoms
      GDotR = dot_product(gCart, sys%atomPositionsCart(:,iA))
        
      sys%exp_iGDotR(iA,ig) = exp(sign_i*ii*cmplx(GDotR, 0.0_dp, kind = dp))
    enddo

    return

  end subroutine getExpiGDotR

!----------------------------------------------------------------------------
  subroutine calcAndWrite2SysMatrixElements(ispSelect, nSpins, braSys, ketSys, pot)

    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    type(crystal) :: braSys, ketSys
       !! The crystal systems to get the
       !! matrix element for

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables 
    integer :: ikLocal, ikGlobal, isp, ibi, ibf
      !! Loop indices

    logical :: calcSpinDepSD, calcSpinDepPC
      !! If spin-dependent subroutines should be called
    logical :: spin1Skipped = .false.
      !! If the first spin channel was skipped
    

    allocate(Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))
  

    do ikLocal = 1, nkPerPool
    
      if(ionode) write(*,'("Beginning k-point loop ", i4, " of ", i4)') ikLocal, nkPerPool
    

      ikGlobal = ikLocal+ikStart_pool-1
        !! Get the global `ik` index from the local one
    

      if(.not. thisKComplete(ikGlobal, ispSelect, nSpins)) then

        !-----------------------------------------------------------------------------------------------
        !> Read projectors

        if(ionode) write(*, '("  Pre-spin-loop: [ ] Read projectors ")') 
        call cpu_time(t1)


        if(braSys%nPWs1kGlobal(ikGlobal) /= ketSys%nPWs1kGlobal(ikGlobal)) &
          call exitError('calAndWrite2SyMatrixElements', 'number of G+k vectors does not match for ik='//trim(int2str(ikGlobal)), 1)

        call distributeItemsInSubgroups(indexInPool, braSys%nPWs1kGlobal(ikGlobal), nProcPerPool, nProcPerPool, nProcPerPool, iGkStart_pool, &
                iGkEnd_pool, nGkVecsLocal)


        allocate(braSys%beta(nGkVecsLocal,braSys%nProj))
        allocate(ketSys%beta(nGkVecsLocal,ketSys%nProj))

        call readProjectors(ikGlobal, nGkVecsLocal, braSys)
        call readProjectors(ikGlobal, nGkVecsLocal, ketSys)


        call cpu_time(t2)
        if(ionode) write(*, '("  Pre-spin-loop: [X] Read projectors (",f10.2," secs)")') t2-t1

      
        !-----------------------------------------------------------------------------------------------
        !> Allocate arrays that are potentially spin-polarized and start spin loop

        allocate(braSys%wfc(nGkVecsLocal))
        allocate(ketSys%wfc(nGkVecsLocal))
      
        allocate(braSys%crossProjection(braSys%nProj))
        allocate(ketSys%crossProjection(ketSys%nProj))

        allocate(braSys%projection(braSys%nProj))
        allocate(ketSys%projection(ketSys%nProj))

        allocate(braSys%pawK(nGVecsLocal))
        allocate(ketSys%pawK(nGVecsLocal))

      
        do isp = 1, nSpins

          if(indexInPool == 0) then
            if(isp == 1 .and. (ispSelect == 2 .or. overlapFileExists(ikGlobal,1))) spin1Skipped = .true.
          endif

          call MPI_BCAST(spin1Skipped, 1, MPI_LOGICAL, root, intraPoolComm, ierr)

          if(loopSpins .or. isp == ispSelect) then

            Ufi(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)

            if(ionode) write(*,'("  Beginning spin ", i2)') isp

            calcSpinDepPC = isp == 1 .or. ketSys%nSpins == 2 .or. spin1Skipped
            calcSpinDepSD = isp == 1 .or. braSys%nSpins == 2 .or. spin1Skipped
              ! If either of the systems only has a single spin 
              ! channel, some inputs and calculations do not need
              ! to be redone for both spin channels. However, we
              ! still need to make sure to calculate these values
              ! for the second spin channel if the first spin channel
              ! was not done for some reason (e.g., the file already
              ! existed or only the second spin channel was selected).

            !-----------------------------------------------------------------------------------------------

            if(.not. overlapFileExists(ikGlobal, isp)) then
      
              do ibi = iBandIinit, iBandIfinal
                !if(calcSpinDepPC) then
                  call readWfc(ibi, ikGlobal, min(isp,ketSys%nSpins), nGkVecsLocal, ketSys)
                  call calculateCrossProjection(ketSys, braSys)
                    ! Get new cross projection with new `ketSys%wfc`
                  call readProjections(ibi, ikGlobal, min(isp,ketSys%nSpins), ketSys)
                  call pawCorrectionK(nGVecsLocal, pot, Ylm, ketSys)
                !endif

                do ibf = iBandFinit, iBandFfinal
                  !if(calcSpinDepSD) then
                    call readWfc(ibf, ikGlobal, min(isp,braSys%nSpins), nGkVecsLocal, braSys)
                    call calculateCrossProjection(braSys, ketSys)
                      ! Get new cross projection with new `braSys%wfc`
                    call readProjections(ibf, ikGlobal, min(isp,braSys%nSpins), braSys)
                    call pawCorrectionK(nGVecsLocal, pot, Ylm, braSys)
                  !endif

                  Ufi(ibf, ibi) = dot_product(braSys%wfc(:),ketSys%wfc(:))
                  if(indexInPool == 0) call pawCorrectionWfc(ketSys, pot)
                  if(indexInPool == 1) call pawCorrectionWfc(braSys, pot)

                  call MPI_BCAST(ketSys%pawWfc, 1, MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
                  call MPI_BCAST(braSys%pawWfc, 1, MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)

                  Ufi(ibf,ibi) = Ufi(ibf,ibi) + (16.0_dp*pi*pi/omega)*dot_product(conjg(braSys%pawK(:)),ketSys%pawK(:))

                  call MPI_ALLREDUCE(MPI_IN_PLACE, Ufi(ibf,ibi), 1, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, intraPoolComm, ierr)

                  Ufi(ibf,ibi) = Ufi(ibf,ibi) + braSys%pawWfc + ketSys%pawWfc

                enddo ! Loop over final bands
              enddo ! Loop over initial bands
        

              !-----------------------------------------------------------------------------------------------
              
              ! Subtract baseline if applicable and write out results
              if(indexInPool == 0) then 
                if(order == 1 .and. subtractBaseline) &
                  call readAndSubtractBaseline(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikLocal, isp, Ufi)
        
                call writeResults(ikLocal,isp, Ufi)
              endif

  
            endif ! If this spin channel file doesn't exist
          endif ! If this spin channel shouldn't be skipped
        enddo ! Spin loop

        deallocate(braSys%wfc, ketSys%wfc)
        deallocate(braSys%beta, ketSys%beta)
        deallocate(braSys%crossProjection, ketSys%crossProjection)
        deallocate(braSys%projection, ketSys%projection)
        deallocate(braSys%pawK, ketSys%pawK)
      endif ! If both spin channels exist
    enddo ! k-point loop

    return

  end subroutine calcAndWrite2SysMatrixElements

!----------------------------------------------------------------------------
  function thisKComplete(ikGlobal, ispSelect, nSpins)

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Global k-point index
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    ! Output variables:
    logical :: thisKComplete
      !! If needed `allElecOverlap.isp.ik` files exist
      !! at the current k-point


    if(indexInPool == 0) then
      if(nSpins == 1 .or. ispSelect == 1) then
        thisKComplete = overlapFileExists(ikGlobal, 1)
      else if(ispSelect == 2) then
        thisKComplete = overlapFileExists(ikGlobal, 2)
      else if(nSpins == 2) then
        thisKComplete = overlapFileExists(ikGlobal, 1) .and. overlapFileExists(ikGlobal, 2)
      endif
    endif

    call MPI_BCAST(thisKComplete, 1, MPI_LOGICAL, root, intraPoolComm, ierr)

  end function thisKComplete    
  
!----------------------------------------------------------------------------
  function overlapFileExists(ikGlobal, isp) result(fileExists)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Output variables:
    character(len=300) :: fName
      !! File name for overlap

    logical :: fileExists
      !! If the overlap file exists for the given 
      !! k-point and spin channel


    fName = trim(getMatrixElementFNameWPath(ikGlobal, isp, outputDir))

    inquire(file=fName, exist=fileExists)

    if(fileExists) write(*,'("Overlap file ", a, " exists and will not be recalculated.")') trim(fName)
    
  end function overlapFileExists

!----------------------------------------------------------------------------
  subroutine readProjectors(ikGlobal, nGkVecsLocal, sys)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor

    ! Output variables:
    type(crystal) :: sys
       !! The crystal system
    
    ! Local variables:
    integer :: reclen
      !! Record length for projectors file
    integer :: igkLocal, igkGlobal, ipr
      !! Loop indices
    
    character(len=300) :: fNameExport
      !! File names
    

    fNameExport = trim(sys%exportDir)//"/projectors."//trim(int2str(ikGlobal)) 


    inquire(iolength=reclen) sys%beta(1,:)
      !! Get the record length needed to write a double complex
      !! array of length `nProj`

    open(unit=72, file=trim(fNameExport), access='direct', recl=reclen, iostat=ierr, status='old', SHARED)


    do igkLocal = 1, nGkVecsLocal

      igkGlobal = igkLocal+iGkStart_pool-1

      read(72,rec=igkGlobal+1) (sys%beta(igkLocal,ipr), ipr=1,sys%nProj)

    enddo

    close(72)
    
    return 

  end subroutine readProjectors

!----------------------------------------------------------------------------
  subroutine readWfc(ib, ikGlobal, isp, nGkVecsLocal, sys)
    !! Read wave function at band ib for given system
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ib
      !! Band index
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor

    ! Output variables
    type(crystal) :: sys
       !! The crystal system

    ! Local variables:
    integer :: displacement(nProcPerPool)
      !! Offset from beginning of array for
      !! scattering coefficients to each process
    integer :: reclen
      !! Record length for projectors file
    integer :: sendCount(nProcPerPool)
      !! Number of items to send to each process
      !! in the pool
    integer :: igk
      !! Loop indices

    complex*8 :: wfcAllPWs(sys%nPWs1KGlobal(ikGlobal))
      !! Wave function read from file

    character(len=300) :: fNameExport
      !! File names


    fNameExport = trim(sys%exportDir)//'/wfc.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal))
    
    sys%wfc(:) = cmplx(0.0_dp, 0.0_dp)

    sendCount = 0
    sendCount(indexInPool+1) = nGkVecsLocal
    call mpiSumIntV(sendCount, intraPoolComm)
      !! * Put the number of G+k vectors on each process
      !!   in a single array per pool

    displacement = 0
    displacement(indexInPool+1) = iGkStart_pool-1
    call mpiSumIntV(displacement, intraPoolComm)
      !! * Put the displacement from the beginning of the array
      !!   for each process in a single array per pool
    

    inquire(iolength=reclen) wfcAllPWs(:)
      !! Get the record length needed to write a complex
      !! array of length nPWs1k

    if(indexInPool == 0) then
      open(unit=72, file=trim(fNameExport), access='direct', recl=reclen, iostat=ierr, status='old', SHARED)

      read(72,rec=ib) (wfcAllPWs(igk), igk=1,sys%nPWs1kGlobal(ikGlobal))
    endif

    call MPI_SCATTERV(wfcAllPWs(:), sendCount, displacement, MPI_COMPLEX, sys%wfc(1:nGkVecsLocal), nGkVecsLocal, &
      MPI_COMPLEX, 0, intraPoolComm, ierr)

    if(indexInPool == 0) close(72)
    
    return
    
  end subroutine readWfc
  
!----------------------------------------------------------------------------
  subroutine readProjections(ib, ikGlobal, isp, sys)
    
    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Input variables:
    integer :: ib
      !! Band index
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Output variables:
    type(crystal) :: sys
       !! The crystal system

    ! Local variables:
    integer :: reclen
      !! Record length for projections files

    character(len=300) :: fNameExport
      !! Export file name
    
    
    sys%projection(:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
    
    if(indexInPool == 0) then

      inquire(iolength=reclen) sys%projection(:)
        !! Get the record length needed to write a complex
        !! array of length nProj

      ! Open the projections file for the given crystal type
      fNameExport = trim(sys%exportDir)//"/projections."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))


      open(72, file=trim(fNameExport), access='direct', form='unformatted', recl=reclen)
    
    
      ! Read the projections
      read(72,rec=ib) sys%projection(:)
    
      close(72)

    endif

    call MPI_BCAST(sys%projection, size(sys%projection), MPI_DOUBLE_COMPLEX, root, intraPoolComm, ierr)
      ! Broadcast entire array to all processes
    
    return
    
  end subroutine readProjections
  
!----------------------------------------------------------------------------
  subroutine calculateCrossProjection(sysWfc, sysBeta)
    !! Calculate the cross projection of one crystal's projectors
    !! on the other crystal's wave function coefficients, distributing
    !! the result to all processors
    
    implicit none

    ! Input variables:
    type(crystal) :: sysWfc
       !! The crystal system to get the wave functions from

    ! Output variables:
    type(crystal) :: sysBeta
       !! The crystal system to get the projectors from
    
    ! Local variables:
    integer :: ipr
      !! Loop indices

    complex(kind=dp) :: crossProjectionLocal
      !! Local version of cross projection to
      !! be summed across processors in pool
    

    sysBeta%crossProjection(:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do ipr = 1, sysBeta%nProj

      crossProjectionLocal = dot_product(sysBeta%beta(:,ipr),sysWfc%wfc(:))
        ! `dot_product` automatically conjugates first argument for 
        ! complex variables.

      call MPI_ALLREDUCE(crossProjectionLocal, sysBeta%crossProjection(ipr), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, intraPoolComm, ierr)

    enddo
    
    return
    
  end subroutine calculateCrossProjection
  
!----------------------------------------------------------------------------
  subroutine pawCorrectionWfc(sys, pot)
    ! calculates the augmentation part of the transition matrix element
    
    implicit none

    ! Input variables:
    type(crystal), intent(inout) :: sys
      !! The crystal system

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables:
    integer :: ia, iT, L, M, LP, MP
      !! Loop indices

    complex(kind = dp) :: projectionI, projectionF
      !! Local storage of projection elements for speed

    integer :: LL, LLP, LMBASE, LM, LMP
    real(kind = dp) :: atomicOverlap
    

    sys%pawWfc = 0.0_dp
    
    LMBASE = 0
    
    do ia = 1, sys%nAtoms
      
      iT = sys%iType(ia)

      LM = 0
      do LL = 1, pot%atom(iT)%nChannels

        L = pot%atom(iT)%angMom(LL)

        do M = -L, L
          LM = LM + 1 !1st index for CPROJ
          
          LMP = 0
          do LLP = 1, pot%atom(iT)%nChannels

            LP = pot%atom(iT)%angMom(LLP)

            do MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) then 
                if(trim(sys%sysType) == 'bra') then

                  atomicOverlap = sum(pot%atom(iT)%F1bra(:,LL, LLP))

                  projectionI = sys%crossProjection(LMP + LMBASE)
                  
                  projectionF = conjg(sys%projection(LM + LMBASE))
                    
                  sys%pawWfc = sys%pawWfc + projectionF*atomicOverlap*projectionI

                else if(trim(sys%sysType) == 'ket') then
                  atomicOverlap = sum(pot%atom(iT)%F1ket(:,LL, LLP))

                  projectionI = sys%projection(LMP + LMBASE)
                  
                  projectionF = conjg(sys%crossProjection(LM + LMBASE))
                    
                  sys%pawWfc = sys%pawWfc + projectionF*atomicOverlap*projectionI
                      
                endif

              endif ! If l = l' and m = m'
            enddo ! Loop over m'
          enddo ! Loop over l'
        enddo ! Loop over m quantum number
      enddo ! Loop over l quantum number

      LMBASE = LMBASE + pot%atom(iT)%lmMax

    enddo ! Loop over atoms
    
    return
    
  end subroutine pawCorrectionWfc

!----------------------------------------------------------------------------
  subroutine pawCorrectionK(nGVecsLocal, pot, Ylm, sys)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    complex(kind=dp), intent(in) :: Ylm((pot%maxAngMom+1)**2,nGVecsLocal)
      !! Spherical harmonics

    ! Output variables:
    type(crystal), intent(inout) :: sys
       !! The crystal system

    ! Local variables:
    integer :: ig, iT, iA, ind, iL
      !! Loop indices
    integer :: LMBASE, LM, L, M
      !! L and M quantum number trackers
    integer :: sign_i
      !! Sign in front of i
    
    complex(kind = dp) :: VifQ_aug
    

    if(trim(sys%sysType) == 'bra') then
      sign_i = 1
    else if(trim(sys%sysType) == 'ket') then
      sign_i = -1
    endif

    sys%pawK(:) = 0.0_dp
    
    do ig = 1, nGVecsLocal
      
      LMBASE = 0
      
      do iA = 1, sys%nAtoms
        
        iT = sys%iType(iA)
        LM = 0
        do iL = 1, pot%atom(iT)%nChannels
          L = pot%atom(iT)%angMom(iL)

          DO M = -L, L
            LM = LM + 1 ! 1st index for projection
            
            ind = L*(L + 1) + M + 1 ! index for spherical harmonics

            if(trim(sys%sysType) == 'ket') then

              VifQ_aug = sys%exp_iGDotR(iA,ig)*Ylm(ind,ig)*iToTheInt(L,sign_i)*pot%atom(iT)%FI(iL,ig)

              sys%pawK(ig) = sys%pawK(ig) + VifQ_aug*sys%projection(LM + LMBASE)

            else if(trim(sys%sysType) == 'bra') then

              VifQ_aug = sys%exp_iGDotR(iA,ig)*conjg(Ylm(ind,ig))*iToTheInt(L,sign_i)*pot%atom(iT)%FI(iL,ig)
                
              sys%pawK(ig) = sys%pawK(ig) + VifQ_aug*conjg(sys%projection(LM + LMBASE))
                
            endif
          enddo
        enddo

        LMBASE = LMBASE + pot%atom(iT)%lmMax
      enddo
      
    enddo
    
    return
    
  end subroutine pawCorrectionK

!----------------------------------------------------------------------------
  function iToTheInt(pow, sign_i)

    implicit none

    ! Input variables:
    integer, intent(in) :: pow
      !! Power for i
    integer, intent(in) :: sign_i
      !! Prefactor in front of i (+1 or -1)

    ! Output variables:
    complex(kind=dp) :: iToTheInt
      !! Result of (+-i)^pow

    
    if(mod(pow,4) == 0) then
      iToTheInt = 1.0_dp
    else if(mod(pow,4) == 1) then
      iToTheInt = sign_i*cmplx(0.0,1.0,kind=dp)
    else if(mod(pow,4) == 2) then
      iToTheInt = -1.0
    else if(mod(pow,4) == 3) then
      iToTheInt = -sign_i*cmplx(0.0,1.0,kind=dp)
    endif

  end function iToTheInt

!----------------------------------------------------------------------------
  subroutine readAndSubtractBaseline(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikLocal, isp, Ufi)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Output variables:
    complex(kind=dp), intent(inout) :: Ufi(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! All-electron overlap

    ! Local variables:
    integer :: iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      !! Energy band bounds for initial and final state from the energy table file
    integer :: ibi, ibf
      !! Loop indices
    integer :: iDum
      !! Dummy integer to ignore input
    integer :: ikGlobal
      !! Current global k-point
    
    real(kind = dp) :: rDum
      !! Dummy real to ignore input

    complex(kind = dp):: baselineOverlap
      !! Input complex overlap 

    character(len=300) :: baselineFName
      !! Name of baseline overlap file


    ikGlobal = ikLocal+ikStart_pool-1

    baselineFName = trim(getMatrixElementFNameWPath(ikGlobal, isp, baselineDir)) 
    open(17, file=trim(baselineFName), status='unknown')

    read(17,*) 
    read(17,*) 
    read(17,*) 
    read(17,*) 
    read(17,'(5i10)') iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_

    ! Check the input band bounds against those in the baseline overlap file
    if(iBandIinit < iBandIinit_ .or. iBandIfinal > iBandIfinal_ .or. iBandFinit < iBandFinit_ .or. iBandFfinal > iBandFfinal_) &
      call exitError('readAndSubtractBaseline', 'given band bounds are outside those in baseline overlap file '//trim(baselineFName), 1)
    
    read(17,*) 

    if(order == 1) read(17,*)
      ! Ignore additional line for phonon mode 

    
    do ibf = iBandFinit_, iBandFfinal_
      do ibi = iBandIinit_, iBandIfinal_
      
        read(17, 1001) iDum, iDum, baselineOverlap, rDum, rDum

        if(ibi >= iBandIinit .and. ibi <= iBandIfinal .and. ibf >= iBandFinit .and. ibf <= iBandFfinal) &
          Ufi(ibf,ibi) = Ufi(ibf,ibi) - baselineOverlap
          
      enddo
    enddo
    
    close(17)
    
 1001 format(2i7,4ES24.15E3)
    
    return
    
  end subroutine readAndSubtractBaseline
  
!----------------------------------------------------------------------------
  subroutine writeResults(ikLocal, isp, Ufi)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    complex(kind=dp), intent(in) :: Ufi(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal)
      !! All-electron overlap

    ! Local variables:
    integer :: ibi, ibf
      !! Loop indices
    integer :: ikGlobal
      !! Current global k-point
    integer :: totalNumberOfElements
      !! Total number of overlaps to output

    real(kind=dp) :: dE(iBandFinit:iBandFfinal,iBandIinit:iBandIFinal,3)
      !! Energy difference to be combined with
      !! overlap for matrix element
    
    character(len = 300) :: text
      !! Text for header


    ikGlobal = ikLocal+ikStart_pool-1
    
    open(17, file=trim(getMatrixElementFNameWPath(ikGlobal, isp, outputDir)), status='unknown')
    
    write(17, '("# Total number of k-points, k-point index, spin index Format : ''(3i10)''")')
    write(17,'(3i10)') nKPoints, ikGlobal, isp

    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') omega
    
    text = "# Total number of <f|i> elements, Initial States (bandI, bandF), Final States (bandI, bandF)"
    write(17,'(a, " Format : ''(5i10)''")') trim(text)
    
    totalNumberOfElements = (iBandIfinal - iBandIinit + 1)*(iBandFfinal - iBandFinit + 1)
    write(17,'(5i10)') totalNumberOfElements, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    

    if(order == 0) then

      text = "# Final Band, Initial Band, Complex <f|i>, |<f|i>|^2, |dE*<f|i>|^2 (Hartree^2)" 
      write(17, '(a, " Format : ''(2i10,3ES24.15E3)''")') trim(text)

    else if(order == 1) then

      write(17,'("# Phonon mode j, dq_j (Bohr*sqrt(elec. mass)). Format: ''(a78, i7, ES24.15E3)'' ", i7, ES24.15E3)') phononModeJ, dq_j
    
      if(subtractBaseline) then
        text = "# Final Band, Initial Band, Complex <f|i>-baseline, |<f|i>|^2, |dE*<f|i>/dq_j|^2 (Hartree^2/(Bohr*sqrt(elec. mass))^2)" 
      else
        text = "# Final Band, Initial Band, Complex <f|i>, |<f|i>|^2, |dE*<f|i>/dq_j|^2 (Hartree^2/(Bohr*sqrt(elec. mass))^2)" 
      endif
      write(17, '(a, " Format : ''(2i7,4ES24.15E3)''")') trim(text)

    endif


    call readCaptureEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, energyTableDir, dE)

    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal
        
        if(order == 0) then
          write(17, 1001) ibf, ibi, Ufi(ibf,ibi), abs(Ufi(ibf,ibi))**2, abs(dE(ibf,ibi,2)*Ufi(ibf,ibi))**2
        else if(order == 1) then
          write(17, 1001) ibf, ibi, Ufi(ibf,ibi), abs(Ufi(ibf,ibi))**2, abs(dE(ibf,ibi,3)*Ufi(ibf,ibi)/dq_j)**2
        endif
            
      enddo
    enddo

    close(17)
    
    call cpu_time(t2)
    write(*, '("    Ufi(:,:) of k-point ", i4, " and spin ", i1, " written.")') ikGlobal, isp
    
 1001 format(2i7,4ES24.15E3)
    
    return
    
  end subroutine writeResults
   
!----------------------------------------------------------------------------
  subroutine bessel_j (x, lmax, jl)
    
    ! x is the argument of j, jl(0:lmax) is the output values.
    implicit none
    integer, intent(in) :: lmax
    real(kind = dp), intent(in) :: x
    real(kind = dp), intent(out) :: jl(0:lmax)
    integer :: l
    
    if (x <= 0.0_dp) then
      jl = 0.0_dp
      jl(0) = 1.0_dp
      return
    end if
    
    jl(0) = sin(x)/x
    if (lmax <= 0) return
    jl(1) = (jl(0)-cos(x))/x
    if (lmax == 1) return
    
    do l = 2, lmax
      jl(l) = dble(2*l-1)*jl(l-1)/x - jl(l-2)
    enddo
    
    return
    
  end subroutine bessel_j
  
  
!----------------------------------------------------------------------------
  subroutine getYlm(v_in,lmax,y)
  !
  ! lmax   : spherical harmonics are calculated for l = 0 to lmax
  ! v      : vector, argument of the spherical harmonics (we calculate
  ! Ylm(v/norm(v))
  ! y      : array containing Ylm(v) for several l,m
  !
  ! !DESCRIPTION:
  !   1.  PURPOSE
  !        The spherical harmonics (Condon and Shortley convention)
  !          Y(0,0),Y(1,-1),Y(1,0),Y(1,1),Y(2,-2) ... Y(LMAX,LMAX)
  !        for vector V (given in Cartesian coordinates)
  !        are calculated. In the Condon Shortley convention the
  !        spherical harmonics are defined as
  !        $$ Y(l,m) = (-1)^m \sqrt{\frac{1}{\pi}} P_{lm}(\cos{\theta})
  !        \rm
  !        e^{\rm i m \phi} $$
  !                        
  !        where  $P_{lm}(\cos{\theta})$ is the normalized Associated
  !        Legendre
  !                  
  !        function. Thus,
  !                                             
  !                     $$  Y(l,-m) = (-1)^m Y^*(l,m) $$
  !
  !   2.  USAGE
  !        DOUBLE PRECISION V(3), Y(5*5)
  !        V(1) = ...
  !        V(2) = ...
  !        V(3) = ...
  !        CALL YLM(V,4,Y)
  !
  !       ARGUMENT-DESCRIPTION
  !          V      - DOUBLE PRECISION vector, dimension 3        (input)
  !                   Must be given in Cartesian coordinates.
  !                   Conversion of V to polar coordinates gives the
  !                   angles Theta and Phi necessary for the calculation
  !                   of the spherical harmonics.
  !          LMAX   - INTEGER value                               (input)
  !                   upper bound of L for which spherical harmonics
  !                   will be calculated
  !                   constraint:
  !                      LMAX >= 0
  !          Y      - COMPLEX*16 array, dimension (LMAX+1)**2    (output)
  !                   contains the calculated spherical harmonics
  !                   Y(1)                   for L .EQ. 0 (M = 0)
  !                   Y(2), ..., Y(4)        for L .EQ. 1 (M = -1, 0, 1)
  !                   ...
  !                   Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))
  !                                          for L .EQ. LMAX
  !                                              (M = -L,...,L)
  !                   constraint:
  !                      Dimension of Y .GE. (LMAX+1)**2 (not checked)
  !        USED SUBROUTINES (DIRECTLY CALLED)
  !           none
  !
  !        INDIRECTLY CALLED SUBROUTINES
  !           none
  !
  !        UTILITY-SUBROUTINES (USE BEFOREHAND OR AFTERWARDS)
  !           none
  !
  !        INPUT/OUTPUT (READ/WRITE)
  !           none
  !
  !        MACHINENDEPENDENT PROGRAMPARTS
  !           Type COMPLEX*16 is used which does not conform to the
  !           FORTRAN 77 standard.
  !           Also the non-standard type conversion function DCMPLX()
  !           is used which combines two double precision values into
  !           one double complex value.
  !
  !   3.     METHOD
  !           The basic algorithm used to calculate the spherical
  !           harmonics for vector V is as follows:
  !
  !           Y(0,0)
  !           Y(1,0)
  !           Y(1,1)
  !           Y(1,-1) = -Y(1,1)
  !           DO L = 2, LMAX
  !              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
  !              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
  !              DO M = L-2, 0, -1
  !                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
  !                 Y(L,-M)= (-1)**M*Y(L,M)
  !              ENDDO
  !           ENDDO
  !
  !           In the following the necessary recursion formulas and
  !           starting values are given:
  !
  !        Start:
  !%                        +------+
  !%                        |   1     
  !%           Y(0,0) =  -+ | -----  
  !%                       \| 4(Pi)  
  !%
  !%                                   +------+
  !%                                   |   3     
  !%           Y(1,0) =  cos(Theta) -+ | -----  
  !%                                  \| 4(Pi)  
  !%
  !%                                     +------+
  !%                                     |   3    i(Phi)
  !%           Y(1,1) =  - sin(Theta) -+ | ----- e
  !%                                    \| 8(Pi)  
  !%
  !%        Formula 1:
  !%
  !%           Y(l,l) =
  !%                           +--------+
  !%                           | (2l+1)   i(Phi)
  !%            -sin(Theta) -+ | ------  e       Y(l-1,l-1)
  !%                          \|   2l  
  !%
  !%        Formula 2:
  !%                                  +---------------+  
  !%                                  |  (2l-1)(2l+1)   
  !%           Y(l,m) = cos(Theta) -+ | -------------- Y(l-1,m)  -
  !%                                 \|   (l-m)(l+m)       
  !%
  !%                                    +--------------------+  
  !%                                    |(l-1+m)(l-1-m)(2l+1)
  !%                              -  -+ |-------------------- Y(l-2,m)
  !%                                   \|  (2l-3)(l-m)(l+m)                 
  !%
  !%        Formula 3: (not used in the algorithm because of the division
  !%                    by sin(Theta) which may be zero)
  !%
  !%                                    +--------------+  
  !%                      cos(Theta)    |  4(m+1)(m+1)   -i(Phi)
  !%           Y(l,m) = - ---------- -+ | ------------  e       Y(l,m+1) -
  !%                      sin(Theta)   \| (l+m+1)(l-m)       
  !%
  !%                                    +--------------+  
  !%                                    |(l-m-1)(l+m+2)  -2i(Phi)
  !%                              -  -+ |-------------- e        Y(l,m+2)
  !%                                   \| (l-m)(l+m+1)                         
  !%                                  
  !%
  ! !REVISION HISTORY:
  !   26. April 1994                                   Version 1.2
  !   Taken 8 1 98 from SRC_lapw2 to SRC_telnes
  !   Updated November 2004 (Kevin Jorissen)
  !   cosmetics March 2005 (Kevin Jorissen)
  !
      implicit none
  
  !   In/Output :
  
      integer, intent(in) :: LMAX
      real(kind = dp), intent(in) :: V_in(3)
      complex(kind = dp), intent(out) :: Y(*)
  !   Local variables :
      real(kind = dp), parameter :: pi = 3.1415926535897932384626433_dp
  
      INTEGER         ::  I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      real(kind = dp) ::  A, B, C, AB, ABC, ABMAX, ABCMAX, V(3)
      real(kind = dp) ::  D4LL1C, D2L13
      real(kind = dp) ::  COSTH, SINTH, COSPH, SINPH
      real(kind = dp) ::  TEMP1, TEMP2, TEMP3
      real(kind = dp) ::  YLLR, YLL1R, YL1L1R, YLMR
      real(kind = dp) ::  YLLI, YLL1I, YL1L1I, YLMI
      
      ! Y(0,0)
      
      do INDEX = 1,3
        V(INDEX) = dble(V_in(INDEX))
      enddo
      YLLR = 1.0_dp/sqrt(4.0_dp*PI)
      YLLI = 0.0_dp
      Y(1) = CMPLX(YLLR, YLLI, kind = dp)
      
      ! continue only if spherical harmonics for (L .GT. 0) are desired
      
      IF (LMAX .LE. 0) GOTO 999
      
      ! calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
      ! Theta, Phi ... polar angles of vector V
      
      ABMAX  = MAX(ABS(V(1)),ABS(V(2)))
      IF (ABMAX .GT. 0.0_dp) THEN
        A = V(1)/ABMAX
        B = V(2)/ABMAX
        AB = SQRT(A*A+B*B)
        COSPH = A/AB
        SINPH = B/AB
      ELSE
        COSPH = 1.0_dp
        SINPH = 0.0_dp
      ENDIF
      ABCMAX = MAX(ABMAX,ABS(V(3)))
      IF (ABCMAX .GT. dble(0)) THEN
        A = V(1)/ABCMAX
        B = V(2)/ABCMAX
        C = V(3)/ABCMAX
        AB = A*A + B*B
        ABC = SQRT(AB + C*C)
        COSTH = C/ABC
        SINTH = SQRT(AB)/ABC
      ELSE
        COSTH = 1.0_dp
        SINTH = 0.0_dp
      ENDIF
      
      ! Y(1,0)
      
      Y(3) = CMPLX(sqrt(3.0_dp)*YLLR*COSTH, 0.0_dp, kind = dp)
      
      ! Y(1,1) ( = -DCONJG(Y(1,-1)))
      
      TEMP1 = -SQRT(1.5_dp)*YLLR*SINTH
      Y(4) = CMPLX(TEMP1*COSPH,TEMP1*SINPH, kind = dp)
      Y(2) = -CONJG(Y(4))
      
      DO L = 2, LMAX
        INDEX  = L*L + 1
        INDEX2 = INDEX + 2*L
        MSIGN  = 1 - 2*MOD(L,2)
        
        ! YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
        
        YL1L1R = DBLE(Y(INDEX-1))
        YL1L1I = DIMAG(Y(INDEX-1))
        TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))*SINTH
        YLLR = TEMP1*(COSPH*YL1L1R - SINPH*YL1L1I)
        YLLI = TEMP1*(COSPH*YL1L1I + SINPH*YL1L1R)
        Y(INDEX2) = CMPLX(YLLR,YLLI, kind = dp)
        Y(INDEX)  = cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
        !Y(INDEX)  = dble(MSIGN)*CONJG(Y(INDEX2))
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        
        ! YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
        ! (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
        
        TEMP2 = SQRT(DBLE(2*L+1))*COSTH
        YLL1R = TEMP2*YL1L1R
        YLL1I = TEMP2*YL1L1I
        Y(INDEX2) = CMPLX(YLL1R,YLL1I, kind = dp)
        Y(INDEX)  = -cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
  !      Y(INDEX)  = -dble(MSIGN)*CONJG(Y(INDEX2))
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        
        I4L2 = INDEX2 - 4*L + 2
        I2L  = INDEX2 - 2*L
        D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
        D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
        
        DO M = L - 2, 0, -1
          
          ! YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
          
          TEMP1 = 1.0_dp/SQRT(DBLE((L+M)*(L-M)))
          TEMP2 = D4LL1C*TEMP1
          TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
          YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
          YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
          Y(INDEX2) = CMPLX(YLMR,YLMI, kind = dp)
          Y(INDEX)  = cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
    !      Y(INDEX)  = dble(MSIGN)*CONJG(Y(INDEX2))
          
          MSIGN  = -MSIGN
          INDEX2 = INDEX2 - 1
          INDEX  = INDEX  + 1
          I4L2   = I4L2   - 1
          I2L    = I2L    - 1
        ENDDO
      ENDDO
      
  999 RETURN
  END subroutine getYlm
  
!----------------------------------------------------------------------------
  subroutine finalizeCalculation(nSys, crystalSystem, pot)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nSys
      !! Number of crystal systems

    type(crystal) :: crystalSystem(nSys)
       !! The crystal systems to get the
       !! matrix element for

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables:
    integer :: iT, isys
      !! Loop index


    deallocate(Ylm)

    do isys = 1, nSys
      deallocate(crystalSystem(isys)%nPWs1kGlobal)
      deallocate(crystalSystem(isys)%atomPositionsCart)
      deallocate(crystalSystem(isys)%iType)
      deallocate(crystalSystem(isys)%exp_iGDotR)
    enddo

    do iT = 1, pot%maxNAtomTypes
      deallocate(pot%atom(iT)%radGrid)
      deallocate(pot%atom(iT)%angMom)
      deallocate(pot%atom(iT)%F)
      deallocate(pot%atom(iT)%F1bra)
      deallocate(pot%atom(iT)%F1ket)
      deallocate(pot%atom(iT)%F2)
      deallocate(pot%atom(iT)%FI)
    enddo

    deallocate(pot%atom)


    call MPI_Barrier(worldComm, ierr)
    
    if(ionode) then

      write(*,'("-----------------------------------------------------------------")')
    
      call cpu_time(tf)
      write(*, '(" Total time needed:                         ", f10.2, " secs.")') tf-t0

    endif

    call MPI_Barrier(worldComm, ierr)
    
    return
    
  end subroutine finalizeCalculation
  
!----------------------------------------------------------------------------
  function getMatrixElementFNameWPath(ikGlobal, isp, path) result(fName)

    use miscUtilities, only: int2str, int2strLeadZero

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    character(*), intent(in) :: path
      !! Path to matrix element file

    ! Output variables:
    character(len=300) :: fName
      !! Matrix element file name


    fName = trim(path)//"/allElecOverlap."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))

  end function getMatrixElementFNameWPath
  
!----------------------------------------------------------------------------
  function getMatrixElementFName(ikGlobal, isp) result(fName)

    use miscUtilities, only: int2str, int2strLeadZero

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Output variables:
    character(len=300) :: fName
      !! Matrix element file name


    fName = "allElecOverlap."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))

  end function getMatrixElementFName
  
end module TMEmod
