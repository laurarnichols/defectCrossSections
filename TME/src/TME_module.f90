module TMEmod
  use constants, only: dp, pi, eVToHartree, ii
  use miscUtilities, only: int2str, int2strLeadZero, getUniqueInts
  use energyTabulatorMod, only: energyTableDir, readCaptureEnergyTable, readScatterEnergyTable
  use base, only: nKPoints, nSpins, order, ispSelect
  use cell, only: volume, recipLattVec

  use errorsAndMPI
  use mpi
  
  implicit none


  integer :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
    !! Optional band bounds
  integer, allocatable :: ibBra(:), ikBra(:), ikKet(:), ibKet(:)
    !! State indices for different systems
  integer :: ibShift_braket
    !! Shift of band indexing of bra system relative
    !! to ket system
  integer, allocatable :: mill_local(:,:)
    !! Local Miller indices
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors
  integer :: nPairs
    !! Number of pairs of bands to get overlaps for
  integer :: nSys
    !! Number of systems
  integer :: phononModeJ
    !! Index of phonon mode for the calculation
    !! of \(M_j\) (only for order=1)

  real(kind=dp) :: dq_j
    !! \(\delta q_j) for displaced wave functions
    !! (only order = 1)

  complex(kind=dp), allocatable :: Ylm(:,:)
    !! Spherical harmonics
  
  character(len=300) :: baselineDir
    !! Directory for baseline overlap to optionally
    !! be subtracted for the first-order term
  character(len=300) :: braExportDir, ketExportDir
    !! Path to export dirs
  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: optimalPairsDir
    !! Path to store or read optimalPairs.out file
  character(len=300) :: outputDir
    !! Path to where matrix elements should be output

  logical :: capture
    !! If considering capture as opposed to scattering 
  logical :: dqOnly
    !! If first-order matrix elements should only be
    !! divided by dq
  logical :: intraK
    !! If overlaps should be calculated across different
    !! k-points (true) or just between single k-points (false)
  logical :: lineUpBands
    !! If calculating preliminary overlaps to line up bands
    !! from different systems
  logical :: overlapOnly
    !! If only the wave function overlap should be
    !! calculated
  logical :: readOptimalPairs
    !! If optimal pairs should be read and states reordered
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
    integer :: iGkStart_pool
      !! Start and end G+k vector for process in pool
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
    integer :: nGkVecsLocal
      !! Local number of G+k vectors on this processor
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
    real(kind=dp) :: recipLattVec(3,3)
      !! Reciprocal-space lattice vectors
    real(kind=dp) :: volume
      !! Supercell volume

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

  type(crystal), allocatable :: crystalSystem(:)
    !! Array containing all crystal systems
  
  real(kind = dp) t0, tf
  
  
contains

!----------------------------------------------------------------------------
  subroutine readInputParams(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibBra, ikBra, ibKet, ikKet, ibShift_braket, &
          ispSelect, nPairs, order, phononModeJ, baselineDir, braExportDir, ketExportDir, dqFName, energyTableDir, &
          optimalPairsDir, outputDir, capture, dqOnly, intraK, lineUpBands, overlapOnly, readOptimalPairs, subtractBaseline)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
    integer, intent(out) :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
      !! Optional band bounds
    integer, allocatable, intent(out) :: ibBra(:), ikBra(:), ibKet(:), ikKet(:)
      !! State indices for different systems
    integer, intent(out) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: nPairs
      !! Number of pairs of bands to get overlaps for
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
    character(len=300), intent(out) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file
    character(len=300), intent(out) :: outputDir
      !! Path to where matrix elements should be output

    logical, intent(out) :: capture
      !! If considering capture as opposed to scattering 
    logical, intent(out) :: dqOnly
      !! If first-order matrix elements should only be
      !! divided by dq
    logical, intent(out) :: intraK
      !! If overlaps should be calculated across different
      !! k-points (true) or just between single k-points (false)
    logical, intent(out) :: lineUpBands
      !! If calculating preliminary overlaps to line up bands
      !! from different systems
    logical, intent(out) :: overlapOnly
      !! If only the wave function overlap should be
      !! calculated
    logical, intent(out) :: readOptimalPairs
      !! If optimal pairs should be read and states reordered
    logical, intent(out) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative

    ! Local variables:
    character(len=300) :: braBands, ketBands
      !! Strings to take input band pairs for overlap-only


    namelist /TME_Input/ ketExportDir, braExportDir, outputDir, energyTableDir, &
                         order, dqFName, phononModeJ, subtractBaseline, baselineDir, &
                         ispSelect, nPairs, braBands, ketBands, lineUpBands, overlapOnly, &
                         dqOnly, iBandLBra, iBandHBra, iBandLKet, iBandHKet, capture, &
                         intraK, ibShift_braket, optimalPairsDir, readOptimalPairs
    

    if(ionode) then
    
      call initialize(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibShift_braket, ispSelect, nPairs, order, &
              phononModeJ, baselineDir, braBands, ketBands, braExportDir, ketExportDir, dqFName, energyTableDir, &
              optimalPairsDir, outputDir, capture, dqOnly, intraK, lineUpBands, overlapOnly, readOptimalPairs, &
              subtractBaseline)
    
      read(5, TME_Input, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputParams', 'reading TME_Input namelist', abs(ierr))
    
      call checkInitialization(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibShift_braket, ispSelect, nPairs, order, phononModeJ, &
              baselineDir, braBands, ketBands, braExportDir, ketExportDir, dqFName, energyTableDir, optimalPairsDir, outputDir, &
              capture, dqOnly, intraK, lineUpBands, overlapOnly, readOptimalPairs, subtractBaseline, ibBra, ikBra, ibKet, ikKet)

    endif

    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(ibShift_braket, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ispSelect, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(phononModeJ, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(capture, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(dqOnly, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(intraK, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(lineUpBands, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(overlapOnly, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(readOptimalPairs, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(subtractBaseline, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(braExportDir, len(braExportDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(ketExportDir, len(ketExportDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(baselineDir, len(baselineDir), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(optimalPairsDir, len(optimalPairsDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(nPairs, 1, MPI_INTEGER, root, worldComm, ierr)
    if(.not. ionode) then
      ! Every case allows for multiple ket band states
      allocate(ibKet(nPairs))

      ! Capture only allows a single bra band state
      if(capture) then
        allocate(ibBra(1))
      else  
        allocate(ibBra(nPairs))
      endif

      ! For capture or non-intra-k overlaps, we do not use the
      ! ikBra or ikKet arrays because we just loop over k-points.
      ! Otherwise, the k-points for each pair/transition are 
      ! specified from the energy table file.
      if(capture .or. (overlapOnly .and. .not. intraK)) then
        allocate(ikBra(1), ikKet(1))
      else
        allocate(ikBra(nPairs), ikKet(nPairs))
      endif
    endif

    call MPI_BCAST(ibBra, size(ibBra), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ibKet, size(ibKet), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikBra, size(ikBra), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikKet, size(ikKet), MPI_INTEGER, root, worldComm, ierr)

    return
    
  end subroutine readInputParams
  
!----------------------------------------------------------------------------
  subroutine initialize(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibShift_braket, ispSelect, nPairs, order, &
          phononModeJ, baselineDir, braBands, ketBands, braExportDir, ketExportDir, dqFName, energyTableDir, &
          optimalPairsDir, outputDir, capture, dqOnly, intraK, lineUpBands, overlapOnly, readOptimalPairs, &
          subtractBaseline)
    
    implicit none

    ! Output variables:
    integer, intent(out) :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
      !! Optional band bounds
    integer, intent(out) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: nPairs
      !! Number of pairs of bands to get overlaps for
    integer, intent(out) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(out) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(out) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term
    character(len=300), intent(out) :: braBands, ketBands
      !! Strings to take input band pairs for overlap-only
    character(len=300), intent(out) :: braExportDir, ketExportDir
      !! Path to export dirs
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(out) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file
    character(len=300), intent(out) :: outputDir
      !! Path to where matrix elements should be output

    logical, intent(out) :: capture
      !! If considering capture as opposed to scattering 
    logical, intent(out) :: dqOnly
      !! If first-order matrix elements should only be
      !! divided by dq
    logical, intent(out) :: intraK
      !! If overlaps should be calculated across different
      !! k-points (true) or just between single k-points (false)
    logical, intent(out) :: lineUpBands
      !! If calculating preliminary overlaps to line up bands
      !! from different systems
    logical, intent(out) :: overlapOnly
      !! If only the wave function overlap should be
      !! calculated
    logical, intent(out) :: readOptimalPairs
      !! If optimal pairs should be read and states reordered
    logical, intent(out) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative
    

    iBandLBra = -1
    iBandHBra = -1
    iBandLKet = -1
    iBandHKet = -1
    ibShift_braket = 0
    order = -1
    nPairs = -1
    ispSelect = -1
    phononModeJ = -1

    braBands = ''
    ketBands = ''
    braExportDir = ''
    ketExportDir = ''
    energyTableDir = ''
    dqFName = ''
    optimalPairsDir = './optimalPairs'
    outputDir = './TMEs'
    baselineDir = ''
    
    capture = .true.
    dqOnly = .false.
    intraK = .false.
    lineUpBands = .false.
    overlapOnly = .false.
    readOptimalPairs = .false.
    subtractBaseline = .false.
    
    return
    
  end subroutine initialize
  
!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibShift_braket, ispSelect, nPairs, order, &
          phononModeJ, baselineDir, braBands, ketBands, braExportDir, ketExportDir, dqFName, energyTableDir, optimalPairsDir, &
          outputDir, capture, dqOnly, intraK, lineUpBands, overlapOnly, readOptimalPairs, subtractBaseline, ibBra, ikBra, &
          ibKet, ikKet)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
      !! Optional band bounds
    integer, intent(in) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(inout) :: nPairs
      !! Number of pairs of bands to get overlaps for
    integer, intent(in) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(in) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    character(len=300), intent(in) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term
    character(len=300), intent(in) :: braBands, ketBands
      !! Strings to take input band pairs for overlap-only
    character(len=300), intent(in) :: braExportDir, ketExportDir
      !! Path to export dirs
    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(in) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file
    character(len=300), intent(in) :: outputDir
      !! Path to where matrix elements should be output

    logical, intent(inout) :: capture
      !! If considering capture as opposed to scattering 
    logical, intent(in) :: dqOnly
      !! If first-order matrix elements should only be
      !! divided by dq
    logical, intent(inout) :: intraK
      !! If overlaps should be calculated across different
      !! k-points (true) or just between single k-points (false)
    logical, intent(inout) :: lineUpBands
      !! If calculating preliminary overlaps to line up bands
      !! from different systems
    logical, intent(inout) :: overlapOnly
      !! If only the wave function overlap should be
      !! calculated
    logical, intent(in) :: readOptimalPairs
      !! If optimal pairs should be read and states reordered
    logical, intent(in) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative

    ! Output variables:
    integer, allocatable, intent(out) :: ibBra(:), ikBra(:), ibKet(:), ikKet(:)
      !! State indices for different systems
    
    ! Local variables
    logical :: abortExecution
      !! If program should stop
    logical :: bandBoundsGiven 
      !! If optional band bounds given
    logical :: energyTableGiven
      !! If energy table given
    logical :: explicitBandStringsGiven
      !! If strings with pairs of bands given;
      !! supercedes other choices of bands
    logical :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel
    
    
    write(*,'("Inputs: ")')
    
    abortExecution = checkDirInitialization('braExportDir', braExportDir, 'input')
    abortExecution = checkDirInitialization('ketExportDir', ketExportDir, 'input') .or. abortExecution
    abortExecution = checkStringInitialization('outputDir', outputDir) .or. abortExecution


    if(ispSelect < 1 .or. ispSelect > 2) then
      write(*,*) "No valid choice for spin channel selection given. Looping over spin."
      loopSpins = .true.
    else
      write(*,'("Only exporting spin channel ", i2)') ispSelect
      loopSpins = .false.
    endif


    write(*,'("readOptimalPairs = ",L1)') readOptimalPairs
    if(readOptimalPairs .and. lineUpBands) then
      write(*,'("readOptimalPairs option selected. Overriding lineUpBands and setting to false!!")')
      lineUpBands = .false.
    endif

    write(*,'("lineUpBands = ",L)') lineUpBands


    if(lineUpBands .and. .not. overlapOnly) then
      write(*,'("lineUpBands option selected. Overriding overlapOnly and setting to true!!")')
      overlapOnly = .true.
    endif

    write(*,'("overlapOnly = ",L)') overlapOnly


    if(overlapOnly .and. capture) then
      write(*,'("Both overlapOnly and capture are true. Overlap only takes precedence!")')
      capture = .false.
    endif
    write(*,'("capture = ",L)') capture


    if(readOptimalPairs .or. lineUpBands) then
      abortExecution = checkStringInitialization('optimalPairsDir', optimalPairsDir) .or. abortExecution

      if(capture) then
        write(*,'("Capture into level in the gap is not compatible with optimizing band matching!")')
        abortExecution = .true.
      endif
    endif



    write(*,'("ibShift_braket = ",i10)') ibShift_braket
    write(*,'("intraK = ",L)') intraK
    
    if(.not. overlapOnly) then
      ! Require that states be read from the energy table for capture and scattering

      abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution

      if(order == 1) then
        abortExecution = checkFileInitialization('dqFName', dqFName) .or. abortExecution
        abortExecution = checkIntInitialization('phononModeJ', phononModeJ, 1, int(1e9)) .or. abortExecution

        write(*,'("dqOnly = ",L1)') dqOnly
        write(*,'("subtractBaseline = ",L1)') subtractBaseline

        if(subtractBaseline) then
          if(ispSelect == 2) then
            abortExecution = checkDirInitialization('baselineDir', baselineDir, 'allElecOverlap.2.1') .or. abortExecution
          else
            abortExecution = checkDirInitialization('baselineDir', baselineDir, 'allElecOverlap.1.1') .or. abortExecution
          endif
        endif
      endif

      call getBandBoundsFromEnergyTable(ispSelect, capture, intraK, overlapOnly, loopSpins, energyTableDir, &
          ibBra, ikBra, ibKet, ikKet, nPairs, abortExecution)

    else
      ! For overlap-only, can get the bands from explicit strings,
      ! band bounds, or the energy table

      explicitBandStringsGiven = trim(braBands) /= '' .and. trim(ketBands) /= '' 
      bandBoundsGiven = iBandLBra > 0 .and. iBandHBra > 0 .and. &
                        iBandLKet > 0 .and. iBandHKet > 0
      energyTableGiven = trim(energyTableDir) /= ''


      if(lineUpBands) then
        if(.not. bandBoundsGiven) then
          write(*,'("Must give a range of bands with lineUpBands option!")')
          abortExecution = .true.
        else if(iBandHBra - iBandLBra /= iBandHKet - iBandLKet) then
          write(*,'("Must have an equal number of bra and ket states to lineUpBands!")')
          abortExecution = .true.
        endif
        
        if(intraK) then
          write(*,'("intraK options is not compatible with lineUpBands option!")')
          write(*,'("Overriding and setting intraK = .false.!")')
          intraK = .false.
        endif
      endif


      if(explicitBandStringsGiven .and. .not. intraK) then

        write(*,'("Reading states from band strings braBands and ketBands.")')
        if(bandBoundsGiven) write(*,'("Band bounds detected; will be ignored.")')
        if(energyTableGiven) write(*,'("Input detected for energy table; will be ignored.")')

        call getBandBoundsFromBandStrings(nPairs, braBands, ketBands, ibBra, ibKet, abortExecution)
        allocate(ikBra(1),ikKet(1))
          ! Allocate the k states to avoid issues with broadcasting

      else if(bandBoundsGiven .and. .not. intraK) then

        write(*,'("Reading states from band ranges.")')
        call getBandBoundsFromRanges(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibBra, ibKet, nPairs, abortExecution)
        allocate(ikBra(1),ikKet(1))
          ! Allocate the k states to avoid issues with broadcasting

      else if(energyTableGiven) then

        write(*,'("Reading states from energy table.")')
        call getBandBoundsFromEnergyTable(ispSelect, capture, intraK, overlapOnly, loopSpins, energyTableDir, &
            ibBra, ikBra, ibKet, ikKet, nPairs, abortExecution)
          ! Reads energy table based on `capture` variable

      else if((explicitBandStringsGiven .or. bandBoundsGiven) .and. intraK) then
        write(*,'("Explicit band strings or ranges not currently supported for intraK")')
        abortExecution = .true.
      else
        write(*,'("No input detected for reading band states!")')
        abortExecution = .true.
      endif

    endif


    call system('mkdir -p '//trim(outputDir))

    
    if(abortExecution) then
      write(*,'(" Program stops!")')
      stop
    endif
    
    return
    
  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine getBandBoundsFromEnergyTable(ispSelect, capture, intraK, overlapOnly, loopSpins, energyTableDir, &
            ibBra, ikBra, ibKet, ikKet, nPairs, abortExecution)
    ! Read nTransitions and band indices from energy table, assuming that
    ! they are the same for the two spin channels and all k-points

    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user

    logical, intent(in) :: capture
      !! If considering capture as opposed to scattering 
    logical, intent(in) :: intraK
      !! If overlaps should be calculated across different
      !! k-points (true) or just between single k-points (false)
    logical, intent(in) :: overlapOnly
      !! If only the wave function overlap should be
      !! calculated
    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables

    ! Output variables:
    integer, allocatable, intent(out) :: ibBra(:), ikBra(:), ibKet(:), ikKet(:)
      !! State indices for different systems
    integer, intent(out) :: nPairs
      !! Number of pairs of bands to get overlaps for

    logical, intent(inout) :: abortExecution
      !! If program should stop

    ! Local variables:
    integer :: ibBra1
      !! Single band index for reading from capture table
    integer :: isp
      !! Spin channel index

    real(kind=dp), allocatable :: rDum(:,:)
      !! Dummy real to ignore energy when just getting 
      !! band bounds and number of transitions


    if(loopSpins) then
      isp = 1
    else 
      isp = ispSelect
    endif


    if(capture .or. (overlapOnly .and. .not. intraK)) then

      abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.'//trim(int2str(isp))//'.1') &
                              .or. abortExecution

      if(.not. abortExecution) then

        allocate(ikBra(1), ikKet(1))

        call readCaptureEnergyTable(1, isp, energyTableDir, ibKet, ibBra1, nPairs, rDum)
          ! nPairs is read from the energy table

        if(capture) then
          allocate(ibBra(1))
        else
          allocate(ibBra(nPairs))
        endif

        ibBra(:) = ibBra1

      endif
    else
      abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.'//trim(int2str(isp))) &
                              .or. abortExecution

      if(.not. abortExecution) call readScatterEnergyTable(isp, .true., energyTableDir, ibKet, ibBra, ikKet, ikBra, nPairs, rDum)
        ! nPairs is read from the energy table. 
        ! Pass true to get matrix-element and delta-function energies.
    endif

    return

  end subroutine getBandBoundsFromEnergyTable

!----------------------------------------------------------------------------
  subroutine getBandBoundsFromBandStrings(nPairs, braBands, ketBands, ibBra, ibKet, abortExecution)

    implicit none

    ! Input variables:
    integer, intent(in) :: nPairs
      !! Number of pairs of bands to get overlaps for

    character(len=300), intent(in) :: braBands, ketBands
      !! Strings to take input band pairs for overlap-only

    ! Output variables:
    integer, allocatable, intent(out) :: ibBra(:), ibKet(:)
      !! Band indices for different systems

    logical, intent(inout) :: abortExecution
      !! If program should stop

    ! Local variables:
    integer :: ip
      !! Loop index


    abortExecution = checkStringInitialization('braBands', braBands) .or. abortExecution
    abortExecution = checkStringInitialization('ketBands', ketBands) .or. abortExecution
    abortExecution = checkIntInitialization('nPairs', nPairs, 1, 100) .or. abortExecution
      ! nPairs must be given explicitly


    if(.not. abortExecution) then
      allocate(ibBra(nPairs), ibKet(nPairs))

      read(braBands,*,iostat=ierr) (ibBra(ip), ip=1,nPairs)
      if(ierr /= 0) call exitError('getBandBoundsFromStrings', 'reading braBands', abs(ierr))

      read(ketBands,*) (ibKet(ip), ip=1,nPairs)
      if(ierr /= 0) call exitError('getBandBoundsFromStrings', 'reading ketBands', abs(ierr))
    endif

    return

  end subroutine getBandBoundsFromBandStrings

!----------------------------------------------------------------------------
  subroutine getBandBoundsFromRanges(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ibBra, ibKet, nPairs, abortExecution)
     
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
      !! Optional band bounds

    ! Output variables:
    integer, allocatable, intent(out) :: ibBra(:), ibKet(:)
      !! Band indices for different systems
    integer, intent(out) :: nPairs
      !! Number of pairs of bands to get overlaps for

    logical, intent(inout) :: abortExecution
      !! If program should stop

    ! Local variables:
    integer :: ip, ibB, ibK
      !! Loop indices



    abortExecution = checkIntInitialization('iBandLBra', iBandLBra, 1, int(1e9))
    abortExecution = checkIntInitialization('iBandHBra', iBandHBra, iBandLBra, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandLKet', iBandLKet, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandHKet', iBandHKet, iBandLKet, int(1e9)) .or. abortExecution 

    if(.not. abortExecution) then
      nPairs = (iBandHBra - iBandLBra + 1)*(iBandHKet - iBandLKet + 1)
        ! Calculate nPairs from input

      allocate(ibBra(nPairs), ibKet(nPairs))

      ip = 0
      do ibB = iBandLBra, iBandHBra
        do ibK = iBandLKet, iBandHKet
          ip = ip + 1
          ibBra(ip) = ibB
          ibKet(ip) = ibK
        enddo
      enddo
    endif

    return

  end subroutine getBandBoundsFromRanges

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
  subroutine completePreliminarySetup(nSys, order, phononModeJ, capture, intraK, dqFName, mill_local, nGVecsGlobal, &
        nGVecsLocal, nKPoints, nSpins, dq_j, recipLattVec, volume, Ylm, crystalSystem, pot)

    implicit none

    ! Input variables:
    integer, intent(in) :: nSys
      !! Number of systems to read files for
    integer, intent(in) :: order
      !! Order of matrix element (0 or 1)
    integer, intent(in) :: phononModeJ
      !! Index of phonon mode for the calculation
      !! of \(M_j\) (only for order=1)

    logical, intent(in) :: capture
      !! If considering capture as opposed to scattering 
    logical, intent(in) :: intraK
      !! If overlaps should be calculated across different
      !! k-points (true) or just between single k-points (false)

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms

    ! Output variables:
    integer, allocatable, intent(out) :: mill_local(:,:)
      !! Local Miller indices
    integer, intent(out) :: nGVecsGlobal
      !! Number of global G-vectors (tested to be consistent
      !! across all systems)
    integer, intent(out) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(out) :: nKPoints
      !! Number of k-points (tested to be consistent
      !! across all systems)
    integer, intent(out) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp), intent(out) :: dq_j
      !! \(\delta q_j) for displaced wave functions
      !! (only order = 1)
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal-space lattice vectors
    real(kind=dp), intent(out) :: volume
      !! Cell volume (tested to be consistent
      !! across all systems)

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

    real(kind = dp) :: t1, t2
      !! For timing different processes


    if(ionode) write(*, '("Pre-k-loop: [ ] Read inputs  [ ] Set up tables ")')
    call cpu_time(t1)

    ! Initialize global variables to be tracked to make sure they
    ! are consistent across all of the systems
    nKPoints = -1
    nGVecsGlobal = -1
    nSpins = 0
    volume = -1.0

    do isys = 1, nSys

      call readInputFileSkipPseudo(intraK, nGVecsGlobal, nKPoints, nSpins, volume, crystalSystem(isys))

    enddo

    recipLattVec = crystalSystem(1)%recipLattVec
      ! Could, in principle, check that all of the systems are consistent,
      ! but I don't want to do that right now. Right now, the code just 
      ! assumes that they are consistent.


    call getGlobalPseudo(nSys, crystalSystem, pot)
    
    pot%maxAngMom = 2*pot%maxAngMom + 1


    if(order == 1) then
      call readDqFile(phononModeJ, dqFName, dq_j)
      call MPI_BCAST(dq_j, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    endif

    
    ! Distribute k-points in pools
    ! No parallelization over k-points for scattering and intra-k overlap only
    if(capture .or. (overlapOnly .and. .not. intraK)) then
      call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    endif

    ! Distribute G-vectors across processes in pool
    call distributeItemsInSubgroups(indexInPool, nGVecsGlobal, nProcPerPool, nProcPerPool, nProcPerPool, iGStart_pool, &
        iGEnd_pool, nGVecsLocal)


    allocate(mill_local(3,nGVecsLocal))

    call getFullPWGrid(nGVecsLocal, nGVecsGlobal, crystalSystem(1), mill_local)
      ! We have already verified at this point that the number of 
      ! G-vectors is the same in all of the systems, so it doesn't
      ! matter which system we read the PW grid from.


    call cpu_time(t2)
    if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [ ] Set up tables (",f10.2," secs)")') t2-t1
    call cpu_time(t1)


    allocate(Ylm((pot%maxAngMom+1)**2,nGVecsLocal))

    call setUpTables(nGVecsLocal, mill_local, nSys, recipLattVec, pot, Ylm, crystalSystem)

    deallocate(mill_local)

  
    call cpu_time(t2)
    if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [X] Set up tables (",f10.2," secs)")') t2-t1
    call cpu_time(t1)

    return

  end subroutine completePreliminarySetup
  
!----------------------------------------------------------------------------
  subroutine readInputFileSkipPseudo(intraK, nGVecsGlobal, nKPoints, nSpins, volume, sys)

    use miscUtilities, only: getFirstLineWithKeyword, ignoreNextNLinesFromFile
    
    implicit none

    ! Input variables:
    logical, intent(in) :: intraK
      !! If overlaps should be calculated across different
      !! k-points (true) or just between single k-points (false)

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

    real(kind=dp) :: volume
      !! Cell volume (tested to be consistent
      !! across all systems)

    type(crystal) :: sys
      !! The crystal system

    ! Local variables:
    integer :: iDum
      !! Dummy integer
    integer :: ik, iA, ix, iT
      !! Loop indices
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
      read(50,*) sys%volume


      if(volume < 0) then
        volume = sys%volume
      else if(abs(volume - sys%volume) > 1e-8) then
        call exitError('readInput', 'volumes don''t match', 1)
      endif

    endif

    call MPI_BCAST(sys%volume, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(volume, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)


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


      if(.not. intraK) then
        if(nKPoints < 0) then
          nKPoints = sys%nKPoints
        else if(sys%nKPoints /= nKPoints) then
          call exitError('readInput', 'number of k-points in system '//trim(sys%ID)//' does not match', 1)
        end if
      endif
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
      read(50,*) ! nBands

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
              call exitError('getGlobalPseudo', 'Elements must be consistent across all systems!', 1)
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
    real(kind = dp) :: t1, t2
      !! For timing different processes

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
  subroutine getAndWriteInterKOnlyOverlaps(iBandLBra, iBandHBra, iBandLKet, iBandHKet, nPairs, ibShift_braket, &
          ibBra, ibKet, ispSelect, nGVecsLocal, nSpins, volume, lineUpBands, readOptimalPairs, optimalPairsDir, &
          braSys, ketSys, pot)

    use optimalBandMatching, only: findOptimalPairsAndOutput, readAllOptimalPairs

    implicit none

    ! Input variables:
    integer :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
      !! Optional band bounds to be input for lineUpBands option
      !! or read from optimalPairs file for readOptimalPairs
      !! Removing the intent because what we do with these variables 
      !! changes and I don't want it to yell at me
    integer, intent(in) :: nPairs
      !! Number of pairs of bands to get overlaps for
    integer, intent(in) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(in) :: ibBra(nPairs), ibKet(nPairs)
      !! Band indices for different systems
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    logical, intent(in) :: lineUpBands
      !! If calculating preliminary overlaps to line up bands
      !! from different systems
    logical, intent(in) :: readOptimalPairs
      !! If optimal pairs should be read and states reordered

    character(len=300), intent(in) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file

    type(crystal) :: braSys, ketSys
       !! The crystal systems to get the
       !! matrix element for

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables:
    integer, allocatable :: ibBra_optimal(:,:)
      !! Optimal indices from the bra system corresponding 
      !! to the input indices from the ket system
    integer :: ikLocal, ikGlobal, ip, isp
      !! Loop index

    real(kind = dp) :: t1, t2
      !! For timing different processes

    complex(kind=dp) :: Ufi(nPairs,nSpins), Ufi_ip(nSpins)
      !! All-electron overlap

    logical :: spin1Skipped, spin2Skipped
      !! If spin channels skipped


    Ufi(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)

    if(ionode .and. lineUpBands) call system('mkdir -p '//trim(optimalPairsDir))

    do ikLocal = 1, nkPerPool
    
      if(ionode) write(*,'("Beginning k-point loop ", i4, " of ", i4)') ikLocal, nkPerPool
      ikGlobal = ikLocal+ikStart_pool-1
        !! Get the global `ik` index from the local one
    

      if(.not. thisKComplete(ikGlobal, ispSelect, nSpins)) then

        if(ionode) write(*,'("  Spin-independent setup")')
        call cpu_time(t1)

        if(indexInPool == 0) then
          if(braSys%nPWs1kGlobal(ikGlobal) /= ketSys%nPWs1kGlobal(ikGlobal)) &
            call exitError('getAndWriteInterKOnlyOverlaps', 'number of G+k vectors does not match for ik='//trim(int2str(ikGlobal)), 1)
        endif

        call getSpinSkipped(ikGlobal, ispSelect, spin1Skipped, spin2Skipped)

        call spinAndBandIndependentSetup(ikGlobal, nGVecsLocal, braSys)
        call spinAndBandIndependentSetup(ikGlobal, nGVecsLocal, ketSys)

        call cpu_time(t2)
        if(ionode) write(*, '("  Spin independent setup complete! (",f10.2," secs)")') t2-t1

        
        if(readOptimalPairs) then
          if(indexInPool == 0) &
            call readAllOptimalPairs(ikGlobal, nSpins, spin1Skipped, spin2Skipped, optimalPairsDir, iBandLKet, iBandHKet, &
                  ibBra_optimal)

          call MPI_BCAST(iBandLKet, 1, MPI_INTEGER, root, intraPoolComm, ierr)
          call MPI_BCAST(iBandHKet, 1, MPI_INTEGER, root, intraPoolComm, ierr)

          if(indexInPool /= 0) allocate(ibBra_optimal(nSpins,iBandLKet:iBandHKet))
          call MPI_BCAST(ibBra_optimal, size(ibBra_optimal), MPI_INTEGER, root, intraPoolComm, ierr)
        else
          allocate(ibBra_optimal(1,1))
        endif


        do ip = 1, nPairs

          if(ionode) write(*,'("  Beginning overlap <", i5, "|",i5">")') ibBra(ip), ibKet(ip)
          call cpu_time(t1)

          if(indexInPool == 0 .and. readOptimalPairs .and. (ibBra(ip) < iBandLKet .or. ibBra(ip) > iBandHKet)) &
            call exitError('getAndWriteInterKOnlyOverlaps',&
              'Index '//trim(int2str(ibBra(ip)))//' not included in optimalPairs file for ik = '//trim(int2str(ikGlobal)),1)

          call calculateBandPairOverlap(ibShift_braket, ibBra(ip), ibKet(ip), ikGlobal, ikGlobal, nSpins, &
                ibBra_optimal(:,ibBra(ip)), nGVecsLocal, volume, readOptimalPairs, spin1Skipped, spin2Skipped, &
                braSys, ketSys, pot, Ufi_ip)

          Ufi(ip,:) = Ufi_ip

          call cpu_time(t2)
          if(ionode) write(*, '("  Overlap <",i5,"|",i5,"> complete! (",f10.2," secs)")') ibBra(ip), ibKet(ip), t2-t1
        enddo


        if(indexInPool == 0) then 
          do isp = 1, nSpins
            if((isp == 1 .and. .not. spin1Skipped) .or. (isp == 2 .and. .not. spin2Skipped)) then
              if(lineUpBands) then
                call findOptimalPairsAndOutput(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ikLocal, isp, nPairs, Ufi(:,isp), &
                    optimalPairsDir)
              endif

              call writeInterKOverlaps(nPairs, ibBra, ibKet, ikLocal, isp, volume, Ufi(:,isp))
            endif
          enddo
        endif


        call deallocateSysArrays(braSys)
        call deallocateSysArrays(ketSys)
        deallocate(ibBra_optimal)

      endif ! If both spin channels exist
    enddo ! k-point loop

    return

  end subroutine getAndWriteInterKOnlyOverlaps

!----------------------------------------------------------------------------
  subroutine getAndWriteCaptureMatrixElements(nTransitions, ibShift_braket, ibi, ibf, ispSelect, nGVecsLocal, nSpins, dq_j, &
          volume, dqOnly, braSys, ketSys, pot)

    implicit none

    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(in) :: ibi(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: ibf
      !! Final-state index
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp), intent(in) :: dq_j
      !! \(\delta q_j) for displaced wave functions
      !! (only order = 1)
    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    logical, intent(in) :: dqOnly
      !! If first-order matrix elements should only be
      !! divided by dq

    type(crystal) :: braSys, ketSys
       !! The crystal systems to get the
       !! matrix element for

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables 
    integer :: iDum(nSpins)
      !! Dummy integer array to pass to calculateBandPairOverlap
      !! for ibBra_optimal
    integer :: ikLocal, ikGlobal, isp, iE
      !! Loop indices

    real(kind = dp) :: t1, t2
      !! For timing different processes

    complex(kind=dp) :: Ufi(nTransitions,nSpins), Ufi_iE(nSpins)
      !! All-electron overlap

    logical :: spin1Skipped, spin2Skipped
      !! If spin channels skipped


    Ufi(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)

    do ikLocal = 1, nkPerPool
    
      if(ionode) write(*,'("Beginning k-point loop ", i4, " of ", i4)') ikLocal, nkPerPool
      ikGlobal = ikLocal+ikStart_pool-1
        !! Get the global `ik` index from the local one


      if(.not. thisKComplete(ikGlobal, ispSelect, nSpins)) then

        if(ionode) write(*,'("  Spin-independent setup")')
        call cpu_time(t1)

        if(indexInPool == 0) then
          if(braSys%nPWs1kGlobal(ikGlobal) /= ketSys%nPWs1kGlobal(ikGlobal)) &
            call exitError('getAndWriteCaptureMatrixElements', 'number of G+k vectors does not match for ik='//trim(int2str(ikGlobal)), 1)
        endif

        call getSpinSkipped(ikGlobal, ispSelect, spin1Skipped, spin2Skipped)

        call spinAndBandIndependentSetup(ikGlobal, nGVecsLocal, braSys)
        call spinAndBandIndependentSetup(ikGlobal, nGVecsLocal, ketSys)


        call cpu_time(t2)
        if(ionode) write(*, '("  Spin independent setup complete! (",f10.2," secs)")') t2-t1
        

        do iE = 1, nTransitions
          if(ionode) write(*,'("  Beginning transition ", i5, " -> ",i5)') ibi(iE), ibf
          call cpu_time(t1)

          call calculateBandPairOverlap(ibShift_braket, ibf, ibi(iE), ikGlobal, ikGlobal, nSpins, iDum, nGVecsLocal, volume, &
                .false., spin1Skipped, spin2Skipped, braSys, ketSys, pot, Ufi_iE)
            ! Pass false for readOptimalPairs and dummy integer for ibBra_optimal
            ! because that option isn't allowed for capture

          Ufi(iE,:) = Ufi_iE

          call cpu_time(t2)
          if(ionode) write(*, '("  Transition ",i5," -> ",i5," complete! (",f10.2," secs)")') ibi(iE), ibf, t2-t1
        enddo 


        ! Subtract baseline if applicable and write out results
        if(indexInPool == 0) then 
          do isp = 1, nSpins
            if((isp == 1 .and. .not. spin1Skipped) .or. (isp == 2 .and. .not. spin2Skipped)) then
                if(order == 1 .and. subtractBaseline) &
                  call readAndSubtractBaseline(ikGlobal, isp, nTransitions, Ufi(:,isp))
        
                call writeCaptureMatrixElements(nTransitions, ibi, ibf, ikLocal, isp, dq_j, volume, Ufi(:,isp), dqOnly)
            endif 
          enddo
        endif

        call deallocateSysArrays(braSys)
        call deallocateSysArrays(ketSys)

      endif ! If both spin channels exist
    enddo ! k-point loop

    return

  end subroutine getAndWriteCaptureMatrixElements

!----------------------------------------------------------------------------
  subroutine getAndWriteScatterMatrixElementsOrOverlaps(nTransitions, ibShift_braket, ibi, iki, ibf, ikf, ispSelect, &
            nGVecsLocal, nSpins, volume, overlapOnly, readOptimalPairs, optimalPairsDir, braSys, ketSys, pot)

    use optimalBandMatching, only: readAllOptimalPairs

    implicit none

    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(in) :: ibi(nTransitions), iki(nTransitions), ibf(nTransitions), ikf(nTransitions)
      !! State indices
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    logical, intent(in) :: overlapOnly
      !! If only the wave function overlap should be
      !! calculated
    logical, intent(in) :: readOptimalPairs
      !! If optimal pairs should be read and states reordered

    character(len=300), intent(in) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file

    type(crystal) :: braSys, ketSys
       !! The crystal systems to get the
       !! matrix element for

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Local variables 
    integer :: iBandLKet, iBandHKet
      !! Band bounds from optimalPairs file
    integer, allocatable :: ibBra_optimal(:,:)
      !! Optimal index from the bra system corresponding 
      !! to the input index from the ket system
    integer :: isp, iE, iU_iki, iU_ikf
      !! Loop indices
    integer, allocatable :: ikiUnique(:), ikfUnique(:)
      !! Unique initial and final k-points
    integer :: nUnique_iki, nUnique_ikf
      !! Number of unique initial and final k-points
    integer :: updateFrequency
      !! Number of transitions to update after

    real(kind = dp) :: t1, t2
      !! For timing different processes

    complex(kind=dp) :: Ufi(nTransitions,nSpins), Ufi_iE(nSpins)
      !! All-electron overlap

    logical :: spin1Skipped, spin2Skipped
      !! If spin channels skipped


    Ufi(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)

    call getSpinSkipped(-1, ispSelect, spin1Skipped, spin2Skipped)
      ! Pass in ikGlobal = -1 to trigger scattering format
      
    if(.not. (spin1Skipped .and. spin2Skipped)) then

      updateFrequency = ceiling(nTransitions/10.0)

      ! First, get the unique k-points for each system
      if(ionode) call getUniqueInts(nTransitions, ikf, nUnique_ikf, ikfUnique)
      call MPI_BCAST(nUnique_ikf, 1, MPI_INTEGER, root, worldComm, ierr)
      if(.not. ionode) allocate(ikfUnique(nUnique_ikf))
      call MPI_BCAST(ikfUnique, nUnique_ikf, MPI_INTEGER, root, worldComm, ierr)

      if(ionode) call getUniqueInts(nTransitions, iki, nUnique_iki, ikiUnique)
      call MPI_BCAST(nUnique_iki, 1, MPI_INTEGER, root, worldComm, ierr)
      if(.not. ionode) allocate(ikiUnique(nUnique_iki))
      call MPI_BCAST(ikiUnique, nUnique_iki, MPI_INTEGER, root, worldComm, ierr)


      ! Have an outside loop over the unique k-points. The EnergyTabulator
      ! code uses a range of k-points and bands and loops over the initial
      ! k-points/bands then the final k-points/bands. Looping over the unique
      ! initial k-points then the unique final k-points and calculating only
      ! the transitions corresponding to those will keep us from duplicating
      ! calculations unnecessarily, and it will match the order from the 
      ! energy table.
      do iU_iki = 1, nUnique_iki


        if(ionode) write(*,'("Spin-independent setup for ket sys, iki =",i5)') ikiUnique(iU_iki)
        call cpu_time(t1)

        call spinAndBandIndependentSetup(ikiUnique(iU_iki), nGVecsLocal, ketSys)

        call cpu_time(t2)
        if(ionode) write(*,'("Spin independent setup for ket sys, iki =",i5," complete! (",f10.2," secs)")') ikiUnique(iU_iki), t2-t1


        do iU_ikf = 1, nUnique_ikf


          if(ionode) write(*,'("  Spin-independent setup for bra sys, ikf =",i5)') ikfUnique(iU_ikf)
          call cpu_time(t1)

          call spinAndBandIndependentSetup(ikfUnique(iU_ikf), nGVecsLocal, braSys)

          call cpu_time(t2)
          if(ionode) write(*,'("  Spin independent setup for bra sys, iki =",i5," complete! (",f10.2," secs)")') ikfUnique(iU_ikf), t2-t1


          if(ionode) write(*,'("    Beginning transitions iki ", i5," -> ",i5)') ikiUnique(iU_iki), ikfUnique(iU_ikf) 
          call cpu_time(t1)

        
          if(readOptimalPairs) then
            if(ionode) &
              call readAllOptimalPairs(ikfUnique(iU_ikf), nSpins, spin1Skipped, spin2Skipped, optimalPairsDir, iBandLKet, iBandHKet, &
                    ibBra_optimal)

            call MPI_BCAST(iBandLKet, 1, MPI_INTEGER, root, worldComm, ierr)
            call MPI_BCAST(iBandHKet, 1, MPI_INTEGER, root, worldComm, ierr)

            if(.not. ionode) allocate(ibBra_optimal(nSpins,iBandLKet:iBandHKet))
            call MPI_BCAST(ibBra_optimal, size(ibBra_optimal), MPI_INTEGER, root, worldComm, ierr)
          else
            allocate(ibBra_optimal(1,1))
          endif


          do iE = 1, nTransitions

            if(ionode .and. mod(iE,updateFrequency) == 0) then
              call cpu_time(t2)
              write(*,'("    ", i2,"% complete with loop over transitions. Time in loop: ",f10.2," secs")') &
                      int((iE*100.0)/nTransitions), t2-t1
            endif

            if((iki(iE) == ikiUnique(iU_iki)) .and. (ikf(iE) == ikfUnique(iU_ikf))) then

              if(ionode .and. readOptimalPairs .and. (ibf(iE) < iBandLKet .or. ibf(iE) > iBandHKet)) &
                call exitError('getAndWriteScatterMatrixElementsOrOverlaps',&
                  'Index '//trim(int2str(ibf(iE)))//' not included in optimalPairs file for ik = '//trim(int2str(ikf(iE))),1)

              call calculateBandPairOverlap(ibShift_braket, ibf(iE), ibi(iE), ikf(iE), iki(iE), nSpins, &
                    ibBra_optimal(:,ibf(iE)), nGVecsLocal, volume, readOptimalPairs, spin1Skipped, &
                    spin2Skipped, braSys, ketSys, pot, Ufi_iE)

              Ufi(iE,:) = Ufi_iE

            endif
          enddo ! Loop over all transitions

          call deallocateSysArrays(braSys)

        enddo ! Loop over unique final k-points

        call deallocateSysArrays(ketSys)

      enddo ! Loop over unique initial k-points


      ! Subtract baseline if applicable and write out results
      if(ionode) then 
        do isp = 1, nSpins
          if((isp == 1 .and. .not. spin1Skipped) .or. (isp == 2 .and. .not. spin2Skipped)) then
            if(order == 1 .and. subtractBaseline) &
              call readAndSubtractBaseline(-1, isp, nTransitions, Ufi(:,isp))
                ! Pass in ikGlobal = -1 to trigger scattering format
        
            call writeScatterMatrixElementsOrOverlaps(nTransitions, ibi, ibf, iki, ikf, isp, volume, Ufi(:,isp), overlapOnly)
          endif 
        enddo
      endif

    endif

    return

  end subroutine getAndWriteScatterMatrixElementsOrOverlaps

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

!-----------------------------------------------------------------------------------------------
  subroutine getSpinSkipped(ikGlobal, ispSelect, spin1Skipped, spin2Skipped)
    ! Passing -1 for ikGlobal will result in checking the
    ! scattering file pattern with only the spin index

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user

    ! Output variables:
    logical, intent(out) :: spin1Skipped, spin2Skipped
      !! If spin channels skipped


    spin1Skipped = .false.
    spin2Skipped = .false.
    if(indexInPool == 0) then
      if(ispSelect == 2 .or. overlapFileExists(ikGlobal,1)) spin1Skipped = .true.
      if(ispSelect == 1 .or. overlapFileExists(ikGlobal,2)) spin2Skipped = .true.
    endif

    call MPI_BCAST(spin1Skipped, 1, MPI_LOGICAL, root, intraPoolComm, ierr)
    call MPI_BCAST(spin2Skipped, 1, MPI_LOGICAL, root, intraPoolComm, ierr)


    return

  end subroutine getSpinSkipped

!-----------------------------------------------------------------------------------------------
  subroutine spinAndBandIndependentSetup(ikGlobal, nGVecsLocal, sys)

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Global k-point index
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors

    ! Output variables:

    type(crystal) :: sys
      !! Crystal system to get the matrix element for

    ! Local variables:
    integer :: iDum
      !! Dummy integer


    call distributeItemsInSubgroups(indexInPool, sys%nPWs1kGlobal(ikGlobal), nProcPerPool, nProcPerPool, nProcPerPool, &
            sys%iGkStart_pool, iDum, sys%nGkVecsLocal)
      ! Ignore the iGkEnd_pool variable because it is never used

    call allocateSysArrays(nGVecsLocal, sys)

    call readProjectors(ikGlobal, sys)

    return

  end subroutine spinAndBandIndependentSetup

!-----------------------------------------------------------------------------------------------
  subroutine allocateSysArrays(nGVecsLocal, sys)

    implicit none
    
    ! Input variables:
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors

    ! Output variables:
    type(crystal) :: sys
       !! The crystal system

        
    allocate(sys%beta(sys%nGkVecsLocal,sys%nProj))
    allocate(sys%wfc(sys%nGkVecsLocal))
    allocate(sys%crossProjection(sys%nProj))
    allocate(sys%projection(sys%nProj))
    allocate(sys%pawK(nGVecsLocal))

    return

  end subroutine allocateSysArrays

!-----------------------------------------------------------------------------------------------
  subroutine deallocateSysArrays(sys)

    implicit none
    
    ! Output variables:
    type(crystal) :: sys
       !! The crystal system

        
    deallocate(sys%wfc)
    deallocate(sys%beta)
    deallocate(sys%crossProjection)
    deallocate(sys%projection)
    deallocate(sys%pawK)

    return

  end subroutine deallocateSysArrays
  
!----------------------------------------------------------------------------
  function overlapFileExists(ikGlobal, isp) result(fileExists)
    ! Passing -1 for ikGlobal will result in checking the
    ! scattering file pattern with only the spin index
    
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
      ! Logic for scattering vs capture pattern is in this function

    inquire(file=fName, exist=fileExists)

    if(fileExists) write(*,'("Overlap file ", a, " exists and will not be recalculated.")') trim(fName)
    
  end function overlapFileExists

!----------------------------------------------------------------------------
  subroutine readProjectors(ikGlobal, sys)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current k point

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


    do igkLocal = 1, sys%nGkVecsLocal

      igkGlobal = igkLocal+sys%iGkStart_pool-1

      read(72,rec=igkGlobal+1) (sys%beta(igkLocal,ipr), ipr=1,sys%nProj)

    enddo

    close(72)
    
    return 

  end subroutine readProjectors
  
!----------------------------------------------------------------------------
  subroutine calculateBandPairOverlap(ibShift_braket, ibBra, ibKet, ikBra, ikKet, nSpins, ibBra_optimal, nGVecsLocal, volume, &
          readOptimalPairs, spin1Skipped, spin2Skipped, braSys, ketSys, pot, Ufi)

    implicit none

    ! Input variables:
    integer, intent(in) :: ibShift_braket
      !! Shift of band indexing of bra system relative
      !! to ket system
    integer, intent(in) :: ibBra, ibKet
      !! Band indices for bra and ket systems
    integer, intent(in) :: ikBra, ikKet
      !! Current k-point for bra and ket systems
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)
    integer, intent(in) :: ibBra_optimal(nSpins)
      !! Optimal index from the bra system corresponding 
      !! to the input index from the ket system
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors

    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    logical, intent(in) :: readOptimalPairs
      !! If optimal pairs should be read and states reordered
    logical, intent(in) :: spin1Skipped, spin2Skipped
      !! If spin channels skipped

    type(crystal) :: braSys, ketSys
       !! The crystal systems to get the
       !! matrix element for

    type(potcar) :: pot
      !! Structure containing all pseudopotential-related
      !! information

    ! Output variables:
    complex(kind=dp), intent(out) :: Ufi(nSpins)
      !! All-electron overlap for this k-point and band pair


    ! Local variables:
    integer :: isp
      !! Loop index

    logical :: calcSpinDepBra, calcSpinDepKet
      !! If spin-dependent subroutines should be called


    do isp = 1, nSpins
      if((isp == 1 .and. .not. spin1Skipped) .or. (isp == 2 .and. .not. spin2Skipped)) then

        calcSpinDepKet = isp == 1 .or. ketSys%nSpins == 2 .or. spin1Skipped
        calcSpinDepBra = isp == 1 .or. braSys%nSpins == 2 .or. spin1Skipped .or. readOptimalPairs
          ! If either of the systems only has a single spin channel, some inputs and
          ! calculations do not need to be redone for both spin channels. However, we
          ! still need to make sure to calculate these values for the second spin
          ! channel if the first spin channel was not done (e.g., the file already 
          ! existed or only the second spin channel was selected).
          !
          ! I added the condition to always calculate the spin-dependent pieces if the
          ! user selected readOptimalPairs because the way the bands line up could be 
          ! different for the different spin channels. It is unlikely that the perfect
          ! crystal (ket) would have multiple spins while the defect system (bra) would
          ! be non-spin-polarized, but I am adding it just in case. 


        if(calcSpinDepBra) then
          if(readOptimalPairs) then
            call calcSpinDep(ibBra_optimal(isp), ikBra, isp, nGVecsLocal, braSys, ketSys)
          else
            call calcSpinDep(ibBra, ikBra, isp, nGVecsLocal, braSys, ketSys)
          endif
        endif
      
        if(calcSpinDepKet) call calcSpinDep(ibKet+ibShift_braket, ikKet, isp, nGVecsLocal, ketSys, braSys)


        Ufi(isp) = dot_product(braSys%wfc(:),ketSys%wfc(:))
        if(indexInPool == 0) call pawCorrectionWfc(ketSys, pot)
        if(indexInPool == 1) call pawCorrectionWfc(braSys, pot)

        call MPI_BCAST(ketSys%pawWfc, 1, MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
        call MPI_BCAST(braSys%pawWfc, 1, MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)

        Ufi(isp) = Ufi(isp) + (16.0_dp*pi*pi/volume)*dot_product(conjg(braSys%pawK(:)),ketSys%pawK(:))

        call MPI_ALLREDUCE(MPI_IN_PLACE, Ufi(isp), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, intraPoolComm, ierr)

        Ufi(isp) = Ufi(isp) + braSys%pawWfc + ketSys%pawWfc

      endif
    enddo 

    return

  end subroutine calculateBandPairOverlap

!----------------------------------------------------------------------------
  subroutine calcSpinDep(ib, ikSysCalc, isp, nGVecsLocal, sysCalc, sysCrossProj)

    implicit none

    ! Input variables:
    integer, intent(in) :: ib
      !! Band index
    integer, intent(in) :: ikSysCalc
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors

    ! Output variables:
    type(crystal) :: sysCalc
      !! The crystal system to calculate spin-dependent 
      !! pieces for
    type(crystal) :: sysCrossProj
      !! The crystal system used for the cross projection


    call readWfc(ib, ikSysCalc, min(isp,sysCalc%nSpins), sysCalc)

    call calculateCrossProjection(sysCalc, sysCrossProj)
      ! Get new cross projection with new `sysCalc%wfc`

    call readProjections(ib, ikSysCalc, min(isp,sysCalc%nSpins), sysCalc)

    call pawCorrectionK(nGVecsLocal, pot, Ylm, sysCalc)

    return

  end subroutine calcSpinDep

!----------------------------------------------------------------------------
  subroutine readWfc(ib, ikGlobal, isp, sys)
    !! Read wave function at band ib for given system
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ib
      !! Band index
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel

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
    sendCount(indexInPool+1) = sys%nGkVecsLocal
    call mpiSumIntV(sendCount, intraPoolComm)
      !! * Put the number of G+k vectors on each process
      !!   in a single array per pool

    displacement = 0
    displacement(indexInPool+1) = sys%iGkStart_pool-1
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

    call MPI_SCATTERV(wfcAllPWs(:), sendCount, displacement, MPI_COMPLEX, sys%wfc(1:sys%nGkVecsLocal), sys%nGkVecsLocal, &
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
  subroutine readAndSubtractBaseline(ikGlobal, isp, nTransitions, Ufi)
    ! Passing -1 for ikGlobal will result in checking the
    ! scattering file pattern with only the spin index

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nTransitions
      !! Total number of transitions 

    ! Output variables:
    complex(kind=dp), intent(inout) :: Ufi(nTransitions)
      !! All-electron overlap

    ! Local variables:
    integer :: iE
      !! Loop index
    integer :: iDum
      !! Dummy integer to ignore input

    complex(kind = dp):: baselineOverlap
      !! Input complex overlap 

    character(len=300) :: baselineFName
      !! Name of baseline overlap file


    baselineFName = trim(getMatrixElementFNameWPath(ikGlobal, isp, baselineDir)) 
    open(17, file=trim(baselineFName), status='unknown')

    call ignoreNextNLinesFromFile(17,8)
      ! Assume that baseline was not run with order = 1 because
      ! that is irrelevant for the baseline

    
    do iE = 1, nTransitions 
      ! Read with format unspecified so that output from both
      ! overlapOnly = .true. and .false. can be used
      ! Use ikGlobal < 0 as a switch for scattering format
      if(ikGlobal > 0) then
        read(17,*) iDum, baselineOverlap
      else
        read(17,*) iDum, iDum, iDum, baselineOverlap
      endif

      Ufi(iE) = Ufi(iE) - baselineOverlap
          
    enddo
    
    close(17)
    
    return
    
  end subroutine readAndSubtractBaseline
  
!----------------------------------------------------------------------------
  subroutine writeCaptureMatrixElements(nTransitions, ibi, ibf, ikLocal, isp, dq_j, volume, Ufi, dqOnly)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: ibf
      !! Final-state index
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    real(kind=dp), intent(in) :: dq_j
      !! \(\delta q_j) for displaced wave functions
      !! (only order = 1)
    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    complex(kind=dp), intent(in) :: Ufi(nTransitions)
      !! All-electron overlap

    logical, intent(in) :: dqOnly
      !! If first-order matrix elements should only be
      !! divided by dq

    ! Local variables:
    integer, allocatable :: iDum1D(:)
      !! Initial-state indices
    integer :: iDum0D_1, iDum0D_2
      !! Final-state index
    integer :: iE
      !! Loop index
    integer :: ikGlobal
      !! Current global k-point

    real(kind=dp), allocatable :: dE(:,:)
      !! Energy difference to be combined with
      !! overlap for matrix element
    
    character(len = 300) :: text
      !! Text for header


    ikGlobal = ikLocal+ikStart_pool-1


    if(.not. (order == 1 .and. dqOnly)) &
      call readCaptureEnergyTable(ikGlobal, isp, energyTableDir, iDum1D, iDum0D_1, iDum0D_2, dE)
    

    open(17, file=trim(getMatrixElementFNameWPath(ikGlobal, isp, outputDir)), status='unknown')
    
    write(17, '("# Total number of k-points, k-point index, spin index Format : ''(3i10)''")')
    write(17,'(3i10)') nKPoints, ikGlobal, isp

    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') volume
    
    ! Include in the output file that these are capture matrix elements
    write(17,'("# Capture matrix elements? Alternative is scattering or overlap-only.)")')
    write(17,'(L4)') .true.

    text = "# Total number of transitions, Initial States (bandI, bandF), Final State (band)"
    write(17,'(a, " Format : ''(4i10)''")') trim(text)   
  
    write(17,'(4i10)') nTransitions, ibi(1), ibi(nTransitions), ibf


    if(order == 1) &
      write(17,'("# Phonon mode j, dq_j (Bohr*sqrt(elec. mass)). Format: ''(a78, i7, ES24.15E3)'' ", i7, ES24.15E3)') phononModeJ, dq_j
    

    text = "# Initial Band, Complex <f|i>"

    if(order == 1 .and. subtractBaseline) text = trim(text)//"-baseline"

    text = trim(text)//", |<f|i>|^2,"

    if(order == 0) then
      text = trim(text)//" |dE*<f|i>|^2 (Hartree^2)"
    else if(order == 1 .and. dqOnly) then
      text = trim(text)//" |<f|i>/dq_j|^2 (1/(Bohr*sqrt(elec. mass))^2)"
    else if(order == 1) then
      text = trim(text)//" |dE*<f|i>/dq_j|^2 (Hartree^2/(Bohr*sqrt(elec. mass))^2)"
    endif

    write(17, '(a, " Format : ''(i10,4ES24.15E3)''")') trim(text)


    do iE = 1, nTransitions
        
        if(order == 0) then
          write(17,'(i10,4ES24.15E3)') ibi(iE), Ufi(iE), abs(Ufi(iE))**2, abs(dE(2,iE)*Ufi(iE))**2
        else if(order == 1) then
          if(dqOnly) then
            write(17,'(i10,4ES24.15E3)') ibi(iE), Ufi(iE), abs(Ufi(iE))**2, abs(Ufi(iE)/dq_j)**2
          else
            write(17,'(i10,4ES24.15E3)') ibi(iE), Ufi(iE), abs(Ufi(iE))**2, abs(dE(3,iE)*Ufi(iE)/dq_j)**2
          endif
        endif
            
    enddo

    close(17)
    
    write(*, '("    Ufi(:) of k-point ", i4, " and spin ", i1, " written.")') ikGlobal, isp
    
    return
    
  end subroutine writeCaptureMatrixElements
  
!----------------------------------------------------------------------------
  subroutine writeScatterMatrixElementsOrOverlaps(nTransitions, ibi, ibf, iki, ikf, isp, volume, Ufi, overlapOnly)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: iki(nTransitions), ikf(nTransitions), ibi(nTransitions), ibf(nTransitions)
      !! State indices
    integer, intent(in) :: isp
      !! Current spin channel

    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    complex(kind=dp), intent(in) :: Ufi(nTransitions)
      !! All-electron overlap
      
    logical, intent(in) :: overlapOnly
      !! If only the wave function overlap should be
      !! calculated

    ! Local variables:
    integer, allocatable :: iDum1D_1(:), iDum1D_2(:), iDum1D_3(:), iDum1D_4(:)
      !! Ignore state bounds already read
    integer :: iDum0D
      !! Dummy integer to ignore nTransitions
    integer :: iE
      !! Loop index

    real(kind=dp), allocatable :: dE(:,:)
      !! Energy difference to be combined with
      !! overlap for matrix element
    
    character(len = 300) :: text
      !! Text for header


    if(.not. overlapOnly) &
      call readScatterEnergyTable(isp, .true., energyTableDir, iDum1D_1, iDum1D_2, iDum1D_3, iDum1D_4, iDum0D, dE)
        ! Pass true to get matrix-element and delta-function energies.
    

    open(17, file=trim(getMatrixElementFNameWPath(-1, isp, outputDir)), status='unknown')
      ! Passing -1 for ikGlobal switches to reading the scatter
      ! matrix element file with only isp
    
    write(17, '("# Spin index Format : ''(3i10)''")')
    write(17,'(3i10)') isp

    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') volume
    
    ! Include in the output file that these are capture matrix elements
    write(17,'("# Capture matrix elements? Alternative is scattering or overlap-only.)")')
    write(17,'(L4)') .false.

    text = "# Total number of transitions, Initial States (kI, kF, bandI, bandF), Final States (kI, kF, bandI, bandF)"
    write(17,'(a, " Format : ''(9i10)''")') trim(text)   
  
    write(17,'(9i10)') nTransitions, iki(1), iki(nTransitions), ibi(1), ibi(nTransitions), &
                                     ikf(1), ikf(nTransitions), ibf(1), ibf(nTransitions)
    


    if(overlapOnly) then
      text = "# iki, ibi, ikf, ibf, Complex <f|i>, |<f|i>|^2" 
    else if(order == 0) then
      text = "# iki, ibi, ikf, ibf, Complex <f|i>, |<f|i>|^2, |dE*<f|i>|^2 (Hartree^2)" 
    else if(order == 1) then

      write(17,'("# Phonon mode j, dq_j (Bohr*sqrt(elec. mass)). Format: ''(a78, i7, ES24.15E3)'' ", i7, ES24.15E3)') phononModeJ, dq_j
    
      if(subtractBaseline) then
        text = "# iki, ibi, ikf, ibf, Complex <f|i>-baseline, |<f|i>|^2, |dE*<f|i>/dq_j|^2 (Hartree^2/(Bohr*sqrt(elec. mass))^2)" 
      else
        text = "# iki, ibi, ikf, ibf, Complex <f|i>, |<f|i>|^2, |dE*<f|i>/dq_j|^2 (Hartree^2/(Bohr*sqrt(elec. mass))^2)" 
      endif
    endif
    
    write(17, '(a, " Format : ''(4i10,4ES24.15E3)''")') trim(text)


    do iE = 1, nTransitions
        
        if(overlapOnly) then
          write(17,'(4i10,3ES24.15E3)') iki(iE), ibi(iE), ikf(iE), ibf(iE), Ufi(iE), abs(Ufi(iE))**2
        else if(order == 0) then
          write(17,'(4i10,4ES24.15E3)') iki(iE), ibi(iE), ikf(iE), ibf(iE), Ufi(iE), abs(Ufi(iE))**2, abs(dE(2,iE)*Ufi(iE))**2
        else if(order == 1) then
          write(17,'(4i10,4ES24.15E3)') iki(iE), ibi(iE), ikf(iE), ibf(iE), Ufi(iE), abs(Ufi(iE))**2, abs(dE(3,iE)*Ufi(iE)/dq_j)**2
        endif
            
    enddo

    close(17)
    
    write(*, '("    Ufi(:) of spin ", i1, " written.")') isp
    
    return
    
  end subroutine writeScatterMatrixElementsOrOverlaps
  
!----------------------------------------------------------------------------
  subroutine writeInterKOverlaps(nPairs, ibBra, ibKet, ikLocal, isp, volume, Ufi)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: nPairs
      !! Number of pairs of bands to get overlaps for
    integer, intent(in) :: ibBra(nPairs), ibKet(nPairs)
      !! Band indices for different systems
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    real(kind=dp), intent(in) :: volume
      !! Volume of unit cell

    complex(kind=dp), intent(in) :: Ufi(nPairs)
      !! All-electron overlap

    ! Local variables:
    integer :: ip
      !! Loop index
    integer :: ikGlobal
      !! Current global k-point


    ikGlobal = ikLocal+ikStart_pool-1
    

    open(17, file=trim(getMatrixElementFNameWPath(ikGlobal, isp, outputDir)), status='unknown')
    
    write(17, '("# Total number of k-points, k-point index, spin index Format : ''(3i10)''")')
    write(17,'(3i10)') nKPoints, ikGlobal, isp

    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') volume
    
    ! Include in the output file that these are capture matrix elements
    write(17,'("# Capture matrix elements? Alternative is scattering or overlap-only.)")')
    write(17,'(L4)') .false.

    write(17,'("# Total number of pairs Format : ''(i10)''")')
    write(17,'(i10)') nPairs
    

    write(17, '("# ibBra, ibKet, Complex <f|i>, |<f|i>|^2 Format : ''(2i10,3ES24.15E3)''")')

    do ip = 1, nPairs
        
      write(17,'(2i10,3ES24.15E3)') ibBra(ip), ibKet(ip), Ufi(ip), abs(Ufi(ip))**2
            
    enddo

    close(17)
    
    write(*, '("    Ufi(:) of k-point ", i4, " and spin ", i1, " written.")') ikGlobal, isp
    
    return
    
  end subroutine writeInterKOverlaps
   
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


    if(ikGlobal > 0) then
      fName = trim(path)//"/allElecOverlap."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))
    else
      fName = trim(path)//"/allElecOverlap."//trim(int2str(isp))
    endif

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


    if(ikGlobal > 0) then
      fName = "allElecOverlap."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))
    else
      fName = "allElecOverlap."//trim(int2str(isp))
    endif

  end function getMatrixElementFName

!----------------------------------------------------------------------------
  subroutine readSingleKMatrixElements(ikGlobal, iSpin, nTransitions, ibi, nModes, jReSort, mDim, order, suffixLength, &
            dE, captured, newEnergyTable, oldFormat, rereadDq, reSortMEs, dqInput, matrixElementDir, MjBaseDir, prefix, &
            matrixElement, volumeLine)
    ! For zeroth-order, this will read a single file. For first-order,
    ! this will read all of the modes. Scattering and capture are
    ! handled outside this subroutine by what is passed in ikGlobal. If
    ! ikGlobal < 0, the scattering file-name pattern will be used without
    ! the k-point index.

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Global k-point index for capture and -1 for scattering
    integer, intent(in) :: iSpin
      !! Spin channel to use
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(in) :: jReSort(nModes)
      !! Indices to optionally resort matrix elements
    integer, intent(in) :: mDim
      !! Size of first dimension for matrix element
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)
    integer, intent(in) :: suffixLength
      !! Length of shifted POSCAR file suffix

    real(kind=dp), intent(in) :: dE(3,nTransitions)
      !! All energy differences from energy table

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: newEnergyTable
      !! If this code and TME are being run with a different
      !! energy table
    logical, intent(in) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used
    logical, intent(in) :: rereadDq
      !! If dq should be read from matrix element file
      !! (.false.) or from the dq.txt file (.true.)
    logical, intent(in) :: reSortMEs
      !! If matrix elements should be resorted

    character(len=300), intent(in) :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(in) :: matrixElementDir
      !! Path to matrix element file `allElecOverlap.isp.ik`. 
      !! For first-order term, the path is just within each 
      !! subdirectory.
    character(len=300), intent(in) :: MjBaseDir
      !! Path to the base directory for the first-order
      !! matrix element calculations
    character(len=300), intent(in) :: prefix
      !! Prefix of directories for first-order matrix
      !! elements

    ! Output variables:
    real(kind=dp), intent(out) :: matrixElement(mDim,nTransitions)
      !! Electronic matrix element

    character(len=300), intent(out) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: j, jStore
      !! Loop indices

    real(kind=dp), allocatable :: dENew(:)
      !! New energy to update matrix element; needed
      !! not to create a temporary array when passing
      !! a slice of dE
    real(kind=dp), allocatable :: ME_tmp(:)
      !! Temporary storage of matrix element

    character(len=300) :: fName
      !! File name to read


    allocate(dENew(nTransitions))
    allocate(ME_tmp(nTransitions))


    ! The zeroth-order term reads only one matrix element, with or without the new energy table
    if(order == 0) then
      fName = getMatrixElementFNameWPath(ikGlobal, iSpin, matrixElementDir)

      dENew = dE(2,:)
        ! The second index holds the zeroth-order energy

      call callSingleMESubroutineWithProperArguments(nTransitions, ibi, -1, order, dENew, captured, newEnergyTable, &
              oldFormat, rereadDq, dqInput, fName, ME_tmp, volumeLine)
        ! Pass -1 for jStore as it is not relevant for the zeroth-order

      matrixElement(1,:) = ME_tmp


    ! The first-order term reads a matrix element for each mode. 
    !
    ! The volume line will get overwritten each time through the loop, but that's okay because
    ! the volume doesn't change between the files. 
    else if(order == 1) then

      dENew = dE(3,:)
        ! The third index holds the first-order energy
    
      do j = 1, nModes

        fName = trim(MjBaseDir)//'/'//trim(prefix)//trim(int2strLeadZero(j,suffixLength))//'/'&
                //trim(getMatrixElementFNameWPath(ikGlobal,iSpin,matrixElementDir))

        ! If resorting the matrix element files based on a different PhononPP output
        ! order, make sure to pass the resorted mode index if re-reading dq's
        if(reSortMEs) then
          jStore = jReSort(j)
        else 
          jStore = j
        endif

        call callSingleMESubroutineWithProperArguments(nTransitions, ibi, jStore, order, dENew, captured, newEnergyTable, &
                oldFormat, rereadDq, dqInput, fName, ME_tmp, volumeLine)


        ! Store the matrix element for this mode
        matrixElement(jStore,:) = ME_tmp

      enddo
    endif

    deallocate(dENew)
    deallocate(ME_tmp)

    return

  end subroutine readSingleKMatrixElements

!----------------------------------------------------------------------------
  subroutine callSingleMESubroutineWithProperArguments(nTransitions, ibi, jStore, order, dENew, captured, newEnergyTable, &
            oldFormat, rereadDq, dqInput, fName, ME_tmp, volumeLine)
    ! The format read from the matrix element file and what factors are used
    ! depends on the arguments passed. Handle that logic here to call the
    ! subroutine to read a single matrix element file with the correct
    ! arguments.

    implicit none

    ! Input variables:
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: ibi(nTransitions)
      !! Initial-state indices
    integer, intent(in) :: jStore
      !! Index to read the new delta q from,
      !! if applicable
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dENew(nTransitions)
      !! New energy to update matrix element; needed
      !! not to create a temporary array when passing
      !! a slice of dE

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: newEnergyTable
      !! If this code and TME are being run with a different
      !! energy table
    logical, intent(in) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used
    logical, intent(in) :: rereadDq
      !! If dq should be read from matrix element file
      !! (.false.) or from the dq.txt file (.true.)

    character(len=300), intent(in) :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(in) :: fName
      !! File name to read

    ! Output variables:
    real(kind=dp), intent(out) :: ME_tmp(nTransitions)
      !! Temporary storage of matrix element

    character(len=300), intent(out) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file


    ! For new energy table and new dqs, must pass the band bounds, mode index
    ! to read (jStore) and dqInput. This is only an option for the first-order
    ! term.
    if(newEnergyTable .and. order == 1 .and. rereadDq) then
      call readSingleMatrixElementFile(minval(ibi), maxval(ibi), nTransitions, order, dENew, captured, newEnergyTable, oldFormat, &
            fName, ME_tmp, volumeLine, jStore, dqInput)

    ! For just new energy table, only pass band bounds. For order = 0 ignore
    ! the value in rereadDq.
    else if(newEnergyTable .and. (order == 0 .or. .not. rereadDq)) then
      call readSingleMatrixElementFile(minval(ibi), maxval(ibi), nTransitions, order, dENew, captured, newEnergyTable, oldFormat, &
            fName, ME_tmp, volumeLine)

    ! For just new dq, only pass jStore and dqInput. Again only an option
    ! for the first-order term.
    else if((.not. newEnergyTable) .and. order == 1 .and. rereadDq) then
      call readSingleMatrixElementFile(-1, -1, nTransitions, order, dENew, captured, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine, &
        jStore, dqInput)
          ! dENew will be ignored here

    ! For neither, don't pass anything and just read matrix elements as-is. Again
    ! ignore rereadDq for zeroth-order.
    else if((.not. newEnergyTable) .and. (order == 0 .or. .not. rereadDq)) then
      call readSingleMatrixElementFile(-1, -1, nTransitions, order, dENew, captured, newEnergyTable, oldFormat, fName, ME_tmp, volumeLine)
        ! dENew will be ignored here
    endif

    return

  end subroutine callSingleMESubroutineWithProperArguments

!----------------------------------------------------------------------------
  subroutine readSingleMatrixElementFile(ibL, ibH, nTransitions, order, dE, capture, newEnergy, oldFormat, fName, matrixElement, volumeLine, &
        phononModeJ, dqInput)

    use constants, only: HartreeToJ

    implicit none

    ! Input variables:
    integer, intent(in) :: ibL, ibH
      !! Band bounds for initial state; pass in <= 0
      !! to ignore and read all
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)
    integer, optional :: phononModeJ
      !! Phonon mode index for rereading dq file

    real(kind=dp), intent(in) :: dE(nTransitions)
      !! Optional new energy to use

    logical, intent(in) :: capture
      !! If matrix elements expected for capture or scattering
    logical, intent(in) :: newEnergy
      !! If we should use a new input energy in dE
    logical, intent(in) :: oldFormat
      !! If the old format of the matrix element files
      !! should be used

    character(len=300), optional :: dqInput
      !! Input file for dq.txt if rereading
    character(len=300), intent(in) :: fName
      !! Path to matrix element file `allElecOverlap.isp.ik`

    ! Output variables:
    real(kind=dp), intent(out) :: matrixElement(nTransitions)
      !! Electronic matrix element

    character(len=300), intent(out) :: volumeLine
      !! Volume line from overlap file to be
      !! output exactly in transition rate file

    ! Local variables:
    integer :: iDum
      !! Dummy integer to ignore input
    integer :: iE, ibi, iE_
      !! Loop index
    integer :: nTransitions_
      !! Store the number of transitions from the file

    real(kind=dp) :: dq_j
      !! dq if order = 1
    real(kind=dp) :: rDum
      !! Dummy real
    real(kind=dp) :: normSqOverlap, overlapWithFactors
      !! Store the values from the input file

    logical :: capture_
      !! If matrix elements were written for capture or scattering

    character(len=300) :: line


    ! Test the input band bounds and transitions if given.
    if(ibL > 0) then
      if(ibH < ibL) call exitError('readSingleMatrixElementFile', 'High band bound is lower than low band bound!', 1)
      if(ibH - ibL + 1 /= nTransitions) call exitError('readSingleMatrixElementFile', 'Number of transitions input does not match band bounds!', 1)
      if(.not. capture) call exitError('readSingleMatrixElementFile','New state indices not currently implemented for scattering', 1)
        ! This would involve testing iki, ikf, ibi, and ibf. I don't want to do that right now.
    endif


    open(12,file=trim(fName))

    read(12,*)
    read(12,*)

    read(12,'(a)') volumeLine


    ! There was an old format that I used with the Si files, so I put this
    ! here to be able to read those files. This test and the oldFormat 
    ! specifier could really be removed.
    if(.not. oldFormat) then
      ! Test if the matrix element file comes from capture so that 
      ! we know what format to read the file in.
      read(12,*)
      read(12,'(L4)') capture_
      if(capture_ /= capture) then
        if(capture) then
          call exitError('readSingleMatrixElementFile', 'This matrix element was not calculated for capture!', 1)
        else
          call exitError('readSingleMatrixElementFile', 'This matrix element was calculated for capture!', 1)
        endif
      endif
    endif


    read(12,*)
    read(12,*) nTransitions_


    ! For the first-order matrix elements, there is an additional line
    ! for the mode index j and delta q_j. The user can either read the
    ! matrix element file as-is, where the dq line will be ignored, or
    ! they can use new energies and/or delta q_j. If reading a new energy
    ! but not a new delta q_j, read the dq_j originally output in the 
    ! matrix element file. If reading a new delta q_j (i.e., dqInput
    ! is given), ignore this line and instead read from the dq.txt file.
    if(order == 1) then
      if(newEnergy .and. .not. present(dqInput)) then
        read(12,'(a78, i7, ES24.15E3)') line, iDum, dq_j
      else
        ! Ignore additional line for phonon mode 
        read(12,*)
      endif

      if(present(dqInput)) then
        call readDqFile(phononModeJ, dqInput, dq_j)
      endif
    endif

    ! Ignore header line
    read(12,*)

    if(ibL <= 0 .and. nTransitions /= nTransitions_) &
      call exitError('readSingleMatrixElementFile', 'Number of transitions to read and from file do not match, but no bounds given!', 1)


    ! iE is the index to store if new bounds given. iE_ is the index from 
    ! the matrix element file and is used for storing if new bounds not given.
    iE = 0
    do iE_ = 1, nTransitions_

      ! For capture, we only need to read ibi, the norm square of the overlap (in
      ! case the user wants to use new energy or dq_j) and the overlap with the
      ! original energy and dq_j factors included. 
      if(capture) then
        if(oldFormat) then
          read(12,'(2i7,4ES24.15E3)') iDum, ibi, rDum, rDum, normSqOverlap, overlapWithFactors
        else
          read(12,'(i10,4ES24.15E3)') ibi, rDum, rDum, normSqOverlap, overlapWithFactors
        endif
      else
        read(12,'(4i10,4ES24.15E3)') iDum, iDum, iDum, iDum, rDum, rDum, normSqOverlap, overlapWithFactors ! in Hartree^2
      endif


      ! If using band bounds from a new energy table, check the band
      ! bounds for each state to make sure they line up.
      ! This option is not currently implemented for scattering.
      if(ibL > 0) then
        if(ibi >= ibL .and. ibi <= ibH) then
          iE = iE + 1
          call storeSingleElement(iE, nTransitions, order, dE, dq_j, normSqOverlap, overlapWithFactors, newEnergy, matrixElement)
        endif
      else
        call storeSingleElement(iE_, nTransitions, order, dE, dq_j, normSqOverlap, overlapWithFactors, newEnergy, matrixElement)
      endif
        
    enddo


    close(12)

    return

  end subroutine readSingleMatrixElementFile

!----------------------------------------------------------------------------
  subroutine storeSingleElement(iE, nTransitions, order, dE, dq_j, normSqOverlap, overlapWithFactors, newEnergy, matrixElement)

    implicit none

    ! Input variables:
    integer, intent(in) :: iE
      !! Index to store the element in
    integer, intent(in) :: nTransitions
      !! Total number of transitions 
    integer, intent(in) :: order
      !! Order to calculate (0 or 1)

    real(kind=dp), intent(in) :: dE(nTransitions)
      !! Optional new energy to use
    real(kind=dp), intent(in) :: dq_j
      !! dq if order = 1
    real(kind=dp), intent(in) :: normSqOverlap, overlapWithFactors
      !! Stores the values from the input file

    logical, intent(in) :: newEnergy
      !! If we should use a new input energy in dE

    ! Output variables:
    real(kind=dp), intent(out) :: matrixElement(nTransitions)
      !! Electronic matrix element


    if(newEnergy) then
      matrixElement(iE) = normSqOverlap*dE(iE)**2

      if(order == 1) matrixElement(iE) = matrixElement(iE)/(dq_j**2)
    else
      matrixElement(iE) = overlapWithFactors
    endif

    return

  end subroutine storeSingleElement
  
end module TMEmod
