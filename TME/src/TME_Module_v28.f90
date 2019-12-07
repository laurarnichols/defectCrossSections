module TMEModule
  !
  !! Declare all global variables
  !! and house all subroutines
  !!
  implicit none
  !
  ! Declare integer parameters
  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision
  integer, parameter :: iostd = 16
    !! Unit number for output file
  integer, parameter :: root  = 0
    !! ID of the root process
  !
  ! Declare real parameters
  real(kind = dp), parameter :: evToHartree = 0.03674932538878_dp
    !! Conversion factor from eV to Hartree
  real(kind = dp), parameter :: HartreeToEv = 27.21138386_dp
    !! Conversion factor from Hartree to eV
  real(kind = dp), parameter :: pi = 3.141592653589793_dp
    !! Pi
  real(kind = dp), parameter :: sq4pi = 3.544907701811032_dp
    !! \(\sqrt{4\pi}\)
  !
  ! Declare complex parameter
  complex(kind = dp), parameter ::    ii = cmplx(0.0_dp, 1.0_dp, kind = dp)
    !! Complex \(i\)
  !
  ! Declare character parameter 
  character(len = 6), parameter ::      output = 'output'
    !! Name of the output file;
    !! used in [[TMEModule(module):readInput(subroutine)]]
    !! @todo Change I/O from file to console so that usage matches that of QE @endtodo
  !
  ! 
  ! Declare scalar integers
  integer :: gx
  integer :: gy
  integer :: gz
  integer :: i
  integer :: ibf
  integer :: ibi
  integer :: id
  integer :: ierr
    !! Error code returned from MPI
  integer :: ig
  integer :: ik
  integer :: ind2
  integer :: ios
    !! Status returned from I/O commands
  integer :: iPn
  integer :: iTypes
  integer :: j
  integer :: JMAX
    !! \(2*L_{\text{max}} + 1\)
  integer :: kf
  integer :: ki
  integer :: maxL
    !! Maximum angular momentum of projector from any atom type
  integer :: myid
    !! ID for each MPI process
  integer :: n
  integer :: n1
  integer :: n2
  integer :: n3
  integer :: n4
  integer :: nF
  integer :: nGf
  integer :: nGi
  integer :: nGvsF
  integer :: nGvsI
  integer :: nI
  integer :: np
  integer :: nPP
  integer :: npw
  integer :: npwMf
  integer :: npwMi
  integer :: npwNf
  integer :: npwNi
  integer :: nSquareProcs
  integer :: numOfPWs
  integer :: numOfUsedGvecsPP
  integer :: numprocs
    !! Number of processes in the MPI pool
  !
  ! Declare scalar reals
  real(kind = dp) :: eBin
  real(kind = dp) t0
    !! Start time for program
  real(kind = dp) tf
    !! End time for program
  real(kind = dp) :: threej
  !
  ! Declare scalar complex numbers
  complex(kind = dp) :: paw
  complex(kind = dp) :: paw2
  complex(kind = dp) :: pseudo1
  complex(kind = dp) :: pseudo2
  !
  ! Define scalar logicals
  logical :: calculateVfis
  logical :: coulomb
  logical :: gamma_only
  logical :: master
  logical :: tmes_file_exists
  !
  ! Declare scalar characters
  character(len = 300) :: elementsPath
  character(len = 320) :: mkdir
    !! Command for creating the elements path directory
  character(len = 300) :: textDum
    !! Dummy variable to hold unneeded lines from input file
  character(len = 200) :: VfisOutput
    !! Output file for ??
  !
  !
  ! Declare matrix/vector integers
  integer, allocatable :: counts(:)
  !integer, allocatable :: displmnt(:)
  integer, allocatable :: igvs(:,:,:)
  integer, allocatable :: iqs(:)
  integer, allocatable :: nFs(:,:)
  integer, allocatable :: ngs(:,:)
  integer, allocatable :: nIs(:,:)
  integer, allocatable :: nPWsI(:)
  integer, allocatable :: nPWsF(:)
  integer, allocatable :: pwGvecs(:,:)
  integer, allocatable :: pwGs(:,:)
  !
  ! Declare matrix/vector reals
  real(kind = dp), allocatable :: absVfi2(:,:)
  real(kind = dp), allocatable :: DE(:,:)
  real(kind = dp), allocatable :: eigvF(:)
  real(kind = dp), allocatable :: eigvI(:)
  real(kind = dp), allocatable :: gvecs(:,:)
  !
  ! Declare matrix/vector complex numbers
  complex(kind = dp), allocatable :: paw_id(:,:)
  complex(kind = dp), allocatable :: paw_fi(:,:)
  complex(kind = dp), allocatable :: pawPsiPC(:,:)
  complex(kind = dp), allocatable :: pawSDPhi(:,:)
  complex(kind = dp), allocatable :: paw_SDKKPC(:,:)
  complex(kind = dp), allocatable :: Ufi(:,:,:)
  !
  !
  !
  type :: atom
    !! Define a new type to represent an atom in the structure. 
    !! Each different type of atom in the structure will be another
    !! variable with the type `atom`. 
    !! @todo Consider changing `atom` type to `element` since it holds more than one atom @endtodo
    !
    ! Define scalar integers
    integer :: iRAugMax
      !! Maximum radius of beta projector (outer radius to integrate);
      !! for PAW augmentation charge may extend a bit further; I think this
      !! is the max index for the augmentation sphere, so I'm changing the 
      !! name; last name was `iRc`
    integer :: numOfAtoms
      !! Number of atoms of a specific type in the structure
    integer :: numProjs
      !! Number of projectors
    integer :: lmMax
      !! Number of channels
    integer :: nMax
      !! Number of radial mesh points
    ! 
    ! Define scalar character
    character(len = 2) :: symbol
      !! Element name for the given atom type
    !
    ! Define matrix/vector integer
    integer, allocatable :: projAngMom(:)
      !! Angular momentum of each projector
    !
    ! Define matrix/vector reals
    real(kind = dp), allocatable :: bes_J_qr(:,:)
    real(kind = dp), allocatable :: F(:,:)
    real(kind = dp), allocatable :: F1(:,:,:)
    real(kind = dp), allocatable :: F2(:,:,:)
    real(kind = dp), allocatable :: r(:)
      !! Radial mesh
    real(kind = dp), allocatable :: rab(:)
      !! Derivative of radial mesh
    real(kind = dp), allocatable :: wae(:,:)
      !! All electron wavefunction
    real(kind = dp), allocatable :: wps(:,:)
      !! Psuedowavefunction
    !
  end type atom
  !
  !
  type :: crystal
    integer :: iBandL
      !! Lower band
    integer :: iBandH
      !! Higher band
    integer :: nKpts
      !! Number of k points
    integer :: numOfPWs
      !! Total number of plane waves
    integer :: nIons
      !! Total number of atoms in system
    integer :: numOfTypes
      !! Number of different types of atoms
    integer :: nProjs
      !! Number of projectors
    integer :: numOfGvecs
      !! Number of G vectors
    !integer :: fftxMax, fftxMin, fftyMax, fftyMin, fftzMax, fftzMin
      !! FFT grid was read from `input` file but not used, so removed
    integer :: nBands
      !! Number of bands
    integer :: nSpins
      !! Number of spins
    integer, allocatable :: npws(:)
      !! Number of plane waves per k point
    integer, allocatable :: atomTypeIndex(:)
      !! Index of the given atom type
    !integer, allocatable :: groundState(:)
      ! Was read from `input` file but not used, so removed
    !
    real(kind = dp) :: omega
      !! Cell volume
    !real(kind = dp) :: at(3,3)
      ! Was read from `input` file but not used, so removed
    real(kind = dp) :: bg(3,3)
    real(kind = dp), allocatable :: wk(:)
    real(kind = dp), allocatable :: xk(:, :)
    real(kind = dp), allocatable :: posIon(:,:)
    !
    complex(kind = dp), allocatable :: wfc(:,:)
    complex(kind = dp), allocatable :: beta(:,:)
    complex(kind = dp), allocatable :: cProj(:,:,:)
    complex(kind = dp), allocatable :: cCrossProj(:,:,:)
    complex(kind = dp), allocatable :: paw_Wfc(:,:)
    complex(kind = dp), allocatable :: pawK(:,:,:)
    !
    character(len = 2) crystalType
      !! 'PC' for pristine crystal and 'SD' for solid defect
    character(len = 200) :: exportDir
      !! Export directory from [[pw_export_for_tme(program)]]
    !
    TYPE(atom), allocatable :: atoms(:)
    !
!    integer :: Jmax, maxL, iTypes, nn, nm
!    integer :: i, j, n1, n2, n3, n4, n, id
!    !
!    !
!    real(kind = dp), allocatable :: eigvI(:), eigvF(:)
!    real(kind = dp), allocatable :: DE(:,:), absVfi2(:,:)
!    !
!    complex(kind = dp), allocatable :: Ufi(:,:,:)
!    !
!    integer, allocatable :: igvs(:,:,:), pwGvecs(:,:), iqs(:)
!    integer, allocatable :: pwGs(:,:), nIs(:,:), nFs(:,:), ngs(:,:)
!
  end type crystal
  !
  TYPE(crystal) :: perfectCrystal
    !! Structure that holds all of the information on the perfect crystal
  !
  TYPE(crystal) :: solidDefect
    !! Structure that holds all of the information on the solid defect
  !
  type :: vec
    !
    integer :: ind
    integer, allocatable :: igN(:)
    integer, allocatable :: igM(:)
  end type vec
  !
  ! Define vectors of vecs
  TYPE(vec), allocatable :: vecs(:)
  TYPE(vec), allocatable :: newVecs(:)
  !
  !
!=====================================================================================================
contains
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine initializeCalculation(solidDefect, perfectCrystal, elementsPath, VFisOutput, ki, kf, eBin, &
                                   calculateVFis, t0)
    !! Initialize the calculation by starting timer,
    !! setting start values for variables to be read from
    !! `.in` file, removing any existing output in the output directory,
    !! and opening a clean output file
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(out) :: ki, kf
    !
    real(kind = dp), intent(out) :: eBin, t0
    !
    character(len = 200), intent(out) :: VfisOutput
    character(len = 300), intent(out) :: elementsPath
    !
    logical, intent(out) :: calculateVfis
    logical :: fileExists
      !! Whether or not the output file already exists
    TYPE(crystal), intent(inout) :: solidDefect, perfectCrystal
    !
    solidDefect%exportDir = ''
    perfectCrystal%exportDir = ''
    elementsPath = ''
    VfisOutput = ''
    !
    ki = -1
    kf = -1
    !
    eBin = -1.0_dp
    !
    solidDefect%iBandL  = -1
    solidDefect%iBandH = -1
    perfectCrystal%iBandL  = -1
    perfectCrystal%iBandH = -1
    !
    calculateVfis = .false.
    !
    perfectCrystal%crystalType = 'PC'
    solidDefect%crystalType = 'SD'
    !
    call cpu_time(t0)
        !! * Start a timer
    !
    inquire(file = output, exist = fileExists)
        !! * Check if file output exists,
    if ( fileExists ) then
        !! and delete it if it does
      open (unit = 11, file = output, status = "old")
      close(unit = 11, status = "delete")
    endif
    !
    open (iostd, file = output, status='new')
        !! * Open new output file
    !
    return
    !
  end subroutine initializeCalculation
  !
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine readInput(perfectCrystal, solidDefect, elementsPath, ki, kf, calculateVfis, VfisOutput)
    !! Delete any previous output, initialize input variables,
    !! start a timer, and read in the input files
    !!
    implicit none
    !
    integer, intent(inout) :: ki, kf
    integer :: iBandL_SD
    integer :: iBandH_SD
    integer :: iBandL_PC
    integer :: iBandH_PC
    !
    character(len = 300), intent(inout) :: elementsPath
    character(len = 200), intent(inout) :: VfisOutput
    character(len = 200) :: exportDirSD
    character(len = 200) :: exportDirPC
    !
    logical, intent(inout) :: calculateVfis
    !
    TYPE(crystal), intent(inout) :: perfectCrystal
      !! Holds all of the information on the perfect crystal
    TYPE(crystal), intent(inout) :: solidDefect
      !! Holds all of the information on the defective crystal
    !
    NAMELIST /TME_Input/ exportDirSD, exportDirPC, elementsPath, &
                       iBandL_SD, iBandH_SD, iBandL_PC, iBandH_PC, &
                       ki, kf, calculateVfis, VfisOutput, eBin
                       !! Used to group the variables read in from the .in file
    !
    !
    READ (5, TME_Input, iostat = ios)
        !! * Read input from command line (or input file if use `< TME_Input.md`)
    solidDefect%exportDir = exportDirSD
    perfectCrystal%exportDir = exportDirPC
    !
    solidDefect%iBandL = iBandL_SD
    solidDefect%iBandH = iBandH_SD
    perfectCrystal%iBandL = iBandL_PC
    perfectCrystal%iBandH = iBandH_PC
    !
    call checkInitialization()
        !! * Check that all required variables were input and have values that make sense
    !
    call readQEExport(perfectCrystal)
        !! * Read perfect crystal inputs
    call readQEExport(solidDefect)
        !! * Read solid defect inputs
    !
    numOfPWs = max( perfectCrystal%numOfPWs, solidDefect%numOfPWs )
        !! * Calculate the number of plane waves as the maximum of the number of PC and SD plane waves
    !
    return
    !
  end subroutine readInput
  !
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine checkInitialization()
    !! Check to see if variables from .in file still
    !! have the values set in [[TMEModule(module):initializeCalculation(subroutine)]]
    !! or if they have values that aren't allowed
    !!
    !! <h2>Walkthrough</h2>
    !!
    !! @todo Change `checkInitialization()` to have arguments to make clear that these variables are getting changed @endtodo
    !!
    implicit none
    !
    logical :: fileExists
      !! Whether or not the exported directory from [[pw_export_for_TME(program)]]
      !! exists
    logical:: abortExecution
    !
    abortExecution = .false.
      !! * Set the default value of abort execution so that the program
      !! will only abort if there is an issue with the inputs
    !
    write(iostd, '(" Inputs : ")')
      !! * Write out a header to the output file
    !
    if ( wasRead(LEN(trim(solidDefect%exportDir))-1, 'exportDirSD', 'exportDirSD = ''./Export/''', abortExecution) ) then
      !! * If the SD export directory variable was read
      !!    * Check if the SD export directory exists
      !!    * If the SD export directory doesn't exist
      !!       * Output an error message and set `abortExecution` to true
      !!    * Output the given SD export directory
      !
      inquire(file= trim(solidDefect%exportDir), exist = fileExists)
      !
      if ( fileExists .eqv. .false. ) then
        !
        write(iostd, '(" exportDirSD :", a, " does not exist !")') trim(solidDefect%exportDir)
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
        !
      endif
      !
      write(iostd, '("exportDirSD = ''", a, "''")') trim(solidDefect%exportDir)
      !
    endif
    !
    !
    if ( wasRead(LEN(trim(perfectCrystal%exportDir))-1, 'exportDirPC', 'exportDirPC = ''./Export/''', abortExecution) ) then
      !! * If the PC export directory variable was read
      !!    * Check if the PC export directory exists
      !!    * If the PC export directory doesn't exist
      !!       * Output an error message and set `abortExecution` to true
      !!    * Output the given PC export directory
      !
      inquire(file= trim(perfectCrystal%exportDir), exist = fileExists)
      !
      if ( fileExists .eqv. .false. ) then
        !
        write(iostd, '(" exportDirPC :", a, " does not exist !")') trim(perfectCrystal%exportDir)
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
        !
      endif
      !
      write(iostd, '("exportDirPC = ''", a, "''")') trim(perfectCrystal%exportDir)
      !
    endif
    !
    if( .not. wasRead(LEN(elementsPath)-1, 'elementsPath', 'elementsPath = ''./''') ) then
      !! * If the elements path was not read, set the default value to `./`
      !
      write(iostd, '(" The current directory will be used as elementsPath.")')
      elementsPath = './'
      !
    endif
    !
    inquire(file= trim(elementsPath), exist = fileExists)
      !! * Check if the elements path folder exists already
    !
    if ( .not. fileExists ) then
      !! * If the elements path folder doesn't already exist
      !!    * Write the `mkdir` command to a string
      !!    * Execute the command to create the directory
      !
      write(mkDir, '("mkdir -p ", a)') trim(elementsPath) 
      !
      call system(mkDir)
      !
    endif
    !
    write(iostd, '("elementsPath = ''", a, "''")') trim(elementsPath)
      !! * Output the elements path
    !
    if( wasRead(perfectCrystal%iBandL, 'iBandL_PC', 'iBandL_PC = 10', abortExecution) ) then
      !! * If `iBandL_PC` was read, output its value
      !
      write(iostd, '("iBandL_PC = ", i4)') perfectCrystal%iBandL
      !
    endif
    !
    if( wasRead(perfectCrystal%iBandH, 'iBandH_PC', 'iBandH_PC = 20', abortExecution) ) then
      !! * If `iBandH_PC` was read, output its value
      !
      write(iostd, '("iBandH_PC = ", i4)') perfectCrystal%iBandH
      !
    endif
    !
    if( wasRead(solidDefect%iBandL, 'iBandL_SD', 'iBandL_SD = 9', abortExecution) ) then
      !! * If `iBandL_SD` was read, output its value
      !
      write(iostd, '("iBandL_SD = ", i4)') solidDefect%iBandL
      !
    endif
    !
    if( wasRead(solidDefect%iBandH, 'iBandH_SD', 'iBandH_SD = 9', abortExecution) ) then
      !! * If `iBandH_SD` was read, output its value
      !
      write(iostd, '("iBandH_SD = ", i4)') solidDefect%iBandH
      !
    endif
    !
    !> * If `calculateVfis` is true and `iBandL_SD` and `iBandH_SD` are not equal
    !>    * Output an error message and set `abortExecution` to true
    if ( ( calculateVfis ) .and. ( solidDefect%iBandL /= solidDefect%iBandH ) ) then
      !
      write(iostd, *)
      write(iostd, '(" Vfis can be calculated only if the final state is one and only one!")')
      write(iostd, '(" ''iBandL_SD'' = ", i10)') solidDefect%iBandL
      write(iostd, '(" ''iBandH_SD'' = ", i10)') solidDefect%iBandH
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    endif
    !
    write(iostd, '("calculateVfis = ", l5 )') calculateVfis
      !! * Output the value of `calculateVfis`
    !
    !> * If the `VfisOutput` file name is blank
    !>    * Output a warning message and set the default value to `VfisVsE`
    if ( trim(VfisOutput) == '' ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""VfisOutput"" is not defined!")')
      write(iostd, '(" usage : VfisOutput = ''VfisVsE''")')
      write(iostd, '(" The default value ''VfisVsE'' will be used.")')
      VfisOutput = 'VfisVsE'
      !
    endif
    !
    write(iostd, '("VfisOutput = ''", a, "''")') trim(VfisOutput)
      !! * Output the value of `VfisOutput`
    !> @todo Remove everything with `ki` and `kf` because never used @endtodo
    !
    !if( .not. wasRead(ki, 'ki', 'ki = 1') ) then
    !  !! * If `ki` wasn't read, set the default value to 1
    !  !
    !  write(iostd, '(" ki = 1 will be used.")')
    !  ki = 1
    !  !
    !endif
    !
    !if( .not. wasRead(kf, 'kf', 'kf = 1') ) then
    !  !! * If `kf` wasn't read, set the default value to the total 
    !  !!   number of k points (actually done in [[TMEModeul(module):readQEInput(subroutine)]]
    !  !!   where the total number of k points is read
    !  !
    !  write(iostd, '(" kf = total number of k-points will be used.")')
    !  !
    !endif
    !
    !if ( ki /= kf ) then
    !  write(iostd, *)
    !  write(iostd, '(" Initial k-point index ''ki'', should be equal to the Final k-point index ''kf'' !")')
    !  write(iostd, '(" Calculation of transition matrix elements with momentum transfer is not implemented!")')
    !  write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
    !  abortExecution = .true.
    !endif
    !
    if ( .not. wasRead(INT(eBin), 'eBin', 'eBin = 0.01') ) then
      !! * If the value of `eBin` was not read
      !!    * Output a warning message and set the default value to 0.01 eV
      !
      write(iostd,'(" A default value of 0.01 eV will be used !")')
      eBin = 0.01_dp ! eV
      !
    endif
    !
    write(iostd, '("eBin = ", f8.4, " (eV)")') eBin
      !! * Output the value of eBin
    !
    eBin = eBin*evToHartree
      !! * Convert `eBin` from eV to Hartree
    !
    if ( abortExecution ) then
      !! * If `abortExecution` was ever set to true
      !!    * Output an error message and stop the program
      write(iostd, '(" Program stops!")')
      stop
    endif
    !
    flush(iostd)
      !! * Make the output file available for other processes
    !
    return
    !
  end subroutine checkInitialization
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine readQEExport(system)
    !! Read input files in the Export directory created by
    !! [[pw_export_for_tme(program)]]
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    !integer, intent(in) :: id
    !
    TYPE(crystal), intent(inout) :: system
      !! Holds the structure for the system you are working on
      !! (either `perfectCrystal` or `solidDefect`
    !
    integer :: i, ik, iType, ni
      !! Loop index
    integer:: iRAugMax
      !! Local value of `iRAugMax` for each atom so don't have to keep accessing in loop
    integer :: l
      !! Angular momentum of each projector read from input file
    integer :: ind
      !! Index of each projector read from input file
    integer :: iDum
      !! Dummy variable to hold trash from input file
    !
    real(kind = dp) :: t1
      !! Local start time
    real(kind = dp) :: t2
      !! Local end time
    !
    character(len = 300) :: textDum
      !! Dummy variable to hold trash from input file
    character(len = 300) :: input
      !! The input file path
    !
    logical :: fileExists
      !! Whether or not the `input` file exists in the given 
      !! Export directory
    !
    call cpu_time(t1)
      !! * Start a local timer
    !
    !> * Output header to output file based on the input crystal type
    !> @note
    !> The program will end if a crystal type other than `PC` or `SD` is used.
    !> @endnote
    write(iostd, *)
    if ( system%crystalType == 'PC' ) then
      !
      write(iostd, '(" Reading perfect crystal inputs.")')
      !
    else if ( system%crystalType == 'SD' ) then
      !
      write(iostd, '(" Reading solid defect inputs.")')
      !
    else
      !
      write(iostd, '("Unknown crystal type", a, ".")') system%crystalType
      write(iostd, '("Please only use PC for pristine crystal or SD for solid defect.")')
      write(iostd, '(" Program stops!")')
      flush(iostd)
      stop
      !
    endif
    !
    write(iostd, *)
    !
    input = trim(trim(system%exportDir)//'/input')
      !! * Set the path for the input file from the PC export directory
    !
    inquire(file =trim(input), exist = fileExists)
      !! * Check if the input file from the PC export directory exists
    !
    !> * If the input file doesn't exist
    !>    * Output an error message and end the program
    if ( fileExists .eqv. .false. ) then
      !
      write(iostd, '(" File : ", a, " , does not exist!")') trim(input)
      write(iostd, '(" Please make sure that folder : ", a, " has been created successfully !")') trim(system%exportDir)
      write(iostd, '(" Program stops!")')
      flush(iostd)
      stop
      !
    endif
    !
    !...............................................................................................
    !> * Open and read the [input](../page/inputOutput/exportedInput.html) file
    !
    open(50, file=trim(input), status = 'old')
    !
    read(50, '(a)') textDum
    !
    read(50, '(ES24.15E3)' ) system%omega
    !
    read(50, '(a)') textDum
    read(50, '(i10)') system%nKpts
    !if ( kf < 0 ) kf = system%nKpts
    !
    read(50, '(a)') textDum
    ! 
    allocate ( system%npws(system%nKpts), system%wk(system%nKpts), system%xk(3,system%nKpts) )
    ! 
    !allocate( system%groundState(system%nKpts) ) 
      ! Don't allocate space for groundState because it is never used
    !
    do ik = 1, system%nKpts
      !
      !read(50, '(3i10,4ES24.15E3)') iDum, system%groundState(ik), system%npws(ik), system%wk(ik), system%xk(1:3,ik)
        ! Don't read in groundState because it is never used
      read(50, '(3i10,4ES24.15E3)') iDum, iDum, system%npws(ik), system%wk(ik), system%xk(1:3,ik)
      !
    enddo
    !
    read(50, '(a)') textDum
    !
    read(50, * ) system%numOfGvecs
    !
    read(50, '(a)') textDum
    read(50, '(i10)') system%numOfPWs
    !
    read(50, '(a)') textDum     
    !
    !read(50, '(6i10)') fftxMin, fftxMax, fftyMin, fftyMax, fftzMin, fftzMax
      ! Don't read in FFT grid because it is never used
    read(50,  * )
    !
    read(50, '(a)') textDum
    !
    !read(50, '(a5, 3ES24.15E3)') textDum, at(1:3,1)
    !read(50, '(a5, 3ES24.15E3)') textDum, at(1:3,2)
    !read(50, '(a5, 3ES24.15E3)') textDum, at(1:3,3)
      ! Don't read in `at` because it is never used
    read(50,  * )
    read(50,  * )
    read(50,  * )
    !
    read(50, '(a)') textDum
    !
    read(50, '(a5, 3ES24.15E3)') textDum, system%bg(1:3,1)
    read(50, '(a5, 3ES24.15E3)') textDum, system%bg(1:3,2)
    read(50, '(a5, 3ES24.15E3)') textDum, system%bg(1:3,3)
    !
    read(50, '(a)') textDum
    read(50, '(i10)') system%nIons
    !
    read(50, '(a)') textDum
    read(50, '(i10)') system%numOfTypes
    !
    allocate( system%posIon(3,system%nIons), system%atomTypeIndex(system%nIons) )
    !
    read(50, '(a)') textDum
    !
    do ni = 1, system%nIons
      !
      read(50,'(i10, 3ES24.15E3)') system%atomTypeIndex(ni), system%posIon(1:3,ni)
      !
    enddo
    !
    read(50, '(a)') textDum
    !
    read(50, '(i10)') system%nBands
    !
    read(50, '(a)') textDum
    !
    read(50, '(i10)') system%nSpins
    !
    allocate ( system%atoms(system%numOfTypes) )
    !
    system%nProjs = 0
    !
    do iType = 1, system%numOfTypes
      !
      read(50, '(a)') textDum
      read(50, *) system%atoms(iType)%symbol
      !
      read(50, '(a)') textDum
      read(50, '(i10)') system%atoms(iType)%numOfAtoms
      !
      read(50, '(a)') textDum
      read(50, '(i10)') system%atoms(iType)%numProjs              ! number of projectors
      !
      allocate ( system%atoms(iType)%projAngMom( system%atoms(iType)%numProjs ) )
      !
      read(50, '(a)') textDum
      do i = 1, system%atoms(iType)%numProjs 
        !
        read(50, '(2i10)') l, ind
        system%atoms(iType)%projAngMom(ind) = l
        !
      enddo
      !
      read(50, '(a)') textDum
      read(50, '(i10)') system%atoms(iType)%lmMax
      !
      read(50, '(a)') textDum
      read(50, '(2i10)') system%atoms(iType)%nMax, system%atoms(iType)%iRAugMax
      !
      allocate ( system%atoms(iType)%r(system%atoms(iType)%nMax), system%atoms(iType)%rab(system%atoms(iType)%nMax) )
      !
      read(50, '(a)') textDum
      do i = 1, system%atoms(iType)%nMax
        !
        read(50, '(2ES24.15E3)') system%atoms(iType)%r(i), system%atoms(iType)%rab(i)
        !
      enddo
      ! 
      allocate ( system%atoms(iType)%wae(system%atoms(iType)%nMax, system%atoms(iType)%numProjs) )
      allocate ( system%atoms(iType)%wps(system%atoms(iType)%nMax, system%atoms(iType)%numProjs) )
      !
      read(50, '(a)') textDum
      do j = 1, system%atoms(iType)%numProjs
        do i = 1, system%atoms(iType)%nMax
          !
          read(50, '(2ES24.15E3)') system%atoms(iType)%wae(i, j), system%atoms(iType)%wps(i, j) 
          ! write(iostd, '(2i5, ES24.15E3)') j, i, abs(system%atoms(iType)%wae(i, j)-system%atoms(iType)%wps(i, j))
          !
        enddo
      enddo
      !  
      allocate ( system%atoms(iType)%F( system%atoms(iType)%iRAugMax, system%atoms(iType)%numProjs ) ) !, system%atoms(iType)%numProjs) )
      allocate ( system%atoms(iType)%F1(system%atoms(iType)%iRAugMax, system%atoms(iType)%numProjs, system%atoms(iType)%numProjs ) )
      allocate ( system%atoms(iType)%F2(system%atoms(iType)%iRAugMax, system%atoms(iType)%numProjs, system%atoms(iType)%numProjs ) )
      !
      system%atoms(iType)%F = 0.0_dp
      system%atoms(iType)%F1 = 0.0_dp
      system%atoms(iType)%F2 = 0.0_dp
      !
      !> * Calculate `F`, `F1`, and `F2` using the all-electron and psuedowvefunctions
      !> @todo Look more into how AE and PS wavefunctions are combined to further understand this @endtodo
      !> @todo Move this behavior to another subroutine for clarity @endtodo
      do j = 1, system%atoms(iType)%numProjs
        !
        iRAugMax = system%atoms(iType)%iRAugMax
        !
        system%atoms(iType)%F(1:iRAugMax,j)=(system%atoms(iType)%wae(1:iRAugMax,j)-system%atoms(iType)%wps(1:iRAugMax,j))* &
              system%atoms(iType)%r(1:iRAugMax)*system%atoms(iType)%rab(1:iRAugMax)
        !
        do i = 1, system%atoms(iType)%numProjs
          !> @todo Figure out if differences in PC and SD `F1` calculations are intentional @endtodo
          !> @todo Figure out if should be `(wae_i wae_j - wps_i wps_j)r_{ab}` @endtodo
          !> @todo Figure out if first term in each should be conjugated for inner product form @endtodo
          !> @todo Figure out if `rab` plays role of \(dr\) within augmentation sphere @endtodo
          if ( system%crystalType == 'PC' ) then
            !
            system%atoms(iType)%F1(1:iRAugMax,i,j) = (system%atoms(iType)%wps(1:iRAugMax,i)*system%atoms(iType)%wae(1:iRAugMax,j)-&
                                                  system%atoms(iType)%wps(1:iRAugMax,i)*system%atoms(iType)%wps(1:iRAugMax,j))* &
                                                  system%atoms(iType)%rab(1:iRAugMax)
            !
          else if ( system%crystalType == 'SD' ) then
            !
            system%atoms(iType)%F1(1:iRAugMax,i,j) = (system%atoms(iType)%wae(1:iRAugMax,i)*system%atoms(iType)%wps(1:iRAugMax,j)-&
                                                  system%atoms(iType)%wps(1:iRAugMax,i)*system%atoms(iType)%wps(1:iRAugMax,j))* & 
                                                  system%atoms(iType)%rab(1:iRAugMax)
            !
          endif
          !
          system%atoms(iType)%F2(1:iRAugMax,i,j) = ( system%atoms(iType)%wae(1:iRAugMax,i)*system%atoms(iType)%wae(1:iRAugMax,j) - &
                                           system%atoms(iType)%wae(1:iRAugMax,i)*system%atoms(iType)%wps(1:iRAugMax,j) - &
                                           system%atoms(iType)%wps(1:iRAugMax,i)*system%atoms(iType)%wae(1:iRAugMax,j) + &
                                           system%atoms(iType)%wps(1:iRAugMax,i)*system%atoms(iType)%wps(1:iRAugMax,j))* & 
                                           system%atoms(iType)%rab(1:iRAugMax)

        enddo
      enddo
      !
      system%nProjs = system%nProjs + system%atoms(iType)%numOfAtoms*system%atoms(iType)%lmMax
      !
      deallocate ( system%atoms(iType)%wae, system%atoms(iType)%wps )
      !deallocate ( system%groundState )
        ! Don't use because groundState is never used
      !
    enddo
    !
    !...............................................................................................
    !
    close(50)
      !! * Close the input file
    !
    !> * Go through the `projAngMom` values for each projector for each atom
    !> and find the max to store in `JMAX`
    !> @todo Figure out if intentional to only use `JMAX` from SD input @endtodo
    JMAX = 0
    do iType = 1, system%numOfTypes
      !
      do i = 1, system%atoms(iType)%numProjs
        !
        if ( system%atoms(iType)%projAngMom(i) > JMAX ) JMAX = system%atoms(iType)%projAngMom(i)
        !
      enddo
      !
    enddo
    !
    maxL = JMAX
    JMAX = 2*JMAX + 1
    !
    do iType = 1, system%numOfTypes
      !
      allocate ( system%atoms(iType)%bes_J_qr( 0:JMAX, system%atoms(iType)%iRAugMax ) )
      system%atoms(iType)%bes_J_qr(:,:) = 0.0_dp
      !
    enddo
    !
    !> * End the local timer and write out the total time to read the inputs
    !> to the output file
    call cpu_time(t2)
    write(iostd, '(" Reading input files done in:                ", f10.2, " secs.")') t2-t1
    write(iostd, *)
    flush(iostd)
    !
    return
    !
  end subroutine readQEExport
  !
  !
  subroutine readPWsSet()
    !! Read the g vectors in Miller indices from `mgrid` file and convert
    !! using reciprocal lattice vectors
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer :: ig, iDum, iGx, iGy, iGz
    !
    open(72, file=trim(solidDefect%exportDir)//"/mgrid")
      !! * Open the `mgrid` file from Export directory from [[pw_export_for_tme(program)]]
    !
    !> * Ignore the first two lines as they are comments
    read(72, * )
    read(72, * )
    !
    allocate ( gvecs(3, solidDefect%numOfGvecs ) )
      !! * Allocate space for the g vectors
    !
    gvecs(:,:) = 0.0_dp
      !! * Initialize all of the g vectors to zero
    !
    do ig = 1, solidDefect%numOfGvecs
      !! * For each g vector
      !!    * Read in the g vector in terms of Miller indices
      !!    * Calculate the g vector using the reciprocal lattice vectors from input file
      read(72, '(4i10)') iDum, iGx, iGy, iGz
      gvecs(1,ig) = dble(iGx)*solidDefect%bg(1,1) + dble(iGy)*solidDefect%bg(1,2) + dble(iGz)*solidDefect%bg(1,3)
      gvecs(2,ig) = dble(iGx)*solidDefect%bg(2,1) + dble(iGy)*solidDefect%bg(2,2) + dble(iGz)*solidDefect%bg(2,3)
      gvecs(3,ig) = dble(iGx)*solidDefect%bg(3,1) + dble(iGy)*solidDefect%bg(3,2) + dble(iGz)*solidDefect%bg(3,3)
    enddo
    !
    close(72)
      !! Close the `mgrid` file
    !
    return
    !
  end subroutine readPWsSet
  !
  !
  subroutine distributePWsToProcs(nOfPWs, nOfBlocks)
    !! Determine how many g vectors each process should get
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in)  :: nOfPWs
      !! Number of g vectors
    integer, intent(in)  :: nOfBlocks
      !! Number of processes
    integer :: iStep
      !! Number of g vectors per number of processes
    integer :: iModu
      !! Number of remaining g vectors after giving
      !! each process the same number of g vectors
    !
    iStep = int(nOfPWs/nOfBlocks)
      !! * Determine the base number of g vectors to give 
      !!   to each process
    iModu = mod(nOfPWs,nOfBlocks)
      !! * Determine the number of g vectors left over after that
    !
    do i = 0, nOfBlocks - 1
      !! * For each process, give the base amount and an extra
      !!   if there were any still left over
      counts(i) = iStep
      !
      if ( iModu > 0 ) then
        !
        counts(i) = counts(i) + 1
        !
        iModu = iModu - 1
        !
      endif
      !
    enddo
    !
    !displmnt(0) = 0
    !do i = 1, nOfBlocks-1
    !  displmnt(i) = displmnt(i-1) + counts(i)
    !enddo
    !
    return
    !
  end subroutine distributePWsToProcs
  !
  !
  subroutine checkIfCalculated(ik, tmes_file_exists)
    !! Determine if the output file for a given k point already exists
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    logical, intent(out) :: tmes_file_exists
      !! Whether or not the output file exists
    !
    character(len = 300) :: Uelements
      !! Output file name
    character(len = 300) :: ikstr
      !! String version of integer input `ik`
    !
    call int2str(ik, ikstr)
    write(Uelements, '("/TMEs_kptI_",a,"_kptF_",a)') trim(ikstr), trim(ikstr)
      !! * Determine what the file name should be based on the k point index
    !
    inquire(file = trim(elementsPath)//trim(Uelements), exist = tmes_file_exists)
      !! * Check if that file already exists
    !
    return
    !
  end subroutine checkIfCalculated
  !
  !
  subroutine calculatePWsOverlap(ik)
    !! Read the wavefunctions and calculate the overlap
    !! \(\langle\Phi_f|\Psi_i\rangle\)
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    integer :: ibi, ibf
      !! Loop index
    !
    call readWfc(ik, perfectCrystal)
      !! * Read the perfect crystal wavefunction ([[TMEModule(module):readWfc(subroutine)]])
    !
    call readWfc(ik, solidDefect)
      !! * Read the solid defect wavefunction ([[TMEModule(module):readWfc(subroutine)]])
    !
    Ufi(:,:,ik) = cmplx(0.0_dp, 0.0_dp, kind = dp)
      !! * Initialize `Ufi` for the given k point to complex double zero
    !
    do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH 
      !
      do ibf = solidDefect%iBandL, solidDefect%iBandH
        !! * For each initial band, calculate \(\sum \phi_f^*\psi_i\) (overlap??) with each final band
        !!
        Ufi(ibf, ibi, ik) = sum(conjg(solidDefect%wfc(:,ibf))*perfectCrystal%wfc(:,ibi))
          !! * Calculate \(\langle\Phi_f|\Psi_i\rangle\)
          !!
        !if ( ibi == ibf ) write(iostd,'(2i4,3ES24.15E3)') ibf, ibi, Ufi(ibf, ibi, ik), abs(Ufi(ibf, ibi, ik))**2
        flush(iostd)
        !
      enddo
      !
    enddo
    !
    return
    !
  end subroutine calculatePWsOverlap
  !
  !
  subroutine readWfc(ik, system)
    !! Open the `grid.ki` file from [[pw_export_for_tme(program)]]
    !! to get the indices for the wavefunction to be stored in, then
    !! open the `wfc.ki` file and read in the wavefunction for the 
    !! proper bands and store in the proper indices in the system's `wfc`
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    integer :: ib, ig
      !! Loop index
    integer :: iDumV(3)
      !! Dummy vector to ignore g vectors from `grid.ki`
    integer, allocatable :: pwGind(:)
      !! Indices for the wavefunction of a given k point
    !
    complex(kind = dp) :: wfc
      !! Wavefunction
    !
    character(len = 300) :: ikstr
      !! String version of the k point index
    !
    TYPE(crystal), intent(inout) :: system
      !! Holds the structure for the system you are working on
      !! (either `perfectCrystal` or `solidDefect`)
    !
    call int2str(ik, ikstr)
      !! * Convert the k point index to a string
    !
    open(72, file=trim(system%exportDir)//"/grid."//trim(ikstr))
      !! * Open the `grid.ki` file from [[pw_export_for_tme(program)]]
    !
    !> * Ignore the first two lines as they are comments
    read(72, * )
    read(72, * )
    !
    allocate ( pwGind(system%npws(ik)) )
      !! * Allocate space for `pwGind`
    !
    do ig = 1, system%npws(ik)
      !! * For each plane wave for a given k point, 
      !!   read in the indices for the plane waves that 
      !!   are held in `wfc.ki`
      !
      read(72, '(4i10)') pwGind(ig), iDumV(1:3)
      !
    enddo
    !
    close(72)
      !! * Close the `grid.ki` file
    !
    open(72, file=trim(system%exportDir)//"/wfc."//trim(ikstr))
      !! * Open the `wfc.ki` file from [[pw_export_for_tme(program)]]
    !
    !> Ignore the first two lines because they are comments
    read(72, * )
    read(72, * )
    !
    do ib = 1, system%iBandL - 1
      do ig = 1, system%npws(ik)
        !! * For each band before `iBandL`, ignore all of the
        !!   plane waves for the given k point
        read(72, *)
        !
      enddo
    enddo
    !
    system%wfc(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)
      !! * Initialize the wavefunction to complex double zero
    !
    do ib = system%iBandL, system%iBandH
      do ig = 1, system%npws(ik)
        !! * For bands between `iBandL` and `iBandH`,
        !!   read in all of the plane waves for the given k point
        !!   and store them in the proper index of the system's `wfc`
        !
        read(72, '(2ES24.15E3)') wfc
        system%wfc(pwGind(ig), ib) = wfc
        !
      enddo
    enddo
    !
    close(72)
      !! * Close the `wfc.ki` file
    !
    deallocate ( pwGind )
      !! * Deallocate space for `pwGind`
    !
    return
    !
  end subroutine readWfc
  !
  !
  subroutine readProjections(ik, system)
    !! Read in the projection \(\langle\beta|\Psi\rangle\) for each band
    !!
    !! <H2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    integer :: i, j
      !! Loop index
    !
    character(len = 300) :: ikstr
      !! String version of k point index
    TYPE(crystal), intent(inout) :: system
      !! Holds the structure for the system you are working on
      !! (either `perfectCrystal` or `solidDefect`)
    !
    call int2str(ik, ikstr)
      !! * Convert the k point index to a string
    !
    system%cProj(:,:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      !! * Initialize `cProj` to all complex double zero
    !
    open(72, file=trim(system%exportDir)//"/projections."//trim(ikstr))
      !! * Open the `projections.ik` file from [[pw_export_for_tme(program)]]
    !
    read(72, *)
      !! * Ignore the first line as it is a comment
    !
    !write(6,'("Solid defect nBands: ", i3)') solidDefect%nBands
    !write(6,'("Solid defect nSpins: ", i3)') solidDefect%nSpins
    !write(6,'("Perfect crystal nBands: ", i3)') perfectCrystal%nBands
    !write(6,'("Perfect crystal nSpins: ", i3)') perfectCrystal%nSpins
    !! @todo Get actual perfect crystal and solid defect output to test @endtodo
    !! @todo Figure out if loop should be over `solidDefect` or `perfectCrystal` @endtodo
    !! @todo Look into `nSpins` to figure out if it is needed @endtodo
    do j = 1, solidDefect%nBands  ! number of bands 
      do i = 1, system%nProjs ! number of projections
        !! * For each band, read in the projections \(\langle\beta|\Psi\rangle\)
        !
        read(72,'(2ES24.15E3)') system%cProj(i,j,1)
        !
      enddo
    enddo
    !
    close(72)
    !
    return
    !
  end subroutine readProjections
  !
  !
  subroutine projectBeta(ik, betaSystem, projectedSystem)
    !! @todo Figure out what this subroutine really does
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    integer :: ig, i, j
      !! Loop index
    integer :: iDumV(3)
      !! Dummy variable to ignore input from file
    integer, allocatable :: pwGind(:)
      !! Indices for the wavefunction of a given k point
    !
    character(len = 300) :: ikstr
      !! String version of the k point index
    !
    TYPE(crystal), intent(inout) :: betaSystem
      !! Holds the structure for the system you are getting \(\beta\) from
      !! (either `perfectCrystal` or `solidDefect`)
    TYPE(crystal), intent(inout) :: projectedSystem
      !! Holds the structure for the system you are projecting
      !! (either `perfectCrystal` or `solidDefect`)
    !
    call int2str(ik, ikstr)
      !! * Convert the k point index to a string
    !
    ! Reading PC projectors
    !
    open(72, file=trim(betaSystem%exportDir)//"/grid."//trim(ikstr))
      !! * Open the `grid.ki` file from [[pw_export_for_tme(program)]]
    !
    !> * Ignore the next two lines as they are comments
    read(72, * )
    read(72, * )
    !
    allocate ( pwGind(betaSystem%npws(ik)) )
      !! * Allocate space for `pwGind`
    !
    do ig = 1, betaSystem%npws(ik)
      !! * Read in the index for each plane wave
      !
      read(72, '(4i10)') pwGind(ig), iDumV(1:3)
      !
    enddo
    !
    close(72)
      !! * Close the `grid.ki` file
    !
    !
    allocate ( betaSystem%beta(numOfPWs, betaSystem%nProjs) )
      !! * Allocate space for \(|\beta\rangle\)
    !
    betaSystem%beta(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
      !! * Initialize all values of \(|\beta\rangle\) to complex double zero
    !
    open(73, file=trim(betaSystem%exportDir)//"/projectors."//trim(ikstr))
      !! * Open the `projectors.ki` file from [[pw_export_for_tme(program)]]
    !
    read(73, *) 
      !! * Ignore the first line because it is a comment
    read(73, *) 
      !! * Ignore the second line because it is the number of projectors that
      !!   was already calculated in [[TMEModule(module):readQEExport(subroutine)]]
      !!   and the number of plane waves for a given k point that was read in in the
      !!   same subroutine
    !
    do j = 1, betaSystem%nProjs
      do i = 1, betaSystem%npws(ik)
        !! * Read in each \(|\beta\rangle\) and store in the proper index of `beta`
        !!   for the system
        !
        read(73,'(2ES24.15E3)') betaSystem%beta(pwGind(i),j)
        !
      enddo
    enddo
    !
    close(73)
    !
    deallocate ( pwGind )
      !! * Deallocate space for `pwGind`
    !
    do j = projectedSystem%iBandL, projectedSystem%iBandH
      do i = 1, betaSystem%nProjs
        !! * Calculate \(\langle\beta|\Phi\rangle\) between 
        !!   `iBandL` and `iBandH`
        !
        betaSystem%cCrossProj(i,j,1) = sum(conjg(betaSystem%beta(:,i))*projectedSystem%wfc(:,j))
        !
      enddo
    enddo
    !
    deallocate ( betaSystem%beta )
      !! * Deallocate space for \(|\beta\rangle\)
    !
    return
    !
  end subroutine projectBeta
  !
  !
  subroutine pawCorrectionWfc(system)
    !! Calculates the augmentation part of the transition matrix element
    !! @todo Figure out what this subroutine really does @endtodo
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    integer :: iIon
      !! Loop index over atoms
    integer :: iProj, jProj
      !! Loop index of projectors
    integer :: ibi, ibf
      !! Loop index over bands
    integer :: m, mPrime
      !! Loop index for magnetic quantum number for a given projector
    integer :: ispin 
    integer :: LMBASE
    integer :: LM, LMP
      !! Index for cProj
    integer :: l, lPrime
      !! Angular momentum quantum number for a given projector
    integer :: iAtomType
      !! Atom type index for a given ion in the system
    !
    real(kind = dp) :: atomicOverlap
    !
    complex(kind = dp) :: cProjIe, cProjFe
    !
    TYPE(crystal), intent(inout) :: system
      !! Holds the structure for the system you are working on
      !! (either `perfectCrystal` or `solidDefect`)
    !
    ispin = 1
      !! * Set the value of `ispin` to 1
      !! @note
      !! `ispin` never has a value other than one, so I'm not sure
      !!  what its purpose is
      !! @endnote
    !
    system%paw_Wfc(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
      !! * Initialize all values in `paw_Wfc` to complex double zero
    !
    LMBASE = 0
      !! * Initialize the base offset for `cProj`'s first index to zero
    !
    do iIon = 1, system%nIons 
      !! * For each atom in the system
      !!    * Get the index for the atom type
      !!    * Loop over the projectors twice, each time finding the
      !!      angular momentum quantum number (\(l\) and \(l^{\prime}\))
      !!      and magnetic quantum number (\(m\) and \(m^{\prime}\))
      !!    * If \(l = l^{\prime}\) and \(m = m^{\prime}\), loop over the bands to
      !!      calculate `paw_Wfc`
      !!
      !! @todo Figure out the significance of \(l = l^{\prime}\) and \(m = m^{\prime}\) @endtodo
      !
      iAtomType = system%atomTypeIndex(iIon)
      !
      LM = 0
      !
      do iProj = 1, system%atoms(iAtomType)%numProjs
        !
        l = system%atoms(iAtomType)%projAngMom(iProj)
        !
        do m = -l, l
          !
          LM = LM + 1 !1st index for CPROJ
          !
          LMP = 0
          !
          do jProj = 1, system%atoms(iAtomType)%numProjs
            !
            lPrime = system%atoms(iAtomType)%projAngMom(jProj)
            !
            do mPrime = -lPrime, lPrime
              !
              LMP = LMP + 1 ! 2nd index for CPROJ
              !
              atomicOverlap = 0.0_dp
              !
              if ( (l == lPrime).and.(m == mPrime) ) then 
                !
                atomicOverlap = sum(system%atoms(iAtomType)%F1(:,iProj, jProj))
                !
                do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
                  !
                  !> @todo Figure out why the difference between SD and PC @endtodo
                  if ( system%crystalType == 'PC' ) then
                    !
                    cProjIe = system%cProj(LMP + LMBASE, ibi, ISPIN)
                    !
                  else if ( system%crystalType == 'SD' ) then
                    !
                    cProjIe = system%cCrossProj(LMP + LMBASE, ibi, ISPIN)
                    !
                  endif
                  !
                  do ibf = solidDefect%iBandL, solidDefect%iBandH
                    !
                    !> @todo Figure out why the difference between SD and PC @endtodo
                    if ( system%crystalType == 'PC' ) then
                      !
                      cProjFe = conjg(system%cCrossProj(LM + LMBASE, ibf, ISPIN))
                      !
                    else if ( system%crystalType == 'SD' ) then
                      !
                      cProjFe = conjg(system%cProj(LM + LMBASE, ibf, ISPIN))
                      !
                    endif
                    !
                    system%paw_Wfc(ibf, ibi) = system%paw_Wfc(ibf, ibi) + cProjFe*atomicOverlap*cProjIe
                    !
                  enddo
                  !
                enddo
                !
              endif
              !
            enddo
            !
          enddo
          !
        enddo
        !
      enddo
      !
      LMBASE = LMBASE + system%atoms(iAtomType)%lmMax
      !
    enddo
    !
    return
    !
  end subroutine pawCorrectionWfc
  !
  !
  subroutine pawCorrectionK(system)
    !! @todo Figure out what this subroutine really does @endtodo
    !!
    !! <h2>Walkthrough</h2>
    !
    implicit none
    !
    !integer, intent(in) :: ik
    !
    integer :: ibi, ibf
      !! Loop index over bands
    integer :: iPW
      !! Loop index over plane waves for a given process
    integer :: iProj
      !! Loop index over projectors
    integer :: iR
      !! Loop index over radial mesh (up to augmentation sphere)
    integer :: iAtomType
      !! Loop index over atom types
    integer :: iIon
      !! Loop index over ions in system
    integer :: l
      !! Angular momentum quantum number
    integer :: m
      !! Magnetic quantum number
    integer :: ispin
    integer :: LMBASE
    integer :: LM
    integer :: ind
    !
    real(kind = dp) :: qDotR
      !! \(\mathbf{G}\cdot\mathbf{r}\)
    real(kind = dp) :: t1
      !! Start time
    real(kind = dp) :: t2
      !! End time
    real(kind = dp) :: v_in(3)
      !! Unit vector in the direction of \(\mathbf{G}\)
    real(kind = dp) :: JL(0:JMAX)
      !! Spherical bessel functions for a point up to `JMAX`
    real(kind = dp) :: q 
    real(kind = dp) :: FI 
    !
    complex(kind = dp) :: Y( (JMAX+1)**2 )
      !! All spherical harmonics up to some max momentum
    complex(kind = dp) :: ATOMIC_CENTER
      !! \(e^{-i\mathbf{G}\cdot\mathbf{r}}\)
    complex(kind = dp) :: VifQ_aug
    !
    TYPE(crystal), intent(inout) :: system
      !! Holds the structure for the system you are working on
      !! (either `perfectCrystal` or `solidDefect`)
    !
    ispin = 1
      !! * Set the value of `ispin` to 1
      !! @note
      !! `ispin` never has a value other than one, so I'm not sure
      !!  what its purpose is
      !! @endnote
    !
    call cpu_time(t1)
      !! * Start a timer
    !
    system%pawK(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
      !! * Initialize all values in `pawK` to complex double zero
    !
    do iPW = nPWsI(myid), nPWsF(myid) 
      !! * Loop through the plane waves for a given process
      !
      if ( myid == root ) then 
        if ( (iPW == nPWsI(myid) + 1000) .or. (mod(iPW, 25000) == 0) .or. (iPW == nPWsF(myid)) ) then
          !! * If this is the root process, output a status update every 1000 plane waves
          !!   and every multiple of 25000, giving an estimate of the time remaining at each step
          !!   @todo Figure out if this output slows things down significantly @endtodo
          !!   @todo Figure out if formula gives accurate representation of time left @endtodo
          !
          call cpu_time(t2)
          !
          write(iostd, '("        Done ", i10, " of", i10, " k-vecs. ETR : ", f10.2, " secs.")') &
                iPW, nPWsF(myid) - nPWsI(myid) + 1, (t2-t1)*(nPWsF(myid) - nPWsI(myid) + 1 -iPW )/iPW
          !
          flush(iostd)
          !
          !call cpu_time(t1)
          !
        endif
      endif
      !
      q = sqrt(sum(gvecs(:,iPW)*gvecs(:,iPW)))
        !! * Calculate `q` as \(\sqrt{\mathbf{G}\cdot\mathbf{G}}\)
        !!   to get length of \(\mathbf{G}\)
      !
      !> * Define a unit vector in the direction of \(\mathbf{G}\), 
      !>   but only divide by the length if it is bigger than 
      !>   \(1\times10^{-6}\) to avoid dividing by very small numbers
      v_in(:) = gvecs(:,iPW)
      if ( abs(q) > 1.0e-6_dp ) v_in = v_in/q 
      !
      Y = cmplx(0.0_dp, 0.0_dp, kind = dp)
        !! * Initialize the spherical harmonics to complex double zero
      call ylm(v_in, JMAX, Y) 
        !! * Calculate spherical harmonics with argument `v_in` up to 
        !!   \(Y_{J_{\text{max}}}^{\pm J_{\text{max}}}\)
      !
      LMBASE = 0
        !! * Initialize the base offset for `cProj`'s first index to zero
      !
      do iAtomType = 1, system%numOfTypes
        !
        do iR = 1, system%atoms(iAtomType)%iRAugMax 
          !! * For each atom type, loop through the r points
          !!   in the augmentation sphere and calculate the 
          !!   spherical Bessel functions from 0 to `JMAX`
          !!   at each point
          !
          JL = 0.0_dp
          !
          call bessel_j(q*solidDefect%atoms(iAtomType)%r(iR), JMAX, JL) ! returns the spherical bessel at qr point
            !! @todo Figure out if this should be `system` @endtodo
            !! @todo Figure out significance of "qr" point @endtodo
          !
          system%atoms(iAtomType)%bes_J_qr(:,iR) = JL(:)
            !! @todo Test if can just directly store in each atom type's `bes_J_qr` @endtodo
          !
        enddo
        !
      enddo
      !
      do iIon = 1, system%nIons 
        !! * For each atom in the system
        !!    * Calculate \(\mathbf{G}\cdot\mathbf{r}\)
        !!    * Calculate \(e^{-i\mathbf{G}\cdot\mathbf{r}}\)
        !!    * Get the index for the atom type
        !!    * Loop over the projectors, finding \(l, m\) for each
        !!    * For each possible m
        !!       * Calculate \(\text{FI} = j_l\cdot F\) where \(j_l\) is
        !!         the Bessel function and \(F\) is for a given projector
        !!       * Calculate \(\text{VifQ_aug} = e^{-i\mathbf{G}\cdot\mathbf{r}}
        !!         Y_l^m(\mathbf{G}/|\mathbf{G}|)(-i)^l\text{FI}\)
        !!       * Loop over the bands, summing `VifQ_aug*cProj` to get `pawK`
        !
        qDotR = sum(gvecs(:,iPW)*system%posIon(:,iIon))
          !! @todo Figure out if this should be `gDotR` @endtodo
        !
        !> @todo Figure out why this is called `ATOMIC_CENTER` @endtodo
        !> @todo Figure out why the difference between SD and PC @endtodo
        if ( system%crystalType == 'PC' ) then
          !
          ATOMIC_CENTER = exp( -ii*cmplx(qDotR, 0.0_dp, kind = dp) )
          !
        else if ( system%crystalType == 'SD' ) then
          !
          ATOMIC_CENTER = exp( ii*cmplx(qDotR, 0.0_dp, kind = dp) )
          !
        endif
        !
        iAtomType = system%atomTypeIndex(iIon)
        !
        LM = 0
        !
        do iProj = 1, system%atoms(iAtomType)%numProjs
          !
          l = system%atoms(iAtomType)%projAngMom(iProj)
          !
          do m = -l, l
            !
            LM = LM + 1 !1st index for CPROJ
            !
            FI = 0.0_dp
            !
            FI = sum(system%atoms(iAtomType)%bes_J_qr(l,:)*system%atoms(iAtomType)%F(:,iProj)) ! radial part integration F contains rab
            !
            ind = l*(l + 1) + m + 1 ! index for spherical harmonics
            !
            !> @todo Figure out why the difference between SD and PC @endtodo
            if ( system%crystalType == 'PC' ) then
              !
              VifQ_aug = ATOMIC_CENTER*Y(ind)*(-II)**l*FI
              !
            else if ( system%crystalType == 'SD' ) then
              !
              VifQ_aug = ATOMIC_CENTER*conjg(Y(ind))*(II)**l*FI
              !
            endif
            !
            do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
              !
              do ibf = solidDefect%iBandL, solidDefect%iBandH
                !
                !> @todo Figure out why the difference between SD and PC @endtodo
                if ( system%crystalType == 'PC' ) then
                  !
                  system%pawK(ibf, ibi, iPW) = system%pawK(ibf, ibi, iPW) + &
                                                     VifQ_aug*system%cProj(LM + LMBASE, ibi, ISPIN)
                  !
                else if ( system%crystalType == 'SD' ) then
                  !
                  system%pawK(ibf, ibi, iPW) = system%pawK(ibf, ibi, iPW) + &
                                                     VifQ_aug*conjg(system%cProj(LM + LMBASE, ibi, ISPIN))
                  !
                endif
                !
              enddo
              !
            enddo
            !
          ENDDO
        ENDDO
        LMBASE = LMBASE + system%atoms(iAtomType)%lmMax
      ENDDO
      !
    enddo
    !
    !system%pawK(:,:,:) = system%pawK(:,:,:)*4.0_dp*pi/sqrt(solidDefect%omega)
    !
    return
    !
  end subroutine pawCorrectionK
  !
  !
  subroutine ylm(v_in,lmax,y)
  !! Returns the [spherical harmonics](http://mathworld.wolfram.com/SphericalHarmonic.html) 
  !! for a given argument vector up to the maximum value of \(l\) given
  !!
  !! <h2>Description</h2>
  !!   <h3>Purpose</h3>
  !!        The spherical harmonics (Condon and Shortley convention)
  !!          \(Y_0^0,Y_1^{-1},Y_1^0,Y_1^1,Y_2^{-2} ... Y_{l_{\text{max}}}^{\pm l_{\text{max}}}\)
  !!        for vector \(\mathbf{V}\) (given in Cartesian coordinates)
  !!        are calculated. In the Condon Shortley convention the
  !!        spherical harmonics are defined as
  !!        \[ Y_l^m = (-1)^m \sqrt{\frac{1}{\pi}} P_l^m(\cos{\theta})
  !!        e^{im\phi} \]
  !!        where  \(P_l^m(\cos{\theta})\) is the normalized Associated
  !!        Legendre function. Thus,
  !!                     \[  Y_l^{-m} = (-1)^m (Y_l^m)^* \]
  !!
  !!   <h3>Usage</h3>
  !!
  !!
  !!        DOUBLE PRECISION V(3), Y(5*5)
  !!        V(1) = ...
  !!        V(2) = ...
  !!        V(3) = ...
  !!        CALL YLM(V,4,Y)
  !!
  !!   <h3>Argument Description</h3>
  !!     <ul>
  !!          <li> 
  !!                  `V`      - `DOUBLE PRECISION` vector, dimension 3        (input)<br/>
  !!                   Must be given in Cartesian coordinates.
  !!                   Conversion of V to polar coordinates gives the
  !!                   angles \(\theta\) and \(\phi\) necessary for the calculation
  !!                   of the spherical harmonics.  
  !!          </li>
  !!          <li>
  !!                   `LMAX`   - `INTEGER` value                               (input)<br/>
  !!                   upper bound of \(l\) for which spherical harmonics
  !!                   will be calculated<br/>
  !!                   constraint: `LMAX >= 0`  
  !!          </li>
  !!          <li>
  !!                   `Y`      - `COMPLEX*16` array, dimension `(LMAX+1)**2`    (output)<br/>
  !!                   contains the calculated spherical harmonics<br/>
  !!                   `Y(1)` for \(l=0\) (\(m = 0\))<br/>
  !!                   `Y(2), ..., Y(4)` for \(l = 1\) (\(m = -1, 0, 1\))<br/>
  !!                   ...<br/>
  !!                   `Y(LMAX*LMAX+1), ..., Y((LMAX+1)*(LMAX+1))` for \(l = l_{\text{max}}\) 
  !!                            (\(m = -l,...,l\))<br/>
  !!                   constraint: Dimension of `Y` \(\geq (l_{\text{max}} + 1)^2\) (not checked)
  !!          </li>
  !!        </ul>
  !!
  !!   <h3>Used Subroutines (Directly Called)</h3>
  !!           none
  !!
  !!   <h3>Indirectly Called Subroutines</h3>
  !!           none
  !!
  !!   <h3>Input/Output (Read/Write)</h3>
  !!           none
  !!
  !!   <h3>Machine Dependent Program Parts</h3>
  !!           Type `COMPLEX*16` is used which does not conform to the
  !!           FORTRAN 77 standard.
  !!           Also the non-standard type conversion function `DCMPLX()`
  !!           is used which combines two double precision values into
  !!           one double complex value.
  !!
  !!   <h3>Method</h3>
  !!           The basic algorithm used to calculate the spherical
  !!           harmonics for vector \(\mathbf{V}\) is as follows:
  !!
  !!
  !!           Y(0,0)
  !!           Y(1,0)
  !!           Y(1,1)
  !!           Y(1,-1) = -Y(1,1)
  !!           DO L = 2, LMAX
  !!              Y(L,L)   = f(Y(L-1,L-1)) ... Formula 1
  !!              Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
  !!              DO M = L-2, 0, -1
  !!                 Y(L,M) = f(Y(L-1,M),Y(L-2,M)) ... Formula 2
  !!                 Y(L,-M)= (-1)**M*CONJG(Y(L,M))
  !!              ENDDO
  !!           ENDDO
  !!
  !!   <h3>Formulas</h3>
  !!        Starting values:
  !!          \[Y_0^0 = \sqrt{\dfrac{1}{4\pi}}\]
  !!          \[Y_1^0 = \sqrt{\dfrac{3}{4\pi}}\cos\theta\]
  !!          \[Y_1^1 = -\sqrt{\dfrac{3}{8\pi}}\sin\theta e^{i\phi}\]
  !!        Formula 1:
  !!          \[Y_l^l = -\sqrt{\dfrac{2l+1}{2l}}\sin\theta e^{i\phi}Y_{l-1}^{l-1}\]
  !!        Formula 2:
  !!          \[Y_l^m = \sqrt{\dfrac{(2l-1)(2l+1)}{(l-m)(l+m)}}\cos\theta Y_{l-1}^m - 
  !!                    \sqrt{\dfrac{(l-1+m)(l-1-m)(2l+1)}{(2l-3)(l-m)(l+m)}} Y_{l-2}^m\]
  !!        Formula 3: (not used in the algorithm because of the division
  !!                    by \(\sin\theta\) which may be zero)
  !!          \[Y_l^m = -\sqrt{\dfrac{4(m+1)(m+1)}{(l+m+1)(l-m)}}\dfrac{\cos\theta}{\sin\theta}e^{i\phi}Y_1^{m+1} -
  !!                    \sqrt{\dfrac{(l-m-1)(l+m+2)}{(l-m)(l+m+1)}}e^{-2i\phi}Y_l^{m+2}\]
  !!
  ! !REVISION HISTORY:
  !   26. April 1994                                   Version 1.2
  !   Taken 8 1 98 from SRC_lapw2 to SRC_telnes
  !   Updated November 2004 (Kevin Jorissen)
  !   cosmetics March 2005 (Kevin Jorissen)
  !
      implicit none
  !
  !   In/Output :
  !
      integer, intent(in) :: LMAX
        !! Spherical harmonics are calculated for 
        !! \(l = 0, 1, ..., l_{\text{max}}\)
      real(kind = dp), intent(in) :: V_in(3)
        !! Vector, argument of the spherical harmonics (we calculate
        !! \(Y_l^m(\mathbf{v}/|\mathbf{v}|)\))
      complex(kind = dp), intent(out) :: Y(*)
        !! Array containing \(Y_l^m(\mathbf{v})\) for several \(l,m\)
  !
  !   Local variables :
      real(kind = dp), parameter :: pi = 3.1415926535897932384626433_dp
  !
      INTEGER         ::  I2L, I4L2, INDEX, INDEX2, L, M, MSIGN
      real(kind = dp) ::  A, B, C, AB, ABC, ABMAX, ABCMAX, V(3)
      real(kind = dp) ::  D4LL1C, D2L13
      real(kind = dp) ::  COSTH, SINTH, COSPH, SINPH
      real(kind = dp) ::  TEMP1, TEMP2, TEMP3
      real(kind = dp) ::  YLLR, YLL1R, YL1L1R, YLMR
      real(kind = dp) ::  YLLI, YLL1I, YL1L1I, YLMI
      !
      ! Y(0,0)
      !
      do INDEX = 1,3
        V(INDEX) = dble(V_in(INDEX))
      enddo
      YLLR = 1.0_dp/sqrt(4.0_dp*PI)
      YLLI = 0.0_dp
      Y(1) = CMPLX(YLLR, YLLI, kind = dp)
      !
      ! continue only if spherical harmonics for (L .GT. 0) are desired
      !
      IF (LMAX .LE. 0) GOTO 999
      !
      ! calculate sin(Phi), cos(Phi), sin(Theta), cos(Theta)
      ! Theta, Phi ... polar angles of vector V
      !
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
      !
      ! Y(1,0)
      !
      Y(3) = CMPLX(sqrt(3.0_dp)*YLLR*COSTH, 0.0_dp, kind = dp)
      !
      ! Y(1,1) ( = -DCONJG(Y(1,-1)))
      !
      TEMP1 = -SQRT(1.5_dp)*YLLR*SINTH
      Y(4) = CMPLX(TEMP1*COSPH,TEMP1*SINPH, kind = dp)
      Y(2) = -CONJG(Y(4))
      !
      DO L = 2, LMAX
        INDEX  = L*L + 1
        INDEX2 = INDEX + 2*L
        MSIGN  = 1 - 2*MOD(L,2)
        !
        ! YLL = Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
        !
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
        !
        ! YLL1 = Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
        ! (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
        !
        TEMP2 = SQRT(DBLE(2*L+1))*COSTH
        YLL1R = TEMP2*YL1L1R
        YLL1I = TEMP2*YL1L1I
        Y(INDEX2) = CMPLX(YLL1R,YLL1I, kind = dp)
        Y(INDEX)  = -cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
  !      Y(INDEX)  = -dble(MSIGN)*CONJG(Y(INDEX2))
        INDEX2 = INDEX2 - 1
        INDEX  = INDEX  + 1
        !
        I4L2 = INDEX2 - 4*L + 2
        I2L  = INDEX2 - 2*L
        D4LL1C = COSTH*SQRT(DBLE(4*L*L-1))
        D2L13  = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
        !
        DO M = L - 2, 0, -1
          !
          ! YLM = Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
          !
          TEMP1 = 1.0_dp/SQRT(DBLE((L+M)*(L-M)))
          TEMP2 = D4LL1C*TEMP1
          TEMP3 = D2L13*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1
          YLMR = TEMP2*DBLE(Y(I2L))  + TEMP3*DBLE(Y(I4L2))
          YLMI = TEMP2*DIMAG(Y(I2L)) + TEMP3*DIMAG(Y(I4L2))
          Y(INDEX2) = CMPLX(YLMR,YLMI, kind = dp)
          Y(INDEX)  = cmplx(MSIGN,0.0_dp,kind=dp)*CONJG(Y(INDEX2))
    !      Y(INDEX)  = dble(MSIGN)*CONJG(Y(INDEX2))
          !
          MSIGN  = -MSIGN
          INDEX2 = INDEX2 - 1
          INDEX  = INDEX  + 1
          I4L2   = I4L2   - 1
          I2L    = I2L    - 1
        ENDDO
      ENDDO
      !
  999 RETURN
  END subroutine ylm
  !
  !
  subroutine bessel_j (x, lmax, jl)
    !! Generates the 
    !! [spherical bessel function of the first kind](http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html)
    !! for the given argument \(x\) and all possible indices from 0 to `lmax` 
    !!
    !! <h2>Walkthrough</h2>
    !!
    !
    implicit none
    !
    integer, intent(in) :: lmax
    integer :: l
    !
    real(kind = dp), intent(in) :: x
    real(kind = dp), intent(out) :: jl(0:lmax)
    !
    if (x <= 0.0_dp) then
      !! * If \(x\) is less than zero, return 0 for all
      !!   indices but 0 which is 1
      !
      jl = 0.0_dp
      jl(0) = 1.0_dp
      !
      return
      !
    end if
    !
    !> * Explicitly calculate the first 2 functions so can use 
    !>   recursive definition for later terms
    jl(0) = sin(x)/x
    if (lmax <= 0) return
    jl(1) = (jl(0)-cos(x))/x
    if (lmax == 1) return
    !
    do l = 2, lmax
      !! * Define the rest of the functions as
      !!   \[j_l = (2l-1)j_{l-1}/x - j_{l-2}\]
      !
      jl(l) = dble(2*l-1)*jl(l-1)/x - jl(l-2)
      !
    enddo
    !
    return
    !
  end subroutine bessel_j
  !
  !
  subroutine writeResults(ik)
    !! Write out the \(\langle f|U|i\rangle\) matrix and
    !! changes in eigenvalues for a given k point
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    !
    integer :: ibi, ibf
      !! Loop index over bands
    integer :: totalNumberOfElements
      !! Total number of matrix elements
    real(kind = dp) :: t1
      !! Start time
    real(kind = dp) :: t2
      !! End time
    !
    character(len = 300) :: text
      !! String to hold long section header
    character(len = 300) :: Uelements
      !! File name for matrix elements output 
    character(len = 300) :: ikstr
      !! String version of k point index
    !
    call cpu_time(t1)
      !! * Start a timer
    !
    call readEigenvalues(ik)
      !! * Read eigenvalues from [[pw_export_for_tme(program)]]
    !
    write(iostd, '(" Writing Ufi(:,:).")')
      !! * Write out a header to `output` file
    !
    call int2str(ik, ikstr)
      !! * Convert the k point index to a string
    !
    write(Uelements, '("/TMEs_kptI_",a,"_kptF_",a)') trim(ikstr), trim(ikstr)
      !! * Determine what the file name should be based on the k point index
    !
    open(17, file=trim(elementsPath)//trim(Uelements), status='unknown')
      !! * Open the matrix elements output file
    !
    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') solidDefect%omega
      !! * Output cell volume to `Uelements` file
    !
    text = "# Total number of <f|U|i> elements, Initial States (bandI, bandF), Final States (bandI, bandF)"
    write(17,'(a, " Format : ''(5i10)''")') trim(text)
      ! * Output header for next section
    !
    totalNumberOfElements = (perfectCrystal%iBandH - perfectCrystal%iBandL + 1)* &
                            (solidDefect%iBandH - solidDefect%iBandL + 1)
      !! * Calculate the total number of matrix elements
    write(17,'(5i10)') totalNumberOfElements, perfectCrystal%iBandL, perfectCrystal%iBandH, &
                       solidDefect%iBandL, solidDefect%iBandH
      !! * Output the total number of elements and band limits
    !
    write(17, '("# Final Band, Initial Band, Delta energy, Complex <f|U|i>, |<f|U|i>|^2 Format : ''(2i10,4ES24.15E3)''")')
      ! * Output header for next section
    !
    do ibf = solidDefect%iBandL, solidDefect%iBandH
      do ibi = perfectCrystal%iBandL, perfectCrystal%iBandH
        !! * Loop through the bands to output the change 
        !!   in eigenvalues, matrix element, and norm 
        !!   squared matrix element
        !
        write(17, 1001) ibf, ibi, eigvI(ibi) - eigvF(ibf), Ufi(ibf,ibi,ik), abs(Ufi(ibf,ibi,ik))**2
        !    
      enddo
    enddo
    !
    close(17)
      !! * Close the matrix elements output file
    !
    call cpu_time(t2)
      !! * Stop the timer
    !
    write(iostd, '(" Writing Ufi(:,:) done in:                   ", f10.2, " secs.")') t2-t1
      !! * Write to `output` file how long it took to write the matrix 
    !
 1001 format(2i10,4ES24.15E3)
    ! Define format for writing matrix  
    !
    return
    !
  end subroutine writeResults
  !
  !
  subroutine readUfis(ik)
    !! Read in matrix elements for a given k point
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    integer :: ibi, ibf
      !! Loop index over bands
    integer :: totalNumberOfElements
      !! Total number of matrix elements
    integer :: iDum
      !! Dummy variable to ignore input
    integer :: iEl
      !! Loop index over matrix elements
    !
    real(kind = dp) :: rDum
      !! Dummy variable to ignore input
    real(kind = dp) :: t1
      !! Start time
    real(kind = dp) :: t2
      !! End time
    !
    complex(kind = dp):: cUfi
      !! Temporarily store value of each matrix element
    !
    character(len = 300) :: Uelements
      !! File name for matrix elements output 
    character(len = 300) :: ikstr
      !! String version of k point index
    !
    call cpu_time(t1)
      !! * Start a timer
    !
    write(iostd, '(" Reading Ufi(:,:) of k-point: ", i4)') ik
      !! * Write out a header to the `output` file
    !
    call int2str(ik, ikstr)
      !! * Convert the k point index to a string
    !
    write(Uelements, '("/TMEs_kptI_",a,"_kptF_",a)') trim(ikstr), trim(ikstr)
      !! * Determine what the file name should be based on the k point index
    !
    open(17, file=trim(elementsPath)//trim(Uelements), status='unknown')
      !! * Open the matrix elements file
    !
    read(17, *) 
    read(17, *) 
      !! * Ignore the first two lines as they are comments
    !
    read(17,'(5i10)') totalNumberOfElements, iDum, iDum, iDum, iDum
      !! * Read in the total number of matrix elements
    !
    read(17, *) 
      !! * Ignore the next line as it is a comment
    !
    do iEl = 1, totalNumberOfElements
      !! * Read in the indices and value for each
      !!   matrix element, ignoring the change in
      !!   eigenvalues and norm squared matrix
      !!   element
      !
      read(17, 1001) ibf, ibi, rDum, cUfi, rDum
      Ufi(ibf,ibi,ik) = cUfi
      !    
    enddo
    !
    close(17)
    !
    call cpu_time(t2)
      !! * Stop timer
    !
    write(iostd, '(" Reading Ufi(:,:) done in:                   ", f10.2, " secs.")') t2-t1
      !! * Write out time to read `Ufi` to `output` file
    !
 1001 format(2i10,4ES24.15E3)
    ! Define format to read in matrix elements
    !
    return
    !
  end subroutine readUfis
  !
  !
  subroutine calculateVfiElements()
    !! @todo Figure out what the purpose of this function is. For plotting? @endtodo
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: ik
      !! Loop index over k points
    integer :: ib
      !! Loop index over bands
    integer :: nOfEnergies, iE
    !
    real(kind = dp) :: eMin, eMax, E, av, sd, x, epsilon_if, DHifMin
    !
    real(kind = dp), allocatable :: sumWk(:), sAbsVfiOfE2(:), absVfiOfE2(:)
    integer, allocatable :: nKsInEbin(:)
    !
    character (len = 300) :: text
    !
    allocate( DE(perfectCrystal%iBandL:perfectCrystal%iBandH, perfectCrystal%nKpts))
    allocate( absVfi2(perfectCrystal%iBandL:perfectCrystal%iBandH, perfectCrystal%nKpts) )
      !! * Allocate space for `DE` and `absVfi2`
    ! 
    DE(:,:) = 0.0_dp
    absVfi2(:,:) = 0.0_dp 
      !! * Initialize `DE` and `absVfi2` to double zero
    !
    do ik = 1, perfectCrystal%nKpts
      !! * For each k point
      !!    * Read in the initial and final eigenvalues
      !!    * For each band between `iBandL_PC` and `iBandH_PC`
      !!       * Calculate \(\epsilon_{if} = \epsilon_f - \epsilon_i\)
      !!         as defined in paper
      !!       * Calculate 
      !!         \[|\Delta H_{if}|^2 = \dfrac{|\langle\Phi_f|\Psi_i\rangle|^2 - |\langle\Phi_f|\Psi_i\rangle|^4}
      !!                                {(1 - 2|\langle\Phi_f|\Psi_i\rangle|^2)^2}\epsilon_{if}\]
      !!         (A8 in paper)
      !!         @note 
      !!         Only numerator is calculated because the denominator
      !!         is approximately zero, assuming \(|\langle\Phi_f|\Psi_i\rangle| \ll 1\)
      !!         @endnote
      !!       * Calculate `DE`\( = E_i - E_f = \sqrt{\epsilon_{if}^2 + 4|\Delta H_{if}|^2}\)
      !
      eigvI(:) = 0.0_dp
      eigvF(:) = 0.0_dp
      !
      call readEigenvalues(ik)
      !
      do ib = perfectCrystal%iBandL, perfectCrystal%iBandH
        !
        epsilon_if = eigvF(solidDefect%iBandL) - eigvI(ib)
        absVfi2(ib,ik) = epsilon_if**2*( abs(Ufi(solidDefect%iBandL,ib,ik))**2 - &
                         abs(Ufi(solidDefect%iBandL,ib,ik))**4 )
        !
        DE(ib, ik) = sqrt(epsilon_if**2 - 4.0_dp*absVfi2(ib,ik))
        !
      enddo
      !
    enddo
    !
    eMin = minval( DE(:,:) )
    eMax = maxval( DE(:,:) )
      !! * Find the max and min of `DE`
    !
    nOfEnergies = int((eMax-eMin)/eBin) + 1
      !! * Use the min and max values to calculate how many
      !!   energy bins result between the energies using
      !!   a bin size of `eBin`
    !
    allocate ( absVfiOfE2(0:nOfEnergies), nKsInEbin(0:nOfEnergies), sumWk(0:nOfEnergies) )
      !! * Allocate space for `absVfiOfE2`, `nKsInEbin`, and `sumWk`
    !
    absVfiOfE2(:) = 0.0_dp
    nKsInEbin(:) = 0
    sumWk(:) = 0.0_dp
      !! * Initialize all to double zero
    !
    do ik = 1, perfectCrystal%nKpts
      !
      do ib = perfectCrystal%iBandL, perfectCrystal%iBandH
        !! * For each k point and band
        !!    * Store `absVfi2` (\(|\Delta H_{if}|^2\)) if have min `DE`
        !!    * Find the "index" of the particular `DE`
        !!    * If \(|\Delta H_{if}|^2 > 0\) for given k point, add the 
        !!      value to a cumulative sum, weighted by the k point weight `wk`,
        !!      and keep track of the total weight and number of k points in
        !!      each bin
        !
        if ( abs( eMin - DE(ib,ik)) < 1.0e-3_dp ) DHifMin = absVfi2(ib, ik)
          !! @todo Figure out why `DHifMin` is needed @endtodo
          !! @todo Figure out why used difference rather than `==
        !
        iE = int((DE(ib, ik)-eMin)/eBin)
        !
        if ( absVfi2(ib, ik) > 0.0_dp ) then
          !! @todo Figure out why this test is here. All of these should be positive, right? @endtodo
          !
          absVfiOfE2(iE) = absVfiOfE2(iE) + perfectCrystal%wk(ik)*absVfi2(ib, ik)
            !! @note I think `absVfiOfE2` is the sum of the numerators weighted by k point @endnote
          !
          sumWk(iE) = sumWk(iE) + perfectCrystal%wk(ik)
          !
          nKsInEbin(iE) = nKsInEbin(iE) + 1
          !
        else
          write(iostd,*) 'lalala', absVfi2(ib, ik)
        endif
        !
      enddo
      !
    enddo
    !
    allocate ( sAbsVfiOfE2(0:nOfEnergies) )
      !! * Allocate space for `sAbsVfiOfE2`
    !
    sAbsVfiOfE2 = 0.0_dp
      !! * Initialize to zero
    !
    open(11, file=trim(VfisOutput)//'ofKpt', status='unknown')
      !! * Open `VfisVsEofKpt` file
    !
    write(11, '("# |<f|V|i>|^2 versus energy for all the k-points.")')
      !! * Output file header
    !
    write(text, '("# Energy (eV) shifted by half eBin, |<f|V|i>|^2 (Hartree)^2,")')
    write(11, '(a, " k-point index. Format : ''(2ES24.15E3,i10)''")') trim(text)
      !! * Output section header
    !
    do ik = 1, perfectCrystal%nKpts
      !
      do ib = perfectCrystal%iBandL, perfectCrystal%iBandH
        !! * For each k point and band, 
        !!    * Figure out the "index" for `DE`
        !!    * Calculate \(|\Delta H_{if}|^2\) average over k points
        !!    * Store \(|\Delta H_{if}|^2\) for a given k point and band
        !!    * Write out the energy shifted by half a bin, \(|\Delta H_{if}|^2\),
        !!      and the k point index to `VfisVsEofKpt`
        !!    * Write out unshifted energy, \(|\Delta H_{if}|^2\), and the k point
        !!      index to some other file??
        !!    * Calculate the standard deviation of \(|\Delta H_{if}|^2\)
        !
        iE = int((DE(ib,ik)-eMin)/eBin)
        !
        av = absVfiOfE2(iE)/sumWk(iE)
        !
        x = absVfi2(ib,ik)
        !
        write(11, '(2ES24.15E3,i10)') (eMin + (iE+0.5_dp)*eBin)*HartreeToEv, x, ik
        write(12, '(2ES24.15E3,i10)') DE(ib,ik)*HartreeToEv, absVfi2(ib, ik), ik
          !! @todo Figure out where unit 12 file is opened and what it is @endtodo
        !write(11, '(2ES24.15E3,i10)') (eMin + iE*eBin + eBin/2.0_dp), x, ik
        !
        sAbsVfiOfE2(iE) = sAbsVfiOfE2(iE) + perfectCrystal%wk(ik)*(x - av)**2/sumWk(iE)
        !
      enddo
      !
    enddo
    !
    close(11)
      !! * Close `VfisVsEofKpt` file
    !
    open(63, file=trim(VfisOutput), status='unknown')
      !! * Open `VfisVsE` file
    !
    write(63, '("# Averaged |<f|V|i>|^2 over K-points versus energy.")')
      !! * Output file header
    !
    write(63, '("#                 Cell volume : ", ES24.15E3, " (a.u.)^3,   Format : ''(ES24.15E3)''")') solidDefect%omega
    write(63, '("#   Minimun transition energy : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eMin
    write(63, '("# |DHif|^2 at minimum Tr. En. : ", ES24.15E3, " (Hartree^2),Format : ''(ES24.15E3)''")') DHifMin
    write(63, '("#                  Energy bin : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eBin
      !! * Output cell volumn, min `DE`, \(|\Delta H_{if}|^2\) at min `DE`, and bin size
    !
    write(text, '("# Energy (Hartree), averaged |<f|V|i>|^2 over K-points (Hartree)^2,")')
    write(63, '(a, " standard deviation (Hartree)^2. Format : ''(3ES24.15E3)''")') trim(text)
      !! * Output section header
    !
    do iE = 0, nOfEnergies
      !! * For each energy 
      !!    * Calculate the average and standard deviation
      !!      of \(|\Delta H_{if}|^2\)
      !!    * Output the energy, average, and standard deviation
      !
      E = iE*eBin
        !! @todo Figure out why use `eMin + iE*eBin` rather than `DE` @endtodo
      !
      av = 0.0_dp
      !
      sd = 0.0_dp
      !
      if (nKsInEbin(iE) > 0) then
        !
        av = absVfiOfE2(iE)/sumWk(iE)
        !
        sd = sqrt(sAbsVfiOfE2(iE))
        !
      endif
      !
      write(63,'(3ES24.15E3)') eMin + E, av, sd
      !
    enddo
    !
    close(63)
      !! * Close `VfisVsE` file
    !
    return
    !
  end subroutine calculateVfiElements
  !
  !
  subroutine readEigenvalues(ik)
    !! Read in the eigenvalues output from [[pw_export_for_tme(program)]]
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(in) :: ik
      !! K point index
    integer :: ib
      !! Loop index over bands
    !
    character(len = 300) :: ikstr
      !! String version of k point index
    !
    call int2str(ik, ikstr)
      !! * Convert k point index to string
    !
    open(72, file=trim(perfectCrystal%exportDir)//"/eigenvalues."//trim(ikstr))
      !! * Open the perfectCrystal `eigenvalues.ik` file from [[pw_export_for_tme(program)]]
    !
    read(72, * )
    read(72, * )
      !! * Ignore the first two lines as they are comments
    !
    do ib = 1, perfectCrystal%iBandL - 1
      !! * Ignore eigenvalues for bands that are before `iBandL_PC`
      !
      read(72, *)
      !
    enddo
    !
    do ib = perfectCrystal%iBandL, perfectCrystal%iBandH
      !! * Read in the eigenvalues from `iBandL_PC` to `iBandH_PC`
      !
      read(72, '(ES24.15E3)') eigvI(ib)
      !
    enddo
    !
    close(72)
      !! * Close the perfect crystal `eigenvalues.ik` file
    !
    open(72, file=trim(solidDefect%exportDir)//"/eigenvalues."//trim(ikstr))
      !! * Open the solid defect `eigenvalues.ik` file from [[pw_export_for_tme(program)]]
    !
    read(72, * )
    read(72, * ) 
      !! * Ignore the first two lines as they are comments
    !
    do ib = 1, solidDefect%iBandL - 1
      !! * Ignore eigenvalues for bands that are before `iBandL_SD`
      !
      read(72, *)
      !
    enddo
    !
    do ib = solidDefect%iBandL, solidDefect%iBandH
      !! * Read in the eigenvalues from `iBandL_SD` to `iBandH_SD`
      !
      read(72, '(ES24.15E3)') eigvF(ib)
      !
    enddo
    !
    close(72)
      !! * Close the solid defect `eigenvalues.ik` file
    !
    return
    !
  end subroutine readEigenvalues
  !
  !
  subroutine finalizeCalculation()
    !! Stop timer, write out total time taken, and close the output file
    !
    implicit none
    !
    write(iostd,'("-----------------------------------------------------------------")')
    !
    call cpu_time(tf)
    write(iostd, '(" Total time needed:                         ", f10.2, " secs.")') tf-t0
    !
    close(iostd)
    !
    return
    !
  end subroutine finalizeCalculation
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
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine int2str(integ, string)
    !! Write a given integer to a string, using only as many digits as needed
    !
    implicit none
    integer :: integ
    character(len = 300) :: string
    !
    if ( integ < 10 ) then
      write(string, '(i1)') integ
    else if ( integ < 100 ) then
      write(string, '(i2)') integ
    else if ( integ < 1000 ) then
      write(string, '(i3)') integ
    else if ( integ < 10000 ) then
      write(string, '(i4)') integ
    endif
    !
    string = trim(string)
    !
    return
    !
  end subroutine int2str
  !
end module TMEModule
