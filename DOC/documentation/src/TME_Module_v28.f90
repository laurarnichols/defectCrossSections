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
  integer :: fftxMin
    !! x min for fft grid; not read in from PC input
  integer :: fftxMax
    !! x max for fft grid; not read in from PC input
  integer :: fftyMin
    !! y min for fft grid; not read in from PC input
  integer :: fftyMax
    !! y max for fft grid; not read in from PC input
  integer :: fftzMin
    !! z min for fft grid; not read in from PC input
  integer :: fftzMax
    !! z max for fft grid; not read in from PC input
  integer :: gx
  integer :: gy
  integer :: gz
  integer :: i
  integer :: iBandFfinal
  integer :: iBandFinit
  integer :: iBandIfinal
  integer :: iBandIinit
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
  integer :: nBands
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
  integer :: nSpins
  integer :: numOfGvecs
    !! Number of G vectors; not read in from PC input
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
  character(len = 200) :: exportDirSD
  character(len = 200) :: exportDirPC
  !
  !
  ! Declare matrix/vector integers
  integer, allocatable :: counts(:)
  integer, allocatable :: displmnt(:)
  integer, allocatable :: igvs(:,:,:)
  integer, allocatable :: iqs(:)
  integer, allocatable :: nFs(:,:)
  integer, allocatable :: ngs(:,:)
  integer, allocatable :: nIs(:,:)
  integer, allocatable :: nPWsI(:)
  integer, allocatable :: nPWsF(:)
  integer, allocatable :: pwGindPC(:)
  integer, allocatable :: pwGindSD(:)
  integer, allocatable :: pwGvecs(:,:)
  integer, allocatable :: pwGs(:,:)
  !
  ! Declare matrix/vector reals
  real(kind = dp) at(3,3)
  real(kind = dp) bg(3,3)
  real(kind = dp), allocatable :: absVfi2(:,:)
  real(kind = dp), allocatable :: DE(:,:)
  real(kind = dp), allocatable :: eigvF(:)
  real(kind = dp), allocatable :: eigvI(:)
  real(kind = dp), allocatable :: gvecs(:,:)
  !
  ! Declare matrix/vector complex numbers
  complex(kind = dp), allocatable :: betaPC(:,:)
  complex(kind = dp), allocatable :: betaSD(:,:)
  complex(kind = dp), allocatable :: cProjBetaPCPsiSD(:,:,:)
  complex(kind = dp), allocatable :: cProjBetaSDPhiPC(:,:,:)
  complex(kind = dp), allocatable :: cProjPC(:,:,:)
  complex(kind = dp), allocatable :: cProjSD(:,:,:)
  complex(kind = dp), allocatable :: paw_id(:,:)
  complex(kind = dp), allocatable :: paw_fi(:,:)
  complex(kind = dp), allocatable :: pawKPC(:,:,:)
  complex(kind = dp), allocatable :: paw_PsiPC(:,:)
  complex(kind = dp), allocatable :: pawPsiPC(:,:)
  complex(kind = dp), allocatable :: pawSDK(:,:,:)
  complex(kind = dp), allocatable :: paw_SDPhi(:,:)
  complex(kind = dp), allocatable :: pawSDPhi(:,:)
  complex(kind = dp), allocatable :: paw_SDKKPC(:,:)
  complex(kind = dp), allocatable :: Ufi(:,:,:)
  complex(kind = dp), allocatable :: wfcPC(:,:)
  complex(kind = dp), allocatable :: wfcSD(:,:)
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
    integer :: iRc
      !! Maximum radius of beta projector (outer radius to integrate);
      !! for PAW augmentation charge may extend a bit further
    integer :: numOfAtoms
      !! Number of atoms of a specific type in the structure
    integer :: lMax
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
    integer, allocatable :: lps(:)
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
  !> @todo Use crystal type instead of all of the explicit variables @endtodo
  type :: crystal
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
    integer, allocatable :: npws(:)
      !! Number of plane waves per k point
    integer, allocatable :: atomTypeIndex(:)
      !! Index of the given atom type
    !integer, allocatable :: groundState(:)
      ! Was read from input file but not used, so removed
    !
    real(kind = dp) :: omega
      !! Cell volume
    real(kind = dp), allocatable :: wk(:)
    real(kind = dp), allocatable :: xk(:, :)
    real(kind = dp), allocatable :: posIon(:,:)
    !
    character(len = 2) crystalType
      !! 'PC' for pristine crystal and 'SD' for solid defect
    character(len = 200) :: exportDir
      !! Export directory from [[pw_export_for_tme(program)]]
    !
    TYPE(atom), allocatable :: atoms(:)
    !
!    integer :: Jmax, maxL, iTypes, nn, nm
!    integer :: nBands, nSpins
!    integer :: i, j, n1, n2, n3, n4, n, id
!    !
!    real(kind = dp) :: at(3,3), bg(3,3)
!    !
!    real(kind = dp), allocatable :: eigvI(:), eigvF(:)
!    real(kind = dp), allocatable :: DE(:,:), absVfi2(:,:)
!    !
!    complex(kind = dp), allocatable :: wfc(:,:), wfcSD(:,:), Ufi(:,:,:)
!    complex(kind = dp), allocatable :: cProjPC(:,:,:), cProjSD(:,:,:)
!    !
!    integer, allocatable :: igvs(:,:,:), pwGvecs(:,:), iqs(:)
!    integer, allocatable :: pwGindI(:), pwGindF(:), pwGs(:,:), nIs(:,:), nFs(:,:), ngs(:,:)
!    integer, allocatable :: npwsPC(:)
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
  NAMELIST /TME_Input/ exportDirSD, exportDirPC, elementsPath, &
                       iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, &
                       ki, kf, calculateVfis, VfisOutput, eBin
                       !! Used to group the variables read in from the .in file
  !
  !
!=====================================================================================================
contains
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine initializeCalculation(solidDefect, pristineCrystal, elementsPath, VFisOutput, ki, kf, eBin, &
                                   iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, calculateVFis, t0)
    !! Initialize the calculation by starting timer,
    !! setting start values for variables to be read from
    !! `.in` file, removing any existing output in the output directory,
    !! and opening a clean output file
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer, intent(out) :: ki, kf, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !
    real(kind = dp), intent(out) :: eBin, t0
    !
    character(len = 200), intent(out) :: VfisOutput
    character(len = 300), intent(out) :: elementsPath
    !
    logical, intent(out) :: calculateVfis
    logical :: fileExists
      !! Whether or not the output file already exists
    TYPE(crystal), intent(inout) :: solidDefect, pristineCrystal
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
    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
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
  subroutine readInput()
    !! Delete any previous output, initialize input variables,
    !! start a timer, and read in the input files
    !!
    !! @todo Change `readInput()` to have arguments to make clear that these variables are getting changed @endtodo
    !!
    implicit none
    !
    READ (5, TME_Input, iostat = ios)
        !! * Read input from command line (or input file if use `< TME_Input.md`)
    solidDefect%exportDir = exportDirSD
    perfectCrystal%exportDir = exportDirPC
    !
    call checkInitialization()
        !! * Check that all required variables were input and have values that make sense
    !
    !> @todo Figure out what the difference in PC and SD is @endtodo
    write(6,'("Calling subroutine with crystal type ", a)') 'PC'
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
    logical :: file_exists
      !! Whether or not the exported directory from [[pw_export_for_TME(program)]]
      !! exists
    !>
    !> @todo Change `file_exists` to `fileExists` in `checkInitialization()` @endtodo
    logical:: abortExecution
    !
    abortExecution = .false.
      !! * Set the default value of abort execution so that the program
      !! will only abort if there is an issue with the inputs
    !
    write(iostd, '(" Inputs : ")')
      !! * Write out a header to the output file
    !
    !> * If the SD export directory variable is blank
    !>    * Output an error message and set `abortExecution` to true
    !> * Otherwise
    !>    * Check if the SD export directory exists
    !>    * If the SD export directory doesn't exist
    !>       * Output an error message and set `abortExecution` to true
    if ( trim(solidDefect%exportDir) == '' ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable: ""exportDirSD"" is not defined!")')
      write(iostd, '(" usage : exportDirSD = ''./Export/''")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    else
      !
      inquire(file= trim(solidDefect%exportDir), exist = file_exists)
      !
      if ( file_exists .eqv. .false. ) then
        !
        write(iostd, '(" exportDirSD :", a, " does not exist !")') trim(solidDefect%exportDir)
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
        !
      endif
      !
    endif
    !
    write(iostd, '("exportDirSD = ''", a, "''")') trim(solidDefect%exportDir)
      !! * Output the given SD export directory
    !
    !> * If the PC export directory variable is blank
    !>    * Output an error message and set `abortExecution` to true
    !> * Otherwise
    !>    * Check if the PC export directory exists
    !>    * If the PC export directory doesn't exist
    !>       * Output an error message and set `abortExecution` to true
    if ( trim(perfectCrystal%exportDir) == '' ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""exportDirPC"" is not defined!")')
      write(iostd, '(" usage : exportDirPC = ''./Export/''")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    else
      !
      inquire(file= trim(perfectCrystal%exportDir), exist = file_exists)
      !
      if ( file_exists .eqv. .false. ) then
        !
        write(iostd, '(" exportDirPC :", a, " does not exist !")') trim(perfectCrystal%exportDir)
        write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
        !
      endif
      !
    endif
    !
    write(iostd, '("exportDirPC = ''", a, "''")') trim(perfectCrystal%exportDir)
      !! * Output the given PC export directory
    !
    !> * If the elements path is blank
    !>    * Output a warning message and set the default value to `./`
    if ( trim(elementsPath) == '' ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""elementsPath"" is not defined!")')
      write(iostd, '(" usage : elementsPath = ''./''")')
      write(iostd, '(" The current directory will be used as elementsPath.")')
      elementsPath = './'
      !
    endif
    !
    inquire(file= trim(elementsPath), exist = file_exists)
      !! * Check if the elements path folder exists already
    !
    !> * If the elements path folder doesn't already exist
    !>    * Create the directory by writing the `mkdir` command to a string
    !>    * Then execute the command
    if ( .not.file_exists ) then
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
    !...............................................................................................
      !! * If `iBandIinit`, `iBandIfinal`, `iBandFinit`, or `iBandFfinal` is still less than zero
      !!    * Output an error message and set `abortExecution` to true
      !! * Then output each of their values
    !
    if ( iBandIinit < 0 ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandIinit"" is not defined!")')
      write(iostd, '(" usage : iBandIinit = 10")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    endif
    !
    write(iostd, '("iBandIinit = ", i4)') iBandIinit
    !
    if ( iBandIfinal < 0 ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandIfinal"" is not defined!")')
      write(iostd, '(" usage : iBandIfinal = 20")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    endif
    !
    write(iostd, '("iBandIfinal = ", i4)') iBandIfinal
    !
    if ( iBandFinit < 0 ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandFinit"" is not defined!")')
      write(iostd, '(" usage : iBandFinit = 9")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    endif
    !
    write(iostd, '("iBandFinit = ", i4)') iBandFinit
    !
    if ( iBandFfinal < 0 ) then
      !
      write(iostd, *)
      write(iostd, '(" Variable : ""iBandFfinal"" is not defined!")')
      write(iostd, '(" usage : iBandFfinal = 9")')
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    endif
    !
    write(iostd, '("iBandFfinal = ", i4)') iBandFfinal
    !
    !...............................................................................................
    !
    !> * If `calculateVfis` is true and `iBandFinit` and `iBandFfinal` are not equal
    !>    * Output an error message and set `abortExecution` to true
    if ( ( calculateVfis ) .and. ( iBandFinit /= iBandFfinal ) ) then
      !
      write(iostd, *)
      write(iostd, '(" Vfis can be calculated only if the final state is one and only one!")')
      write(iostd, '(" ''iBandFInit'' = ", i10)') iBandFinit
      write(iostd, '(" ''iBandFfinal'' = ", i10)') iBandFfinal
      write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
      abortExecution = .true.
      !
    endif
    !
    write(iostd, '("calculateVfis = ", l )') calculateVfis
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
    !> @todo Check if there is any kind of check on `ki` and `kf`. Why was this commented out? @endtodo
    !
    !if ( ki < 0 ) then
    !  write(iostd, *)
    !  write(iostd, '(" Initial k-point index is not defined!")')
    !  write(iostd, '(" Variable : ""ki"" is not defined!")')
    !  write(iostd, '(" usage : ki = 1")')
    !  !write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
    !  !abortExecution = .true.
    !  write(iostd, '(" ki = 1 will be used.")')
    !  ki = 1
    !endif
    !
    !if ( kf < 0 ) then
    !  write(iostd, *)
    !  write(iostd, '(" Final k-point index is not defined!")')
    !  write(iostd, '(" Variable : ""kf"" is not defined!")')
    !  write(iostd, '(" usage : kf = 1")')
    !  !write(iostd, '(" This variable is mandatory and thus the program will not be executed!")')
    !  write(iostd, '(" kf = total number of k-points will be used.")')
    !  !abortExecution = .true.
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
    !> * If the value of `eBin` is still less than zero
    !>    * Output a warning message and set the default value to 0.01 eV
    if ( eBin < 0.0_dp ) then
      !
      write(iostd,'(" Variable : ""eBin"" is not defined!")')
      write(iostd,'(" usage : eBin = 0.01")')
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
    !> * If `abortExecution` was ever set to true
    !>    * Output an error message and stop the program
    if ( abortExecution ) then
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
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine readQEExport(system)
    !! Read input files in the Export directory created by
    !! [[pw_export_for_tme(program)]]
    !!
    !! <h2>Walkthrough</h2>
    !!
    !! @todo Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed @endtodo
    !!
    !
    implicit none
    !
    !integer, intent(in) :: id
    !
    TYPE(crystal), intent(inout) :: system
    !
    integer :: i, j, ik, iType, ni
      !! Loop indices
    integer:: irc
      !! Local value of iRc for each atom so don't have to keep accessing in loop
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
    real(kind = dp) :: ef
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
    !> @todo Add information about these variables to top @endtodo
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
    !> @note
    !> Only read the number of G vectors from the solid
    !> defect crystal because??
    !> @endnote
    if ( system%crystalType == 'PC' ) then
      !
      read(50, * ) 
      !
    else if ( system%crystalType == 'SD' ) then
      !
      read(50, * ) numOfGvecs
      !
    endif
    !
    read(50, '(a)') textDum
    read(50, '(i10)') system%numOfPWs
    !
    read(50, '(a)') textDum     
    !
    !> @note
    !> Only read the fft grid from the solid defect crystal
    !> because??
    !> @endnote
    if ( system%crystalType == 'PC' ) then
      !
      read(50, * ) 
      !
    else if ( system%crystalType == 'SD' ) then
      !
      read(50, '(6i10)') fftxMin, fftxMax, fftyMin, fftyMax, fftzMin, fftzMax
      !
    endif
    !
    read(50, '(a)') textDum
    !
    if ( system%crystalType == 'PC' ) then
      !
      read(50,  * )
      read(50,  * )
      read(50,  * )
      !
    else if ( system%crystalType == 'SD' ) then
      !
      read(50, '(a5, 3ES24.15E3)') textDum, at(1:3,1)
      read(50, '(a5, 3ES24.15E3)') textDum, at(1:3,2)
      read(50, '(a5, 3ES24.15E3)') textDum, at(1:3,3)
      !
    endif
    !
    read(50, '(a)') textDum
    !
    if ( system%crystalType == 'PC' ) then
      !
      read(50,  * )
      read(50,  * )
      read(50,  * )
      !
    else if ( system%crystalType == 'SD' ) then
      !
      read(50, '(a5, 3ES24.15E3)') textDum, bg(1:3,1)
      read(50, '(a5, 3ES24.15E3)') textDum, bg(1:3,2)
      read(50, '(a5, 3ES24.15E3)') textDum, bg(1:3,3)
      !
    endif
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
      read(50,'(i10, 3ES24.15E3)') system%atomTypeIndex(ni), (system%posIon(j,ni) , j = 1,3)
        !! @todo Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:3,ni)` in `readQEExport()` for clarity @endtodo
      !
    enddo
    !
    read(50, '(a)') textDum
    !
    if ( system%crystalType == 'PC' ) then
      !
      read(50, * )
      !
    else if ( system%crystalType == 'SD' ) then
      !
      read(50, '(i10)') nBands
      !
    endif
    !
    read(50, '(a)') textDum
    !
    if ( system%crystalType == 'PC' ) then
      !
      read(50, * )
      !
    else if ( system%crystalType == 'SD' ) then
      !
      read(50, '(i10)') nSpins
      !
    endif
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
      read(50, '(i10)') system%atoms(iType)%lMax              ! number of projectors
      !
      allocate ( system%atoms(iType)%lps( system%atoms(iType)%lMax ) )
      !
      read(50, '(a)') textDum
      do i = 1, system%atoms(iType)%lMax 
        !
        read(50, '(2i10)') l, ind
        system%atoms(iType)%lps(ind) = l
        !
      enddo
      !
      read(50, '(a)') textDum
      read(50, '(i10)') system%atoms(iType)%lmMax
      !
      read(50, '(a)') textDum
      read(50, '(2i10)') system%atoms(iType)%nMax, system%atoms(iType)%iRc
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
      allocate ( system%atoms(iType)%wae(system%atoms(iType)%nMax, system%atoms(iType)%lMax) )
      allocate ( system%atoms(iType)%wps(system%atoms(iType)%nMax, system%atoms(iType)%lMax) )
      !
      read(50, '(a)') textDum
      do j = 1, system%atoms(iType)%lMax
        do i = 1, system%atoms(iType)%nMax
          !
          read(50, '(2ES24.15E3)') system%atoms(iType)%wae(i, j), system%atoms(iType)%wps(i, j) 
          ! write(iostd, '(2i5, ES24.15E3)') j, i, abs(system%atoms(iType)%wae(i, j)-system%atoms(iType)%wps(i, j))
          !
        enddo
      enddo
      !  
      allocate ( system%atoms(iType)%F( system%atoms(iType)%iRc, system%atoms(iType)%lMax ) ) !, system%atoms(iType)%lMax) )
      allocate ( system%atoms(iType)%F1(system%atoms(iType)%iRc, system%atoms(iType)%lMax, system%atoms(iType)%lMax ) )
      allocate ( system%atoms(iType)%F2(system%atoms(iType)%iRc, system%atoms(iType)%lMax, system%atoms(iType)%lMax ) )
      !
      system%atoms(iType)%F = 0.0_dp
      system%atoms(iType)%F1 = 0.0_dp
      system%atoms(iType)%F2 = 0.0_dp
      !
      !> * Calculate `F`, `F1`, and `F2` using the all-electron and psuedowvefunctions
      !> @todo Look more into how AE and PS wavefunctions are combined to further understand this @endtodo
      !> @todo Move this behavior to another subroutine for clarity @endtodo
      do j = 1, system%atoms(iType)%lMax
        !
        irc = system%atoms(iType)%iRc
        !
        system%atoms(iType)%F(1:irc,j)=(system%atoms(iType)%wae(1:irc,j)-system%atoms(iType)%wps(1:irc,j))* &
              system%atoms(iType)%r(1:irc)*system%atoms(iType)%rab(1:irc)
        !
        do i = 1, system%atoms(iType)%lMax
          !> @todo Figure out if differences in PC and SD `F1` calculations are intentional @endtodo
          !> @todo Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}` @endtodo
          if ( system%crystalType == 'PC' ) then
            !
            system%atoms(iType)%F1(1:irc,i,j) = ( system%atoms(iType)%wps(1:irc,i)*system%atoms(iType)%wae(1:irc,j) - &
                                                  system%atoms(iType)%wps(1:irc,i)*system%atoms(iType)%wps(1:irc,j))* &
                                                  system%atoms(iType)%rab(1:irc)
            !
          else if ( system%crystalType == 'SD' ) then
            !
            system%atoms(iType)%F1(1:irc,i,j) = ( system%atoms(iType)%wae(1:irc,i)*system%atoms(iType)%wps(1:irc,j) - &
                                                  system%atoms(iType)%wps(1:irc,i)*system%atoms(iType)%wps(1:irc,j))* & 
                                                  system%atoms(iType)%rab(1:irc)
            !
          endif
          !
          system%atoms(iType)%F2(1:irc,i,j) = ( system%atoms(iType)%wae(1:irc,i)*system%atoms(iType)%wae(1:irc,j) - &
                                           system%atoms(iType)%wae(1:irc,i)*system%atoms(iType)%wps(1:irc,j) - &
                                           system%atoms(iType)%wps(1:irc,i)*system%atoms(iType)%wae(1:irc,j) + &
                                           system%atoms(iType)%wps(1:irc,i)*system%atoms(iType)%wps(1:irc,j))* & 
                                           system%atoms(iType)%rab(1:irc)

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
    !> * Go through the `lps` values for each projector for each atom
    !> and find the max to store in `JMAX`
    !> @todo Figure out if intentional to only use `JMAX` from SD input @endtodo
    JMAX = 0
    do iType = 1, system%numOfTypes
      !
      do i = 1, system%atoms(iType)%lMax
        !
        if ( system%atoms(iType)%lps(i) > JMAX ) JMAX = system%atoms(iType)%lps(i)
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
      allocate ( system%atoms(iType)%bes_J_qr( 0:JMAX, system%atoms(iType)%iRc ) )
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
  subroutine distributePWsToProcs(nOfPWs, nOfBlocks)
    !
    implicit none
    !
    integer, intent(in)  :: nOfPWs, nOfBlocks
    !
    integer :: iStep, iModu
    !
    iStep = int(nOfPWs/nOfBlocks)
    iModu = mod(nOfPWs,nOfBlocks)
    !
    do i = 0, nOfBlocks - 1
      counts(i) = iStep
      if ( iModu > 0 ) then
        counts(i) = counts(i) + 1
        iModu = iModu - 1
      endif
    enddo
    !
    displmnt(0) = 0
    do i = 1, nOfBlocks-1
      displmnt(i) = displmnt(i-1) + counts(i)
    enddo
    !
    return
    !
  end subroutine distributePWsToProcs
  !
  !
  subroutine int2str(integ, string)
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
  !
  subroutine finalizeCalculation()
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
  !
  subroutine readPWsSet()
    !
    implicit none
    !
    integer :: ig, iDum, iGx, iGy, iGz
    !
    open(72, file=trim(solidDefect%exportDir)//"/mgrid")
    !
    read(72, * )
    read(72, * )
    !
    allocate ( gvecs(3, numOfGvecs ) )
    !
    gvecs(:,:) = 0.0_dp
    !
    do ig = 1, numOfGvecs
      read(72, '(4i10)') iDum, iGx, iGy, iGz
      gvecs(1,ig) = dble(iGx)*bg(1,1) + dble(iGy)*bg(1,2) + dble(iGz)*bg(1,3)
      gvecs(2,ig) = dble(iGx)*bg(2,1) + dble(iGy)*bg(2,2) + dble(iGz)*bg(2,3)
      gvecs(3,ig) = dble(iGx)*bg(3,1) + dble(iGy)*bg(3,2) + dble(iGz)*bg(3,3)
    enddo
    !
    close(72)
    !
    return
    !
  end subroutine readPWsSet
  !
  !
  subroutine readWfcPC(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: ib, ig, iDumV(3)
    !
    complex(kind = dp) :: wfc
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    open(72, file=trim(perfectCrystal%exportDir)//"/grid."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    allocate ( pwGindPC(perfectCrystal%npws(ik)) )
    !
    do ig = 1, perfectCrystal%npws(ik)
      read(72, '(4i10)') pwGindPC(ig), iDumV(1:3)
    enddo
    !
    close(72)
    !
    open(72, file=trim(perfectCrystal%exportDir)//"/wfc."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    do ib = 1, iBandIinit - 1
      do ig = 1, perfectCrystal%npws(ik)
        read(72, *)
      enddo
    enddo
    !
    wfcPC(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)
    !
    do ib = iBandIinit, iBandIfinal
      do ig = 1, perfectCrystal%npws(ik)
        read(72, '(2ES24.15E3)') wfc
        wfcPC(pwGindPC(ig), ib) = wfc
      enddo
    enddo
    !
    close(72)
    !
    deallocate ( pwGindPC )
    !
    return
    !
  end subroutine readWfcPC
  !
  !
  subroutine projectBetaPCwfcSD(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: ig, iDumV(3)
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    ! Reading PC projectors
    !
    open(72, file=trim(perfectCrystal%exportDir)//"/grid."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    allocate ( pwGindPC(perfectCrystal%npws(ik)) )
    !
    do ig = 1, perfectCrystal%npws(ik)
      read(72, '(4i10)') pwGindPC(ig), iDumV(1:3)
    enddo
    !
    close(72)
    !
    allocate ( betaPC(numOfPWs, perfectCrystal%nProjs) )
    !
    betaPC(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    open(73, file=trim(perfectCrystal%exportDir)//"/projectors."//trim(iks))
    !
    read(73, '(a)') textDum
    read(73, '(2i10)') perfectCrystal%nProjs, npw
    !
    do j = 1, perfectCrystal%nProjs
      do i = 1, npw
        read(73,'(2ES24.15E3)') betaPC(pwGindPC(i),j)
      enddo
    enddo
    !
    close(73)
    !
    deallocate ( pwGindPC )
    !
    do j = iBandFinit, iBandFfinal
      do i = 1, perfectCrystal%nProjs
        cProjBetaPCPsiSD(i,j,1) = sum(conjg(betaPC(:,i))*wfcSD(:,j))
        !write(65,'(2f17.12)') cProjPC(i,j,1) - cProjBetaPCPsiSD(i,j,1)
      enddo
    enddo
    !
    !close(65)
    !
    deallocate ( betaPC )
    !
    return
    !
  end subroutine projectBetaPCwfcSD
  !
  !
  subroutine readWfcSD(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: ib, ig, iDumV(3)
    !
    complex(kind = dp) :: wfc
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    open(72, file=trim(solidDefect%exportDir)//"/grid."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    allocate ( pwGindSD(solidDefect%npws(ik)) )
    !
    do ig = 1, solidDefect%npws(ik)
      read(72, '(4i10)') pwGindSD(ig), iDumV(1:3)
    enddo
    !
    close(72)
    !
    open(72, file=trim(solidDefect%exportDir)//"/wfc."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    do ib = 1, iBandFinit - 1
      do ig = 1, solidDefect%npws(ik)
        read(72, *)
      enddo
    enddo
    !
    wfcSD(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)
    !
    do ib = iBandFinit, iBandFfinal
      do ig = 1, solidDefect%npws(ik)
        read(72, '(2ES24.15E3)') wfc
        wfcSD(pwGindSD(ig), ib) = wfc
      enddo
    enddo
    !
    close(72)
    !
    deallocate ( pwGindSD )
    !
    return
    !
  end subroutine readWfcSD
  !
  !
  subroutine calculatePWsOverlap(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: ibi, ibf
    !
    call readWfcPC(ik)
    !
    call readWfcSD(ik)
    !
    Ufi(:,:,ik) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    do ibi = iBandIinit, iBandIfinal 
      !
      do ibf = iBandFinit, iBandFfinal
        Ufi(ibf, ibi, ik) = sum(conjg(wfcSD(:,ibf))*wfcPC(:,ibi))
        !if ( ibi == ibf ) write(iostd,'(2i4,3ES24.15E3)') ibf, ibi, Ufi(ibf, ibi, ik), abs(Ufi(ibf, ibi, ik))**2
        flush(iostd)
      enddo
      !
    enddo
    !
    return
    !
  end subroutine calculatePWsOverlap
  !
  !
  subroutine readProjectionsPC(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: i, j
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    cProjPC(:,:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
    !
    ! Reading projections
    !
    open(72, file=trim(perfectCrystal%exportDir)//"/projections."//trim(iks))
    !
    read(72, *)
    !
    do j = 1, nBands  ! number of bands 
      do i = 1, perfectCrystal%nProjs ! number of projections
        read(72,'(2ES24.15E3)') cProjPC(i,j,1)
      enddo
    enddo
    !
    close(72)
    !
    return
    !
  end subroutine readProjectionsPC
  !
  !
  subroutine readProjectionsSD(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: i, j
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    cProjSD(:,:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
    !
    ! Reading projections
    !
    open(72, file=trim(solidDefect%exportDir)//"/projections."//trim(iks))
    !
    read(72, *)
    !
    do j = 1, nBands  ! number of bands 
      do i = 1, solidDefect%nProjs ! number of projections
        read(72,'(2ES24.15E3)') cProjSD(i,j,1)
      enddo
    enddo
    !
    close(72)
    !
    return
    !
  end subroutine readProjectionsSD
  !
  !
  subroutine projectBetaSDwfcPC(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: ig, iDumV(3)
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    ! Reading SD projectors
    !
    open(72, file=trim(solidDefect%exportDir)//"/grid."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    allocate ( pwGindSD(solidDefect%npws(ik)) )
    !
    do ig = 1, solidDefect%npws(ik)
      read(72, '(4i10)') pwGindSD(ig), iDumV(1:3)
    enddo
    !
    close(72)
    !
    allocate ( betaSD(numOfPWs, solidDefect%nProjs) )
    !
    betaSD(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    open(73, file=trim(solidDefect%exportDir)//"/projectors."//trim(iks))
    !
    read(73, '(a)') textDum
    read(73, '(2i10)') solidDefect%nProjs, npw
    !
    do j = 1, solidDefect%nProjs
      do i = 1, npw
        read(73,'(2ES24.15E3)') betaSD(pwGindSD(i),j)
      enddo
    enddo
    !
    close(73)
    !
    deallocate ( pwGindSD )
    !
    do j = iBandIinit, iBandIfinal
      do i = 1, solidDefect%nProjs
        cProjBetaSDPhiPC(i,j,1) = sum(conjg(betaSD(:,i))*wfcPC(:,j))
        !write(66,'(2f17.12)') cProjSD(i,j,1) - cProjBetaSDPhiPC(i,j,1)
      enddo
    enddo
    !
    !close(66)
    !
    deallocate ( betaSD )
    !
    return
    !
  end subroutine projectBetaSDwfcPC
  !
  !
  subroutine pawCorrectionKPC()
    !
    implicit none
    !
    !integer, intent(in) :: ik
    !
    integer :: ibi, ibf, ispin, ig
    integer :: LL, I, NI, LMBASE, LM
    integer :: L, M, ind, iT
    real(kind = dp) :: q, qDotR, FI, t1, t2
    !
    real(kind = dp) :: JL(0:JMAX), v_in(3)
    complex(kind = dp) :: Y( (JMAX+1)**2 )
    complex(kind = dp) :: VifQ_aug, ATOMIC_CENTER
    !
    ispin = 1
    !
    call cpu_time(t1)
    !
    pawKPC(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    do ig = nPWsI(myid), nPWsF(myid) ! 1, numOfGvecs
      !
      if ( myid == root ) then 
        if ( (ig == nPWsI(myid) + 1000) .or. (mod(ig, 25000) == 0) .or. (ig == nPWsF(myid)) ) then
          call cpu_time(t2)
          write(iostd, '("        Done ", i10, " of", i10, " k-vecs. ETR : ", f10.2, " secs.")') &
                ig, nPWsF(myid) - nPWsI(myid) + 1, (t2-t1)*(nPWsF(myid) - nPWsI(myid) + 1 -ig )/ig
          flush(iostd)
          !call cpu_time(t1)
        endif
      endif
      !
      q = sqrt(sum(gvecs(:,ig)*gvecs(:,ig)))
      !
      v_in(:) = gvecs(:,ig)
      if ( abs(q) > 1.0e-6_dp ) v_in = v_in/q ! i have to determine v_in = q
      Y = cmplx(0.0_dp, 0.0_dp, kind = dp)
      CALL ylm(v_in, JMAX, Y) ! calculates all the needed spherical harmonics once
      !
      LMBASE = 0
      !
      do iT = 1, perfectCrystal%numOfTypes
        !
        DO I = 1, perfectCrystal%atoms(iT)%iRc ! nMax - 1
          !
          JL = 0.0_dp
          CALL bessel_j(q*solidDefect%atoms(iT)%r(I), JMAX, JL) ! returns the spherical bessel at qr point
          perfectCrystal%atoms(iT)%bes_J_qr(:,I) = JL(:)
          !
        ENDDO
        !
      enddo
      !
      do ni = 1, perfectCrystal%nIons ! LOOP OVER THE IONS
        !
        qDotR = sum(gvecs(:,ig)*perfectCrystal%posIon(:,ni))
        !
        ATOMIC_CENTER = exp( -ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        !
        iT = perfectCrystal%atomTypeIndex(ni)
        LM = 0
        DO LL = 1, perfectCrystal%atoms(iT)%lMax
          L = perfectCrystal%atoms(iT)%LPS(LL)
          DO M = -L, L
            LM = LM + 1 !1st index for CPROJ
            !
            FI = 0.0_dp
            !
            FI = sum(perfectCrystal%atoms(iT)%bes_J_qr(L,:)*perfectCrystal%atoms(iT)%F(:,LL)) ! radial part integration F contains rab
            !
            ind = L*(L + 1) + M + 1 ! index for spherical harmonics
            VifQ_aug = ATOMIC_CENTER*Y(ind)*(-II)**L*FI
            !
            do ibi = iBandIinit, iBandIfinal
              !
              do ibf = iBandFinit, iBandFfinal
                !
                pawKPC(ibf, ibi, ig) = pawKPC(ibf, ibi, ig) + VifQ_aug*cProjPC(LM + LMBASE, ibi, ISPIN)
                !
              enddo
              !
            enddo
            !
          ENDDO
        ENDDO
        LMBASE = LMBASE + perfectCrystal%atoms(iT)%lmMax
      ENDDO
      !
    enddo
    !
    !pawKPC(:,:,:) = pawKPC(:,:,:)*4.0_dp*pi/sqrt(solidDefect%omega)
    !
    return
    !
  end subroutine pawCorrectionKPC
  !
  !
  subroutine pawCorrectionSDK()
    !
    implicit none
    !
    !integer, intent(in) :: ik
    !
    integer :: ibi, ibf, ispin, ig
    integer :: LL, I, NI, LMBASE, LM
    integer :: L, M, ind, iT
    real(kind = dp) :: q, qDotR, FI, t1, t2
    !
    real(kind = dp) :: JL(0:JMAX), v_in(3)
    complex(kind = dp) :: Y( (JMAX+1)**2 )
    complex(kind = dp) :: VifQ_aug, ATOMIC_CENTER
    !
    ispin = 1
    !
    call cpu_time(t1)
    !
    pawSDK(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    do ig = nPWsI(myid), nPWsF(myid) ! 1, numOfGvecs
      !
      if ( myid == root ) then
        if ( (ig == nPWsI(myid) + 1000) .or. (mod(ig, 25000) == 0) .or. (ig == nPWsF(myid)) ) then
          call cpu_time(t2)
          write(iostd, '("        Done ", i10, " of", i10, " k-vecs. ETR : ", f10.2, " secs.")') &
                ig, nPWsF(myid) - nPWsI(myid) + 1, (t2-t1)*(nPWsF(myid) - nPWsI(myid) + 1 -ig )/ig
          flush(iostd)
          call cpu_time(t1)
        endif
      endif
      q = sqrt(sum(gvecs(:,ig)*gvecs(:,ig)))
      !
      v_in(:) = gvecs(:,ig)
      if ( abs(q) > 1.0e-6_dp ) v_in = v_in/q ! i have to determine v_in = q
      Y = cmplx(0.0_dp, 0.0_dp, kind = dp)
      CALL ylm(v_in, JMAX, Y) ! calculates all the needed spherical harmonics once
      !
      LMBASE = 0
      !
      do iT = 1, solidDefect%numOfTypes
        !
        DO I = 1, solidDefect%atoms(iT)%iRc ! nMax - 1
          !
          JL = 0.0_dp
          CALL bessel_j(q*solidDefect%atoms(iT)%r(I), JMAX, JL) ! returns the spherical bessel at qr point
          solidDefect%atoms(iT)%bes_J_qr(:,I) = JL(:)
          !
        ENDDO
      enddo
      !
      do ni = 1, solidDefect%nIons ! LOOP OVER THE IONS
        !
        qDotR = sum(gvecs(:,ig)*solidDefect%posIon(:,ni))
        !
        ATOMIC_CENTER = exp( ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        !
        iT = solidDefect%atomTypeIndex(ni)
        LM = 0
        DO LL = 1, solidDefect%atoms(iT)%lMax
          L = solidDefect%atoms(iT)%LPS(LL)
          DO M = -L, L
            LM = LM + 1 !1st index for CPROJ
            !
            FI = 0.0_dp
            !
            FI = sum(solidDefect%atoms(iT)%bes_J_qr(L,:)*solidDefect%atoms(iT)%F(:,LL)) ! radial part integration F contains rab
            !
            ind = L*(L + 1) + M + 1 ! index for spherical harmonics
            VifQ_aug = ATOMIC_CENTER*conjg(Y(ind))*(II)**L*FI
            !
            do ibi = iBandIinit, iBandIfinal
              !
              do ibf = iBandFinit, iBandFfinal
                !
                pawSDK(ibf, ibi, ig) = pawSDK(ibf, ibi, ig) + VifQ_aug*conjg(cProjSD(LM + LMBASE, ibf, ISPIN))
                !
              enddo
              !
            enddo
            !
          ENDDO
        ENDDO
        LMBASE = LMBASE + solidDefect%atoms(iT)%lmMax
      ENDDO
      !
    enddo
    !
    !pawSDK(:,:,:) = pawSDK(:,:,:)*4.0_dp*pi/sqrt(solidDefect%omega)
    !
    return
    !
  end subroutine pawCorrectionSDK
  !
  !
  subroutine pawCorrectionPsiPC()
    !
    ! calculates the augmentation part of the transition matrix element
    !
    implicit none
    integer :: ibi, ibf, niPC, ispin 
    integer :: LL, LLP, LMBASE, LM, LMP
    integer :: L, M, LP, MP, iT
    real(kind = dp) :: atomicOverlap
    !
    complex(kind = dp) :: cProjIe, cProjFe
    !
    ispin = 1
    !
    paw_PsiPC(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    LMBASE = 0
    !
    do niPC = 1, perfectCrystal%nIons ! LOOP OVER THE IONS
      !
      iT = perfectCrystal%atomTypeIndex(niPC)
      LM = 0
      DO LL = 1, perfectCrystal%atoms(iT)%lMax
        L = perfectCrystal%atoms(iT)%LPS(LL)
        DO M = -L, L
          LM = LM + 1 !1st index for CPROJ
          !
          LMP = 0
          DO LLP = 1, perfectCrystal%atoms(iT)%lMax
            LP = perfectCrystal%atoms(iT)%LPS(LLP)
            DO MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              !
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) then 
                atomicOverlap = sum(perfectCrystal%atoms(iT)%F1(:,LL, LLP))
                !
                do ibi = iBandIinit, iBandIfinal
                  cProjIe = cProjPC(LMP + LMBASE, ibi, ISPIN)
                  !
                  do ibf = iBandFinit, iBandFfinal
                    cProjFe = conjg(cProjBetaPCPsiSD(LM + LMBASE, ibf, ISPIN))
                    !
                    paw_PsiPC(ibf, ibi) = paw_PsiPC(ibf, ibi) + cProjFe*atomicOverlap*cProjIe
                    !write(iostd,*) LL, LLP, L, M, atomicOverlap
                    flush(iostd)
                    !
                  enddo
                  !
                enddo
                !
              endif
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      LMBASE = LMBASE + perfectCrystal%atoms(iT)%lmMax
    ENDDO
    !
    return
    !
  end subroutine pawCorrectionPsiPC
  !
  !
  subroutine pawCorrectionSDPhi()
    !
    ! calculates the augmentation part of the transition matrix element
    !
    implicit none
    integer :: ibi, ibf, ni, ispin 
    integer :: LL, LLP, LMBASE, LM, LMP
    integer :: L, M, LP, MP, iT
    real(kind = dp) :: atomicOverlap
    !
    complex(kind = dp) :: cProjIe, cProjFe
    !
    ispin = 1
    !
    paw_SDPhi(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    LMBASE = 0
    !
    do ni = 1, solidDefect%nIons ! LOOP OVER THE IONS
      !
      iT = solidDefect%atomTypeIndex(ni)
      LM = 0
      DO LL = 1, solidDefect%atoms(iT)%lMax
        L = solidDefect%atoms(iT)%LPS(LL)
        DO M = -L, L
          LM = LM + 1 !1st index for CPROJ
          !
          LMP = 0
          DO LLP = 1, solidDefect%atoms(iT)%lMax
            LP = solidDefect%atoms(iT)%LPS(LLP)
            DO MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              !
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) then
                atomicOverlap = sum(solidDefect%atoms(iT)%F1(:,LL,LLP))
                !
                do ibi = iBandIinit, iBandIfinal
                  cProjIe = cProjBetaSDPhiPC(LMP + LMBASE, ibi, ISPIN)
                  !
                  do ibf = iBandFinit, iBandFfinal
                    cProjFe = conjg(cProjSD(LM + LMBASE, ibf, ISPIN))
                    !
                    paw_SDPhi(ibf, ibi) = paw_SDPhi(ibf, ibi) + cProjFe*atomicOverlap*cProjIe
                    !
                  enddo
                  !
                enddo
                !
              endif
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      LMBASE = LMBASE + solidDefect%atoms(iT)%lmMax
    ENDDO
    !
    return
    !
  end subroutine pawCorrectionSDPhi
  !
  !
  subroutine pawCorrection()
    !
    ! calculates the augmentation part of the transition matrix element
    !
    implicit none
    integer :: ibi, ibf, niPC, ispin 
    integer :: LL, LLP, LMBASE, LM, LMP
    integer :: L, M, LP, MP, iT
    real(kind = dp) :: atomicOverlap
    !
    !real(kind = dp), allocatable :: Qij(:,:)
    !
    complex(kind = dp) :: cProjIe, cProjFe
    !
    ispin = 1
    !
    !open(52, file=trim(perfectCrystal%exportDir)//"/Qij")
    !read(52,*)
    !allocate ( Qij (8,8) )
    !do LL = 1, 8
    !  do LLP = 1, 8
    !    read(52,'(3i3,ES24.15E3)') LMBASE, LMBASE, LMBASE, Qij(LL, LLP)
    !  enddo
    !enddo
    !close(52)
    !
    paw_fi(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    !
    LMBASE = 0
    !
    do niPC = 1, perfectCrystal%nIons ! LOOP OVER THE IONS
      !
      iT = perfectCrystal%atomTypeIndex(niPC)
      LM = 0
      DO LL = 1, perfectCrystal%atoms(iT)%lMax
        L = perfectCrystal%atoms(iT)%LPS(LL)
        DO M = -L, L
          LM = LM + 1 !1st index for CPROJ
          !
          LMP = 0
          DO LLP = 1, perfectCrystal%atoms(iT)%lMax
            LP = perfectCrystal%atoms(iT)%LPS(LLP)
            DO MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              !
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) atomicOverlap = sum(perfectCrystal%atoms(iT)%F2(:,LL,LLP))
              !
              do ibi = iBandIinit, iBandIfinal
                cProjIe = cProjPC(LMP + LMBASE, ibi, ISPIN)
                !
                do ibf = iBandFinit, iBandFfinal
                  cProjFe = conjg(cProjPC(LM + LMBASE, ibf, ISPIN))
                  !
                  paw_fi(ibf, ibi) = paw_fi(ibf, ibi) + cProjFe*atomicOverlap*cProjIe
                  !
                enddo
                !
              enddo
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      LMBASE = LMBASE + perfectCrystal%atoms(iT)%lmMax
    ENDDO
    !
    return
    !
  end subroutine pawCorrection
  !
  !
  subroutine readEigenvalues(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    integer :: ib
    !
    character(len = 300) :: iks
    !
    call int2str(ik, iks)
    !
    open(72, file=trim(solidDefect%exportDir)//"/eigenvalues."//trim(iks))
    !
    read(72, * )
    read(72, * )
    !
    do ib = 1, iBandIinit - 1
      read(72, *)
    enddo
    !
    do ib = iBandIinit, iBandIfinal
      read(72, '(ES24.15E3)') eigvI(ib)
    enddo
    !
    close(72)
    !
    open(72, file=trim(solidDefect%exportDir)//"/eigenvalues."//trim(iks))
    !
    read(72, * )
    read(72, * ) 
    !
    do ib = 1, iBandFinit - 1
      read(72, *)
    enddo
    !
    do ib = iBandFinit, iBandFfinal
      read(72, '(ES24.15E3)') eigvF(ib)
    enddo
    !
    close(72)
    !
    return
    !
  end subroutine readEigenvalues
  !
  !
  subroutine calculateVfiElements()
    !
    implicit none
    !
    integer :: ik, ib, nOfEnergies, iE
    !
    real(kind = dp) :: eMin, eMax, E, av, sd, x, EiMinusEf, A, DHifMin
    !
    real(kind = dp), allocatable :: sumWk(:), sAbsVfiOfE2(:), absVfiOfE2(:)
    integer, allocatable :: nKsInEbin(:)
    !
    character (len = 300) :: text
    !
    allocate( DE(iBandIinit:iBandIfinal, perfectCrystal%nKpts), absVfi2(iBandIinit:iBandIfinal, perfectCrystal%nKpts) )
    ! 
    DE(:,:) = 0.0_dp
    absVfi2(:,:) = 0.0_dp 
    !
    do ik = 1, perfectCrystal%nKpts
      !
      eigvI(:) = 0.0_dp
      eigvF(:) = 0.0_dp
      !
      call readEigenvalues(ik)
      !
      do ib = iBandIinit, iBandIfinal
        !
        EiMinusEf = eigvI(ib) - eigvF(iBandFinit)
        absVfi2(ib,ik) = EiMinusEf**2*( abs(Ufi(iBandFinit,ib,ik))**2 - abs(Ufi(iBandFinit,ib,ik))**4 )
        !
        DE(ib, ik) = sqrt(EiMinusEf**2 - 4.0_dp*absVfi2(ib,ik))
        !
      enddo
      !
    enddo
    !
    eMin = minval( DE(:,:) )
    eMax = maxval( DE(:,:) )
    !
    nOfEnergies = int((eMax-eMin)/eBin) + 1
    !
    allocate ( absVfiOfE2(0:nOfEnergies), nKsInEbin(0:nOfEnergies), sumWk(0:nOfEnergies) )
    !
    absVfiOfE2(:) = 0.0_dp
    nKsInEbin(:) = 0
    sumWk(:) = 0.0_dp
    !
    do ik = 1, perfectCrystal%nKpts
      !
      do ib = iBandIinit, iBandIfinal
        !
        if ( abs( eMin - DE(ib,ik)) < 1.0e-3_dp ) DHifMin = absVfi2(ib, ik)
        iE = int((DE(ib, ik)-eMin)/eBin)
        if ( absVfi2(ib, ik) > 0.0_dp ) then
          absVfiOfE2(iE) = absVfiOfE2(iE) + perfectCrystal%wk(ik)*absVfi2(ib, ik)
          sumWk(iE) = sumWk(iE) + perfectCrystal%wk(ik)
          nKsInEbin(iE) = nKsInEbin(iE) + 1
        else
          write(iostd,*) 'lalala', absVfi2(ib, ik)
        endif
        !
      enddo
      !
    enddo
    !
    allocate ( sAbsVfiOfE2(0:nOfEnergies) )
    !
    sAbsVfiOfE2 = 0.0_dp
    !
    open(11, file=trim(VfisOutput)//'ofKpt', status='unknown')
    !
    write(11, '("# |<f|V|i>|^2 versus energy for all the k-points.")')
    write(text, '("# Energy (eV) shifted by half eBin, |<f|V|i>|^2 (Hartree)^2,")')
    write(11, '(a, " k-point index. Format : ''(2ES24.15E3,i10)''")') trim(text)
    !
    do ik = 1, perfectCrystal%nKpts
      !
      do ib = iBandIinit, iBandIfinal
        !
        iE = int((DE(ib,ik)-eMin)/eBin)
        av = absVfiOfE2(iE)/sumWk(iE)
        x = absVfi2(ib,ik)
        write(11, '(2ES24.15E3,i10)') (eMin + (iE+0.5_dp)*eBin)*HartreeToEv, x, ik
        write(12, '(2ES24.15E3,i10)') DE(ib,ik)*HartreeToEv, absVfi2(ib, ik), ik
        !write(11, '(2ES24.15E3,i10)') (eMin + iE*eBin + eBin/2.0_dp), x, ik
        sAbsVfiOfE2(iE) = sAbsVfiOfE2(iE) + perfectCrystal%wk(ik)*(x - av)**2/sumWk(iE)
        !
      enddo
      !
    enddo
    !
    close(11)
    !
    open(63, file=trim(VfisOutput), status='unknown')
    !
    write(63, '("# Averaged |<f|V|i>|^2 over K-points versus energy.")')
    write(63, '("#                 Cell volume : ", ES24.15E3, " (a.u.)^3,   Format : ''(ES24.15E3)''")') solidDefect%omega
    write(63, '("#   Minimun transition energy : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eMin
    write(63, '("# |DHif|^2 at minimum Tr. En. : ", ES24.15E3, " (Hartree^2),Format : ''(ES24.15E3)''")') DHifMin
    write(63, '("#                  Energy bin : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eBin
    write(text, '("# Energy (Hartree), averaged |<f|V|i>|^2 over K-points (Hartree)^2,")')
    write(63, '(a, " standard deviation (Hartree)^2. Format : ''(3ES24.15E3)''")') trim(text)
    !
    do iE = 0, nOfEnergies
      E = iE*eBin
      av = 0.0_dp
      sd = 0.0_dp
      if (nKsInEbin(iE) > 0) then
        av = absVfiOfE2(iE)/sumWk(iE)
        sd = sqrt(sAbsVfiOfE2(iE))
      endif
      write(63,'(3ES24.15E3)') eMin + E, av, sd
    enddo
    !
    close(63)
    !
    return
    !
  end subroutine calculateVfiElements
  !
  !
  subroutine checkIfCalculated(ik, tmes_file_exists)
    !
    implicit none
    !
    integer, intent(in) :: ik
    logical, intent(out) :: tmes_file_exists
    !
    character(len = 300) :: Uelements
    !
    if ( ik < 10 ) then
      write(Uelements, '("/TMEs_kptI_",i1,"_kptF_",i1)') ik, ik
    else if ( ik < 100 ) then
      write(Uelements, '("/TMEs_kptI_",i2,"_kptF_",i2)') ik, ik
    else if ( ik < 1000 ) then
      write(Uelements, '("/TMEs_kptI_",i3,"_kptF_",i3)') ik, ik
    else if ( ik < 10000 ) then
      write(Uelements, '("/TMEs_kptI_",i4,"_kptF_",i4)') ik, ik
    else if ( ik < 10000 ) then
      write(Uelements, '("/TMEs_kptI_",i5,"_kptF_",i5)') ik, ik
    endif
    !
    inquire(file = trim(elementsPath)//trim(Uelements), exist = tmes_file_exists)
    !
    return
    !
  end subroutine checkIfCalculated
  !
  !
  subroutine readUfis(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    !
    integer :: ibi, ibf, totalNumberOfElements, iDum, i
    real(kind = dp) :: rDum, t1, t2
    complex(kind = dp):: cUfi
    !
    character(len = 300) :: Uelements
    !
    call cpu_time(t1)
    write(iostd, '(" Reading Ufi(:,:) of k-point: ", i4)') ik
    !
    if ( ik < 10 ) then
      write(Uelements, '("/TMEs_kptI_",i1,"_kptF_",i1)') ik, ik
    else if ( ik < 100 ) then
      write(Uelements, '("/TMEs_kptI_",i2,"_kptF_",i2)') ik, ik
    else if ( ik < 1000 ) then
      write(Uelements, '("/TMEs_kptI_",i3,"_kptF_",i3)') ik, ik
    else if ( ik < 10000 ) then
      write(Uelements, '("/TMEs_kptI_",i4,"_kptF_",i4)') ik, ik
    else if ( ik < 10000 ) then
      write(Uelements, '("/TMEs_kptI_",i5,"_kptF_",i5)') ik, ik
    endif
    !
    open(17, file=trim(elementsPath)//trim(Uelements), status='unknown')
    !
    read(17, *) 
    read(17, *) 
    read(17,'(5i10)') totalNumberOfElements, iDum, iDum, iDum, iDum
    read(17, *) 
    !
    do i = 1, totalNumberOfElements
      !
      read(17, 1001) ibf, ibi, rDum, cUfi, rDum
      Ufi(ibf,ibi,ik) = cUfi
      !    
    enddo
    !
    close(17)
    !
    call cpu_time(t2)
    write(iostd, '(" Reading Ufi(:,:) done in:                   ", f10.2, " secs.")') t2-t1
    !
 1001 format(2i10,4ES24.15E3)
    !
    return
    !
  end subroutine readUfis
  !
  !
  subroutine writeResults(ik)
    !
    implicit none
    !
    integer, intent(in) :: ik
    !
    integer :: ibi, ibf, totalNumberOfElements
    real(kind = dp) :: t1, t2
    !
    character(len = 300) :: text, Uelements
    !
    call cpu_time(t1)
    !
    call readEigenvalues(ik)
    !
    write(iostd, '(" Writing Ufi(:,:).")')
    !
    if ( ik < 10 ) then
      write(Uelements, '("/TMEs_kptI_",i1,"_kptF_",i1)') ik, ik
    else if ( ik < 100 ) then
      write(Uelements, '("/TMEs_kptI_",i2,"_kptF_",i2)') ik, ik
    else if ( ik < 1000 ) then
      write(Uelements, '("/TMEs_kptI_",i3,"_kptF_",i3)') ik, ik
    else if ( ik < 10000 ) then
      write(Uelements, '("/TMEs_kptI_",i4,"_kptF_",i4)') ik, ik
    else if ( ik < 10000 ) then
      write(Uelements, '("/TMEs_kptI_",i5,"_kptF_",i5)') ik, ik
    endif
    !
    open(17, file=trim(elementsPath)//trim(Uelements), status='unknown')
    !
    write(17, '("# Cell volume (a.u.)^3. Format: ''(a51, ES24.15E3)'' ", ES24.15E3)') solidDefect%omega
    !
    text = "# Total number of <f|U|i> elements, Initial States (bandI, bandF), Final States (bandI, bandF)"
    write(17,'(a, " Format : ''(5i10)''")') trim(text)
    !
    totalNumberOfElements = (iBandIfinal - iBandIinit + 1)*(iBandFfinal - iBandFinit + 1)
    write(17,'(5i10)') totalNumberOfElements, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !
    write(17, '("# Final Band, Initial Band, Delta energy, Complex <f|U|i>, |<f|U|i>|^2 Format : ''(2i10,4ES24.15E3)''")')
    !
    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal
        !
        write(17, 1001) ibf, ibi, eigvI(ibi) - eigvF(ibf), Ufi(ibf,ibi,ik), abs(Ufi(ibf,ibi,ik))**2
        !    
      enddo
    enddo
    !
    close(17)
    !
    call cpu_time(t2)
    write(iostd, '(" Writing Ufi(:,:) done in:                   ", f10.2, " secs.")') t2-t1
    !
 1001 format(2i10,4ES24.15E3)
    !
    return
    !
  end subroutine writeResults
   !
   !
  subroutine bessel_j (x, lmax, jl)
    !
    ! x is the argument of j, jl(0:lmax) is the output values.
    implicit none
    integer, intent(in) :: lmax
    real(kind = dp), intent(in) :: x
    real(kind = dp), intent(out) :: jl(0:lmax)
    integer :: l
    !
    if (x <= 0.0_dp) then
      jl = 0.0_dp
      jl(0) = 1.0_dp
      return
    end if
    !
    jl(0) = sin(x)/x
    if (lmax <= 0) return
    jl(1) = (jl(0)-cos(x))/x
    if (lmax == 1) return
    !
    do l = 2, lmax
      jl(l) = dble(2*l-1)*jl(l-1)/x - jl(l-2)
    enddo
    !
    return
    !
  end subroutine bessel_j
  !
  !
  !
  subroutine ylm(v_in,lmax,y)
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
  !
  !   In/Output :
  !
      integer, intent(in) :: LMAX
      real(kind = dp), intent(in) :: V_in(3)
      complex(kind = dp), intent(out) :: Y(*)
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
end module TMEModule
