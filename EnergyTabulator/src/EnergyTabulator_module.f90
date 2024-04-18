module energyTabulatorMod
  
  use constants, only: dp, eVToHartree
  use miscUtilities, only: int2str
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nKPoints, nSpins, loopSpins, ispSelect
  use errorsAndMPI
  use mpi

  implicit none

  ! Global variables not passed as arguments:
  integer :: ikStart, ikEnd
    !! Start and end k-points for each process
  integer :: nkPerProc
    !! Number of k-points on each process


  ! Variables that should be passed as arguments
  integer :: ibShift_eig
    !! Optional shift of eigenvalue bands relative to
    !! total-energy bands
  integer :: ikIinit, ikIfinal, ikFinit, ikFfinal
    !! K-point bounds for initial and final state
  integer :: refBand
    !! Band of WZP reference carrier

  real(kind=dp) :: eCorrectTot
    !! Total-energy correction, if any
  real(kind=dp) :: eCorrectEigRef
    !! Correction to eigenvalue difference with reference carrier, if any

  character(len=300) :: allStatesBaseDir
    !! Base dir for sets of relaxed files if not captured
  character(len=300) :: energyTableDir
    !! Path to energy tables
  character(len=300) :: exportDirEigs
    !! Path to export for system to get eigenvalues
  character(len=300) :: exportDirInitInit
    !! Path to export for initial charge state
    !! in the initial positions
  character(len=300) :: exportDirFinalInit
    !! Path to export for final charge state
    !! in the initial positions
  character(len=300) :: exportDirFinalFinal
    !! Path to export for final charge state
    !! in the final positions
  character(len=300) :: singleStateExportDir
    !! Export dir name within each subfolder

  logical :: captured
    !! If carrier is captured as opposed to scattered
  logical :: elecCarrier
    !! If carrier is electron as opposed to hole


  contains

!----------------------------------------------------------------------------
  subroutine readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, ikFinit, ikFfinal, &
        ispSelect, refBand, eCorrectTot, eCorrectEigRef, allStatesBaseDir, energyTableDir, exportDirEigs, &
        exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, singleStateExportDir, captured, elecCarrier, loopSpins)

    implicit none

    ! Output variables:
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ibShift_eig
      !! Optional shift of eigenvalue bands relative to
      !! total-energy bands
    integer, intent(out) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(out) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(out) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

    character(len=300), intent(out) :: allStatesBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(out) :: exportDirEigs
      !! Path to export for system to get eigenvalues
    character(len=300), intent(out) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions
    character(len=300), intent(out) :: exportDirFinalInit
      !! Path to export for final charge state
       !! in the initial positions
    character(len=300), intent(out) :: exportDirFinalFinal
      !! Path to export for final charge state
      !! in the final positions
    character(len=300), intent(out) :: singleStateExportDir
      !! Export dir name within each subfolder

    logical, intent(out) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(out) :: elecCarrier
      !! If carrier is electron as opposed to hole
    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel
  
    namelist /inputParams/ exportDirEigs, exportDirFinalFinal, exportDirFinalInit, exportDirInitInit, energyTableDir, &
                           eCorrectTot, eCorrectEigRef, captured, elecCarrier, ispSelect, allStatesBaseDir, singleStateExportDir, &
                           iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, refBand, ikIinit, ikIfinal, ikFinit, ikFfinal, &
                           ibShift_eig


    if(ionode) then

      ! Set default values for input variables and start timers
      call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, ikFinit, ikFfinal, ispSelect, &
            refBand, eCorrectTot, eCorrectEigRef, allStatesBaseDir, energyTableDir, exportDirEigs, exportDirInitInit, &
            exportDirFinalInit, exportDirFinalFinal, singleStateExportDir, captured, elecCarrier)
    
      ! Read input variables
      read(5, inputParams, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputs', 'reading inputParams namelist', abs(ierr))

      ! Check that all variables were properly set
      call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, ikFinit, ikFfinal, &
            ispSelect, refBand, eCorrectTot, eCorrectEigRef, allStatesBaseDir, energyTableDir, exportDirEigs, &
            exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, singleStateExportDir, captured, elecCarrier, loopSpins)

    endif


    ! Broadcast all input variables
    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ibShift_eig, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikFfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(refBand, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(ispSelect, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(loopSpins, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(eCorrectTot, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eCorrectEigRef, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(allStatesBaseDir, len(allStatesBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirEigs, len(exportDirEigs), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirInitInit, len(exportDirInitInit), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirFinalInit, len(exportDirFinalInit), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirFinalFinal, len(exportDirFinalFinal), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(singleStateExportDir, len(singleStateExportDir), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(captured, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(elecCarrier, 1, MPI_LOGICAL, root, worldComm, ierr)

    return

  end subroutine readInputs

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, ikFinit, ikFfinal, &
        ispSelect, refBand, eCorrectTot, eCorrectEigRef, allStatesBaseDir, energyTableDir, exportDirEigs, exportDirInitInit, &
        exportDirFinalInit, exportDirFinalFinal, singleStateExportDir, captured, elecCarrier)
    !! Set the default values for input variables and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    ! Input variables:
    !integer, intent(in) :: nProcs
      ! Number of processes


    ! Output variables:
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ibShift_eig
      !! Optional shift of eigenvalue bands relative to
      !! total-energy bands
    integer, intent(out) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(out) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(out) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

    character(len=300), intent(out) :: allStatesBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(out) :: exportDirEigs
      !! Path to export for system to get eigenvalues
    character(len=300), intent(out) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions
    character(len=300), intent(out) :: exportDirFinalInit
      !! Path to export for final charge state
      !! in the initial positions
    character(len=300), intent(out) :: exportDirFinalFinal
      !! Path to export for final charge state
      !! in the final positions
    character(len=300), intent(out) :: singleStateExportDir
      !! Export dir name within each subfolder

    logical, intent(out) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(out) :: elecCarrier
      !! If carrier is electron as opposed to hole


    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    ibShift_eig = 0
    ikIinit = -1
    ikIfinal = -1
    ikFinit = -1
    ikFfinal = -1
    refBand = -1

    ispSelect = -1

    eCorrectTot = 0.0_dp
    eCorrectEigRef = 0.0_dp

    allStatesBaseDir = ''
    energyTableDir = './'
    exportDirEigs = ''
    exportDirInitInit = ''
    exportDirFinalInit = ''
    exportDirFinalFinal = ''
    singleStateExportDir = ''

    captured = .true.
    elecCarrier = .true.

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, ikFinit, ikFfinal, &
        ispSelect, refBand, eCorrectTot, eCorrectEigRef, allStatesBaseDir, energyTableDir, exportDirEigs, &
        exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, singleStateExportDir, captured, elecCarrier, loopSpins)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ibShift_eig
      !! Optional shift of eigenvalue bands relative to
      !! total-energy bands
    integer, intent(in) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(inout) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(inout) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

    character(len=300), intent(in) :: allStatesBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues
    character(len=300), intent(in) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions
    character(len=300), intent(in) :: exportDirFinalInit
      !! Path to export for final charge state
      !! in the initial positions
    character(len=300), intent(in) :: exportDirFinalFinal
      !! Path to export for final charge state
      !! in the final positions
    character(len=300), intent(in) :: singleStateExportDir
      !! Export dir name within each subfolder

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: elecCarrier
      !! If carrier is electron as opposed to hole

    ! Output variables:
    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9))
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 
    abortExecution = checkIntInitialization('refBand', refBand, 1, int(1e9)) .or. abortExecution

    write(*,'("ibShift_eig = ", i4)') ibShift_eig


    if(ispSelect < 1 .or. ispSelect > 2) then
      write(*,*) "No valid choice for spin channel selection given. Looping over spin."
      loopSpins = .true.
    else
      write(*,'("Only exporting spin channel ", i2)') ispSelect
      loopSpins = .false.
    endif


    write(*,'("eCorrectEigRef = ", f8.4, " (eV)")') eCorrectEigRef


    write(*,'("captured = ",L)') captured
    if(captured) then
      if(iBandFinit /= iBandFfinal) then
        write(*,'("Capture only expected for a single final-state band.")')
        abortExecution = .true.
      endif

      abortExecution = checkDirInitialization('exportDirInitInit', exportDirInitInit, 'input') .or. abortExecution
      abortExecution = checkDirInitialization('exportDirFinalInit', exportDirFinalInit, 'input') .or. abortExecution
      abortExecution = checkDirInitialization('exportDirFinalFinal', exportDirFinalFinal, 'input') .or. abortExecution

      write(*,'("eCorrectTot = ", f8.4, " (eV)")') eCorrectTot

    else
      abortExecution = checkIntInitialization('ikIinit', ikIinit, 1, int(1e9))
      abortExecution = checkIntInitialization('ikIfinal', ikIfinal, ikIinit, int(1e9)) .or. abortExecution
      abortExecution = checkIntInitialization('ikFinit', ikFinit, 1, int(1e9)) .or. abortExecution
      abortExecution = checkIntInitialization('ikFfinal', ikFfinal, ikFinit, int(1e9)) .or. abortExecution 

      write(*,'("allStatesBaseDir = ",a)') trim(allStatesBaseDir)
      write(*,'("singleStateExportDir = ",a)') trim(singleStateExportDir)
    endif

    write(*,'("elecCarrier = ",L)') elecCarrier


    eCorrectTot = eCorrectTot*eVToHartree
    eCorrectEigRef = eCorrectEigRef*eVToHartree


    abortExecution = checkDirInitialization('exportDirEigs', exportDirEigs, 'input') .or. abortExecution

    call system('mkdir -p '//trim(energyTableDir))


    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif
    
    return

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine getnSpinsAndnKPoints(exportDirEigs, nKPoints, nSpins)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Input variables:
    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues
     
    ! Output variables:
    integer, intent(out) :: nKPoints
      !! Number of k-points
    integer, intent(out) :: nSpins
      !! Number of spin channels
    
    
    if(ionode) then
    
      open(50, file=trim(exportDirEigs)//'/input', status = 'old')
    
      ! Ignore cell volume, number of G-vectors, and 
      ! FFT grid size
      call ignoreNextNLinesFromFile(50,6)
      read(50,*) ! nSpins comment
      read(50, '(i10)') nSpins
      read(50,*) ! nKPoints comment
      read(50, '(i10)') nKPoints

    endif

    call MPI_BCAST(nSpins, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
   
    return

  end subroutine getnSpinsAndnKPoints

!----------------------------------------------------------------------------
  subroutine calcAndWriteCaptureEnergies(iBandIinit, iBandIfinal, iBandFinit, ispSelect, nSpins, refBand, eCorrectTot, &
        eCorrectEigRef, elecCarrier, loopSpins, energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nSpins
      !! Number of spin channels
    integer, intent(in) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(in) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(in) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

    logical, intent(in) :: elecCarrier
      !! If carrier is electron as opposed to hole
    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues
    character(len=300), intent(in) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions
    character(len=300), intent(in) :: exportDirFinalInit
      !! Path to export for final charge state
      !! in the initial positions
    character(len=300), intent(in) :: exportDirFinalFinal
      !! Path to export for final charge state
      !! in the final positions

    ! Local variables:
    integer :: isp, ikLocal, ikGlobal, ibi
      !! Loop indices
    integer :: nTransitions
      !! Total number of overlaps to output

    real(kind=dp) :: dEDelta
      !! Energy to be used in delta function
    real(kind=dp) :: dEEigRef
      !! Eigenvalue difference from initial to reference
    real(kind=dp) :: dEEigRefDefect
      !! Eigenvalue difference between the reference
      !! eigenvalue and the defect level
    real(kind=dp) :: dEFirst
      !! Energy to be used in first-order matrix element
    real(kind=dp) :: dETotElecOnly
      !! Total energy difference between charge states
      !! with no change in atomic positions to get the
      !! electronic-only energy to be used in the 
      !! zeroth-order matrix element
    real(kind=dp) :: dETotWRelax
      !! Total energy difference between relaxed
      !! charge states to be used in delta function
    real(kind=dp) :: dEZeroth
      !! Energy to be used in zeroth-order matrix element
    real(kind=dp) :: eigvI(iBandIinit:iBandIfinal)
      !! Initial-state eigenvalues
    real(kind=dp) :: eTotInitInit
      !! Total energy of the relaxed initial charge
      !! state (initial positions)
    real(kind=dp) :: eTotFinalInit
      !! Total energy of the unrelaxed final charge
      !! state (initial positions)
    real(kind=dp) :: eTotFinalFinal
      !! Total energy of the relaxed final charge
      !! state (final positions)
    real(kind=dp) :: refEig
      !! Eigenvalue of WZP reference carrier
    real(kind=dp) :: t1, t2
      !! Timers

    logical :: fileExists
      !! If the input file exists in the given exportDir
    
    character(len = 300) :: text
      !! Text for header


    if(ionode) then
      ! Get total energies from exports of all different structures
      call getTotalEnergy(exportDirInitInit, eTotInitInit, fileExists)
      if(.not. fileExists) call exitError('calcAndWriteCaptureEnergies', 'exportDirInitInit/input must exist!', 1)

      call getTotalEnergy(exportDirFinalInit, eTotFinalInit, fileExists)
      if(.not. fileExists) call exitError('calcAndWriteCaptureEnergies', 'exportDirFinalInit/input must exist!', 1)

      call getTotalEnergy(exportDirFinalFinal, eTotFinalFinal, fileExists)
      if(.not. fileExists) call exitError('calcAndWriteCaptureEnergies', 'exportDirFinalFinal/input must exist!', 1)
    endif

    call MPI_BCAST(eTotInitInit, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eTotFinalInit, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eTotFinalFinal, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)


    do isp = 1, nSpins
      if(loopSpins .or. isp == ispSelect) then

        ! Get reference eigenvalue from the Gamma point
        if(ionode) call getSingleEig(refBand, 1, isp, exportDirEigs, refEig)
        call MPI_BCAST(refBand, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

        call getRefToDefectEigDiff(iBandFinit, isp, refBand, exportDirInitInit, elecCarrier, dEEigRefDefect)

        do ikLocal = 1, nkPerProc

          ikGlobal = ikLocal+ikStart-1
    

          ! Update status to user
          call cpu_time(t1)
          write(*, '(" Writing energy table of k-point ", i2, " and spin ", i1, ".")') ikGlobal, isp
    

          ! Open file and write header (same for capture and scattering)
          open(17, file=trim(energyTableDir)//"/energyTable."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)), status='unknown')


          ! Include in the output file that these energies are tabulated for capture
          write(17,'("# Energies tabulated for capture? Alternative is scattering.)")')
          write(17,'(L4)') .true.

    
          text = "# Total number of transitions, Initial States (bandI, bandF), Final State (band)"
          write(17,'(a, " Format : ''(5i10)''")') trim(text)
    
          nTransitions = (iBandIfinal - iBandIinit + 1)
          write(17,'(5i10)') nTransitions, iBandIinit, iBandIfinal, iBandFinit
    

          ! Get the total energy difference between the two charge states, 
          ! including the atomic relaxation energy (with a potential energy
          ! correction defined by the user). This dE is used in the delta 
          ! function.
          dETotWRelax = eTotFinalFinal - eTotInitInit + eCorrectTot

          write(17,'("# Total-energy difference (Hartree). Format: ''(ES24.15E3)''")') 
          write(17,'("# With relaxation (for delta function)")')
          write(17,'(ES24.15E3)') dETotWRelax


          ! Get the total energy difference between the two charge states, 
          ! not including atomic relaxation (with a potential energy correction
          ! defined by the user). This dE represents the total electronic-only
          ! energy difference between the two charge states and goes in the
          ! zeroth-order matrix element.
          dETotElecOnly = eTotFinalInit - eTotInitInit + eCorrectTot

          write(17,'("# Electronic only without relaxation (for zeroth-order matrix element)")')
          write(17,'(ES24.15E3)') dETotElecOnly
    

          ! Output header for main data and loop over bands
          text = "# Final Band, Initial Band, Delta Function (Hartree), Zeroth-order (Hartree), First-order (Hartree), Plotting (eV)" 
          write(17, '(a, " Format : ''(2i10,4ES24.15E3)''")') trim(text)

          
          open(27,file="dEPlot."//trim(int2str(isp))//trim(int2str(ikGlobal)))
          write(27,'("# ib, Total elec. energy diff. from ref. (eV)")')


          call readEigenvalues(iBandIinit, iBandIfinal, ikGlobal, isp, exportDirEigs, eigvI)

          do ibi = iBandIinit, iBandIfinal

            ! All of the energies require an eigenvalue difference from a reference band.
            ! Calculate this once to be used for all of the energies.
            !
            ! The energy correction `eCorrectEigRef` should be zero if the reference state
            ! and initial state are both in the conduction band or both in the valence band,
            ! since eigenvalue differences within the bands are okay at the PBE level and 
            ! do not need to be corrected. One example where this correction would be needed 
            ! is setting up the negative charge state of the triply-hydrogenated Si defect. 
            ! The simplest treatment of that charge state has the reference carrier in the 
            ! valence band top, so the distance between the valence band and the conduction 
            ! band would need to be corrected from PBE to HSE.
            !
            ! Switch the order of the eigenvalue subtraction for hole vs electron capture
            ! to represent that the actual energy we need is that of the electron.
            if(elecCarrier) then
              dEEigRef = eigvI(ibi) - refEig + eCorrectEigRef
            else
              dEEigRef = refEig - eigvI(ibi) + eCorrectEigRef
            endif

    
            ! To get the total energy that needs to be conserved (what goes into the delta
            ! function), add the total energy difference between the two relaxed charge states
            ! and the additional eigenvalue energy difference between the initial state and
            ! the WZP reference-carrier state. 
            dEDelta = dETotWRelax - dEEigRef


            ! The zeroth-order matrix element contains the electronic-only energy difference.
            ! We get that from a total energy difference between the two charge states in the
            ! initial positions. Like in the energy for the delta function, the additional 
            ! carrier energy must also be included with a potential correction.
            dEZeroth = dETotElecOnly - dEEigRef


            ! The first-order term contains only the unperturbed eigenvalue difference. The
            ! perturbative expansion has \(\varepsilon_i - \varepsilon_f\), in terms of the 
            ! actual electron (rather than hole). The absolute value is needed for the hole 
            ! case.
            ! 
            ! For capture, the defect level should not have dispersion, and the level of the 
            ! defect changes between charge states, so we use a single reference energy and
            ! distance to the defect, taken from the initial-charge-state, Gamma-only HSE 
            ! calculation. 
            dEFirst = dEEigRef + dEEigRefDefect
        

            write(27,'(2i7,2f12.4)') ibi, abs(dEZeroth)/eVToHartree
            write(17,'(2i7,3ES24.15E3)') iBandFinit, ibi, dEDelta, dEZeroth, dEFirst
            
          enddo ! Loop over initial bands

    
          close(27)
          close(17)
    
          call cpu_time(t2)
          write(*, '(" Writing energy table of k-point ", i4, "and spin ", i1, " done in:                   ", f10.2, " secs.")') &
            ikGlobal, isp, t2-t1
    
 1001     format(2i7,3ES24.15E3)

        enddo ! Loop over local k-points
      endif ! If we should calculate this spin
    enddo ! Loop over spins

    return

  end subroutine calcAndWriteCaptureEnergies

!----------------------------------------------------------------------------
  subroutine calcAndWriteScatterEnergies(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ibShift_eig, ikIinit, ikIfinal, &
        ikFinit, ikFfinal, ispSelect, nSpins, refBand, eCorrectEigRef, elecCarrier, loopSpins, allStatesBaseDir, energyTableDir, &
        exportDirEigs, singleStateExportDir)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ibShift_eig
      !! Optional shift of eigenvalue bands relative to
      !! total-energy bands
    integer, intent(in) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nSpins
      !! Number of spin channels
    integer, intent(in) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(in) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

    logical, intent(in) :: elecCarrier
      !! If carrier is electron as opposed to hole
    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    character(len=300), intent(in) :: allStatesBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues
    character(len=300), intent(in) :: singleStateExportDir
      !! Export dir name within each subfolder

    ! Local variables:
    integer :: ikMin, ikMax, ibMin, ibMax
      !! Overall bounds on k-points and bands for initial 
      !! and final states
    integer :: isp, ik, ib, iki, ikf, ibi, ibf, iE
      !! Loop indices
    integer :: ib_minETot, ik_minETot
      !! Indices of state with minimum total energy
    integer :: nTransitions
      !! Total number of energies to output

    real(kind=dp) :: dEDelta
      !! Energy to be used in delta function
    real(kind=dp) :: dEFirst
      !! Energy to be used in first-order matrix element
    real(kind=dp) :: dEZeroth
      !! Energy to be used in zeroth-order matrix element
    real(kind=dp), allocatable :: eigv(:,:)
      !! Eigenvalues
    real(kind=dp), allocatable :: eTot(:,:)
      !! Total energies for all states
    real(kind=dp) :: minETot
      !! Minimum total energy (total potential)
    real(kind=dp) :: refEig
      !! Eigenvalue of WZP reference carrier

    logical :: fileExists
      !! If the input file exists in the given exportDir
    logical :: inInitkRange, inFinalkRange, inInitBandRange, inFinalBandRange
      !! If in k-point or band range
    logical, allocatable :: skipState(:,:)
      !! If a state should be skipped
    
    character(len=300) :: baseFName
      !! Base file name for energy table
    character(len=300) :: path
      !! Path to the export for each band state
    character(len=300) :: text
      !! Text for header


    ikMin = min(ikIinit, ikFinit)
    ikMax = max(ikIfinal, ikFfinal)

    ibMin = min(iBandIinit, iBandFinit)
    ibMax = max(iBandIfinal, iBandFfinal)

    allocate(eTot(ibMin:ibMax,ikMin:ikMax))
    allocate(skipState(ibMin:ibMax,ikMin:ikMax))
    allocate(eigv(ibMin+ibShift_eig:ibMax+ibShift_eig,ikMin:ikMax))
      ! The bands in the eigenvalue system may not line up
      ! with the bands in the total-energy (defect) system,
      ! so give the user the option to shift the bands so
      ! that they line up.


    ! Get total energies from exports of all different structures
    eTot = 0.0_dp
    minETot = 0.0_dp
    skipState = .true.
    do ik = ikMin, ikMax

      ! Test if in range of either k bounds
      inInitkRange = ik >= ikIinit .or. ik <= ikIfinal
      inFinalkRange = ik >= ikFinit .or. ik <= ikFfinal

      if(inInitkRange .or. inFinalkRange) then
        do ib = ibMin, ibMax

          ! Test if in range of either band bounds
          inInitBandRange = ib >= iBandIinit .or. ib <= iBandIfinal
          inFinalBandRange = ib >= iBandFinit .or. ib <= iBandFfinal

          ! Only consider if band and k bounds line up for initial/final states
          if((inInitkRange .and. inInitBandRange) .or. (inFinalkRange .and. inFinalBandRange)) then

            ! Assume subfolders have pattern <allStatesBaseDir>/k<ik>_b<ib>/<singleStateExportDir>
            ! e.g., ../VASP/k1_b1616/export/
            path = trim(allStatesBaseDir)//'/k'//trim(int2str(ik))//'_b'//trim(int2str(ib))//'/'//trim(singleStateExportDir)
            call getTotalEnergy(path, eTot(ib,ik), fileExists)

            ! Skip the consideration of this state if the necessary file
            ! doesn't exist. Otherwise, track the minimum total energy to 
            ! be output.
            if(fileExists) then
              skipState(ib,ik) = .false.

              if(eTot(ib,ik) < minETot) then
                minETot = eTot(ib,ik)
                ib_minETot = ib
                ik_minETot = ik
              endif
            else
              write(*,'("Skipping state ik,ib = ",2i7,"! Input file does not exist in path ",a)') ik, ib, trim(path)
            endif

          endif
        enddo
      endif
    enddo


    do isp = 1, nSpins
      if(loopSpins .or. isp == ispSelect) then

        ! Update status to user
        write(*, '(" Writing energy table of spin ", i1, ".")') isp

        ! Get reference eigenvalue from the Gamma point
        call getSingleEig(refBand+ibShift_eig, 1, isp, exportDirEigs, refEig)
    
        ! Read in all eigenvalues needed and output different plotting energies
        open(27,file="dEPlot."//trim(int2str(isp)))
        write(27,'("# Min. total energy of ",f12.4," eV at ik,ib = ",2i7)') minETot/eVToHartree, ik_minETot, ib_minETot
        write(27,'("# ik, ib, Eig. diff. from ref. (eV), Total energy diff. from min. (eV)")')

        eigv = 0.0_dp
        do ik = ikMin, ikMax

          inInitkRange = ik >= ikIinit .or. ik <= ikIfinal
          inFinalkRange = ik >= ikFinit .or. ik <= ikFfinal

          if(inInitkRange .or. inFinalkRange) &
            call readEigenvalues(ibMin+ibShift_eig, ibMax+ibShift_eig, ik, isp, exportDirEigs, eigv(:,ik))
              ! Go ahead and read in all possible bands at all k-points because it
              ! is simple to include more bands even if they aren't all needed.

          do ib = ibMin, ibMax
            if(.not. skipState(ib,ik)) write(27,'(2i7,2f12.4)') &
              ik, ib, (eigv(ib+ibShift_eig,ik)-refEig+eCorrectEigRef)/eVToHartree, (eTot(ib,ik)-minETot)/eVToHartree
          enddo
        enddo

        close(27)


        ! Open temporary file without header
        baseFName = trim(energyTableDir)//"/energyTable."//trim(int2str(isp))
        open(17, file=trim(baseFName)//".tmp", status='unknown')


        iE = 0

        write(*,'("States skiped:")')
        write(*,'("   iki  ibi  ikf  ibf  Reason Skipped")')
        write(*,'("--------------------------------------")')

        do iki = ikIinit, ikIfinal
          do ibi = iBandIinit, iBandIfinal

            if(skipState(ibi,iki)) cycle

            do ikf = ikFinit, ikFfinal
              do ibf = iBandFinit, iBandFfinal

                if(skipState(ibf,ikf) .or. (ikf == iki .and. ibf == ibi)) cycle

                iE = iE + 1

    
                ! For scattering, we do separate calculations for each band state, so we only need
                ! to use the total energy difference
                !
                ! Switch the order for hole vs electron capture to represent that the actual energy 
                ! we need is that of the electron. This is needed because it is expected that the
                ! user will index initial/final bands based on the carrier (electron/hole) and not
                ! just the actual electron that is transitioning.
                if(elecCarrier) then
                  dEDelta = eTot(ibf,ikf) - eTot(ibi,iki)
                else
                  dEDelta = eTot(ibi,iki) - eTot(ibf,ikf)
                endif

                if(dEDelta >= 0.0_dp) then
                  ! If this delta-function energy is positive, the phonon energy will have to be
                  ! negative, and we don't want to consider cooling processes.

                  write(*,'("   ", 4i5, " Zero or neg. energy transfer")') iki, ibi, ikf, ibf

                  iE = iE - 1

                  cycle
                endif


                ! The zeroth-order matrix element contains the electronic-only energy difference.
                ! For scattering, this is just an eigenvalue difference. The eigenvalues must be
                ! read from the same system to make the difference meaningful. Switch the order
                ! for different carriers
                if(elecCarrier) then
                  dEZeroth = eigv(ibf+ibShift_eig,ikf) - eigv(ibi+ibShift_eig,iki)
                else
                  dEZeroth = eigv(ibi+ibShift_eig,iki) - eigv(ibf+ibShift_eig,ikf)
                endif

                
                if(abs(dEZeroth) < 1e-6*eVToHartree) then
                  ! If there is zero electronic energy difference, the matrix elements will be zero.

                  write(*,'("   ", 4i5, " Zero elec. energy diff.")') iki, ibi, ikf, ibf

                  iE = iE - 1

                  cycle
                endif

                
                ! The first-order term contains only the unperturbed eigenvalue difference. The
                ! perturbative expansion has \(\varepsilon_i - \varepsilon_f\), in terms of the 
                ! actual electron (rather than hole). We cannot take the absolute value here as
                ! in capture because we do not assume that the carrier always loses energy.
                !
                ! Because the zeroth-order term only has an eigenvalue difference, we can just
                ! take the negative of that value. The sign doesn't actually matter because it
                ! is squared in the matrix element, but we calculate it correctly for clarity.
                dEFirst = -dEZeroth


                write(17,'(4i7,3ES24.15E3)') iki, ibi, ikf, ibf, dEDelta, dEZeroth, dEFirst
        
              enddo ! Loop over final bands 
            enddo ! Loop over final k-points
          enddo ! Loop over initial bands
        enddo ! Loop over initial k-points

        nTransitions = iE

        close(17)


        ! Open file and write header
        open(17, file=trim(baseFName), status='unknown')


        ! Include in the output file that these energies are tabulated for capture
        write(17,'("# Energies tabulated for capture? Alternative is scattering.)")')
        write(17,'(L4)') .false.
    
        text = "# Total number of transitions, Initial States (kI, kF, bandI, bandF), Final States (kI, kF, bandI, bandF)"
        write(17,'(a, " Format : ''(9i10)''")') trim(text)
    
        write(17,'(9i10)') nTransitions, ikIinit, ikIfinal, iBandIinit, iBandIfinal, &
                                                  ikFinit, ikFfinal, iBandFinit, iBandFfinal
    

        ! Output header for main data and loop over bands
        text = "# iki, ibi, ikf, ibf, Delta Function (Hartree), Zeroth-order (Hartree), First-order (Hartree)" 
        write(17, '(a, " Format : ''(2i10,3ES24.15E3)''")') trim(text)

        close(17)

        call system('cat '//trim(baseFName)//'.tmp >> '//trim(baseFName))
        call system('rm '//trim(baseFName)//'.tmp')

      endif
    enddo

    return

  end subroutine calcAndWriteScatterEnergies

!----------------------------------------------------------------------------
  subroutine getTotalEnergy(exportDir, eTot, fileExists)

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    character(len=300), intent(in) :: exportDir
      !! Path to export directory

    ! Output variables:
    real(kind=dp), intent(out) :: eTot
      !! Total energy read from the `input` file

    logical, intent(out) :: fileExists
      !! If the input file exists in the given exportDir

    ! Local variables:
    character(len=300) :: line
      !! Line from file


    inquire(file=trim(exportDir)//'/input', exist=fileExists)

    if(.not. fileExists) return

    open(30,file=trim(exportDir)//'/input')
    line = getFirstLineWithKeyword(30, 'Total Energy')
    read(30,*) eTot
    close(30)

    return

  end subroutine getTotalEnergy

!----------------------------------------------------------------------------
  subroutine getSingleEig(ib, ik, isp, exportDir, eig_ib)

    use miscUtilities, only: ignoreNextNLinesFromFile, int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: ib
      !! Band index to get eigenvalue for
    integer, intent(in) :: ik
      !! K-point index
    integer, intent(in) :: isp
      !! Current spin channel

    character(len=300), intent(in) :: exportDir
      !! Path to export for system to get eigenvalue

    ! Output variables:
    real(kind=dp), intent(out) :: eig_ib
      !! Eigenvalue at band index

    ! Local variables:
    integer :: iDum
      !! Dummy integer to ignore band index

    real(kind=dp) :: rDum
      !! Dummy real to ignore occupation

    character(len=300) :: fName
      !! File name


    fName = trim(exportDir)//"/eigenvalues."//trim(int2str(isp))//"."//trim(int2str(ik))

    open(72, file=fName)

    call ignoreNextNLinesFromFile(72, 2 + (ib-1))
      ! Ignore header and all bands before reference band
    
    read(72, '(i10, 2ES24.15E3)') iDum, eig_ib, rDum
  
    close(72)

    return

  end subroutine getSingleEig

!----------------------------------------------------------------------------
  subroutine getRefToDefectEigDiff(iBandFinit, isp, refBand, exportDirInitInit, elecCarrier, dEEigRefDefect)

    use miscUtilities, only: ignoreNextNLinesFromFile, int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandFinit
      !! Energy band for final state
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: refBand
      !! Band of WZP reference carrier

    character(len=300), intent(in) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions

    logical, intent(in) :: elecCarrier
      !! If carrier is electron as opposed to hole

    ! Output variables:
    real(kind=dp), intent(out) :: dEEigRefDefect
      !! Eigenvalue difference between the reference
      !! eigenvalue and the defect level

    ! Local variables:
    real(kind=dp) :: defectEig
      !! Eigenvalue of the defect for capture
    real(kind=dp) :: refEig
      !! Eigenvalue of WZP reference carrier (not the same
      !! as `refEig` used in other places in the code; this
      !! is just to determine the correct distance between
      !! the defect and the band states


    if(ionode) then

      ! Get eigenvalues for the reference band and defect level from 
      ! the relaxed initial charge state because that correctly reflects
      ! the binding energy.
      call getSingleEig(refBand, 1, isp, exportDirInitInit, refEig)
      call getSingleEig(iBandFinit, 1, isp, exportDirInitInit, defectEig)


      ! Switch the order of the eigenvalue subtraction for hole vs electron
      ! capture to represent that the actual energy we need is that of the 
      ! electron
      if(elecCarrier) then
        dEEigRefDefect = refEig - defectEig
      else
        dEEigRefDefect = defectEig - refEig
      endif

    endif

    call MPI_BCAST(dEEigRefDefect, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine getRefToDefectEigDiff
  
!----------------------------------------------------------------------------
  subroutine readEigenvalues(iBandInit, iBandFinal, ikGlobal, isp, exportDirEigs, eigv)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none
    
    ! Input variables
    integer, intent(in) :: iBandInit, iBandFinal
      !! Energy band bounds
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues

    ! Output variables:
    real(kind=dp), intent(out) :: eigv(iBandInit:iBandFinal)
      !! Eigenvalues

    ! Local variables:
    integer :: ib
      !! Loop index
    integer :: iDum
      !! Dummy integer to ignore band index

    real(kind=dp) :: rDum
      !! Dummy real to ignore occupation

    character(len=300) :: fName
      !! File name

    
    fName = trim(exportDirEigs)//"/eigenvalues."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))

    open(72, file=fName)

    ! Ignore header and all bands before lowest band
    call ignoreNextNLinesFromFile(72, 2 + (iBandInit-1))
    
    do ib = iBandInit, iBandFinal
      read(72, '(i10, 2ES24.15E3)') iDum, eigv(ib), rDum
    enddo
    
    close(72)
    
    return
    
  end subroutine readEigenvalues
  
!----------------------------------------------------------------------------
  subroutine readCaptureEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, energyTableDir, dE)
    !! Read all energies from energy table and store in dE
    
    use miscUtilities, only: ignoreNextNLinesFromFile, int2str
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table

    ! Output variables:
    real(kind=dp), intent(out) :: dE(iBandFinit:iBandFfinal,iBandIinit:iBandIFinal,3)
      !! All energy differences from energy table

    ! Local variables:
    integer :: iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      !! Energy band bounds for initial and final state from the energy table file
    integer :: ibi, ibf
      !! Loop indices
    integer :: iDum
      !! Dummy integer to ignore input

    real(kind=dp) :: dE_(3)
      !! Input energy

    logical :: captured
      !! If energy table was tabulated for capture

    character(len=300) :: fName
      !! Energy table file name
    
    
    fName = trim(energyTableDir)//"/energyTable."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))

    open(27, file=trim(fName), status='unknown')


    read(27,*)
    read(27,'(L4)') captured

    if(.not. captured) call exitError('readCaptureEnergyTable', 'This energy table was not tabulated for capture!', 1)


    read(27,*)
    read(27,'(5i10)') iDum, iBandIinit_, iBandIfinal_, iBandFinit_

    ! Check the input band bounds against those in the energy file
    if(iBandIinit < iBandIinit_ .or. iBandIfinal > iBandIfinal_ .or. iBandFinit /= iBandFinit_ .or. iBandFfinal /= iBandFinit_) &
      call exitError('readCaptureEnergyTable', 'given band bounds are outside those in energy table '//trim(fName), 1)
    
    call ignoreNextNLinesFromFile(27,6)
    

    do ibi = iBandIinit_, iBandIfinal_
      
      read(27,*) iDum, iDum, dE_

      if(ibi >= iBandIinit .and. ibi <= iBandIfinal) dE(iBandFinit_,ibi,:) = dE_(:)
          
    enddo
    
    close(27)
    
    return
    
  end subroutine readCaptureEnergyTable
  
!----------------------------------------------------------------------------
  subroutine readScatterEnergyTable(isp, energyTableDir, ibi, ibf, iki, ikf, nTransitions, dE)
    !! Read all energies from energy table and store in dE
    
    use miscUtilities, only: ignoreNextNLinesFromFile, int2str
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: isp
      !! Current spin channel

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table

    ! Output variables:
    integer, allocatable, intent(out) :: iki(:), ikf(:), ibi(:), ibf(:)
      !! State indices
    integer, intent(out) :: nTransitions
      !! Total number of transitions 

    real(kind=dp), allocatable, intent(out) :: dE(:,:)
      !! All energy differences from energy table

    ! Local variables:
    integer :: iE
      !! Loop indices

    logical :: captured
      !! If energy table was tabulated for capture

    character(len=300) :: fName
      !! Energy table file name
    
    
    fName = trim(energyTableDir)//"/energyTable."//trim(int2str(isp))

    open(27, file=trim(fName), status='unknown')


    read(27,*)
    read(27,'(L4)') captured

    if(captured) call exitError('readScatterEnergyTable', 'This energy table was tabulated for capture!', 1)


    read(27,*)
    read(27,*) nTransitions

    allocate(ibi(nTransitions), ibf(nTransitions), iki(nTransitions), ikf(nTransitions))
    allocate(dE(3,nTransitions))
    
    read(27,*)
    
    do iE = 1, nTransitions
      
      read(27,'(4i7,3ES24.15E3)') iki(iE), ibi(iE), ikf(iE), ibf(iE), dE(:,iE)

    enddo
    
    close(27)
    
    return
    
  end subroutine readScatterEnergyTable

end module energyTabulatorMod
