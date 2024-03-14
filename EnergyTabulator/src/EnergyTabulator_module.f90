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
  integer :: refBand
    !! Band of WZP reference carrier

  real(kind=dp) :: eCorrectTot
    !! Total-energy correction, if any
  real(kind=dp) :: eCorrectEigRef
    !! Correction to eigenvalue difference with reference carrier, if any

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

  logical :: captured
    !! If carrier is captured as opposed to scattered
  logical :: elecCarrier
    !! If carrier is electron as opposed to hole

  namelist /inputParams/ exportDirEigs, exportDirFinalFinal, exportDirFinalInit, exportDirInitInit, energyTableDir, &
                         eCorrectTot, eCorrectEigRef, captured, elecCarrier, ispSelect, &
                         iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, refBand


  contains

!----------------------------------------------------------------------------
  subroutine readInputs(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef, &
        energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier)

    implicit none

    ! Output variables:
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(out) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(out) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

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

    logical, intent(out) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(out) :: elecCarrier
      !! If carrier is electron as opposed to hole


    if(ionode) then

      ! Set default values for input variables and start timers
      call initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef, energyTableDir, &
            exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier)
    
      ! Read input variables
      read(5, inputParams, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputs', 'reading inputParams namelist', abs(ierr))

      ! Check that all variables were properly set
      call checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef,&
            energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier, loopSpins)

    endif


    ! Broadcast all input variables
    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(refBand, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(ispSelect, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(loopSpins, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(eCorrectTot, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eCorrectEigRef, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirEigs, len(exportDirEigs), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirInitInit, len(exportDirInitInit), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirFinalInit, len(exportDirFinalInit), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirFinalFinal, len(exportDirFinalFinal), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(captured, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(elecCarrier, 1, MPI_LOGICAL, root, worldComm, ierr)

    return

  end subroutine readInputs

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef, &
        energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier)
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
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(out) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(out) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

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

    logical, intent(out) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(out) :: elecCarrier
      !! If carrier is electron as opposed to hole


    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    refBand = -1

    ispSelect = -1

    eCorrectTot = 0.0_dp
    eCorrectEigRef = 0.0_dp

    energyTableDir = './'
    exportDirEigs = ''
    exportDirInitInit = ''
    exportDirFinalInit = ''
    exportDirFinalFinal = ''

    captured = .true.
    elecCarrier = .true.

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, refBand, eCorrectTot, eCorrectEigRef, &
      energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured, elecCarrier, loopSpins)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: refBand
      !! Band of WZP reference carrier

    real(kind=dp), intent(inout) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(inout) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any

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


    if(ispSelect < 1 .or. ispSelect > 2) then
      write(*,*) "No valid choice for spin channel selection given. Looping over spin."
      loopSpins = .true.
    else
      write(*,'("Only exporting spin channel ", i2)') ispSelect
      loopSpins = .false.
    endif


    write(*,'("eCorrectTot = ", f8.4, " (eV)")') eCorrectTot
    write(*,'("eCorrectEigRef = ", f8.4, " (eV)")') eCorrectEigRef


    if(captured .and. iBandFinit /= iBandFfinal) then
      write(*,'("Capture only expected for a single final-state band.")')
      abortExecution = .true.
    endif
    write(*,'("captured = ",L)') captured
    write(*,'("elecCarrier = ",L)') elecCarrier


    eCorrectTot = eCorrectTot*eVToHartree
    eCorrectEigRef = eCorrectEigRef*eVToHartree


    abortExecution = checkDirInitialization('exportDirEigs', exportDirEigs, 'input') .or. abortExecution
    abortExecution = checkDirInitialization('exportDirInitInit', exportDirInitInit, 'input') .or. abortExecution
    abortExecution = checkDirInitialization('exportDirFinalInit', exportDirFinalInit, 'input') .or. abortExecution
    abortExecution = checkDirInitialization('exportDirFinalFinal', exportDirFinalFinal, 'input') .or. abortExecution

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
  subroutine calcAndWriteCaptureEnergies(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ispSelect, nSpins, refBand, eCorrectTot, &
        eCorrectEigRef, elecCarrier, loopSpins, energyTableDir, exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
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
    integer :: isp, ikLocal, ikGlobal, ibi, ibf
      !! Loop indices
    integer :: totalNumberOfElements
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
    real(kind=dp) :: dEPlot
      !! Energy to be used for plotting
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
    real(kind=dp) :: eigvF(iBandFinit:iBandFfinal)
      !! Final-state eigenvalues
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
    
    character(len = 300) :: text
      !! Text for header


    ! Get total energies from exports of all different structures
    call getTotalEnergy(exportDirInitInit, eTotInitInit)
    call getTotalEnergy(exportDirFinalInit, eTotFinalInit)
    call getTotalEnergy(exportDirFinalFinal, eTotFinalFinal)


    do isp = 1, nSpins
      if(loopSpins .or. isp == ispSelect) then

        ! Get reference eigenvalue
        call getRefEig(isp, refBand, exportDirEigs, refEig)

        call getRefToDefectEigDiff(iBandFinit, isp, refBand, exportDirInitInit, elecCarrier, dEEigRefDefect)

        do ikLocal = 1, nkPerProc

          ikGlobal = ikLocal+ikStart-1
    

          ! Update status to user
          call cpu_time(t1)
          write(*, '(" Writing energy table of k-point ", i2, " and spin ", i1, ".")') ikGlobal, isp
    

          ! Open file and write header (same for capture and scattering)
          open(17, file=trim(energyTableDir)//"/energyTable."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)), status='unknown')
    
          text = "# Total number of <f|i> elements, Initial States (bandI, bandF), Final States (bandI, bandF)"
          write(17,'(a, " Format : ''(5i10)''")') trim(text)
    
          totalNumberOfElements = (iBandIfinal - iBandIinit + 1)*(iBandFfinal - iBandFinit + 1)
          write(17,'(5i10)') totalNumberOfElements, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal


          ! Include in the output file that these energies are tabulated for capture
          write(17,'("# Energies tabulated for capture? Alternative is scattering.)")')
          write(17,'(L4)') .true.
    

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


          call readEigenvalues(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, exportDirEigs, eigvF, eigvI)

          do ibf = iBandFinit, iBandFfinal
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
              ! distance to the defect, taken from the initial-charge-state HSE calculation. 
              ! For the scattering case, the bands do not have those issues, so we can directly 
              ! use the final-state eigenvalue.
              if(captured) then
                dEFirst = dEEigRef + dEEigRefDefect
              else
                dEFirst = abs(eigvI(ibi) - eigvF(ibf))
              endif

          
              ! Energy plotted should be the positive, electronic-only energy difference (same
              ! as in zeroth-order) in eV
              dEPlot = abs(dEZeroth)/eVToHartree
        

              write(17, 1001) ibf, ibi, dEDelta, dEZeroth, dEFirst, dEPlot
            
            enddo ! Loop over initial bands
          enddo ! Loop over final bands

    
          close(17)
    
          call cpu_time(t2)
          write(*, '(" Writing energy table of k-point ", i4, "and spin ", i1, " done in:                   ", f10.2, " secs.")') &
            ikGlobal, isp, t2-t1
    
 1001     format(2i7,4ES24.15E3)

        enddo ! Loop over local k-points
      endif ! If we should calculate this spin
    enddo ! Loop over spins

    return

  end subroutine calcAndWriteCaptureEnergies

!----------------------------------------------------------------------------
  subroutine calcAndWriteScatterEnergies()

    implicit none

    return

  end subroutine calcAndWriteScatterEnergies

!----------------------------------------------------------------------------
  subroutine getTotalEnergy(exportDir, eTot)

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    character(len=300), intent(in) :: exportDir
      !! Path to export directory

    ! Output variables:
    real(kind=dp), intent(out) :: eTot
      !! Total energy read from the `input` file

    ! Local variables:
    character(len=300) :: line
      !! Line from file


    if(ionode) then
      open(30,file=trim(exportDir)//'/input')
      line = getFirstLineWithKeyword(30, 'Total Energy')
      read(30,*) eTot
      close(30)
    endif

    call MPI_BCAST(eTot, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine getTotalEnergy

!----------------------------------------------------------------------------
  subroutine getRefEig(isp, refBand, exportDirEigs, refEig)

    use miscUtilities, only: ignoreNextNLinesFromFile, int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: refBand
      !! Band of WZP reference carrier

    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues

    ! Output variables:
    real(kind=dp), intent(out) :: refEig
      !! Eigenvalue of WZP reference carrier

    ! Local variables:
    integer :: iDum
      !! Dummy integer to ignore band index

    real(kind=dp) :: rDum
      !! Dummy real to ignore occupation

    character(len=300) :: fName
      !! File name


    if(ionode) then

      fName = trim(exportDirEigs)//"/eigenvalues."//trim(int2str(isp))//".1"
        !! Use first k-point as reference and assume that initial
        !! state is not spin polarized

      open(72, file=fName)

      call ignoreNextNLinesFromFile(72, 2 + (refBand-1))
        ! Ignore header and all bands before reference band
    
      read(72, '(i10, 2ES24.15E3)') iDum, refEig, rDum
    
      close(72)

    endif

    call MPI_BCAST(refEig, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine getRefEig

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
    integer :: iDum
      !! Dummy integer to ignore band index

    real(kind=dp) :: defectEig
      !! Eigenvalue of the defect for capture
    real(kind=dp) :: rDum
      !! Dummy real to ignore occupation
    real(kind=dp) :: refEig
      !! Eigenvalue of WZP reference carrier (not the same
      !! as `refEig` used in other places in the code; this
      !! is just to determine the correct distance between
      !! the defect and the band states

    character(len=300) :: fName
      !! File name


    if(ionode) then

      fName = trim(exportDirInitInit)//"/eigenvalues."//trim(int2str(isp))//".1"
        !! Use first k-point as reference and assume that initial
        !! state is not spin polarized

      open(72, file=fName)

      call ignoreNextNLinesFromFile(72, 2 + (iBandFinit-1))
        ! Ignore header and all bands before reference band
    
      read(72, '(i10, 2ES24.15E3)') iDum, defectEig, rDum
    
      close(72)


      open(72, file=fName)

      call ignoreNextNLinesFromFile(72, 2 + (refBand-1))
        ! Ignore header and all bands before reference band
    
      read(72, '(i10, 2ES24.15E3)') iDum, refEig, rDum
    
      close(72)


      !> Switch the order of the eigenvalue subtraction
      !> for hole vs electron capture to represent that 
      !> the actual energy we need is that of the electron
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
  subroutine writeEnergyTable(iBandIInit, iBandIFinal, iBandFInit, iBandFFinal, ikLocal, isp, dEEigRefDefect, eCorrectTot, eCorrectEigRef, &
      eTotInitInit, eTotFinalInit, eTotFinalFinal, refEig, energyTableDir, exportDirEigs, captured, elecCarrier)
  
    implicit none
    
    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    real(kind=dp), intent(in) :: dEEigRefDefect
      !! Eigenvalue difference between the reference
      !! eigenvalue and the defect level
    real(kind=dp), intent(in) :: eCorrectTot
      !! Total-energy correction, if any
    real(kind=dp), intent(in) :: eCorrectEigRef
      !! Correction to eigenvalue difference with reference carrier, if any
    real(kind=dp), intent(in) :: eTotInitInit
      !! Total energy of the relaxed initial charge
      !! state (initial positions)
    real(kind=dp), intent(in) :: eTotFinalInit
      !! Total energy of the unrelaxed final charge
      !! state (initial positions)
    real(kind=dp), intent(in) :: eTotFinalFinal
      !! Total energy of the relaxed final charge
      !! state (final positions)
    real(kind=dp), intent(in) :: refEig
      !! Eigenvalue of WZP reference carrier

    character(len=300), intent(in) :: energyTableDir
      !! Path to energy tables
    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues

    logical, intent(in) :: captured
      !! If carrier is captured as opposed to scattered
    logical, intent(in) :: elecCarrier
      !! If carrier is electron as opposed to hole

    ! Local variables:
    integer :: ibi, ibf
      !! Loop indices
    integer :: ikGlobal
      !! Current global k-point
    integer :: totalNumberOfElements
      !! Total number of overlaps to output

    real(kind=dp) :: dEDelta
      !! Energy to be used in delta function
    real(kind=dp) :: dEEigRef
      !! Eigenvalue difference from initial to reference
    real(kind=dp) :: dEFirst
      !! Energy to be used in first-order matrix element
    real(kind=dp) :: dEPlot
      !! Energy to be used for plotting
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
    real(kind=dp) :: eigvF(iBandFinit:iBandFfinal)
      !! Final-state eigenvalues
    real(kind=dp) :: eigvI(iBandIinit:iBandIfinal)
      !! Initial-state eigenvalues
    
    character(len = 300) :: text
      !! Text for header



    return

  end subroutine writeEnergyTable
  
!----------------------------------------------------------------------------
  subroutine readEigenvalues(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, exportDirEigs, eigvF, eigvI)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none
    
    ! Input variables
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: isp
      !! Current spin channel

    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues

    ! Output variables:
    real(kind=dp), intent(out) :: eigvF(iBandFinit:iBandFfinal)
      !! Final-state eigenvalues
    real(kind=dp), intent(out) :: eigvI(iBandIinit:iBandIfinal)
      !! Initial-state eigenvalues

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

    call ignoreNextNLinesFromFile(72, 2 + (iBandIinit-1))
      ! Ignore header and all bands before lowest initial-state band
    
    do ib = iBandIinit, iBandIfinal
      read(72, '(i10, 2ES24.15E3)') iDum, eigvI(ib), rDum
    enddo
    
    close(72)
    

    open(72, file=fName)

    call ignoreNextNLinesFromFile(72, 2 + (iBandFinit-1))
      ! Ignore header and all bands before lowest final-state band
    
    do ib = iBandFinit, iBandFfinal
      read(72, '(i10, 2ES24.15E3)') iDum, eigvF(ib), rDum
    enddo
    
    close(72)


    return
    
  end subroutine readEigenvalues
  
!----------------------------------------------------------------------------
  subroutine readEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, energyTableDir, dE)
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

    character(len=300) :: energyTableDir
      !! Path to energy table

    ! Output variables:
    real(kind=dp), intent(out) :: dE(iBandFinit:iBandFfinal,iBandIinit:iBandIFinal,4)
      !! All energy differences from energy table

    ! Local variables:
    integer :: iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_
      !! Energy band bounds for initial and final state from the energy table file
    integer :: ibi, ibf
      !! Loop indices
    integer :: iDum
      !! Dummy integer to ignore input

    real(kind=dp) :: dE_(4)
      !! Input energy

    character(len=300) :: fName
      !! Energy table file name
    
    
    fName = trim(energyTableDir)//"/energyTable."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))

    open(27, file=trim(fName), status='unknown')

    read(27,*)
    read(27,'(5i10)') iDum, iBandIinit_, iBandIfinal_, iBandFinit_, iBandFfinal_

    ! Check the input band bounds against those in the energy file
    if(iBandIinit < iBandIinit_ .or. iBandIfinal > iBandIfinal_ .or. iBandFinit < iBandFinit_ .or. iBandFfinal > iBandFfinal_) &
      call exitError('readEnergyTable', 'given band bounds are outside those in energy table '//trim(fName), 1)
    
    call ignoreNextNLinesFromFile(27,6)
    
    do ibf = iBandFinit_, iBandFfinal_
      do ibi = iBandIinit_, iBandIfinal_
      
        read(27,*) iDum, iDum, dE_

        if(ibi >= iBandIinit .and. ibi <= iBandIfinal .and. ibf >= iBandFinit .and. ibf <= iBandFfinal) dE(ibf,ibi,:) = dE_(:)
          
      enddo
    enddo
    
    close(27)
    
    return
    
  end subroutine readEnergyTable

end module energyTabulatorMod
