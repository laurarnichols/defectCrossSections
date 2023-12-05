module energyTabulatorMod
  
  use constants, only: dp, eVToHartree
  use miscUtilities, only: int2str
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nKPoints, nSpins
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

  real(kind=dp) :: dEEigRefDefect
    !! Eigenvalue difference between the reference
    !! eigenvalue and the defect level
  real(kind=dp) :: eCorrectTot
    !! Total-energy correction, if any
  real(kind=dp) :: eCorrectEigRef
    !! Correction to eigenvalue difference with reference carrier, if any
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

  namelist /inputParams/ exportDirEigs, exportDirFinalFinal, exportDirFinalInit, exportDirInitInit, energyTableDir, &
                         eCorrectTot, eCorrectEigRef, captured,  &
                         iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, refBand


  contains

!----------------------------------------------------------------------------
  subroutine initialize(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, refBand, eCorrectTot, eCorrectEigRef, energyTableDir, &
        exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured)
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


    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    refBand = -1

    eCorrectTot = 0.0_dp
    eCorrectEigRef = 0.0_dp

    energyTableDir = './'
    exportDirEigs = ''
    exportDirInitInit = ''
    exportDirFinalInit = ''
    exportDirFinalFinal = ''

    captured = .false.

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, refBand, eCorrectTot, eCorrectEigRef, energyTableDir, &
      exportDirEigs, exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, captured)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
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

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9))
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 
    abortExecution = checkIntInitialization('refBand', refBand, 1, int(1e9)) .or. abortExecution

    write(*,'("eCorrectTot = ", f8.4, " (eV)")') eCorrectTot
    write(*,'("eCorrectEigRef = ", f8.4, " (eV)")') eCorrectEigRef

    if(captured .and. iBandFinit /= iBandFfinal) then
      write(*,'("Capture only expected for a single final-state band.")')
      abortExecution = .true.
    endif
    write(*,'("captured = ",L)') captured

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
    
      read(50,*) 
      read(50,*) 
      read(50,*) 
      read(50, '(i10)') nSpins
      read(50,*) 
      read(50, '(i10)') nKPoints

    endif

    call MPI_BCAST(nSpins, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
   
    return

  end subroutine getnSpinsAndnKPoints

!----------------------------------------------------------------------------
  subroutine getTotalEnergies(exportDirInitInit, exportDirFinalInit, exportDirFinalFinal, eTotInitInit, eTotFinalInit, eTotFinalFinal)

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    character(len=300), intent(in) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions
    character(len=300), intent(in) :: exportDirFinalInit
      !! Path to export for final charge state
      !! in the initial positions
    character(len=300), intent(in) :: exportDirFinalFinal
      !! Path to export for final charge state
      !! in the final positions

    ! Output variables:
    real(kind=dp), intent(out) :: eTotInitInit
      !! Total energy of the relaxed initial charge
      !! state (initial positions)
    real(kind=dp), intent(out) :: eTotFinalInit
      !! Total energy of the unrelaxed final charge
      !! state (initial positions)
    real(kind=dp), intent(out) :: eTotFinalFinal
      !! Total energy of the relaxed final charge
      !! state (final positions)

    ! Local variables:
    character(len=300) :: line
      !! Line from file


    if(ionode) then
      open(30,file=trim(exportDirFinalFinal)//'/input')
      line = getFirstLineWithKeyword(30, 'Total Energy')
      read(30,*) eTotFinalFinal
      close(30)
        !! Get the total energy of the relaxed final charge state
        !! (final positions)

      open(30,file=trim(exportDirFinalInit)//'/input')
      line = getFirstLineWithKeyword(30, 'Total Energy')
      read(30,*) eTotFinalInit
      close(30)
        !! Get the total energy of the unrelaxed final charge state
        !! (initial positions)

      open(30,file=trim(exportDirInitInit)//'/input')
      line = getFirstLineWithKeyword(30, 'Total Energy')
      read(30,*) eTotInitInit
      close(30)
        !! Get the total energy of the unrelaxed final charge state
        !! (initial positions)

    endif

    call MPI_BCAST(eTotFinalFinal, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eTotFinalInit, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eTotInitInit, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine getTotalEnergies

!----------------------------------------------------------------------------
  subroutine getRefEig(exportDirEigs, refEig)

    use miscUtilities, only: ignoreNextNLinesFromFile

    implicit none

    ! Input variables:
    character(len=300), intent(in) :: exportDirEigs
      !! Path to export for system to get eigenvalues

    ! Output variables:
    real(kind=dp), intent(out) :: refEig
      !! Eigenvalue of WZP reference carrier

    ! Local variables:
    character(len=300) :: fName
      !! File name


    if(ionode) then

      fName = trim(exportDirEigs)//"/eigenvalues.1.1"
        !! Use first k-point as reference and assume that initial
        !! state is not spin polarized

      open(72, file=fName)

      call ignoreNextNLinesFromFile(72, 2 + (refBand-1))
        ! Ignore header and all bands before reference band
    
      read(72, '(ES24.15E3)') refEig
    
      close(72)

    endif

    call MPI_BCAST(refEig, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine getRefEig

!----------------------------------------------------------------------------
  subroutine getRefToDefectEigDiff(iBandFinit, isp, refEig, exportDirInitInit, dEEigRefDefect)

    use miscUtilities, only: ignoreNextNLinesFromFile

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandFinit
      !! Energy band for final state
    integer, intent(in) :: isp
      !! Current spin channel

    real(kind=dp), intent(in) :: refEig
      !! Eigenvalue of WZP reference carrier

    character(len=300), intent(in) :: exportDirInitInit
      !! Path to export for initial charge state
      !! in the initial positions

    ! Output variables:
    real(kind=dp), intent(out) :: dEEigRefDefect
      !! Eigenvalue difference between the reference
      !! eigenvalue and the defect level

    ! Local variables:
    real(kind=dp) :: defectEig
      !! Eigenvalue of the defect for capture

    character(len=300) :: fName
      !! File name


    if(ionode) then

      fName = trim(exportDirInitInit)//"/eigenvalues.1.1"
        !! Use first k-point as reference and assume that initial
        !! state is not spin polarized

      open(72, file=fName)

      call ignoreNextNLinesFromFile(72, 2 + (iBandFinit-1))
        ! Ignore header and all bands before reference band
    
      read(72, '(ES24.15E3)') defectEig
    
      close(72)

      dEEigRefDefect = abs(refEig - defectEig)

    endif

    call MPI_BCAST(dEEigRefDefect, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine getRefToDefectEigDiff

!----------------------------------------------------------------------------
  subroutine writeEnergyTable(iBandIInit, iBandIFinal, iBandFInit, iBandFFinal, ikLocal, isp, dEEigRefDefect, eCorrectTot, eCorrectEigRef, &
      eTotInitInit, eTotFinalInit, eTotFinalFinal, refEig, energyTableDir, exportDirEigs, captured)
  
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

    ! Local variables:
    integer :: ibi, ibf
      !! Loop indices
    integer :: ikGlobal
      !! Current global k-point
    integer :: totalNumberOfElements
      !! Total number of overlaps to output

    real(kind=dp) :: dE
      !! Energy difference to be output
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
    real(kind=dp) :: t1, t2
      !! Timers
    
    character(len = 300) :: text
      !! Text for header


    ikGlobal = ikLocal+ikStart-1
    
    call cpu_time(t1)
    
    write(*, '(" Writing energy table of k-point ", i2, " and spin ", i1, ".")') ikGlobal, isp
    
    open(17, file=trim(energyTableDir)//"/energyTable."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)), status='unknown')
    
    text = "# Total number of <f|i> elements, Initial States (bandI, bandF), Final States (bandI, bandF)"
    write(17,'(a, " Format : ''(5i10)''")') trim(text)
    
    totalNumberOfElements = (iBandIfinal - iBandIinit + 1)*(iBandFfinal - iBandFinit + 1)
    write(17,'(5i10)') totalNumberOfElements, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    

    dETotWRelax = eTotFinalFinal - eTotInitInit + eCorrectTot
      !! Get the total energy difference between the two charge states, 
      !! including the atomic relaxation energy (with a potential energy
      !! correction defined by the user). This dE is used in the delta 
      !! function.

    write(17,'("# Total-energy difference (Hartree). Format: ''(ES24.15E3)''")') 
    write(17,'("# With relaxation (for delta function)")')
    write(17,'(ES24.15E3)') dETotWRelax


    dETotElecOnly = eTotFinalInit - eTotInitInit + eCorrectTot
      !! Get the total energy difference between the two charge states, 
      !! not including atomic relaxation (with a potential energy correction
      !! defined by the user). This dE represents the total electronic-only
      !! energy difference between the two charge states and goes in the
      !! zeroth-order matrix element.

    write(17,'("# Electronic only without relaxation (for zeroth-order matrix element)")')
    write(17,'(ES24.15E3)') dETotElecOnly
    

    text = "# Final Band, Initial Band, Delta Function (Hartree), Zeroth-order (Hartree), First-order (Hartree), Plotting (eV)" 
    write(17, '(a, " Format : ''(2i10,4ES24.15E3)''")') trim(text)


    call readEigenvalues(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, exportDirEigs, eigvF, eigvI)

    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal

        dEEigRef = abs(eigvI(ibi) - refEig + eCorrectEigRef)
          !! All of the energies require an eigenvalue difference 
          !! from a reference band. Calculate this once to be used
          !! for all of the energies.
          !!
          !! In both hole and electron capture, the actual electron
          !! energy decreases, so the negative absolute value of the
          !! eigenvalue difference is used.
          !!
          !! The energy correction `eCorrectEigRef` should be zero if
          !! the reference state and initial state are both in the conduction
          !! band or both in the valence band, since eigenvalue differences
          !! within the bands are okay at the PBE level and do not need
          !! to be corrected. One example where this correction would 
          !! be needed is setting up the negative charge state of the
          !! triply-hydrogenated Si defect. The simplest treatment of that
          !! charge state has the reference carrier in the valence band top,
          !! so the distance between the valence band and the conduction band
          !! would need to be corrected from PBE to HSE.

        dEDelta = dETotWRelax - dEEigRef
          !! To get the total energy that needs to be conserved (what
          !! goes into the delta function), add the total energy 
          !! difference between the two relaxed charge states and the
          !! additional eigenvalue energy difference between the initial
          !! state and the WZP reference-carrier state. 

        dEZeroth = dETotElecOnly - dEEigRef
          !! The zeroth-order matrix element contains the electronic-only
          !! energy difference. We get that from a total energy difference
          !! between the two charge states in the initial positions. Like
          !! in the energy for the delta function, the additional carrier
          !! energy must also be included with a potential correction.

        if(captured) then
          dEFirst = dEEigRef + dEEigRefDefect
        else
          dEFirst = abs(eigvI(ibi) - eigvF(ibf))
        endif
          !! The first-order term contains only the unperturbed eigenvalue
          !! difference. The perturbative expansion has 
          !! \(\varepsilon_i - \varepsilon_f\), in terms of the actual 
          !! electron. The absolute value is needed for the hole case.
          !! 
          !! For capture, the defect level should not have dispersion, and
          !! the level of the defect changes between charge states, so we
          !! use a single reference energy and distance to the defect, taken
          !! from the initial-charge-state HSE calculation. For the scattering
          !! case, the bands do not have those issues, so we can directly 
          !! use the final-state eigenvalue.

        dEPlot = abs(dEZeroth)/eVToHartree
          !! Energy plotted should be the positive, electronic-only energy 
          !! difference (same as in zeroth-order) in eV
        
        write(17, 1001) ibf, ibi, dEDelta, dEZeroth, dEFirst, dEPlot
            
      enddo
    enddo

    
    close(17)
    
    call cpu_time(t2)
    write(*, '(" Writing energy table of k-point ", i4, "and spin ", i1, " done in:                   ", f10.2, " secs.")') &
      ikGlobal, isp, t2-t1
    
 1001 format(2i7,4ES24.15E3)

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

    character(len=300) :: fName
      !! File name

    
    fName = trim(exportDirEigs)//"/eigenvalues."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))

    open(72, file=fName)

    call ignoreNextNLinesFromFile(72, 2 + (iBandIinit-1))
      ! Ignore header and all bands before lowest initial-state band
    
    do ib = iBandIinit, iBandIfinal
      read(72, '(ES24.15E3)') eigvI(ib)
    enddo
    
    close(72)
    

    open(72, file=fName)

    call ignoreNextNLinesFromFile(72, 2 + (iBandFinit-1))
      ! Ignore header and all bands before lowest final-state band
    
    do ib = iBandFinit, iBandFfinal
      read(72, '(ES24.15E3)') eigvF(ib)
    enddo
    
    close(72)


    return
    
  end subroutine readEigenvalues
  
!----------------------------------------------------------------------------
  subroutine readEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, order, energyTableDir, dE)
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
    integer, intent(in) :: order
      !! Order of matrix element (0 or 1)

    character(len=300) :: energyTableDir
      !! Path to energy table

    ! Output variables:
    real(kind=dp), intent(out) :: dE(iBandFinit:iBandFfinal,iBandIinit:iBandIFinal,4)
      !! All energy differences from energy table

    ! Local variables:
    integer :: ibi, ibf
      !! Loop indices
    integer :: iDum
      !! Dummy integer to ignore input
    
    
    open(27, file=trim(energyTableDir)//"/energyTable."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)), status='unknown')
    
    call ignoreNextNLinesFromFile(27,8)
    
    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal
      
        read(27,*) iDum, iDum, dE(ibf,ibi,:)
          
      enddo
    enddo
    
    close(27)
    
    return
    
  end subroutine readEnergyTable

end module energyTabulatorMod
