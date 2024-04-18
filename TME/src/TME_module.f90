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


  real(kind=dp) :: dq_j
    !! \(\delta q_j) for displaced wave functions
    !! (only order = 1)
  real(kind=dp), allocatable :: gCart(:,:)
    !! G-vectors in Cartesian coordinates
  real(kind=dp) :: omega
    !! Supercell volume
  real(kind = dp) :: t1, t2
    !! For timing different processes


  integer :: maxAngMom
    !! Maximum angular momentum of the projectors
  integer, allocatable :: mill_local(:,:)
    !! Local Miller indices
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors
  integer :: nGkVecsLocalPC, nGkVecsLocalSD
    !! Local number of G+k vectors on this processor
  integer :: phononModeJ
    !! Index of phonon mode for the calculation
    !! of \(M_j\) (only for order=1)

  complex(kind=dp), allocatable :: betaPC(:,:), betaSD(:,:)
    !! Projectors
  complex(kind=dp), allocatable :: Ufi(:,:,:,:)
    !! All-electron overlap
  complex*8, allocatable :: wfcPC(:,:), wfcSD(:,:)
    !! Wave function coefficients
  complex(kind=dp), allocatable :: Ylm(:,:)
    !! Spherical harmonics
  
  character(len=300) :: baselineDir
    !! Directory for baseline overlap to optionally
    !! be subtracted for the first-order term
  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: exportDirSD, exportDirPC
    !! Paths to exports for SD (left) and PC (right) 
    !! systems to use for wave function overlaps <SD|PC>
  character(len=300) :: outputDir
    !! Path to where matrix elements should be output

  logical :: subtractBaseline
    !! If baseline should be subtracted from first-order
    !! overlap for increased numerical accuracy in the 
    !! derivative

  ! Define a type to match the export code for the variables
  ! that come straight from the POTCAR file
  type potcar
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

    real(kind=dp), allocatable :: dRadGrid(:)
      !! Derivative of radial grid
    real(kind=dp), allocatable :: radGrid(:)
      !! Radial grid points
    real(kind=dp), allocatable :: wae(:,:)
      !! AE wavefunction
    real(kind=dp), allocatable :: wps(:,:)
      !! PS wavefunction

    character(len=2) :: element
  end type potcar

  ! Define a type for new items calculated for PAW
  type :: calcPAW
    real(kind=dp), allocatable :: F(:,:)
    real(kind=dp), allocatable :: F1(:,:,:)
    real(kind=dp), allocatable :: F2(:,:,:)
    real(kind=dp), allocatable :: bes_J_qr(:,:,:)
  end type calcPAW

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

    character(len=300) :: ID
      !! How this system should be identified in output
    character(len=3) :: sysType
      !! Type of system for overlap:
      !! bra or ket

    type(potcar), allocatable :: pot(:)
      !! Pseudopotential information from this system

    type(calcPAW), allocatable :: paw(:)
      !! Values needed for PAW formalism
  end type crystal

  type(crystal) :: braSys, ketSys
    !! Define variables for the system used as the
    !! <bra| and the |ket> in the overlap matrix 
    !! element

  character(len = 300) :: textDum
  character(len = 320) :: mkdir
  
  integer :: ik, ig, ibi, ibf
  integer :: iTypes, iPn
  integer :: numOfUsedGvecsPP, npwNi, npwNf, npwMi, npwMf
  integer :: np, nI, nF, nPP, ind2
  integer :: i, j, n1, n2, n3, n4, n, id, npw
  
  real(kind = dp) t0, tf
  
  
  real(kind = dp), allocatable :: gvecs(:,:)
  real(kind = dp), allocatable :: DE(:,:,:), absVfi2(:,:,:)
  
  complex(kind = dp), allocatable :: paw_SDKKPC(:,:), paw_id(:,:)
  complex(kind = dp), allocatable :: pawKPC(:,:,:), pawSDK(:,:,:), pawPsiPC(:,:), pawSDPhi(:,:)
  complex(kind = dp), allocatable :: cProjPC(:,:), cProjSD(:,:)
  complex(kind = dp), allocatable :: paw_PsiPC(:,:), paw_SDPhi(:,:)
  complex(kind = dp), allocatable :: cProjBetaPCPsiSD(:,:)
  complex(kind = dp), allocatable :: cProjBetaSDPhiPC(:,:)
  
  integer, allocatable :: igvs(:,:,:), pwGvecs(:,:), iqs(:)
  integer, allocatable :: pwGs(:,:), nIs(:,:), nFs(:,:), ngs(:,:)
  
  type :: vec
    integer :: ind
    integer, allocatable :: igN(:), igM(:)
  end type vec
  
  TYPE(vec), allocatable :: vecs(:), newVecs(:)
  

  namelist /TME_Input/ exportDirSD, exportDirPC, outputDir, energyTableDir, &
                       iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, &
                       order, dqFName, phononModeJ, subtractBaseline, baselineDir, &
                       ispSelect
  
  
contains

!----------------------------------------------------------------------------
  subroutine readInput(ispSelect, maxAngMom, nGVecsGlobal, nKPoints, nSpins, omega, baselineDir, loopSpins, subtractBaseline)

    use miscUtilities, only: getFirstLineWithKeyword, ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: maxAngMom
      !! Maximum angular momentum of the projectors
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors
    integer, intent(out) :: nKPoints
      !! Total number of k-points
    integer, intent(inout) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    real(kind=dp) :: omega
      !! Cell volume (tested to be consistent
      !! across all systems)

    character(len=300), intent(out) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term

    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel
    logical, intent(out) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative

    ! Local variables:    
    integer :: iDum
      !! Dummy integer to ignore file input
    

    if(ionode) then
    
      call initialize(ispSelect, baselineDir, subtractBaseline)
    
      read(5, TME_Input, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputParams', 'reading TME_Input namelist', abs(ierr))
    
      call checkInitialization(ispSelect, baselineDir, subtractBaseline, loopSpins)

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

    call MPI_BCAST(exportDirSD, len(exportDirSD), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirPC, len(exportDirPC), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(baselineDir, len(baselineDir), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(outputDir, len(outputDir), MPI_CHARACTER, root, worldComm, ierr)


    ! Initialize global variables to be tracked to make sure they
    ! are consistent across all of the systems
    maxAngMom = 0
    nKPoints = -1
    nGVecsGlobal = -1
    omega = -1.0

    braSys%sysType = 'bra'
    braSys%ID = 'bra'
    call readInputFile(exportDirPC, maxAngMom, nGVecsGlobal, nKPoints, nSpins, omega, braSys)

    ketSys%sysType = 'ket'
    ketSys%ID = 'ket'
    call readInputFile(exportDirSD, maxAngMom, nGVecsGlobal, nKPoints, nSpins, omega, ketSys)

    
    maxAngMom = 2*maxAngMom + 1

    
    if(ionode) then
      
      if(order == 1) then

        open(30,file=trim(dqFName))
        call ignoreNextNLinesFromFile(30, 1+phononModeJ-1)
          ! Ignore header and all modes before phononModeJ
        read(30,*) iDum, dq_j
        close(30)

      endif

    endif

    if(order == 1) then
      call MPI_BCAST(dq_j, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    endif

    nSpins = max(braSys%nSpins,ketSys%nSpins)
    
    return
    
  end subroutine readInput
  
!----------------------------------------------------------------------------
  subroutine initialize(ispSelect, baselineDir, subtractBaseline)
    
    implicit none

    ! Output variables:
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user

    character(len=300), intent(out) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term

    logical, intent(out) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative
    

    exportDirSD = ''
    exportDirPC = ''
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
  subroutine checkInitialization(ispSelect, baselineDir, subtractBaseline, loopSpins)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user

    character(len=300), intent(in) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term

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
    
    abortExecution = checkDirInitialization('exportDirSD', exportDirSD, 'input')
    abortExecution = checkDirInitialization('exportDirPC', exportDirPC, 'input') .or. abortExecution
    abortExecution = checkDirInitialization('energyTableDir', energyTableDir, 'energyTable.1.1') .or. abortExecution
    abortExecution = checkIntInitialization('order', order, 0, 1) .or. abortExecution

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

    if(ispSelect < 1 .or. ispSelect > 2) then
      write(*,*) "No valid choice for spin channel selection given. Looping over spin."
      loopSpins = .true.
    else
      write(*,'("Only exporting spin channel ", i2)') ispSelect
      loopSpins = .false.
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
  subroutine readInputFile(exportDir, maxAngMom, nGVecsGlobal, nKPoints, nSpins, omega, sys)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Input variables:
    character(len=300) :: exportDir
      !! Path to Export dir

    ! Output variables:
    integer, intent(inout) :: maxAngMom
      !! Maximum angular momentum of projectors across
      !! all systems
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
    integer :: ik, iA, ix, iT, ip, ir, ip1, ip2, irm
      !! Loop indices
    integer :: nBands
      !! Number of bands

    real(kind=dp), allocatable :: aepsDiff1(:), aepsDiff2(:)
      !! Difference between wae and wps for different channels
    real(kind=dp) :: t1, t2 
      !! Timers
    real(kind=dp) :: rDum
      !! Dummy real variable

    character(len=300) :: inputFName
      !! File name for the input file 
    character(len=300) :: textDum
      !! Dummy text
    
    
    if(ionode) then
      call cpu_time(t1)
    
      inputFName = trim(trim(exportDir)//'/input')
    
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
    
      read(50,*) 
      read(50,'(i10)') sys%nKPoints


      if(nKPoints < 0) then
        nKPoints = sys%nKPoints
      else if(sys%nKPoints /= nKPoints) then
        call exitError('readInput', 'number of k-points in system '//trim(sys%ID)//' does not match', 1)
      end if
    endif

    call MPI_BCAST(sys%nSpins, 1, MPI_INTEGER, root, worldComm, ierr)
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
    
    allocate(sys%pot(sys%nAtomTypes), sys%paw(sys%nAtomTypes), sys%nAtomsEachType(sys%nAtomTypes))

    sys%nProj = 0
    do iT = 1, sys%nAtomTypes
      
      if(ionode) then

        read(50,*) 
        read(50,*) !sys%pot(iT)%element
          ! Don't use this now, but might, so I want
          ! to easily know where it is.
      
        read(50,*)
        read(50,'(i10)') sys%nAtomsEachType(iT)

        read(50,*)
        read(50,'(i10)') sys%pot(iT)%nChannels              ! number of projectors

      endif

      call MPI_BCAST(sys%pot(iT)%nChannels, 1, MPI_INTEGER, root, worldComm, ierr)
      
      allocate(sys%pot(iT)%angMom(sys%pot(iT)%nChannels))

      if(ionode) then

        read(50,*)

        do ip = 1, sys%pot(iT)%nChannels

          read(50,'(2i10)') sys%pot(iT)%angMom(ip), iDum

        enddo

      endif

      call MPI_BCAST(sys%pot(iT)%angMom, sys%pot(iT)%nChannels, MPI_INTEGER, root, worldComm, ierr)

      if(ionode) then
      
        read(50,*)
        read(50,'(i10)') sys%pot(iT)%lmMax
      
        read(50,*)
        read(50,'(2i10)') sys%pot(iT)%nMax, sys%pot(iT)%iRAugMax

      endif
    
      call MPI_BCAST(sys%pot(iT)%lmMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(sys%pot(iT)%nMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(sys%pot(iT)%iRAugMax, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(sys%pot(iT)%radGrid(sys%pot(iT)%nMax))
      
      if(ionode) then
      
        allocate(sys%pot(iT)%dRadGrid(sys%pot(iT)%nMax))

        read(50,*)

        do ir = 1, sys%pot(iT)%nMax
          read(50,'(2ES24.15E3)') sys%pot(iT)%radGrid(ir), sys%pot(iT)%dRadGrid(ir)
        enddo
       
        allocate(sys%pot(iT)%wae(sys%pot(iT)%nMax, sys%pot(iT)%nChannels))
        allocate(sys%pot(iT)%wps(sys%pot(iT)%nMax, sys%pot(iT)%nChannels))
      
        read(50,*)
        do ip = 1, sys%pot(iT)%nChannels
          do ir = 1, sys%pot(iT)%nMax
            read(50,'(2ES24.15E3)') sys%pot(iT)%wae(ir,ip), sys%pot(iT)%wps(ir,ip) 
          enddo
        enddo
        
      endif

      allocate(sys%paw(iT)%F(sys%pot(iT)%iRAugMax, sys%pot(iT)%nChannels))
      allocate(sys%paw(iT)%F1(sys%pot(iT)%iRAugMax, sys%pot(iT)%nChannels, sys%pot(iT)%nChannels))
      allocate(sys%paw(iT)%F2(sys%pot(iT)%iRAugMax, sys%pot(iT)%nChannels, sys%pot(iT)%nChannels))

      
      if(ionode) then
        
        irm = sys%pot(iT)%iRAugMax
        allocate(aepsDiff1(irm), aepsDiff2(irm))

        sys%paw(iT)%F = 0.0_dp
        sys%paw(iT)%F1 = 0.0_dp
        sys%paw(iT)%F2 = 0.0_dp
      
        do ip1 = 1, sys%pot(iT)%nChannels

          aepsDiff1 = sys%pot(iT)%wae(1:irm,ip1) - sys%pot(iT)%wps(1:irm,ip1)

          sys%paw(iT)%F(1:irm,ip1)= aepsDiff1(:)*sys%pot(iT)%radGrid(1:irm)*sys%pot(iT)%dRadGrid(1:irm)
        
          do ip2 = 1, sys%pot(iT)%nChannels

            aepsDiff2 = sys%pot(iT)%wae(1:irm,ip2) - sys%pot(iT)%wps(1:irm,ip2)

            if(trim(sys%sysType) == 'bra') then
              sys%paw(iT)%F1(1:irm,ip2,ip1) = sys%pot(iT)%wps(1:irm,ip2)*aepsDiff1(:)*sys%pot(iT)%dRadGrid(1:irm)
            else if(trim(sys%sysType) == 'ket') then
              sys%paw(iT)%F1(1:irm,ip2,ip1) = sys%pot(iT)%wps(1:irm,ip1)*aepsDiff2(:)*sys%pot(iT)%dRadGrid(1:irm)
            endif
          
            sys%paw(iT)%F2(1:irm,ip2,ip1) = aepsDiff1(:)*aepsDiff2(:)*sys%pot(iT)%dRadGrid(1:irm)

          enddo
        enddo

        deallocate(sys%pot(iT)%wae)
        deallocate(sys%pot(iT)%wps)
        deallocate(sys%pot(iT)%dRadGrid)
        deallocate(aepsDiff1, aepsDiff2)
      
        sys%nProj = sys%nProj + sys%nAtomsEachType(iT)*sys%pot(iT)%lmMax

      endif

      call MPI_BCAST(sys%pot(iT)%radGrid, size(sys%pot(iT)%radGrid), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(sys%paw(iT)%F, size(sys%paw(iT)%F), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(sys%paw(iT)%F1, size(sys%paw(iT)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(sys%paw(iT)%F2, size(sys%paw(iT)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      
    enddo
  
    call MPI_BCAST(sys%nAtomsEachType, sys%nAtomTypes, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(sys%nProj, 1, MPI_INTEGER, root, worldComm, ierr)

    if(ionode) then
    
      close(50)

      do iT = 1, sys%nAtomTypes
        do ip = 1, sys%pot(iT)%nChannels
          if(sys%pot(iT)%angMom(ip) > maxAngMom) maxAngMom = sys%pot(iT)%angMom(ip)
        enddo
      enddo
    
      call cpu_time(t2)
      write(*,'(" Reading input files done in:                ", f10.2, " secs.")') t2-t1
      write(*,*)

    endif

    call MPI_BCAST(maxAngMom, 1, MPI_INTEGER, root, worldComm, ierr)
    
    return
    
  end subroutine readInputFile

!----------------------------------------------------------------------------
  subroutine getFullPWGrid(iGStart_pool, nGVecsLocal, nGVecsGlobal, mill_local)
    !! Read full PW grid from mgrid file
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iGStart_pool
      !! Start G-vector for each process in pool
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors

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

      open(72, file=trim(exportDirSD)//"/mgrid")
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
  subroutine setUpTables(maxAngMom, nGVecsLocal, mill_local, ketSys, gCart, Ylm)

    implicit none

    ! Input variables:
    integer, intent(in) :: maxAngMom
      !! Max J index from L and M
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors
    integer, intent(in) :: mill_local(3,nGVecsLocal)
      !! Miller indices for local G-vectors

    type(crystal) :: ketSys
       !! The ket crystal system used to set up tables

    ! Output variables:
    real(kind=dp), intent(out) :: gCart(3,nGVecsLocal)
      !! G-vectors in Cartesian coordinates

    complex(kind=dp), intent(out) :: Ylm((maxAngMom+1)**2,nGVecsLocal)
      !! Spherical harmonics

    ! Local variables:
    integer :: ig, iT, iR
      !! Loop indices

    real(kind=dp) :: gUnit(3)
      !! Unit G-vector
    real(kind=dp) :: JL(0:maxAngMom)
      !! Bessel_j temporary variable
    real(kind=dp) :: q
      !! Magnitude of G-vector


    Ylm = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    do iT = 1, ketSys%nAtomTypes

      allocate(ketSys%paw(iT)%bes_J_qr( 0:maxAngMom, ketSys%pot(iT)%iRAugMax, nGVecsLocal))
      ketSys%paw(iT)%bes_J_qr(:,:,:) = 0.0_dp
      
    enddo

    do ig = 1, nGVecsLocal

      gCart(:,ig) = matmul(ketSys%recipLattVec, mill_local(:,ig))

      q = sqrt(dot_product(gCart(:,ig),gCart(:,ig)))

      gUnit(:) = gCart(:,ig)
      if(abs(q) > 1.0e-6_dp) gUnit = gUnit/q
        !! Get unit vector for Ylm calculation

      call getYlm(gUnit, maxAngMom, Ylm(:,ig))
        !! Calculate all the needed spherical harmonics

      do iT = 1, ketSys%nAtomTypes

        do iR = 1, ketSys%pot(iT)%iRAugMax

          JL = 0.0_dp
          call bessel_j(q*ketSys%pot(iT)%radGrid(iR), maxAngMom, JL) ! returns the spherical bessel at qr point
            ! Previously used SD atoms structure here for both PC and SD
          ketSys%paw(iT)%bes_J_qr(:,iR,ig) = JL(:)

        enddo

      enddo

    enddo


    return

  end subroutine setUpTables
  
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
  subroutine readProjectors(crystalType, iGkStart_pool, ikGlobal, nGkVecsLocal, sys, beta)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iGkStart_pool
      !! Start G+k index for this process in pool
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD) for projectors

    type(crystal) :: sys
       !! The crystal system

    ! Output variables:
    complex(kind=dp), intent(out) :: beta(nGkVecsLocal,sys%nProj)
      !! Projector of `projCrystalType` with all projectors
      !! and only local PWs/G+k vectors
    
    ! Local variables:
    integer :: reclen
      !! Record length for projectors file
    integer :: igkLocal, igkGlobal, ipr
      !! Loop indices
    
    character(len=300) :: fNameExport
      !! File names
    

    if(crystalType == 'PC') then
      fNameExport = trim(exportDirPC)//"/projectors."//trim(int2str(ikGlobal)) 
    else
      fNameExport = trim(exportDirSD)//"/projectors."//trim(int2str(ikGlobal))
    endif


    inquire(iolength=reclen) beta(1,:)
      !! Get the record length needed to write a double complex
      !! array of length `nProj`

    open(unit=72, file=trim(fNameExport), access='direct', recl=reclen, iostat=ierr, status='old', SHARED)


    do igkLocal = 1, nGkVecsLocal

      igkGlobal = igkLocal+iGkStart_pool-1

      read(72,rec=igkGlobal+1) (beta(igkLocal,ipr), ipr=1,sys%nProj)

    enddo

    close(72)
    
    return 

  end subroutine readProjectors

!----------------------------------------------------------------------------
  subroutine readWfc(crystalType, iBandinit, iBandfinal, iGkStart_pool, ikGlobal, isp, nGkVecsLocal, npws, wfc)
    !! Read wave function for given `crystalType` from `iBandinit`
    !! to `iBandfinal`
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: iBandinit
      !! Starting band
    integer, intent(in) :: iBandfinal
      !! Ending band
    integer, intent(in) :: iGkStart_pool
      !! Start G+k index for this process in pool
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor
    integer, intent(in) :: npws
      !! Number of G+k vectors less than
      !! the cutoff at this k-point

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    ! Output variables
    complex*8, intent(out) :: wfc(nGkVecsLocal,iBandinit:iBandfinal)
      !! Wave function coefficients for 
      !! local G-vectors

    ! Local variables:
    integer :: displacement(nProcPerPool)
      !! Offset from beginning of array for
      !! scattering coefficients to each process
    integer :: reclen
      !! Record length for projectors file
    integer :: sendCount(nProcPerPool)
      !! Number of items to send to each process
      !! in the pool
    integer :: ib, igk
      !! Loop indices

    complex*8 :: wfcAllPWs(npws)
      !! Wave function read from file

    character(len=300) :: fNameExport
      !! File names


    if(crystalType == 'PC') then
      fNameExport = trim(exportDirPC)//'/wfc.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal))
    else
      fNameExport = trim(exportDirSD)//'/wfc.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal))
    endif
    
    wfc(:,:) = cmplx(0.0_dp, 0.0_dp)

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

    if(indexInPool == 0) open(unit=72, file=trim(fNameExport), access='direct', recl=reclen, iostat=ierr, status='old', SHARED)

    do ib = iBandinit, iBandfinal

      if(indexInPool == 0) read(72,rec=ib) (wfcAllPWs(igk), igk=1,npws)

      call MPI_SCATTERV(wfcAllPWs(:), sendCount, displacement, MPI_COMPLEX, wfc(1:nGkVecsLocal,ib), nGkVecsLocal, &
        MPI_COMPLEX, 0, intraPoolComm, ierr)

    enddo

    if(indexInPool == 0) close(72)
    
    return
    
  end subroutine readWfc
  
!----------------------------------------------------------------------------
  subroutine calculatePWsOverlap(ikLocal,isp)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ibi, ibf

    
    Ufi(:,:,ikLocal,isp) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    do ibi = iBandIinit, iBandIfinal 
      
      do ibf = iBandFinit, iBandFfinal

        Ufi(ibf, ibi, ikLocal,isp) = dot_product(wfcSD(:,ibf),wfcPC(:,ibi))
          !! Calculate local overlap
          ! `dot_product` automatically conjugates first argument for 
          ! complex variables.

      enddo
      
    enddo
    
    return
    
  end subroutine calculatePWsOverlap
  
!----------------------------------------------------------------------------
  subroutine readProjections(crystalType, iBandinit, iBandfinal, ikGlobal, isp, nProjs, cProj)
    
    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandinit
      !! Starting band
    integer, intent(in) :: iBandfinal
      !! Ending band
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nProjs
      !! Number of projectors

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    ! Output variables:
    complex(kind=dp), intent(out) :: cProj(nProjs,iBandinit:iBandfinal)
      !! Projections <beta|wfc>

    ! Local variables:
    integer :: ib
      !! Loop indices
    integer :: reclen
      !! Record length for projections files

    character(len=300) :: fNameExport
      !! Export file name
    
    
    cProj(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
    
    if(indexInPool == 0) then

      inquire(iolength=reclen) cProj(:,iBandinit)
        !! Get the record length needed to write a complex
        !! array of length nProjs

      ! Open the projections file for the given crystal type
      if(crystalType == 'PC') then
        fNameExport = trim(exportDirPC)//"/projections."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))
      else
        fNameExport = trim(exportDirSD)//"/projections."//trim(int2str(isp))//"."//trim(int2str(ikGlobal))
      endif


      open(72, file=trim(fNameExport), access='direct', form='unformatted', recl=reclen)
    
    
      ! Read the projections
      do ib = iBandinit, iBandfinal

        read(72,rec=ib) cProj(:,ib)

      enddo
    
      close(72)

    endif

    call MPI_BCAST(cProj, size(cProj), MPI_DOUBLE_COMPLEX, root, intraPoolComm, ierr)
      ! Broadcast entire array to all processes
    
    return
    
  end subroutine readProjections
  
!----------------------------------------------------------------------------
  subroutine calculateCrossProjection(iBandinit, iBandfinal, nGkVecsLocal1, nGkVecsLocal2, nProjs, beta, wfc, crossProjection)
    !! Calculate the cross projection of one crystal's projectors
    !! on the other crystal's wave function coefficients, distributing
    !! the result to all processors
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iBandinit
      !! Starting band for crystal wfc comes from
      !! (not `projCrystalType`)
    integer, intent(in) :: iBandfinal
      !! Ending band for crystal wfc comes from
      !! (not `projCrystalType`)
    integer, intent(in) :: nGkVecsLocal1
      !! Local number of G+k vectors on this processor
      !! for projectors
    integer, intent(in) :: nGkVecsLocal2
      !! Local number of G+k vectors on this processor
      !! for wave function
    integer, intent(in) :: nProjs
      !! Number of projectors

    complex(kind=dp) :: beta(nGkVecsLocal1,nProjs)
      !! Projector of one crystal type
    complex*8, intent(in) :: wfc(nGkVecsLocal2,iBandinit:iBandfinal)
      !! Wave function coefficients for local G-vectors
      !! for other crystal type

    ! Output variables:
    complex(kind=dp), intent(out) :: crossProjection(nProjs,iBandinit:iBandfinal)
      !! Projections <beta|wfc>
    
    ! Local variables:
    integer :: ib, ipr
      !! Loop indices

    complex(kind=dp) :: crossProjectionLocal
      !! Local version of cross projection to
      !! be summed across processors in pool
    

    crossProjection(:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do ib = iBandinit, iBandfinal
      do ipr = 1, nProjs

        crossProjectionLocal = dot_product(beta(:,ipr),wfc(:,ib))
          ! `dot_product` automatically conjugates first argument for 
          ! complex variables.

        call MPI_ALLREDUCE(crossProjectionLocal, crossProjection(ipr,ib), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, intraPoolComm, ierr)

      enddo
    enddo
    
    return
    
  end subroutine calculateCrossProjection
  
!----------------------------------------------------------------------------
  subroutine pawCorrectionWfc(nProjs, cProjI, cProjF, sys, pawWfc)
    ! calculates the augmentation part of the transition matrix element
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nProjs
      !! First index for `cProjI` and `cProjF`

    complex(kind = dp) :: cProjI(nProjs,iBandIinit:iBandIfinal)
      !! Initial-system (PC) projection
    complex(kind = dp) :: cProjF(nProjs,iBandFinit:iBandFfinal)
      !! Final-system (SD) projection

    type(crystal), intent(in) :: sys
       !! The crystal system

    ! Output variables:
    complex(kind=dp), intent(out) :: pawWfc(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal)
      !! Augmentation part of the transition matrix element

    ! Local variables:
    integer :: ia
      !! Loop index

    complex(kind = dp) :: cProjIe
      !! Single element of initial-system (PC)
      !! projection
    complex(kind = dp) :: cProjFe
      !! Single element of final-system (SD)
      !! projection

    integer :: ibi, ibf
    integer :: LL, LLP, LMBASE, LM, LMP
    integer :: L, M, LP, MP, iT
    real(kind = dp) :: atomicOverlap
    
    
    pawWfc(:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    LMBASE = 0
    
    do ia = 1, sys%nAtoms
      
      iT = sys%iType(ia)

      LM = 0
      DO LL = 1, sys%pot(iT)%nChannels
        L = sys%pot(iT)%angMom(LL)
        DO M = -L, L
          LM = LM + 1 !1st index for CPROJ
          
          LMP = 0
          DO LLP = 1, sys%pot(iT)%nChannels
            LP = sys%pot(iT)%angMom(LLP)
            DO MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) then 
                atomicOverlap = sum(sys%paw(iT)%F1(:,LL, LLP))
                
                do ibi = iBandIinit, iBandIfinal
                  cProjIe = cProjI(LMP + LMBASE, ibi)
                  
                  do ibf = iBandFinit, iBandFfinal
                    cProjFe = conjg(cProjF(LM + LMBASE, ibf))
                    
                    pawWfc(ibf, ibi) = pawWfc(ibf, ibi) + cProjFe*atomicOverlap*cProjIe
                    
                  enddo
                  
                enddo
                
              endif
              
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      LMBASE = LMBASE + sys%pot(iT)%lmMax
    ENDDO
    
    return
    
  end subroutine pawCorrectionWfc

!----------------------------------------------------------------------------
  subroutine pawCorrectionK(crystalType, maxAngMom, nGVecsLocal, gCart, Ylm, sys, ketSys, pawK)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: maxAngMom
      !! Max J index from L and M
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors

    real(kind=dp), intent(in) :: gCart(3,nGVecsLocal)
      !! G-vectors in Cartesian coordinates

    complex(kind=dp), intent(in) :: Ylm((maxAngMom+1)**2,nGVecsLocal)
      !! Spherical harmonics

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    type(crystal), intent(in) :: sys
       !! The crystal system
    type(crystal) :: ketSys
       !! The ket crystal system used to set up tables


    ! Output variables:
    complex(kind=dp), intent(out) :: pawK(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal,nGVecsLocal)

    ! Local variables:
    integer :: ig, iT, ni, ibi, ibf, ind
      !! Loop indices
    integer :: LL, LMBASE, LM, L, M
      !! L and M quantum number trackers
    real(kind = dp) :: qDotR, FI
    
    complex(kind = dp) :: VifQ_aug, ATOMIC_CENTER
    
    
    pawK(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    do ig = 1, nGVecsLocal
      
      LMBASE = 0
      
      do ni = 1, sys%nAtoms ! LOOP OVER THE IONS
        
        qDotR = dot_product(gCart(:,ig), sys%atomPositionsCart(:,ni))
        
        if(crystalType == 'PC') then
          ATOMIC_CENTER = exp( -ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        else
          ATOMIC_CENTER = exp( ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        endif
        
        iT = sys%iType(ni)
        LM = 0
        DO LL = 1, sys%pot(iT)%nChannels
          L = sys%pot(iT)%angMom(LL)
          DO M = -L, L
            LM = LM + 1 !1st index for CPROJ
            
            FI = 0.0_dp
            
            FI = dot_product(ketSys%paw(iT)%bes_J_qr(L,:,ig),sys%paw(iT)%F(:,LL)) 
              ! radial part integration F contains dRadGrid
            
            ind = L*(L + 1) + M + 1 ! index for spherical harmonics

            if(crystalType == 'PC') then

              VifQ_aug = ATOMIC_CENTER*Ylm(ind,ig)*(-II)**L*FI

              do ibi = iBandIinit, iBandIfinal

                pawK(:, ibi, ig) = pawK(:, ibi, ig) + VifQ_aug*cProjPC(LM + LMBASE, ibi)
              
              enddo

            else

              VifQ_aug = ATOMIC_CENTER*conjg(Ylm(ind,ig))*(II)**L*FI

              do ibf = iBandFinit, iBandFfinal
                
                pawK(ibf, :, ig) = pawK(ibf, :, ig) + VifQ_aug*conjg(cProjSD(LM + LMBASE, ibf))
                
              enddo
            endif
          ENDDO
        ENDDO

        LMBASE = LMBASE + sys%pot(iT)%lmMax
      ENDDO
      
    enddo
    
    return
    
  end subroutine pawCorrectionK

!----------------------------------------------------------------------------
  subroutine readAndSubtractBaseline(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikLocal, isp, nSpins, Ufi)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nSpins
      !! Number of spin channels

    ! Output variables:
    complex(kind=dp), intent(inout) :: Ufi(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal,nKPerPool,nSpins)
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
          Ufi(ibf,ibi,ikLocal,isp) = Ufi(ibf,ibi,ikLocal,isp) - baselineOverlap
          
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

    complex(kind=dp), intent(in) :: Ufi(iBandFinit:iBandFfinal,iBandIinit:iBandIfinal,nKPerPool,nSpins)
      !! All-electron overlap

    ! Local variables:
    integer :: ibi, ibf, i
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
          write(17, 1001) ibf, ibi, Ufi(ibf,ibi,ikLocal,isp), abs(Ufi(ibf,ibi,ikLocal,isp))**2, abs(dE(ibf,ibi,2)*Ufi(ibf,ibi,ikLocal,isp))**2
        else if(order == 1) then
          write(17, 1001) ibf, ibi, Ufi(ibf,ibi,ikLocal,isp), abs(Ufi(ibf,ibi,ikLocal,isp))**2, abs(dE(ibf,ibi,3)*Ufi(ibf,ibi,ikLocal,isp)/dq_j)**2
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
  subroutine finalizeCalculation()
    
    implicit none

    integer :: iT
      !! Loop index


    deallocate(braSys%nPWs1kGlobal)
    deallocate(braSys%atomPositionsCart)
    deallocate(braSys%iType)
    deallocate(gCart)
    deallocate(Ylm)

    do iT = 1, braSys%nAtomTypes
      deallocate(braSys%pot(iT)%radGrid)
      deallocate(braSys%pot(iT)%angMom)
      deallocate(braSys%paw(iT)%F)
      deallocate(braSys%paw(iT)%F1)
      deallocate(braSys%paw(iT)%F2)
    enddo

    deallocate(braSys%pot)
    deallocate(braSys%paw)

    deallocate(ketSys%nPWs1kGlobal)
    deallocate(ketSys%atomPositionsCart)
    deallocate(ketSys%iType)

    do iT = 1, ketSys%nAtomTypes
      deallocate(ketSys%pot(iT)%radGrid)
      deallocate(ketSys%pot(iT)%angMom)
      deallocate(ketSys%paw(iT)%F)
      deallocate(ketSys%paw(iT)%F1)
      deallocate(ketSys%paw(iT)%F2)
      deallocate(ketSys%paw(iT)%bes_J_qr)
    enddo

    deallocate(ketSys%pot)
    deallocate(ketSys%paw)

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
