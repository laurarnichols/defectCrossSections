! NOTE: This file has different programming styles because
! I haven't been able to fully clean this file up yet.
! Please excuse the mess.
module TMEmod
  use constants, only: dp, pi, eVToHartree, ii
  use miscUtilities, only: int2str, int2strLeadZero
  use energyTabulatorMod, only: energyTableDir, readEnergyTable
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nKPoints, nSpins, order, ispSelect, loopSpins

  use errorsAndMPI
  use mpi
  
  implicit none


  real(kind=dp) :: dq_j
    !! \(\delta q_j) for displaced wave functions
    !! (only order = 1)
  real(kind=dp), allocatable :: gCart(:,:)
    !! G-vectors in Cartesian coordinates
  real(kind=dp) :: realLattVec(3,3)
    !! Real space lattice vectors
  real(kind=dp) :: recipLattVec(3,3)
    !! Reciprocal lattice vectors
  real(kind = dp) :: t1, t2
    !! For timing different processes


  integer :: maxAngMom
    !! Maximum angular momentum of the projectors
  integer :: maxGIndexGlobal
    !! Maximum G-vector index among all \(G+k\)
    !! and processors for PC and SD
  integer, allocatable :: mill_local(:,:)
    !! Local Miller indices
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors
  integer :: nGkVecsLocalPC, nGkVecsLocalSD
    !! Local number of G+k vectors on this processor
  integer :: nSpinsPC, nSpinsSD
    !! Number of spins for PC/SD system
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

  character(len = 300) :: input, inputPC, textDum
  character(len = 320) :: mkdir
  
  integer :: ik, ig, ibi, ibf
  integer :: iTypes, iPn
  integer :: nIonsSD, nIonsPC, nProjsPC, numOfTypesPC
  integer :: numOfTypes, nProjsSD
  integer :: numOfUsedGvecsPP, npwNi, npwNf, npwMi, npwMf
  integer :: np, nI, nF, nPP, ind2
  integer :: i, j, n1, n2, n3, n4, n, id, npw
  
  real(kind = dp) t0, tf
  
  real(kind = dp) :: omega
  
  real(kind = dp), allocatable :: gvecs(:,:), posIonSD(:,:), posIonPC(:,:)
  real(kind = dp), allocatable :: wk(:), xk(:,:)
  real(kind = dp), allocatable :: DE(:,:,:), absVfi2(:,:,:)
  
  complex(kind = dp), allocatable :: paw_SDKKPC(:,:), paw_id(:,:)
  complex(kind = dp), allocatable :: pawKPC(:,:,:), pawSDK(:,:,:), pawPsiPC(:,:), pawSDPhi(:,:)
  complex(kind = dp), allocatable :: cProjPC(:,:), cProjSD(:,:)
  complex(kind = dp), allocatable :: paw_PsiPC(:,:), paw_SDPhi(:,:)
  complex(kind = dp), allocatable :: cProjBetaPCPsiSD(:,:)
  complex(kind = dp), allocatable :: cProjBetaSDPhiPC(:,:)
  
  integer, allocatable :: TYPNISD(:), TYPNIPC(:), igvs(:,:,:), pwGvecs(:,:), iqs(:)
  integer, allocatable :: npwsSD(:), pwGs(:,:), nIs(:,:), nFs(:,:), ngs(:,:)
  integer, allocatable :: npwsPC(:)
  real(kind = dp), allocatable :: wkPC(:), xkPC(:,:)
  
  type :: atom
    integer :: numOfAtoms, lMax, lmMax, nMax, iRc
    integer, allocatable :: lps(:)
    real(kind = dp), allocatable :: r(:), rab(:), wae(:,:), wps(:,:), F(:,:), F1(:,:,:), F2(:,:,:), bes_J_qr(:,:,:)
  end type atom
  
  TYPE(atom), allocatable :: atoms(:), atomsPC(:)
  
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
  subroutine readInput(ispSelect, maxAngMom, maxGIndexGlobal, nKPoints, nGVecsGlobal, realLattVec, recipLattVec, baselineDir, &
        loopSpins, subtractBaseline)

    use miscUtilities, only: getFirstLineWithKeyword, ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: maxAngMom
      !! Maximum angular momentum of the projectors
    integer, intent(out) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC and SD
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors
    integer, intent(out) :: nKPoints
      !! Total number of k-points

    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors

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
    integer :: maxGIndexGlobalPC
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC
    integer :: maxGIndexGlobalSD
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for SD
    

    if(ionode) then
    
      call initialize(ispSelect, baselineDir, subtractBaseline)
    
      read(5, TME_Input, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputParams', 'reading TME_Input namelist', abs(ierr))
    
      call checkInitialization(ispSelect, baselineDir, loopSpins, subtractBaseline)

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


    maxAngMom = 0

    call readInputPC(maxAngMom, nKPoints, maxGIndexGlobalPC)
    call readInputSD(maxAngMom, nKPoints, maxGIndexGlobalSD, nGVecsGlobal, realLattVec, recipLattVec)

    call MPI_BCAST(maxAngMom, 1, MPI_INTEGER, root, worldComm, ierr)
    
    maxAngMom = 2*maxAngMom + 1

    
    if(ionode) then
      
      if(order == 1) then

        open(30,file=trim(dqFName))
        call ignoreNextNLinesFromFile(30, 1+phononModeJ-1)
          ! Ignore header and all modes before phononModeJ
        read(30,*) iDum, dq_j
        close(30)

      endif

      maxGIndexGlobal = max(maxGIndexGlobalPC, maxGIndexGlobalSD)

      if(maxGIndexGlobal > nGVecsGlobal) call exitError('readInput', &
          'Trying to reference G vecs outside of max grid size. Try switching which grid is read.', 1)

    endif

    if(order == 1) then
      call MPI_BCAST(dq_j, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    endif

    call MPI_BCAST(maxGIndexGlobal, 1, MPI_INTEGER, root, worldComm, ierr)

    nSpins = max(nSpinsPC,nSpinsSD)
    
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
  subroutine checkInitialization(ispSelect, baselineDir, loopSpins, subtractBaseline)
    
    implicit none

    ! Output variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user

    character(len=300), intent(in) :: baselineDir
      !! File name for baseline overlap to optionally
      !! be subtracted for the first-order term

    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel
    logical, intent(in) :: subtractBaseline
      !! If baseline should be subtracted from first-order
      !! overlap for increased numerical accuracy in the 
      !! derivative
    
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

      if(subtractBaseline) abortExecution = checkDirInitialization('baselineDir', baselineDir, 'allElecOverlap.1.1') .or. abortExecution
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
  subroutine readInputPC(maxAngMom, nKPoints, maxGIndexGlobalPC)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
    integer, intent(inout) :: maxAngMom
      !! Maximum angular momentum of the projectors
    integer, intent(out) :: nKPoints
      !! Total number of k-points

    ! Output variables:
    integer, intent(out) :: maxGIndexGlobalPC
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC
    
    ! Local variables:
    integer :: nBands
      !! Number of bands
    integer :: i, j, l, ind, ik, iDum, iType, ni, irc
    
    real(kind = dp) :: t1, t2 
    
    character(len = 300) :: textDum
    
    
    if(ionode) then
      call cpu_time(t1)
    
      write(*,*)
      write(*,'(" Reading perfect crystal inputs.")')
      write(*,*)
    
      inputPC = trim(trim(exportDirPC)//'/input')
    
      open(50, file=trim(inputPC), status = 'old')
    
      read(50, '(a)') textDum
      read(50, * ) 
    
      read(50, '(a)') textDum
      read(50, '(i10)') nSpinsPC
    
      read(50, '(a)') textDum
      read(50, '(i10)') nKPoints

    endif

    call MPI_BCAST(nSpinsPC, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(npwsPC(nKPoints), wkPC(nKPoints), xkPC(3,nKPoints))
    
    if(ionode) then

      read(50, '(a)') textDum
    
      do ik = 1, nKPoints
      
        read(50, '(2i10,4ES24.15E3)') iDum, npwsPC(ik), wkPC(ik), xkPC(1:3,ik)
      
      enddo

    endif
    
    call MPI_BCAST(npwsPC, nKPoints, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(wkPC, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(xkPC, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if(ionode) then

      read(50, '(a)') textDum
      read(50, * ) ! nGVecsGlobal
    
      read(50, '(a)') textDum
      read(50, '(i10)') maxGIndexGlobalPC

      call ignoreNextNLinesFromFile(50, 10)
    
      read(50, '(a)') textDum
      read(50, '(i10)') nIonsPC
    
      read(50, '(a)') textDum
      read(50, '(i10)') numOfTypesPC
    
    endif

    call MPI_BCAST(nIonsPC, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(numOfTypesPC, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(posIonPC(3,nIonsPC), TYPNIPC(nIonsPC))


    if(ionode) then
    
      read(50, '(a)') textDum
      do ni = 1, nIonsPC
        read(50,'(i10, 3ES24.15E3)') TYPNIPC(ni), (posIonPC(j,ni) , j = 1,3)
      enddo
    
      read(50, '(a)') textDum
      read(50,'(i10)') nBands

      if(iBandIfinal > nBands .or. iBandFfinal > nBands) &
        call exitError('readInputPC', 'band limits outside the number of bands in the system '//trim(int2str(nBands)), 1)
        ! Only need to test these bands because we tested in
        ! the `checkInitialization` subroutine to make sure
        ! that the `initial` bands are lower than the `final`
        ! bands

    endif

    call MPI_BCAST(TYPNIPC,  size(TYPNIPC),  MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(posIonPC, size(posIonPC), MPI_DOUBLE_PRECISION,root,worldComm,ierr)
    
    allocate(atomsPC(numOfTypesPC))

    nProjsPC = 0
    do iType = 1, numOfTypesPC
      
      if(ionode) then

        read(50, '(a)') textDum
        read(50, *) 
      
        read(50, '(a)') textDum
        read(50, '(i10)') atomsPC(iType)%numOfAtoms

        read(50, '(a)') textDum
        read(50, '(i10)') atomsPC(iType)%lMax              ! number of projectors

      endif

      call MPI_BCAST(atomsPC(iType)%numOfAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%lMax, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(atomsPC(iType)%lps(atomsPC(iType)%lMax))
      
      if(ionode) then

        read(50, '(a)') textDum

        do i = 1, atomsPC(iType)%lMax 

          read(50, '(2i10)') l, ind
          atomsPC(iType)%lps(ind) = l

        enddo

      endif

      call MPI_BCAST(atomsPC(iType)%lps, atomsPC(iType)%lMax, MPI_INTEGER, root, worldComm, ierr)

      if(ionode) then
      
        read(50, '(a)') textDum
        read(50, '(i10)') atomsPC(iType)%lmMax
      
        read(50, '(a)') textDum
        read(50, '(2i10)') atomsPC(iType)%nMax, atomsPC(iType)%iRc

      endif
    
      call MPI_BCAST(atomsPC(iType)%lmMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%nMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%iRc, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(atomsPC(iType)%r(atomsPC(iType)%nMax))
      
      if(ionode) then
      
        allocate(atomsPC(iType)%rab(atomsPC(iType)%nMax))

        read(50, '(a)') textDum

        do i = 1, atomsPC(iType)%nMax
          read(50, '(2ES24.15E3)') atomsPC(iType)%r(i), atomsPC(iType)%rab(i)
        enddo
       
        allocate(atomsPC(iType)%wae(atomsPC(iType)%nMax, atomsPC(iType)%lMax))
        allocate(atomsPC(iType)%wps(atomsPC(iType)%nMax, atomsPC(iType)%lMax))
      
        read(50, '(a)') textDum
        do j = 1, atomsPC(iType)%lMax
          do i = 1, atomsPC(iType)%nMax
            read(50, '(2ES24.15E3)') atomsPC(iType)%wae(i, j), atomsPC(iType)%wps(i, j) 
          enddo
        enddo
        
      endif

      allocate(atomsPC(iType)%F(atomsPC(iType)%iRc, atomsPC(iType)%lMax))
      allocate(atomsPC(iType)%F1(atomsPC(iType)%iRc, atomsPC(iType)%lMax, atomsPC(iType)%lMax))
      allocate(atomsPC(iType)%F2(atomsPC(iType)%iRc, atomsPC(iType)%lMax, atomsPC(iType)%lMax))
      
      if(ionode) then

        atomsPC(iType)%F = 0.0_dp
        atomsPC(iType)%F1 = 0.0_dp
        atomsPC(iType)%F2 = 0.0_dp
      
        do j = 1, atomsPC(iType)%lMax
        
          irc = atomsPC(iType)%iRc
          atomsPC(iType)%F(1:irc,j)=(atomsPC(iType)%wae(1:irc,j)-atomsPC(iType)%wps(1:irc,j))* &
                atomsPC(iType)%r(1:irc)*atomsPC(iType)%rab(1:irc)
        
          do i = 1, atomsPC(iType)%lMax
            atomsPC(iType)%F1(1:irc,i,j) = ( atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wae(1:irc,j) - &
                                            atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wps(1:irc,j))*atomsPC(iType)%rab(1:irc)
          
            atomsPC(iType)%F2(1:irc,i,j) = ( atomsPC(iType)%wae(1:irc,i)*atomsPC(iType)%wae(1:irc,j) - &
                                             atomsPC(iType)%wae(1:irc,i)*atomsPC(iType)%wps(1:irc,j) - &
                                             atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wae(1:irc,j) + &
                                             atomsPC(iType)%wps(1:irc,i)*atomsPC(iType)%wps(1:irc,j))*atomsPC(iType)%rab(1:irc)

          enddo
        enddo

        deallocate(atomsPC(iType)%wae)
        deallocate(atomsPC(iType)%wps)
        deallocate(atomsPC(iType)%rab)
      
        nProjsPC = nProjsPC + atomsPC(iType)%numOfAtoms*atomsPC(iType)%lmMax

      endif

      call MPI_BCAST(atomsPC(iType)%r, size(atomsPC(iType)%r), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%F, size(atomsPC(iType)%F), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%F1, size(atomsPC(iType)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atomsPC(iType)%F2, size(atomsPC(iType)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      
    enddo
  
    call MPI_BCAST(nProjsPC, 1, MPI_INTEGER, root, worldComm, ierr)

    if(ionode) then
    
      close(50)

      do iType = 1, numOfTypes
        do i = 1, atoms(iType)%lMax
          if(atoms(iType)%lps(i) > maxAngMom) maxAngMom = atoms(iType)%lps(i)
        enddo
      enddo
    
      call cpu_time(t2)
      write(*,'(" Reading input files done in:                ", f10.2, " secs.")') t2-t1
      write(*,*)

    endif
    
    return
    
  end subroutine readInputPC
  
!----------------------------------------------------------------------------
  subroutine readInputSD(maxAngMom, nKPoints, maxGIndexGlobalSD, nGVecsGlobal, realLattVec, recipLattVec)
    !
    implicit none
    
    ! Input variables:
    integer, intent(inout) :: maxAngMom
      !! Maximum angular momentum of the projectors
    integer, intent(in) :: nKPoints
      !! Total number of k-points

    ! Output variables:
    integer, intent(out) :: maxGIndexGlobalSD
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for SD
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors

    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors
    
    ! Local variables:
    integer :: nBands
      !! Number of bands
    integer :: nKpts
      !! Number of k-points read from SD file

    integer :: i, j, l, ind, ik, iDum, iType, ni, irc
    
    real(kind = dp) :: t1, t2
    
    character(len = 300) :: textDum
    

    if(ionode) then
      call cpu_time(t1)
    
      write(*,*)
      write(*,'(" Reading solid defect inputs.")')
      write(*,*)
    
      input = trim(trim(exportDirSD)//'/input')
    
      open(50, file=trim(input), status = 'old')
    
      read(50, '(a)') textDum
      read(50, '(ES24.15E3)' ) omega

    endif

    call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    
    if(ionode) then

      read(50, '(a)') textDum
      read(50, '(i10)') nSpinsSD

      read(50, '(a)') textDum
      read(50, '(i10)') nKpts

      if(nKpts /= nKPoints) call exitError('readInputsSD', 'Number of k-points in systems must match', 1)
    
      read(50, '(a)') textDum

    endif

    call MPI_BCAST(nSpinsSD, 1, MPI_INTEGER, root, worldComm, ierr)

    
    allocate(npwsSD(nKPoints), wk(nKPoints), xk(3,nKPoints))

    if(ionode) then
    
      do ik = 1, nKPoints
      
        read(50, '(2i10,4ES24.15E3)') iDum, npwsSD(ik), wk(ik), xk(1:3,ik)
      
      enddo

    endif
    
    call MPI_BCAST(npwsSD, nKPoints, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(wk, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(xk, nKPoints, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if(ionode) then
    
      read(50, '(a)') textDum
      read(50, '(i10)') nGVecsGlobal
    
      read(50, '(a)') textDum
      read(50, '(i10)') maxGIndexGlobalSD
    
      read(50, '(a)') textDum     
      read(50,*) 
      !read(50, '(6i10)') fftxMin, fftxMax, fftyMin, fftyMax, fftzMin, fftzMax
    
      read(50, '(a)') textDum
      read(50, '(a5, 3ES24.15E3)') textDum, realLattVec(1:3,1)
      read(50, '(a5, 3ES24.15E3)') textDum, realLattVec(1:3,2)
      read(50, '(a5, 3ES24.15E3)') textDum, realLattVec(1:3,3)
    
      read(50, '(a)') textDum
      read(50, '(a5, 3ES24.15E3)') textDum, recipLattVec(1:3,1)
      read(50, '(a5, 3ES24.15E3)') textDum, recipLattVec(1:3,2)
      read(50, '(a5, 3ES24.15E3)') textDum, recipLattVec(1:3,3)
    
      read(50, '(a)') textDum
      read(50, '(i10)') nIonsSD
    
      read(50, '(a)') textDum
      read(50, '(i10)') numOfTypes

    endif

    call MPI_BCAST(nGVecsGlobal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(recipLattVec, size(recipLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(nIonsSD, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(numOfTypes, 1, MPI_INTEGER, root, worldComm, ierr)
    
    allocate(posIonSD(3,nIonsSD), TYPNISD(nIonsSD))

    if(ionode) then
    
      read(50, '(a)') textDum
      do ni = 1, nIonsSD
        read(50,'(i10, 3ES24.15E3)') TYPNISD(ni), (posIonSD(j,ni), j = 1,3)
      enddo
    
      read(50, '(a)') textDum
      read(50,'(i10)') nBands

      if(iBandIfinal > nBands .or. iBandFfinal > nBands) &
        call exitError('readInputSD', 'band limits outside the number of bands in the system '//trim(int2str(nBands)), 1)
        ! Only need to test these bands because we tested in
        ! the `checkInitialization` subroutine to make sure
        ! that the `initial` bands are lower than the `final`
        ! bands

    endif
  
    call MPI_BCAST(TYPNISD, size(TYPNISD), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(posIonSD, size(posIonSD), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


    allocate(atoms(numOfTypes))
    
    nProjsSD = 0
    do iType = 1, numOfTypes
      
      if(ionode) then

        read(50, '(a)') textDum
        read(50, *) 
      
        read(50, '(a)') textDum
        read(50, '(i10)') atoms(iType)%numOfAtoms
      
        read(50, '(a)') textDum
        read(50, '(i10)') atoms(iType)%lMax              ! number of projectors

      endif

      call MPI_BCAST(atoms(iType)%numOfAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%lMax, 1, MPI_INTEGER, root, worldComm, ierr)
      
      allocate(atoms(iType)%lps(atoms(iType)%lMax))
      
      if(ionode) then

        read(50, '(a)') textDum
        do i = 1, atoms(iType)%lMax 
          read(50, '(2i10)') l, ind
          atoms(iType)%lps(ind) = l
        enddo

      endif

      call MPI_BCAST(atoms(iType)%lps, size(atoms(iType)%lps), MPI_INTEGER, root, worldComm, ierr)

      if(ionode) then
      
        read(50, '(a)') textDum
        read(50, '(i10)') atoms(iType)%lmMax
      
        read(50, '(a)') textDum
        read(50, '(2i10)') atoms(iType)%nMax, atoms(iType)%iRc

      endif
    
      call MPI_BCAST(atoms(iType)%lmMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%nMax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%iRc, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(atoms(iType)%r(atoms(iType)%nMax))

      if(ionode) then
      
        allocate(atoms(iType)%rab(atoms(iType)%nMax))
      
        read(50, '(a)') textDum
        do i = 1, atoms(iType)%nMax
          read(50, '(2ES24.15E3)') atoms(iType)%r(i), atoms(iType)%rab(i)
        enddo
       
        allocate(atoms(iType)%wae(atoms(iType)%nMax, atoms(iType)%lMax))
        allocate(atoms(iType)%wps(atoms(iType)%nMax, atoms(iType)%lMax))
      
        read(50, '(a)') textDum
        do j = 1, atoms(iType)%lMax
          do i = 1, atoms(iType)%nMax
            read(50, '(2ES24.15E3)') atoms(iType)%wae(i, j), atoms(iType)%wps(i, j) 
          enddo
        enddo

      endif

        
      allocate(atoms(iType)%F(atoms(iType)%iRc, atoms(iType)%lMax))
      allocate(atoms(iType)%F1(atoms(iType)%iRc, atoms(iType)%lMax, atoms(iType)%lMax))
      allocate(atoms(iType)%F2(atoms(iType)%iRc, atoms(iType)%lMax, atoms(iType)%lMax))
      
      if(ionode) then

        atoms(iType)%F = 0.0_dp
        atoms(iType)%F1 = 0.0_dp
        atoms(iType)%F2 = 0.0_dp
      
        do j = 1, atoms(iType)%lMax
          
          irc = atoms(iType)%iRc
          atoms(iType)%F(1:irc,j)=(atoms(iType)%wae(1:irc,j)-atoms(iType)%wps(1:irc,j))*atoms(iType)%r(1:irc) * &
              atoms(iType)%rab(1:irc)
          
          do i = 1, atoms(iType)%lMax
                  
            atoms(iType)%F1(1:irc,i,j) = ( atoms(iType)%wae(1:irc,i)*atoms(iType)%wps(1:irc,j) - &
                                           atoms(iType)%wps(1:irc,i)*atoms(iType)%wps(1:irc,j))*atoms(iType)%rab(1:irc)
          
            atoms(iType)%F2(1:irc,i,j) = ( atoms(iType)%wae(1:irc,i)*atoms(iType)%wae(1:irc,j) - &
                                           atoms(iType)%wae(1:irc,i)*atoms(iType)%wps(1:irc,j) - &
                                           atoms(iType)%wps(1:irc,i)*atoms(iType)%wae(1:irc,j) + &
                                           atoms(iType)%wps(1:irc,i)*atoms(iType)%wps(1:irc,j))*atoms(iType)%rab(1:irc)
          enddo
        enddo
      
        nProjsSD = nProjsSD + atoms(iType)%numOfAtoms*atoms(iType)%lmMax
      
        deallocate(atoms(iType)%wae)
        deallocate(atoms(iType)%wps)
        deallocate(atoms(iType)%rab)

      endif

      call MPI_BCAST(atoms(iType)%r, size(atoms(iType)%r), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%F, size(atoms(iType)%F), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%F1, size(atoms(iType)%F1), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(atoms(iType)%F2, size(atoms(iType)%F2), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      
    enddo

    call MPI_BCAST(nProjsSD, 1, MPI_INTEGER, root, worldComm, ierr)
    
    if(ionode) then
    
      close(50)

      do iType = 1, numOfTypes
        do i = 1, atoms(iType)%lMax
          if(atoms(iType)%lps(i) > maxAngMom) maxAngMom = atoms(iType)%lps(i)
        enddo
      enddo

      call cpu_time(t2)
      write(*,'(" Reading solid defect inputs done in:                ", f10.2, " secs.")') t2-t1
      write(*,*)

    endif
    
    return
    
  end subroutine readInputSD
  
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
  subroutine setUpTables(maxAngMom, nGVecsLocal, mill_local, numOfTypes, recipLattVec, atomsSD, gCart, Ylm)

    implicit none

    ! Input variables:
    integer, intent(in) :: maxAngMom
      !! Max J index from L and M
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors
    integer, intent(in) :: mill_local(3,nGVecsLocal)
      !! Miller indices for local G-vectors
    integer, intent(in) :: numOfTypes
      !! Number of different atom types

    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal-space lattice vectors

    type(atom), intent(inout) :: atomsSD(numOfTypes)
      !! Structure to hold details for each atom for SD

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
    
    do iT = 1, numOfTypes

      allocate(atomsSD(iT)%bes_J_qr( 0:maxAngMom, atomsSD(iT)%iRc, nGVecsLocal))
      atomsSD(iT)%bes_J_qr(:,:,:) = 0.0_dp
      
    enddo

    do ig = 1, nGVecsLocal

      gCart(:,ig) = matmul(recipLattVec, mill_local(:,ig))

      q = sqrt(dot_product(gCart(:,ig),gCart(:,ig)))

      gUnit(:) = gCart(:,ig)
      if(abs(q) > 1.0e-6_dp) gUnit = gUnit/q
        !! Get unit vector for Ylm calculation

      call getYlm(gUnit, maxAngMom, Ylm(:,ig))
        !! Calculate all the needed spherical harmonics

      do iT = 1, numOfTypes

        do iR = 1, atomsSD(iT)%iRc

          JL = 0.0_dp
          call bessel_j(q*atomsSD(iT)%r(iR), maxAngMom, JL) ! returns the spherical bessel at qr point
            ! Previously used SD atoms structure here for both PC and SD
          atomsSD(iT)%bes_J_qr(:,iR,ig) = JL(:)

        enddo

      enddo

    enddo


    return

  end subroutine setUpTables

!----------------------------------------------------------------------------
  subroutine readProjectors(crystalType, iGkStart_pool, ikGlobal, nGkVecsLocal, nProjs, beta)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iGkStart_pool
      !! Start G+k index for this process in pool
    integer, intent(in) :: ikGlobal
      !! Current k point
    integer, intent(in) :: nGkVecsLocal
      !! Local number of G+k vectors on this processor
    integer, intent(in) :: nProjs
      !! Number of projectors

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD) for projectors

    ! Output variables:
    complex(kind=dp), intent(out) :: beta(nGkVecsLocal,nProjs)
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
      !! array of length `nProjs`

    open(unit=72, file=trim(fNameExport), access='direct', recl=reclen, iostat=ierr, status='old', SHARED)


    do igkLocal = 1, nGkVecsLocal

      igkGlobal = igkLocal+iGkStart_pool-1

      read(72,rec=igkGlobal+1) (beta(igkLocal,ipr), ipr=1,nProjs)

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
  subroutine pawCorrectionWfc(nAtoms, iType, nProjs, numOfTypes, cProjI, cProjF, atoms, pawWfc)
    ! calculates the augmentation part of the transition matrix element
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: nProjs
      !! First index for `cProjI` and `cProjF`
    integer, intent(in) :: numOfTypes
      !! Number of different atom types

    complex(kind = dp) :: cProjI(nProjs,iBandIinit:iBandIfinal)
      !! Initial-system (PC) projection
    complex(kind = dp) :: cProjF(nProjs,iBandFinit:iBandFfinal)
      !! Final-system (SD) projection

    type(atom), intent(in) :: atoms(numOfTypes)
      !! Structure to hold details for each atom

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
    
    do ia = 1, nAtoms
      
      iT = iType(ia)

      LM = 0
      DO LL = 1, atoms(iT)%lMax
        L = atoms(iT)%LPS(LL)
        DO M = -L, L
          LM = LM + 1 !1st index for CPROJ
          
          LMP = 0
          DO LLP = 1, atoms(iT)%lMax
            LP = atoms(iT)%LPS(LLP)
            DO MP = -LP, LP
              LMP = LMP + 1 ! 2nd index for CPROJ
              
              atomicOverlap = 0.0_dp
              if ( (L == LP).and.(M == MP) ) then 
                atomicOverlap = sum(atoms(iT)%F1(:,LL, LLP))
                
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
      LMBASE = LMBASE + atoms(iT)%lmMax
    ENDDO
    
    return
    
  end subroutine pawCorrectionWfc

!----------------------------------------------------------------------------
  subroutine pawCorrectionK(crystalType, nAtoms, iType, maxAngMom, nGVecsLocal, numOfTypes, numOfTypesSD, atomPositions, gCart, Ylm, atoms, atomsSD, pawK)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: maxAngMom
      !! Max J index from L and M
    integer, intent(in) :: nGVecsLocal
      !! Number of local G-vectors
    integer, intent(in) :: numOfTypes
      !! Number of different atom types
    integer, intent(in) :: numOfTypesSD
      !! Number of different atom types for SD

    real(kind=dp), intent(in) :: atomPositions(3,nAtoms)
      !! Atom positions in Cartesian coordinates
    real(kind=dp), intent(in) :: gCart(3,nGVecsLocal)
      !! G-vectors in Cartesian coordinates

    complex(kind=dp), intent(in) :: Ylm((maxAngMom+1)**2,nGVecsLocal)
      !! Spherical harmonics

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD)

    type(atom), intent(inout) :: atoms(numOfTypes)
      !! Structure to hold details for each atom
    type(atom), intent(in) :: atomsSD(numOfTypesSD)
      !! Structure to hold details for each atom
      !! for SD (needed for grid for Bessel functions)


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
      
      do ni = 1, nAtoms ! LOOP OVER THE IONS
        
        qDotR = dot_product(gCart(:,ig), atomPositions(:,ni))
        
        if(crystalType == 'PC') then
          ATOMIC_CENTER = exp( -ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        else
          ATOMIC_CENTER = exp( ii*cmplx(qDotR, 0.0_dp, kind = dp) )
        endif
        
        iT = iType(ni)
        LM = 0
        DO LL = 1, atoms(iT)%lMax
          L = atoms(iT)%LPS(LL)
          DO M = -L, L
            LM = LM + 1 !1st index for CPROJ
            
            FI = 0.0_dp
            
            FI = dot_product(atomsSD(iT)%bes_J_qr(L,:,ig),atoms(iT)%F(:,LL)) ! radial part integration F contains rab
            
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

        LMBASE = LMBASE + atoms(iT)%lmMax
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
    integer :: ibi, ibf
      !! Loop indices
    integer :: ikGlobal
      !! Current global k-point
    integer :: totalNumberOfElements
      !! Total number of overlaps to output

    real(kind=dp) :: dE(iBandFinit:iBandFfinal,iBandIinit:iBandIFinal,4)
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


    call readEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, energyTableDir, dE)

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

    integer :: iType
      !! Loop index


    deallocate(npwsPC)
    deallocate(wkPC)
    deallocate(xkPC)
    deallocate(posIonPC)
    deallocate(TYPNIPC)
    deallocate(gCart)
    deallocate(Ylm)

    do iType = 1, numOfTypesPC
      deallocate(atomsPC(iType)%r)
      deallocate(atomsPC(iType)%lps)
      deallocate(atomsPC(iType)%F)
      deallocate(atomsPC(iType)%F1)
      deallocate(atomsPC(iType)%F2)
    enddo

    deallocate(atomsPC)

    deallocate(npwsSD)
    deallocate(wk)
    deallocate(xk)
    deallocate(posIonSD)
    deallocate(TYPNISD)

    do iType = 1, numOfTypes
      deallocate(atoms(iType)%lps)
      deallocate(atoms(iType)%r)
      deallocate(atoms(iType)%F)
      deallocate(atoms(iType)%F1)
      deallocate(atoms(iType)%F2)
      deallocate(atoms(iType)%bes_J_qr)
    enddo

    deallocate(atoms)

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
