module TMEmod
  use constants, only: dp, pi, sq4pi, eVToHartree, ii
  use miscUtilities, only: int2str, int2strLeadZero
  use energyTabulatorMod, only: energyTableDir, readEnergyTable
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, nKPoints, nSpins, order

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
  
  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: exportDirSD, exportDirPC
    !! Paths to exports for SD (left) and PC (right) 
    !! systems to use for wave function overlaps <SD|PC>


  character(len = 200) :: VfisOutput
  character(len = 300) :: input, inputPC, textDum, elementsPath
  character(len = 320) :: mkdir
  
  integer :: ik, ki, kf, ig, ibi, ibf
  integer :: JMAX, iTypes, iPn
  integer :: nIonsSD, nIonsPC, nProjsPC, numOfTypesPC
  integer :: numOfTypes, nBands, nProjsSD
  integer :: numOfUsedGvecsPP, ios, npwNi, npwNf, npwMi, npwMf
  integer :: gx, gy, gz, nGvsI, nGvsF, nGi, nGf
  integer :: np, nI, nF, nPP, ind2
  integer :: i, j, n1, n2, n3, n4, n, id, npw
  
  integer, allocatable :: counts(:), displmnt(:), nPWsI(:), nPWsF(:)
  
  real(kind = dp) t0, tf
  
  real(kind = dp) :: omega, threej
  
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
  
  real(kind = dp) :: eBin
  complex(kind = dp) :: paw, pseudo1, pseudo2, paw2
  
  logical :: gamma_only, master, calculateVfis, coulomb, tmes_file_exists
  
  namelist /TME_Input/ exportDirSD, exportDirPC, elementsPath, energyTableDir, iBandIinit, &
                       iBandIfinal, iBandFinit, iBandFfinal, ki, kf, calculateVfis, VfisOutput, &
                       eBin, order, dqFName, phononModeJ
  
  
contains

!----------------------------------------------------------------------------
  subroutine readInput(maxGIndexGlobal, nKPoints, nGVecsGlobal, realLattVec, recipLattVec)

    use miscUtilities, only: getFirstLineWithKeyword, ignoreNextNLinesFromFile
    
    implicit none

    ! Output variables:
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

    ! Local variables:    
    integer :: iDum
      !! Dummy integer to ignore file input
    integer :: maxGIndexGlobalPC
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC
    integer :: maxGIndexGlobalSD
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for SD

    character(len=300) :: line
      !! Line from file
    

    if(ionode) then
    
      call initialize()
    
      read(5, TME_Input, iostat=ios)
    
      call checkInitialization()

    endif

    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(order, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(phononModeJ, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(calculateVfis, 1, MPI_LOGICAL, root, worldComm, ierr)

    call MPI_BCAST(eBin, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(exportDirSD, len(exportDirSD), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(exportDirPC, len(exportDirPC), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(elementsPath, len(elementsPath), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(VfisOutput, len(exportDirSD), MPI_CHARACTER, root, worldComm, ierr)

    call readInputPC(nKPoints, maxGIndexGlobalPC)
    call readInputSD(nKPoints, maxGIndexGlobalSD, nGVecsGlobal, realLattVec, recipLattVec)
    
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
  subroutine initialize()
    
    implicit none

    ! Output variables:
    

    exportDirSD = ''
    exportDirPC = ''
    energyTableDir = ''
    dqFName = ''
    elementsPath = './TMEs'
    VfisOutput = ''
    
    ki = -1
    kf = -1

    order = -1

    phononModeJ = -1
    
    eBin = 0.01_dp
    
    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1
    
    calculateVfis = .false.
    
    return
    
  end subroutine initialize
  
!----------------------------------------------------------------------------
  subroutine checkInitialization()
    
    implicit none
    
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
    endif

    abortExecution = checkStringInitialization('elementsPath', elementsPath) .or. abortExecution
    abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
    abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 


    call system('mkdir -p '//trim(elementsPath))

    
    write(*,'("calculateVfis = ", l )') calculateVfis

    if(calculateVfis) then
      if(iBandFinit /= iBandFfinal) then
        write(*,*)
        write(*,'(" Vfis can be calculated only if there is a single final-state band!")')
        write(*,'(" This variable is mandatory and thus the program will not be executed!")')
        abortExecution = .true.
      endif

      abortExecution = checkStringInitialization('VfisOutput', VfisOutput) .or. abortExecution
      abortExecution = checkDoubleInitialization('eBin', eBin, 0.0_dp, 2.0_dp) .or. abortExecution

    endif

    eBin = eBin*eVToHartree
    

    if(abortExecution) then
      write(*,'(" Program stops!")')
      stop
    endif
    
    return
    
  end subroutine checkInitialization
  
!----------------------------------------------------------------------------
  subroutine readInputPC(nKPoints, maxGIndexGlobalPC)

    use miscUtilities, only: ignoreNextNLinesFromFile
    
    implicit none

    ! Outpu variables:
    integer, intent(out) :: nKPoints
      !! Total number of k-points

    ! Output variables:
    integer, intent(out) :: maxGIndexGlobalPC
      !! Maximum G-vector index among all \(G+k\)
      !! and processors for PC
    
    ! Local variables:
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
      read(50, * )

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
    
      call cpu_time(t2)
      write(*,'(" Reading input files done in:                ", f10.2, " secs.")') t2-t1
      write(*,*)

    endif
    
    return
    
  end subroutine readInputPC
  
!----------------------------------------------------------------------------
  subroutine readInputSD(nKPoints, maxGIndexGlobalSD, nGVecsGlobal, realLattVec, recipLattVec)
    !
    implicit none
    
    ! Input variables:
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

    endif
  
    call MPI_BCAST(TYPNISD, size(TYPNISD), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(posIonSD, size(posIonSD), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    if(ionode) then
    
      read(50, '(a)') textDum
      read(50, '(i10)') nBands

    endif

    call MPI_BCAST(nBands, 1, MPI_INTEGER, root, worldComm, ierr)
    
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

      JMAX = 0
      do iType = 1, numOfTypes
        do i = 1, atoms(iType)%lMax
          if ( atoms(iType)%lps(i) > JMAX ) JMAX = atoms(iType)%lps(i)
        enddo
      enddo
    
      JMAX = 2*JMAX + 1

    endif

    call MPI_BCAST(JMAX, 1, MPI_INTEGER, root, worldComm, ierr)
    
    if(ionode) then

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
    logical :: fileExists
      !! If the overlap file exists for the given 
      !! k-point and spin channel


    inquire(file=trim(getMatrixElementFName(ikGlobal, isp, elementsPath)), exist=fileExists)
    
  end function overlapFileExists

!----------------------------------------------------------------------------
  subroutine getFullPWGrid(iGStart_pool, iGEnd_pool, nGVecsLocal, nGVecsGlobal, mill_local)
    !! Read full PW grid from mgrid file
    
    implicit none

    ! Input variables:
    integer, intent(in) :: iGStart_pool, iGEnd_pool
      !! Start and end G-vectors for each process in pool
    integer, intent(in) :: nGVecsLocal
      !! Local number of G-vectors
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors

    ! Output variables:
    integer, intent(out) :: mill_local(3,nGVecsLocal)
      !! Integer coefficients for G-vectors
    
    ! Local variables:
    integer :: gVecMillerIndicesGlobal(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer :: iDum
      !! Ignore dummy integer
    integer :: ig
      !! Loop index
    

    if(ionode) then

      open(72, file=trim(exportDirSD)//"/mgrid")
        !! Read full G-vector grid from defect folder.
        !! This assumes that the grids are the same.
    
      read(72, * )
      read(72, * )
    
      do ig = 1, nGVecsGlobal

        read(72, '(4i10)') iDum, gVecMillerIndicesGlobal(1:3,ig)

      enddo
    
      close(72)

    endif

    call MPI_BCAST(gVecMillerIndicesGlobal, size(gVecMillerIndicesGlobal), MPI_INTEGER, root, worldComm, ierr)

    mill_local(:,:) = gVecMillerIndicesGlobal(:,iGStart_pool:iGEnd_pool)
    
    return
    
  end subroutine getFullPWGrid

!----------------------------------------------------------------------------
  subroutine setUpTables(JMAX, nGVecsLocal, mill_local, numOfTypes, recipLattVec, atomsSD, gCart, Ylm)

    implicit none

    ! Input variables:
    integer, intent(in) :: JMAX
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

    complex(kind=dp), intent(out) :: Ylm((JMAX+1)**2,nGVecsLocal)
      !! Spherical harmonics

    ! Local variables:
    integer :: ig, iT, iR
      !! Loop indices

    real(kind=dp) :: gUnit(3)
      !! Unit G-vector
    real(kind=dp) :: JL(0:JMAX)
      !! Bessel_j temporary variable
    real(kind=dp) :: q
      !! Magnitude of G-vector


    Ylm = cmplx(0.0_dp, 0.0_dp, kind = dp)
    
    do iT = 1, numOfTypes

      allocate(atomsSD(iT)%bes_J_qr( 0:JMAX, atomsSD(iT)%iRc, nGVecsLocal))
      atomsSD(iT)%bes_J_qr(:,:,:) = 0.0_dp
      
    enddo

    do ig = 1, nGVecsLocal

      gCart(:,ig) = matmul(recipLattVec, mill_local(:,ig))

      q = sqrt(dot_product(gCart(:,ig),gCart(:,ig)))

      gUnit(:) = gCart(:,ig)
      if(abs(q) > 1.0e-6_dp) gUnit = gUnit/q
        !! Get unit vector for Ylm calculation

      call getYlm(gUnit, JMAX, Ylm(:,ig))
        !! Calculate all the needed spherical harmonics

      do iT = 1, numOfTypes

        do iR = 1, atomsSD(iT)%iRc

          JL = 0.0_dp
          call bessel_j(q*atomsSD(iT)%r(iR), JMAX, JL) ! returns the spherical bessel at qr point
            ! Previously used SD atoms structure here for both PC and SD
          atomsSD(iT)%bes_J_qr(:,iR,ig) = JL(:)

        enddo

      enddo

    enddo


    return

  end subroutine setUpTables

!----------------------------------------------------------------------------
  subroutine readProjectors(crystalType, iGkStart_pool, ikGlobal, nGkVecsLocal, nProjs, npws, betaLocalPWs)
    
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
    integer, intent(in) :: npws
      !! Number of G+k vectors less than
      !! the cutoff at each k-point

    character(len=2), intent(in) :: crystalType
      !! Crystal type (PC or SD) for projectors

    ! Output variables:
    complex(kind=dp), intent(out) :: betaLocalPWs(nGkVecsLocal,nProjs)
      !! Projector of `projCrystalType` with all projectors
      !! and only local PWs/G+k vectors
    
    ! Local variables:
    integer :: displacement(nProcPerPool)
      !! Offset from beginning of array for
      !! scattering coefficients to each process
    integer :: endingProj(nProcPerPool)
      !! End projector for each process in pool
    integer :: iprStart_pool, iprEnd_pool
      !! Start and end projector for process in pool
    integer :: nProjsLocal
      !! Number of projectors read by this process
    integer :: reclen
      !! Record length for projectors file
    integer :: sendCount(nProcPerPool)
      !! Number of items to send to each process
      !! in the pool
    integer :: ipr, igk, iproc, irec
      !! Loop indices

    complex(kind=dp), allocatable :: betaLocalProjs(:,:)
      !! Projector of `projCrystalType` with all PWs/G+k 
      !! vectors and only local projectors
    
    character(len=300) :: fNameExport
      !! File names
    

    if(crystalType == 'PC') then
      fNameExport = trim(exportDirPC)//"/projectors."//trim(int2str(ikGlobal)) 
    else
      fNameExport = trim(exportDirSD)//"/projectors."//trim(int2str(ikGlobal))
    endif


    call distributeItemsInSubgroups(indexInPool, nProjs, nProcPerPool, nProcPerPool, nProcPerPool, iprStart_pool, iprEnd_pool, nProjsLocal)

    allocate(betaLocalProjs(npws,iprStart_pool:iprEnd_pool))

    betaLocalProjs(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)
    betaLocalPWs(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp)

    call MPI_Barrier(intraPoolComm, ierr)


    inquire(iolength=reclen) betaLocalProjs(:,iprStart_pool)
      !! Get the record length needed to write a double complex
      !! array of length nPWs1k

    open(unit=72, file=trim(fNameExport), access='direct', recl=reclen, iostat=ierr, status='old', SHARED)

    irec = 1

    do ipr = 1, nProjs
      
      irec = irec + 1

      if(ipr >= iprStart_pool .and. ipr <= iprEnd_pool) read(72,rec=irec) (betaLocalProjs(igk,ipr), igk=1,npws)

    enddo

    close(72)
    
    call MPI_Barrier(intraPoolComm, ierr)

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

    endingProj = 0
    endingProj(indexInPool+1) = iprEnd_pool
    call mpiSumIntV(endingProj, intraPoolComm)


    !> Distribute projectors across processors so that PWs are local
    !> instead of projectors
    iproc = 0
    do ipr = 1, nProjs

      if(ipr == endingProj(iproc+1)+1) iproc = iproc + 1

      call MPI_SCATTERV(betaLocalProjs(:,ipr), sendCount, displacement, MPI_DOUBLE_COMPLEX, betaLocalPWs(1:nGkVecsLocal,ipr), nGkVecsLocal, &
          MPI_DOUBLE_COMPLEX, iproc, intraPoolComm, ierr)

    enddo

    deallocate(betaLocalProjs)
    
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
    integer :: ib, igk, iproc
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
    integer :: ipr, ib
      !! Loop indices
    
    
    cProj(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
    
    if(indexInPool == 0) then

      ! Open the projections file for the given crystal type
      if(crystalType == 'PC') then
        open(72, file=trim(exportDirPC)//"/projections."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)))
      else
        open(72, file=trim(exportDirSD)//"/projections."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)))
      endif
    

      call ignoreNextNLinesFromFile(72, 1 + (iBandinit-1)*nProjs)
        ! Ignore header and all bands before initial band
    
      ! Read the projections
      do ib = iBandinit, iBandfinal
        do ipr = 1, nProjs 

          read(72,'(2ES24.15E3)') cProj(ipr,ib)

        enddo
      enddo
    
      close(72)

    endif

    call MPI_BCAST(cProj, size(cProj), MPI_DOUBLE_COMPLEX, root, intraPoolComm, ierr)
      ! Broadcast entire array to all processes
    
    return
    
  end subroutine readProjections
  
!----------------------------------------------------------------------------
  subroutine calculateCrossProjection(iBandinit, iBandfinal, ikGlobal, nGkVecsLocal1, nGkVecsLocal2, nProjs, beta, wfc, crossProjection)
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
    integer, intent(in) :: ikGlobal
      !! Current k point
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
  subroutine pawCorrectionK(crystalType, nAtoms, iType, JMAX, nGVecsLocal, numOfTypes, numOfTypesSD, atomPositions, gCart, Ylm, atoms, atomsSD, pawK)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: JMAX
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

    complex(kind=dp), intent(in) :: Ylm((JMAX+1)**2,nGVecsLocal)
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
    integer :: ig, iT, iR, ni, ibi, ibf, ind
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
  subroutine writeResults(ikLocal, isp)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ibi, ibf
      !! Loop indices
    integer :: ikGlobal
      !! Current global k-point
    integer :: totalNumberOfElements
      !! Total number of overlaps to output

    real(kind=dp) :: dE(4,iBandFinit:iBandFfinal,iBandIinit:iBandIFinal)
      !! Energy difference to be combined with
      !! overlap for matrix element
    
    character(len = 300) :: text
      !! Text for header


    ikGlobal = ikLocal+ikStart_pool-1
    
    open(17, file=trim(getMatrixElementFName(ikGlobal, isp, elementsPath)), status='unknown')
    
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
    
      text = "# Final Band, Initial Band, Complex <f|i>, |<f|i>|^2, |dE*<f|i>/dq_j|^2 (Hartree^2/(Bohr*sqrt(elec. mass))^2)" 
      write(17, '(a, " Format : ''(2i10,3ES24.15E3)''")') trim(text)

    endif


    call readEnergyTable(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikGlobal, isp, order, energyTableDir, dE)

    do ibf = iBandFinit, iBandFfinal
      do ibi = iBandIinit, iBandIfinal
        
        if(order == 0) then
          write(17, 1001) ibf, ibi, Ufi(ibf,ibi,ikLocal,isp), abs(Ufi(ibf,ibi,ikLocal,isp))**2, abs(dE(2,ibf,ibi)*Ufi(ibf,ibi,ikLocal,isp))**2
        else if(order == 1) then
          write(17, 1001) ibf, ibi, Ufi(ibf,ibi,ikLocal,isp), abs(Ufi(ibf,ibi,ikLocal,isp))**2, abs(dE(3,ibf,ibi)*Ufi(ibf,ibi,ikLocal,isp)/dq_j)**2
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
  subroutine readUfis(ikLocal,isp)
    
    implicit none
    
    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel

    ! Local variables:
    integer :: ikGlobal
      !! Current global k-point
    
    integer :: ibi, ibf, totalNumberOfElements, iDum, i
    real(kind = dp) :: rDum
    complex(kind = dp):: cUfi


    ikGlobal = ikLocal+ikStart_pool-1
    
    open(17, file=trim(getMatrixElementFName(ikGlobal, isp, elementsPath)), status='unknown')
    
    read(17, *) 
    read(17, *) 
    read(17,'(5i10)') totalNumberOfElements, iDum, iDum, iDum, iDum
    read(17, *) 
    
    do i = 1, totalNumberOfElements
      
      read(17, 1001) ibf, ibi, cUfi, rDum, rDum
      Ufi(ibf,ibi,ikLocal,isp) = cUfi
          
    enddo
    
    close(17)
    
    write(*, '("    Ufi(:,:) of k-point ", i4, " and spin ", i1, " read from file.")') ikGlobal, isp
    
 1001 format(2i10,4ES24.15E3)
    
    return
    
  end subroutine readUfis
  
!----------------------------------------------------------------------------
  subroutine calculateVfiElements()
    !! THIS SUBROUTINE IS NOT UP TO DATE!! DO NOT USE IT!
    
    implicit none
    
    integer :: ikLocal, ikGlobal, ib, nOfEnergies, iE, isp
    
    real(kind = dp) :: eMin, eMax, E, av, sd, x, EiMinusEf, A, DHifMin
    real(kind=dp) :: eigvI(iBandIinit:iBandIfinal), eigvF(iBandFinit:iBandFfinal)
    
    real(kind = dp), allocatable :: sumWk(:), sAbsVfiOfE2(:), absVfiOfE2(:)
    integer, allocatable :: nKsInEbin(:)
    
    character (len = 300) :: text
    character (len = 300) :: fNameBase
    character (len = 300) :: fNameSK

    call exitError('calculateVfiElements', 'K-point binning is out of date! Do not use!', 1)

    allocate(DE(iBandIinit:iBandIfinal, nKPerPool, nSpins))
    allocate(absVfi2(iBandIinit:iBandIfinal, nKPerPool, nSpins))
     
    DE(:,:,:) = 0.0_dp
    absVfi2(:,:,:) = 0.0_dp 
    
    do isp = 1, nSpins
  
      do ikLocal = 1, nKPerPool

        ikGlobal = ikLocal+ikStart_pool-1
      
        eigvI(:) = 0.0_dp
        eigvF(:) = 0.0_dp
      
        !call readEigenvalues(ikGlobal, isp)
      
        do ib = iBandIinit, iBandIfinal

          EiMinusEf = eigvI(ib) - eigvF(iBandFinit)
          absVfi2(ib,ikLocal,isp) = EiMinusEf**2*( abs(Ufi(iBandFinit,ib,ikLocal,isp))**2 - abs(Ufi(iBandFinit,ib,ikLocal,isp))**4 )
        
          DE(ib,ikLocal,isp) = sqrt(EiMinusEf**2 - 4.0_dp*absVfi2(ib,ikLocal,isp))

        enddo

      enddo
      
    enddo
    
    eMin = minval(DE(:,:,:))
    call MPI_ALLREDUCE(MPI_IN_PLACE, eMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, interPoolComm, ierr)

    eMax = maxval(DE(:,:,:))
    call MPI_ALLREDUCE(MPI_IN_PLACE, eMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, interPoolComm, ierr)

    nOfEnergies = int((eMax-eMin)/eBin) + 1
    
    allocate(absVfiOfE2(0:nOfEnergies), nKsInEbin(0:nOfEnergies), sumWk(0:nOfEnergies))
    
    absVfiOfE2(:) = 0.0_dp
    nKsInEbin(:) = 0
    sumWk(:) = 0.0_dp
    DHifMin = 0.0_dp
    
    do isp = 1, nSpins
      
      do ikLocal = 1, nKPerPool

        ikGlobal = ikLocal+ikStart_pool-1
      
        do ib = iBandIinit, iBandIfinal
        
          if(abs(eMin - DE(ib, ikLocal, isp)) < 1.0e-3_dp) DHifMin = absVfi2(ib, ikLocal,isp)

          iE = int((DE(ib, ikLocal, isp) - eMin)/eBin)

          if(absVfi2(ib, ikLocal, isp) > 0.0_dp) then

            absVfiOfE2(iE) = absVfiOfE2(iE) + wkPC(ikGlobal)*absVfi2(ib, ikLocal, isp)
  
            sumWk(iE) = sumWk(iE) + wkPC(ikGlobal)

            nKsInEbin(iE) = nKsInEbin(iE) + 1

          else
            write(*,*) 'absVfi2', absVfi2(ib, ikLocal, isp)
          endif
        
        enddo
      
      enddo

    enddo

    call MPI_ALLREDUCE(MPI_IN_PLACE, DHifMin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, interPoolComm, ierr)
    call mpiSumDoubleV(absVfiOfE2, interPoolComm)
    call mpiSumDoubleV(sumWk, interPoolComm)
    call mpiSumIntV(nKsInEbin, interPoolComm)

    allocate(sAbsVfiOfE2(0:nOfEnergies))
    
    sAbsVfiOfE2 = 0.0_dp
    
    do isp = 1, nSpins
  
      do ikLocal = 1, nKPerPool

        ikGlobal = ikLocal+ikStart_pool-1
    
        open(11, file=trim(VfisOutput)//'ofKpt.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal)), status='unknown')
      
        do ib = iBandIinit, iBandIfinal
        
          iE = int((DE(ib, ikLocal, isp) - eMin)/eBin)

          av = absVfiOfE2(iE)/sumWk(iE)

          x = absVfi2(ib, ikLocal, isp)

          write(11, '(2ES24.15E3,2i10)') (eMin + iE*eBin), x, isp, ikGlobal
          !write(12, '(2ES24.15E3,i10)') DE(ib,ik), absVfi2(ib, ik), ik

          sAbsVfiOfE2(iE) = sAbsVfiOfE2(iE) + wkPC(ikGlobal)*(x - av)**2/sumWk(iE)
        
        enddo

        close(11)
      
      enddo

    enddo

    call mpiSumDoubleV(sAbsVfiOfE2, interPoolComm)
    
    if(myPoolId == 0) then

      open(11, file=trim(VfisOutput)//'ofKpt', status='unknown')
    
      write(11, '("# |<f|V|i>|^2 versus energy for all the k-points.")')
      write(text, '("# Energy (shifted by eBin/2) (Hartree), |<f|V|i>|^2 (Hartree)^2,")')
      write(11, '(a, " spin index, k-point index. Format : ''(2ES24.15E3,,i2,i10)''")') trim(text)

      close(11)

      do isp = 1, nSpins

        do ikGlobal = 1, nKPoints

          fNameBase = trim(VfisOutput)//'ofKpt'
          fNameSK = trim(fNameBase)//'.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal))

          call execute_command_line('cat '//trim(fNameSK)//' >> '//trim(fNameBase)//'&& rm '//trim(fNameSK))

        enddo

      enddo

      open(63, file=trim(VfisOutput), status='unknown')
    
      write(63, '("# Averaged |<f|V|i>|^2 over K-points versus energy.")')
      write(63, '("#                 Cell volume : ", ES24.15E3, " (a.u.)^3,   Format : ''(ES24.15E3)''")') omega
      write(63, '("#   Minimun transition energy : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eMin
      write(63, '("# |DHif|^2 at minimum Tr. En. : ", ES24.15E3, " (Hartree^2),Format : ''(ES24.15E3)''")') DHifMin
      write(63, '("#                  Energy bin : ", ES24.15E3, " (Hartree),  Format : ''(ES24.15E3)''")') eBin
      write(text, '("# Energy (Hartree), averaged |<f|V|i>|^2 over K-points (Hartree)^2,")')
      write(63, '(a, " standard deviation (Hartree)^2. Format : ''(3ES24.15E3)''")') trim(text)
    
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
    
      close(63)
    
    endif  
    
    return
    
  end subroutine calculateVfiElements
   
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
  function getMatrixElementFName(ikGlobal, isp, path) result(fName)

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

  end function getMatrixElementFName
  
end module TMEmod
