module PhononPPMod
  
  use constants, only: dp, angToBohr, angToM, daltonToElecM, elecMToKg, THzToHz, pi, hbar
  use base, only: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
  use cell, only: nAtoms, omega, realLattVec, cartDispProjOnPhononEigsNorm, readPOSCAR, writePOSCARNewPos
  use miscUtilities, only: int2str, int2strLeadZero
  use generalComputations, only: direct2cart
  use errorsAndMPI
  use mpi

  implicit none

  ! Variables that should not be passed as arguments:
  integer :: iModeStart, iModeEnd
    !! Start and end mode for this process

  real(kind=dp) :: t0, t1, t2
    !! Timers

  ! Variables that should be passed as arguments:
  integer :: disp2AtomInd(2)
    !! Index of atoms to check displacement
    !! between if calcMaxDisp
  integer :: ikIinit, ikIfinal, ikFinit, ikFfinal
    !! K-point bounds for initial and final state
  integer :: nModes
    !! Number of phonon modes
  integer :: nModesLocal
    !! Local number of phonon modes


  real(kind=dp), allocatable :: coordFromPhon(:,:)
    !! Corodinates from phonon file
  real(kind=dp), allocatable :: eigenvector(:,:,:)
    !! Eigenvectors for each atom for each mode
  real(kind=dp) :: freqThresh
    !! Threshold for frequency to determine if the mode 
    !! should be skipped or not
  real(kind=dp), allocatable :: mass(:)
    !! Masses of atoms
  real(kind=dp), allocatable :: omegaFreq(:)
    !! Frequency for each mode
  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors

  character(len=300) :: basePOSCARFName
    !! File name for intial POSCAR to calculate shift from
  character(len=300) :: CONTCARsBaseDir
    !! Base dir for sets of relaxed files if not captured
  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: initPOSCARFName, finalPOSCARFName
    !! File name for POSCAR for relaxed initial and final charge states
  character(len=300) :: phononFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: prefix
    !! Prefix for shifted POSCARs

  logical :: calcDq
    !! If dq output file should be generated
  logical :: calcMaxDisp
    !! If calculating the maximum displacement between a pair of atoms
  logical :: calcSj
    !! If Sj should be calculated
  logical :: singleDisp
    !! If there is just a single displacement to consider
  logical :: generateShiftedPOSCARs
    !! If shifted POSCARs should be generated


  namelist /inputParams/ initPOSCARFName, finalPOSCARFName, phononFName, prefix, shift, dqFName, generateShiftedPOSCARs, singleDisp, &
                         iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, CONTCARsBaseDir, basePOSCARFName, freqThresh, calcSj, &
                         calcDq, calcMaxDisp, disp2AtomInd, ikIinit, ikIfinal, ikFinit, ikFfinal


  contains

!----------------------------------------------------------------------------
  subroutine readInputs(disp2AtomInd, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, freqThresh, &
        shift, basePOSCARFName, CONTCARsBaseDir, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, & 
        calcMaxDisp, calcSj, generateShiftedPOSCARs, singleDisp)

    implicit none

    ! Output variables:
    integer, intent(out) :: disp2AtomInd(2)
      !! Index of atoms to check displacement
      !! between if calcMaxDisp
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state

    real(kind=dp), intent(out) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(out) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(out) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states
    character(len=300), intent(out) :: prefix
      !! Prefix for shifted POSCARs

    logical, intent(out) :: calcDq
      !! If dq output file should be generated
    logical, intent(out) :: calcMaxDisp
      !! If calculating the maximum displacement between a pair of atoms
    logical, intent(out) :: calcSj
      !! If Sj should be calculated
    logical, intent(out) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
    logical, intent(out) :: singleDisp
      !! If there is just a single displacement to consider


    if(ionode) then

      call initialize(disp2AtomInd, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, &
            freqThresh, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, &
            prefix, calcDq, calcMaxDisp, calcSj, generateShiftedPOSCARs, singleDisp)
        ! Set default values for input variables and start timers
    
      read(5, inputParams, iostat=ierr)
        !! * Read input variables
    
      if(ierr /= 0) call exitError('readInputs', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(disp2AtomInd, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, &
          freqThresh, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, &
          calcDq, calcMaxDisp, calcSj, generateShiftedPOSCARs, singleDisp)

    endif

    ! Send to other processes only what they need to know
    call MPI_BCAST(disp2AtomInd, 2, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(iBandFfinal, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(ikIinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikIfinal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikFinit, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(ikFfinal, 1, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(freqThresh, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(shift, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(basePOSCARFName, len(basePOSCARFName), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(CONTCARsBaseDir, len(CONTCARsBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(finalPOSCARFName, len(finalPOSCARFName), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(initPOSCARFName, len(initPOSCARFName), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(calcDq, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(calcSj, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(calcMaxDisp, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(generateShiftedPOSCARs, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(singleDisp, 1, MPI_LOGICAL, root, worldComm, ierr)

    return

  end subroutine readInputs

!----------------------------------------------------------------------------
  subroutine initialize(disp2AtomInd, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, &
        freqThresh, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, &
        prefix, calcDq, calcMaxDisp, calcSj, generateShiftedPOSCARs, singleDisp)
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    ! Output variables:
    integer, intent(out) :: disp2AtomInd(2)
      !! Index of atoms to check displacement
      !! between if calcMaxDisp
    integer, intent(out) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(out) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state

    real(kind=dp), intent(out) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(out) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(out) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(out) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states
    character(len=300), intent(out) :: prefix
      !! Prefix for shifted POSCARs

    logical, intent(out) :: calcDq
      !! If dq output file should be generated
    logical, intent(out) :: calcMaxDisp
      !! If calculating the maximum displacement between a pair of atoms
    logical, intent(out) :: calcSj
      !! If Sj should be calculated
    logical, intent(out) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
    logical, intent(out) :: singleDisp
      !! If there is just a single displacement to consider


    disp2AtomInd = -1
    iBandIinit  = -1
    iBandIfinal = -1
    iBandFinit  = -1
    iBandFfinal = -1

    ikIinit  = -1
    ikIfinal = -1
    ikFinit  = -1
    ikFfinal = -1

    CONTCARsBaseDir = ''
    dqFName = 'dq.txt'
    basePOSCARFName = 'POSCAR_init'
    initPOSCARFName = 'POSCAR_init'
    finalPOSCARFName = 'POSCAR_final'
    phononFName = 'mesh.yaml'
    prefix = 'ph_POSCAR'

    freqThresh = 0.5_dp
    shift = 0.01_dp

    calcDq = .true.
    calcSj = .true.
    calcMaxDisp = .false.
    generateShiftedPOSCARs = .true.
    singleDisp = .true.

    return

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(disp2AtomInd, iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, &
      freqThresh, shift, basePOSCARFName, CONTCARsBaseDir, dqFName, phononFName, finalPOSCARFName, initPOSCARFName, prefix, &
      calcDq, calcMaxDisp, calcSj, generateShiftedPOSCARs, singleDisp)

    implicit none

    ! Input variables:
    integer, intent(in) :: disp2AtomInd(2)
      !! Index of atoms to check displacement
      !! between if calcMaxDisp
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state

    real(kind=dp), intent(in) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors

    character(len=300), intent(inout) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(in) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: phononFName
      !! File name for mesh.yaml phonon file
    character(len=300), intent(in) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states
    character(len=300), intent(in) :: prefix
      !! Prefix for shifted POSCARs

    logical, intent(in) :: calcDq
      !! If dq output file should be generated
    logical, intent(in) :: calcMaxDisp
      !! If calculating the maximum displacement between a pair of atoms
    logical, intent(in) :: calcSj
      !! If Sj should be calculated
    logical, intent(in) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
    logical, intent(in) :: singleDisp
      !! If there is just a single displacement to consider

    ! Local variables:
    integer :: ibi, ibf, iki, ikf
      !! Loop indices

    logical :: abortExecution
      !! Whether or not to abort the execution


    if(.not. (calcDq .or. calcMaxDisp .or. calcSj .or. generateShiftedPOSCARs)) &
      call exitError('checkInitialization','You must choose at least one option to run this program!',1)


    ! Must-haves for reading phonons:
    abortExecution = checkFileInitialization('phononFName', phononFName)
    abortExecution = checkDoubleInitialization('freqThres', freqThresh, 0.0_dp, 1.0_dp) .or. abortExecution


    ! Need for calculating Sj:
    write(*,'("calcSj = ''",L1,"''")') calcSj
    if(calcSj) then

      write(*,'("singleDisp = ''",L1,"''")') singleDisp

      if(singleDisp) then

        abortExecution = checkFileInitialization('initPOSCARFName', initPOSCARFName) .or. abortExecution
        abortExecution = checkFileInitialization('finalPOSCARFName', finalPOSCARFName) .or. abortExecution

      else

        abortExecution = checkIntInitialization('iBandIinit', iBandIinit, 1, int(1e9)) .or. abortExecution
        abortExecution = checkIntInitialization('iBandIfinal', iBandIfinal, iBandIinit, int(1e9)) .or. abortExecution
        abortExecution = checkIntInitialization('iBandFinit', iBandFinit, 1, int(1e9)) .or. abortExecution
        abortExecution = checkIntInitialization('iBandFfinal', iBandFfinal, iBandFinit, int(1e9)) .or. abortExecution 
        
        abortExecution = checkIntInitialization('ikIinit', ikIinit, 0, int(1e9)) .or. abortExecution
        abortExecution = checkIntInitialization('ikIfinal', ikIfinal, ikIinit, int(1e9)) .or. abortExecution
        abortExecution = checkIntInitialization('ikFinit', ikFinit, 0, int(1e9)) .or. abortExecution
        abortExecution = checkIntInitialization('ikFfinal', ikFfinal, ikFinit, int(1e9)) .or. abortExecution 

        do iki = ikIinit, ikIfinal
          do ibi = iBandIinit, iBandIfinal
            abortExecution = checkFileInitialization('CONTCARsBaseDir/k'//trim(int2str(iki))//'_b'//trim(int2str(ibi))//'/CONTCAR', &
                             trim(CONTCARsBaseDir)//'/k'//trim(int2str(iki))//'_b'//trim(int2str(ibi))//'/CONTCAR') .or. abortExecution
          enddo
        enddo

        do ikf = ikFinit, ikFfinal
          do ibf = iBandFinit, iBandFfinal
            ! Only test the existence of files that haven't already been tested
            if((ikf < ikIinit .or. ikf > ikIfinal) .or. (ibf < iBandIinit .or. ibf > iBandIfinal)) & 
              abortExecution = checkFileInitialization('CONTCARsBaseDir/k'//trim(int2str(ikf))//'_b'//trim(int2str(ibf))//'/CONTCAR', &
                             trim(CONTCARsBaseDir)//'/k'//trim(int2str(ikf))//'_b'//trim(int2str(ibf))//'/CONTCAR') .or. abortExecution

          enddo
        enddo

      endif
    endif


    ! Need for dq:
    write(*,'("calcDq = ''",L1,"''")') calcDq
    if(calcDq) abortExecution = checkStringInitialization('dqFName', dqFName) .or. abortExecution

    ! Need for shifted POSCARs:
    write(*,'("generateShiftedPOSCARs = ''",L1,"''")') generateShiftedPOSCARs
    if(generateShiftedPOSCARs) abortExecution = checkStringInitialization('prefix', prefix) .or. abortExecution

    ! Needed for dq and shifted POSCARs and calculating max displacement:
    if(calcDq .or. generateShiftedPOSCARs .or. calcMaxDisp) then
      abortExecution = checkDoubleInitialization('shift', shift, 0.0_dp, 10.0_dp) .or. abortExecution
      abortExecution = checkFileInitialization('basePOSCARFName', basePOSCARFName) .or. abortExecution
    endif

    ! Need for calculating max displacement:
    if(calcMaxDisp) then
      abortExecution = checkIntInitialization('disp2AtomInd(1)', disp2AtomInd(1), 1, int(1e9)) .or. abortExecution
      abortExecution = checkIntInitialization('disp2AtomInd(2)', disp2AtomInd(2), 1, int(1e9)) .or. abortExecution
    endif

    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif
    
    return

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine standardizeCoordinates(nAtoms, atomPositionsDir)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms in system

    ! Output variables:
    real(kind=dp), intent(inout) :: atomPositionsDir(3,nAtoms)
      !! Atom positions in direct coordinates


    where(atomPositionsDir(:,:) < 0)
      atomPositionsDir(:,:) = 1.0_dp + atomPositionsDir(:,:)
    end where

    return

  end subroutine standardizeCoordinates

!----------------------------------------------------------------------------
  subroutine readPhonons(freqThresh, phononFName, nAtoms, nModes, coordFromPhon, eigenvector, mass, omegaFreq)

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not

    character(len=300), intent(in) :: phononFName
      !! File name for mesh.yaml phonon file

    ! Output variables:
    integer, intent(out) :: nAtoms
      !! Number of atoms
    integer, intent(out) :: nModes
      !! Number of modes

    real(kind=dp), allocatable, intent(out) :: coordFromPhon(:,:)
      !! Corodinates from phonon file
    real(kind=dp), allocatable, intent(out) :: eigenvector(:,:,:)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), allocatable, intent(out) :: mass(:)
      !! Masses of atoms
    real(kind=dp), allocatable, intent(out) :: omegaFreq(:)
      !! Frequency for each mode

    ! Local variables:
    integer :: j, ia, ix
      !! Loop index
    integer :: nUsed
      !! Number of phonon modes used based on frequency threshold

    real(kind=dp) :: omegaFreq_j
      !! Frequency for a single mode
    real(kind=dp) :: qPos(3)
      !! Phonon q position

    character(len=300) :: line
      !! Line read from file

    logical :: fileExists
      !! If phonon file exists


    if(ionode) then

      inquire(file=phononFName, exist=fileExists)

      if(.not. fileExists) call exitError('readPhonons', 'Phonon file '//trim(phononFName)//' does not exist', 1)

      open(57, file=phononFName)


      line = getFirstLineWithKeyword(57,'natom')
      read(line(7:len(trim(line))),*) nAtoms

    endif


    call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
    allocate(coordFromPhon(3,nAtoms))
    allocate(mass(nAtoms))


    if(ionode) then

      line = getFirstLineWithKeyword(57,'points')
        !! Ignore everything next until you get to points line

      do ia = 1, nAtoms

        read(57,'(A)') ! Ignore symbol 
        read(57,'(A)') line
        read(line(17:len(trim(line))),*) coordFromPhon(:,ia)
        read(57,'(a7,f)') line, mass(ia)
          !! Read mass

      enddo

    endif

    call MPI_BCAST(coordFromPhon, size(coordFromPhon), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(mass, size(mass), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


    nModes = 3*nAtoms - 3

    allocate(eigenvector(3,nAtoms,nModes))
    allocate(omegaFreq(nModes))


    if(ionode) then

      line = getFirstLineWithKeyword(57,'q-position')
        !! Ignore everything next until you get to q-position line

      read(line(16:len(trim(line))),*) qPos
        !! Read in the q position

      if(qPos(1) > 1e-8 .or. qPos(2) > 1e-8 .or. qPos(3) > 1e-8) &
        call exitError('readPhonons', 'Code assumes phonons have no momentum (at q=0)!',1)

      read(57,'(A)') line
      read(57,'(A)') line
      read(57,'(A)') line
        !! Ignore next 3 lines


      nUsed = 0
      do j = 1, nModes+3

        read(57,'(A)') ! Ignore mode number


        read(57,'(A)') line
        read(line(15:len(trim(line))),*) omegaFreq_j
        if(abs(omegaFreq_j) > freqThresh) then
          nUsed = nUsed + 1
          omegaFreq(nUsed) = omegaFreq_j

          if(omegaFreq_j < 0) write(*,'("WARNING: Negative frequency detected in mode ", i7,&
            "! Re-relax the structure using the POSCAR shifted along this mode''s eigenvector.")') j
        else
          write(*,'("Skipping mode ", i7, " with freqency ", ES24.15E3,".")') j, omegaFreq_j
        endif


        read(57,'(A)') ! Ignore eigenvector section header

        do ia = 1, nAtoms

          read(57,'(A)') line ! Ignore atom number
  
          do ix = 1, 3

            read(57,'(A)') line
            if(abs(omegaFreq_j) > freqThresh)  read(line(10:len(trim(line))),*) eigenvector(ix,ia,nUsed)
              !! Only store the eigenvectors if the mode isn't skipped
              !! because those are the translational modes 

          enddo
        enddo
      enddo

      if(nUsed /= nModes) call exitError('readPhonons', 'Did not skip 3 modes based on threshold set. Check threshold.', 1)

      omegaFreq(:) = omegaFreq(:)*2*pi

      close(57)

    endif

    call MPI_BCAST(eigenvector, size(eigenvector), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omegaFreq, size(omegaFreq), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine readPhonons

!----------------------------------------------------------------------------
  subroutine calculateSj(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikIinit, ikIfinal, ikFinit, ikFfinal, nAtoms, nModes, &
          coordFromPhon, eigenvector, mass, omegaFreq, singleDisp, CONTCARsBaseDir, initPOSCARFName, finalPOSCARFName)

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
      !! Energy band bounds for initial and final state
    integer, intent(in) :: ikIinit, ikIfinal, ikFinit, ikFfinal
      !! K-point bounds for initial and final state
    integer, intent(inout) :: nAtoms
      !! Number of atoms in intial system
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: coordFromPhon(3,nAtoms)
      !! Corodinates from phonon file
    real(kind=dp), intent(in) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms
    real(kind=dp), intent(in) :: omegaFreq(nModes)
      !! Frequency for each mode

    logical, intent(in) :: singleDisp
      !! If there is just a single displacement to consider

    character(len=300), intent(in) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(inout) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states

    ! Local variables:
    integer :: ibi, ibf, iki, ikf
      !! Loop indices

    real(kind=dp), allocatable :: atomPositionsDirInit(:,:)
      !! Atom positions in initial relaxed positions
    real(kind=dp) :: displacement(3,nAtoms)
      !! Displacement vector; not used here, only output 
      !! from check compatibility
    real(kind=dp) :: omega
      !! Volume of supercell
    real(kind=dp) :: realLattVec(3,3)
      !! Real space lattice vectors

    character(len=300) :: SjFName
      !! File name for the Sj output


    if(singleDisp) then
      if(ionode) then
        call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
        call standardizeCoordinates(nAtoms, atomPositionsDirInit)
      endif

      call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)

      if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
      call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

      call getRelaxDispAndCheckCompatibility(nAtoms, coordFromPhon, atomPositionsDirInit, displacement)

      SjFName = 'Sj.out'
      call getSingleDisp(nAtoms, nModes, atomPositionsDirInit, eigenvector, mass, omega, omegaFreq, finalPOSCARFName, SjFName)

      deallocate(atomPositionsDirInit)

    else

      do iki = ikIinit, ikIfinal
        ! I do a loop over all of the bands. Depending on the band selection,
        ! this could result in duplicate pairs (e.g., .1.2 and .2.1), but I 
        ! can't think of a general way to exlude these pairs without a significant
        ! amount of work to track the pairs that have already been calcualted.
        ! Any solution I can think of would not hold for both hole and electron
        ! capture and/or different ranges of bands selected by the user. For
        ! now, I just have the code output all of the duplicate pairs and the
        ! user/LSF code can use whatever they need from that.
        do ibi = iBandIinit, iBandIfinal

          ! Get the initial positions for this band 
          if(ionode) then
            initPOSCARFName = trim(CONTCARsBaseDir)//'/k'//trim(int2str(iki))//'_b'//trim(int2str(ibi))//'/CONTCAR'

            call readPOSCAR(initPOSCARFName, nAtoms, atomPositionsDirInit, omega, realLattVec)
            call standardizeCoordinates(nAtoms, atomPositionsDirInit)
          endif

          call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)

          if(.not. ionode) allocate(atomPositionsDirInit(3,nAtoms))
          call MPI_BCAST(atomPositionsDirInit, size(atomPositionsDirInit), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

          call getRelaxDispAndCheckCompatibility(nAtoms, coordFromPhon, atomPositionsDirInit, displacement)

          do ikf = ikFinit, ikFfinal
            do ibf = iBandFinit, iBandFfinal

              if(iki /= ikf .or. ibf /= ibi) then
                ! Set the final positions and output file names
                finalPOSCARFName = trim(CONTCARsBaseDir)//'/k'//trim(int2str(ikf))//'_b'//trim(int2str(ibf))//'/CONTCAR'
                SjFName = 'Sj.k'//trim(int2str(iki))//'_b'//trim(int2str(ibi))//'.k'&
                                //trim(int2str(ikf))//'_b'//trim(int2str(ibf))//'.out'

                ! Get displacement and output ranked Sj for a single pair of states
                call getSingleDisp(nAtoms, nModes, atomPositionsDirInit, eigenvector, mass, omega, omegaFreq, finalPOSCARFName, SjFName)

              endif
            enddo
          enddo

          ! Positions are allocated in readPOSCAR, so we need to 
          ! deallocate them after every loop
          deallocate(atomPositionsDirInit)
        enddo
      enddo

    endif


    return

  end subroutine calculateSj

!----------------------------------------------------------------------------
  subroutine getSingleDisp(nAtoms, nModes, atomPositionsDirInit, eigenvector, mass, omega, omegaFreq, finalPOSCARFName, SjFName)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms in intial system
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: atomPositionsDirInit(3,nAtoms)
      !! Atom positions in initial relaxed positions
    real(kind=dp), intent(in) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms
    real(kind=dp), intent(in) :: omega
      !! Volume of initial-state supercell
    real(kind=dp), intent(in) :: omegaFreq(nModes)
      !! Frequency for each mode

    character(len=300), intent(in) :: finalPOSCARFName
      !! File name for CONTCAR for relaxed final state
    character(len=300), intent(in) :: SjFName
      !! File name for the Sj output

    ! Local variables:
    integer :: j
      !! Loop index
    integer :: nAtomsFinal
      !! Number of atoms in final system

    real(kind=dp), allocatable :: atomPositionsDirFinal(:,:)
      !! Atom positions in final relaxed positions
    real(kind=dp) :: displacement(3,nAtoms)
      !! Displacement 
    real(kind=dp) :: omegaFinal
      !! Volume of final-state supercell
    real(kind=dp) :: projNorm(nModes)
      !! Generalized norms after displacement
    real(kind=dp) :: realLattVec(3,3)
      !! Real space lattice vectors

    
    if(ionode) then
      call readPOSCAR(finalPOSCARFName, nAtomsFinal, atomPositionsDirFinal, omegaFinal, realLattVec)
        ! I define realLattVec locally here because I don't use the input value
        ! nor do I want to output the value from here. It is assumed for now
        ! that the lattice vectors for the two systems are compatible. This
        ! would be a good test to add at some point in the future.

      call standardizeCoordinates(nAtomsFinal, atomPositionsDirFinal)

      if(nAtoms /= nAtomsFinal) &
        call exitError('getSingleDisp', 'number of atoms does not match: '//trim(int2str(nAtoms))//' '//trim(int2str(nAtomsFinal)), 1)

      if(abs(omega - omegaFinal) > 1e-8) call exitError('getSingleDisp', 'volumes don''t match', 1)

    endif


    if(.not. ionode) allocate(atomPositionsDirFinal(3,nAtoms))
    call MPI_BCAST(atomPositionsDirFinal, size(atomPositionsDirFinal), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
  

    ! Define the displacement for the relaxation, 
    ! project onto the phonon eigenvectors, and
    ! get Sj
    call getRelaxDispAndCheckCompatibility(nAtoms, atomPositionsDirFinal, atomPositionsDirInit, displacement)
  
    projNorm = 0.0_dp
    do j = iModeStart, iModeEnd

      projNorm(j) = cartDispProjOnPhononEigsNorm(nAtoms, displacement, eigenvector(:,:,j), mass, realLattVec)

    enddo

    projNorm = projNorm*angToM*sqrt(daltonToElecM*elecMToKg)
    call calcAndWriteSj(nModes, omegaFreq, projNorm, SjFName)


    deallocate(atomPositionsDirFinal)

    return

  end subroutine getSingleDisp

!----------------------------------------------------------------------------
  subroutine getRelaxDispAndCheckCompatibility(nAtoms, atomPositionsDirFinal, atomPositionsDirInit, displacement)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: atomPositionsDirInit(3,nAtoms), atomPositionsDirFinal(3,nAtoms)
      !! Atom positions in initial and final relaxed positions

    ! Output variables:
    real(kind=dp), intent(out) :: displacement(3,nAtoms)
      !! Displacement vector

    ! Local variables:
    integer :: ia, ix

    real(kind=dp) :: dispInner, dispOuter
      !! Inner (not crossing boundary) and outer
      !! (crossing boundary) distance between two 
      !! points
    real(kind=dp) :: posInit, posFinal
      !! Local storage of single coordinate

    logical :: abortExecution
      !! Whether or not to abort execution


    if(ionode) then

      abortExecution = .false.

      do ia= 1, nAtoms

        do ix = 1, 3

          posInit = atomPositionsDirInit(ix,ia)
          posFinal = atomPositionsDirFinal(ix,ia)

          dispInner = posFinal - posInit

          if(posFinal > posInit) then
            dispOuter = -posInit + posFinal - 1
          else
            dispOuter = 1 - posInit + posFinal
          endif

          if(abs(dispInner) < abs(dispOuter)) then
            displacement(ix,ia) = dispInner
          else
            displacement(ix,ia) = dispOuter
          endif

        enddo

        if(maxval(abs(displacement(:,ia))) > 0.2) then
          abortExecution = .true.
          write(*,'("Large displacement detected. Atom possibly out of order: ", i5, 3f10.7)') ia, displacement(:,ia)
        endif

      enddo

      if(abortExecution) then
        write(*, '(" Program stops!")')
        stop
      endif

    endif

    call MPI_BCAST(displacement, size(displacement), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return 

  end subroutine getRelaxDispAndCheckCompatibility

!----------------------------------------------------------------------------
  function centerOfMass(nAtoms, coords, mass) result(centerOfMassCoords)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: coords(3,nAtoms)
      !! Coordinates of atoms
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms

    ! Output variables:
    real(kind=dp) :: centerOfMassCoords(3,nAtoms)
      !! Coordinates of the center of mass

    ! Local variables:
    integer :: ia
      !! Loop index

    real(kind=dp) :: totalMass
      !! Total mass of atoms
    real(kind=dp) :: totalMassTimesPos(3)
      !! Total of mass times positions


    totalMass = sum(mass(:))

    totalMassTimesPos = 0.0_dp

    do ia = 1, nAtoms
      totalMassTimesPos = totalMassTimesPos + mass(ia) * coords(:, ia)
    end do

    centerOfMassCoords(1,:) = totalMassTimesPos(1)/totalMass
    centerOfMassCoords(2,:) = totalMassTimesPos(2)/totalMass
    centerOfMassCoords(3,:) = totalMassTimesPos(3)/totalMass

  end function centerOfMass
  
!----------------------------------------------------------------------------
  subroutine calcAndWriteSj(nModes, omegaFreq, projNorm, SjFName)

    use miscUtilities, only: hpsort_eps

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: omegaFreq(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(inout) :: projNorm(nModes)
      !! Generalized norms after displacement

    character(len=300), intent(in) :: SjFName
      !! File name for the Sj output

    ! Local variables:
    integer :: j, jSort
      !! Loop index
    integer :: modeIndex(nModes)
      !! Track mode indices after sorting

    real(kind=dp) :: Sj(nModes)
      !! Hold the Sj values to be sorted


    call mpiSumDoubleV(projNorm, worldComm)
      !! * Get the generalized-displacement norms
      !!   from all processes

      if(ionode) then

      do j = 1, nModes

        modeIndex(j) = j

        Sj(j) = projNorm(j)**2*omegaFreq(j)*THzToHz/(2*hbar) 

      enddo

      call hpsort_eps(nModes, Sj, modeIndex, 1e-14_dp)
        ! Sort in ascending order

      open(60, file=trim(SjFName))

      write(60,'(1i7)') nModes

      do j = 1, nModes

        jSort = modeIndex(nModes-(j-1))

        write(60,'(1i7, 2ES24.15E3)') modeIndex(jSort), Sj(nModes-(j-1)), omegaFreq(jSort)
          ! Write out in descending order

      enddo

      close(60)

    endif


    return

  end subroutine calcAndWriteSj

!----------------------------------------------------------------------------
  subroutine calculateShiftAndDq(disp2AtomInd, nAtoms, nModes, coordFromPhon, eigenvector, mass, shift, calcDq, calcMaxDisp, &
        generateShiftedPOSCARs, basePOSCARFName, dqFName, prefix)

    implicit none

    ! Input variables:
    integer, intent(in) :: disp2AtomInd(2)
      !! Index of atoms to check displacement
      !! between if calcMaxDisp
    integer, intent(inout) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: coordFromPhon(3,nAtoms)
      !! Corodinates from phonon file
    real(kind=dp), intent(in) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors

    logical, intent(in) :: calcDq
      !! If dq output file should be generated
    logical, intent(in) :: calcMaxDisp
      !! If calculating the maximum displacement between a pair of atoms
    logical, intent(in) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
      
    character(len=300), intent(in) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: prefix
      !! Prefix for shifted POSCARs

    ! Local variables:
    integer :: j
      !! Loop index
    integer :: suffixLength
      !! Length of shifted POSCAR file suffix

    real(kind=dp), allocatable :: atomPositionsDirBase(:,:)
      !! Atom positions in initial relaxed positions
    real(kind=dp), allocatable :: displacement(:,:)
      !! Atom displacements in angstrom
    real(kind=dp) :: omega
      !! Volume of supercell
    real(kind=dp), allocatable :: projNorm(:)
      !! Generalized norms after displacement
    real(kind=dp) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp) :: relDisp(3)
      !! Relative displacement between selected atoms
      !! if calculating the maximum displacement
    real(kind=dp), allocatable :: relDispMag(:)
      !! Magnitude of relative displacement between 
      !! selected atoms if checking the maximum
      !! displacement
    real(kind=dp), allocatable :: shiftedPositions(:,:)
      !! Positions after shift along eigenvector

    character(len=300) :: memoLine
      !! Memo line to output shift and mode number to shifted POSCARs
    character(len=300) :: shiftedPOSCARFName
      !! File name for shifted POSCAR


    if(generateShiftedPOSCARs) then
      if(nModes < 10) then
        suffixLength = 1
      else if(nModes < 100) then
        suffixLength = 2
      else if(nModes < 1000) then
        suffixLength = 3
      else if(nModes < 10000) then
        suffixLength = 4
      else if(nModes < 100000) then
        suffixLength = 5
      endif
    endif


    if(ionode) then

      call readPOSCAR(basePOSCARFName, nAtoms, atomPositionsDirBase, omega, realLattVec)
      call standardizeCoordinates(nAtoms, atomPositionsDirBase)

    endif

    call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
    if(.not. ionode) allocate(atomPositionsDirBase(3,nAtoms))
    call MPI_BCAST(atomPositionsDirBase, size(atomPositionsDirBase), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    allocate(displacement(3,nAtoms))

    call getRelaxDispAndCheckCompatibility(nAtoms, coordFromPhon, atomPositionsDirBase, displacement)


    allocate(shiftedPositions(3,nAtoms))
    allocate(relDispMag(nModes))
    allocate(projNorm(nModes))
    projNorm = 0.0_dp
    relDispMag = 0.0_dp

    write(memoLine,'("  shift = ", ES9.2E1)') shift


    ! Get the displacement for each mode to 
    ! calculate the derivative of the wave function.
    ! Project onto the phonon eigenvectors (here,
    ! the effect is just to convert back to generalized
    ! coordinates because the displacement is already
    ! a scaled form of the eigenvectors), then (if needed)
    ! write the shifted positions and generalized displacement
    ! norms.
    do j = iModeStart, iModeEnd

      displacement = getShiftDisplacement(nAtoms, eigenvector(:,:,j), realLattVec, mass, shift)

      if(calcMaxDisp) then
        relDisp = displacement(:,disp2AtomInd(1)) - displacement(:,disp2AtomInd(2))

        relDispMag(j) = sqrt(dot_product(relDisp,relDisp))
      endif

      if(calcDq) projNorm(j) = cartDispProjOnPhononEigsNorm(nAtoms, displacement, eigenvector(:,:,j), mass, realLattVec)

      if(generateShiftedPOSCARs) then

        shiftedPositions = direct2cart(nAtoms, atomPositionsDirBase, realLattVec) + displacement

        shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

        call writePOSCARNewPos(nAtoms, shiftedPositions, basePOSCARFName, shiftedPOSCARFName, trim(memoLine)//' j = '//int2str(j), .true.)

      endif

    enddo


    deallocate(atomPositionsDirBase)
    deallocate(shiftedPositions)
    deallocate(displacement)

    if(calcMaxDisp) then
      call mpiSumDoubleV(relDispMag, worldComm)

      if(ionode) write(*,'("Maximum displacement between atoms ", i4," and ",i4," is ",1ES12.3E2," in mode ",i6)') &
            disp2AtomInd(1), disp2AtomInd(2), relDispMag(maxloc(relDispMag)), maxloc(relDispMag)
    endif

    deallocate(relDispMag)


    if(calcDq) then
      projNorm = projNorm*angToBohr*sqrt(daltonToElecM)
      call writeDqs(nModes, projNorm, dqFName)
    endif


    deallocate(projNorm)

    return 

  end subroutine calculateShiftAndDq

!----------------------------------------------------------------------------
  function getShiftDisplacement(nAtoms, eigenvector, realLattVec, mass, shift) result(displacement)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: eigenvector(3,nAtoms)
      !! Eigenvectors for each atom for this mode
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Masses of atoms
    real(kind=dp), intent(in) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors

    ! Output variables:
    real(kind=dp) :: displacement(3,nAtoms)
      !! Displacements for each atom for this mode

    ! Local variables:
    integer :: ia
      !! Loop indices

    real(kind=dp) :: eig(3)
      !! Eigenvector for single mode and atom
    real(kind=dp) :: cartNorm
      !! Norm of eigenvectors in Cartesian coordinates


    !> Convert from generalized to Cartesian coordinates
    !> and scale displacement based on norm of entire
    !> displacement vector
    cartNorm = 0.0_dp
    do ia = 1, nAtoms

      eig = eigenvector(:,ia)/sqrt(mass(ia))

      displacement(:,ia) = eig

      eig = matmul(realLattVec, eig)
        ! Convert to Cartesian coordinates before 
        ! getting norm

      cartNorm = cartNorm + dot_product(eig,eig)

    enddo

    cartNorm = sqrt(cartNorm)

    displacement = displacement*shift/cartNorm
      !! Scale the displacement so that the 
      !! norm of the entire vector is `shift`

    return

  end function getShiftDisplacement

!----------------------------------------------------------------------------
  subroutine writeDqs(nModes, projNorm, dqFName)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(inout) :: projNorm(nModes)
      !! Generalized norms after displacement

    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms

    ! Local variables:
    integer :: j
      !! Loop index


    call mpiSumDoubleV(projNorm, worldComm)
      !! * Get the generalized-displacement norms
      !!   from all processes

    if(ionode) then

      open(60, file=dqFName)

      write(60,'("# Norm of generalized displacement vectors after scaling Cartesian displacement Format: ''(1i7, 1ES24.15E3)''")')

      do j = 1, nModes

        write(60,'(1i7, 1ES24.15E3)') j, projNorm(j)
          !! Output norm to file

      enddo

      close(60)

    endif
      
    return

  end subroutine writeDqs

end module PhononPPMod
