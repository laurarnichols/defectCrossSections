module PhononPPMod
  
  use constants, only: dp, angToBohr, daltonToElecM, THzToHartree, pi, hbar_atomic, kB_atomic
  use energyTabulatorMod, only: energyTableDir, readScatterEnergyTable
  use cell, only: nAtoms, volume, realLattVec
  use miscUtilities, only: int2str, int2strLeadZero
  use hungarianOptimizationMod, only: hungarianOptimalMatching
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
  integer :: nAtomsPrime
    !! Number of atoms in optional final phonon file 
  integer :: nModes, nModesPrime
    !! Number of phonon modes
  integer :: nModesLocal
    !! Local number of phonon modes


  real(kind=dp) :: beta
    !! 1/kB*T
  real(kind=dp), allocatable :: coordFromPhon(:,:), coordFromPhonPrime(:,:)
    !! Coordinates from phonon file
  real(kind=dp), allocatable :: eigenvector(:,:,:), eigenvectorPrime(:,:,:), dqEigenvectors(:,:,:)
    !! Eigenvectors for each atom for each mode 
  real(kind=dp) :: freqThresh
    !! Threshold for frequency to determine if the mode 
    !! should be skipped or not
  real(kind=dp), allocatable :: mass(:)
    !! Masses of atoms
  real(kind=dp),allocatable :: omega(:)
    !! Initial-state frequency for each mode
  real(kind=dp), allocatable :: omegaPrime(:)
    !! Optional final-state frequency for each mode
  real(kind=dp), allocatable :: Sj(:)
    !! Huang-Rhys factor for each mode with omega_j
  real(kind=dp), allocatable :: SjPrime(:)
    !! Huang-Rhys factor for each mode with omega_j'
  real(kind=dp) :: shift
    !! Magnitude of shift along phonon eigenvectors
  real(kind=dp) :: temperature    
  real(kind=dp) :: SjThresh
    !! Threshold to count modes above

  character(len=300) :: basePOSCARFName
    !! File name for intial POSCAR to calculate shift from
  character(len=300) :: CONTCARsBaseDir
    !! Base dir for sets of relaxed files if not captured
  character(len=300) :: dqFName
    !! File name for generalized-coordinate norms
  character(len=300) :: initPOSCARFName, finalPOSCARFName
    !! File name for POSCAR for relaxed initial and final charge states
  character(len=300) :: phononFName, phononPrimeFName
    !! File name for mesh.yaml phonon file
  character(len=300) :: prefix
    !! Prefix for shifted POSCARs

  logical :: calcDq
    !! If dq output file should be generated
  logical :: calcMaxDisp
    !! If calculating the maximum displacement between a pair of atoms
  logical :: calcSj
    !! If Sj should be calculated
  logical :: diffOmega
    !! If initial- and final-state frequencies 
    !! should be treated as different
  logical :: dqEigvecsFinal
    !! If final phonon eigenvectors should be used for
    !! Delta q, if applicable
  logical :: singleDisp
    !! If there is just a single displacement to consider
  logical :: generateShiftedPOSCARs
    !! If shifted POSCARs should be generated


  contains

!----------------------------------------------------------------------------
  subroutine readInputs(disp2AtomInd, freqThresh, shift, SjThresh, temperature, basePOSCARFName, CONTCARsBaseDir, dqFName, &
        energyTableDir, phononFName, phononPrimeFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, calcMaxDisp, &
        calcSj, diffOmega, dqEigvecsFinal, generateShiftedPOSCARs, singleDisp)

    implicit none

    ! Output variables:
    integer, intent(out) :: disp2AtomInd(2)
      !! Index of atoms to check displacement
      !! between if calcMaxDisp

    real(kind=dp), intent(out) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors
    real(kind=dp), intent(out) :: SjThresh
      !! Threshold to count modes above
    real(kind=dp), intent(out) :: temperature

    character(len=300), intent(out) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(out) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to get allowed transitions
    character(len=300), intent(out) :: phononFName, phononPrimeFName
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
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(out) :: dqEigvecsFinal
      !! If final phonon eigenvectors should be used for
      !! Delta q, if applicable
    logical, intent(out) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
    logical, intent(out) :: singleDisp
      !! If there is just a single displacement to consider


    namelist /inputParams/ initPOSCARFName, finalPOSCARFName, phononFName, prefix, shift, dqFName, generateShiftedPOSCARs, &
                           singleDisp, CONTCARsBaseDir, basePOSCARFName, freqThresh, calcSj, calcDq, calcMaxDisp, & 
                           disp2AtomInd, energyTableDir, diffOmega, phononPrimeFName, dqEigvecsFinal, SjThresh, &
                           temperature


    if(ionode) then

      call initialize(disp2AtomInd, freqThresh, shift, SjThresh, temperature, basePOSCARFName, CONTCARsBaseDir, &
          dqFName, energyTableDir, phononFName, phononPrimeFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, &
          calcMaxDisp, calcSj, diffOmega, dqEigvecsFinal, generateShiftedPOSCARs, singleDisp)
        ! Set default values for input variables and start timers
    
      read(5, inputParams, iostat=ierr)
        !! * Read input variables
    
      if(ierr /= 0) call exitError('readInputs', 'reading inputParams namelist', abs(ierr))
        !! * Exit calculation if there's an error

      call checkInitialization(disp2AtomInd, freqThresh, shift, SjThresh, temperature, basePOSCARFName, CONTCARsBaseDir, &
        dqFName, energyTableDir, phononFName, phononPrimeFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, &
        calcMaxDisp, calcSj, diffOmega, dqEigvecsFinal, generateShiftedPOSCARs, singleDisp)

    endif

    ! Send to other processes only what they need to know
    call MPI_BCAST(disp2AtomInd, 2, MPI_INTEGER, root, worldComm, ierr)

    call MPI_BCAST(freqThresh, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(shift, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(SjThresh, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(temperature, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call MPI_BCAST(basePOSCARFName, len(basePOSCARFName), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(CONTCARsBaseDir, len(CONTCARsBaseDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energyTableDir, len(energyTableDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(finalPOSCARFName, len(finalPOSCARFName), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(initPOSCARFName, len(initPOSCARFName), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(prefix, len(prefix), MPI_CHARACTER, root, worldComm, ierr)

    call MPI_BCAST(calcDq, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(calcSj, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(calcMaxDisp, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(diffOmega, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(dqEigvecsFinal, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(generateShiftedPOSCARs, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(singleDisp, 1, MPI_LOGICAL, root, worldComm, ierr)

    return

  end subroutine readInputs

!----------------------------------------------------------------------------
  subroutine initialize(disp2AtomInd, freqThresh, shift, SjThresh, temperature, basePOSCARFName, CONTCARsBaseDir, &
      dqFName, energyTableDir, phononFName, phononPrimeFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, &
      calcMaxDisp, calcSj, diffOmega, dqEigvecsFinal, generateShiftedPOSCARs, singleDisp)
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

    real(kind=dp), intent(out) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not
    real(kind=dp), intent(out) :: shift
      !! Magnitude of shift along phonon eigenvectors
    real(kind=dp), intent(out) :: SjThresh
      !! Threshold to count modes above
    real(kind=dp), intent(out) :: temperature

    character(len=300), intent(out) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(out) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(out) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(out) :: energyTableDir
      !! Path to energy table to get allowed transitions
    character(len=300), intent(out) :: phononFName, phononPrimeFName
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
    logical, intent(out) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(out) :: dqEigvecsFinal
      !! If final phonon eigenvectors should be used for
      !! Delta q, if applicable
    logical, intent(out) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
    logical, intent(out) :: singleDisp
      !! If there is just a single displacement to consider


    disp2AtomInd = -1

    CONTCARsBaseDir = ''
    dqFName = 'dq.txt'
    energyTableDir = ''
    basePOSCARFName = ''
    initPOSCARFName = 'POSCAR_init'
    finalPOSCARFName = 'POSCAR_final'
    phononFName = 'mesh.yaml'
    phononPrimeFName = ''
    prefix = 'ph_POSCAR'

    freqThresh = 0.5_dp
    shift = 0.01_dp
    SjThresh = 1e-1_dp
    temperature = 300_dp

    calcDq = .false.
    calcSj = .false.
    calcMaxDisp = .false.
    diffOmega = .false.
    dqEigvecsFinal = .true.
    generateShiftedPOSCARs = .false.
    singleDisp = .true.

    return

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(disp2AtomInd, freqThresh, shift, SjThresh, temperature, basePOSCARFName, CONTCARsBaseDir, &
      dqFName, energyTableDir, phononFName, phononPrimeFName, finalPOSCARFName, initPOSCARFName, prefix, calcDq, calcMaxDisp, &
      calcSj, diffOmega, dqEigvecsFinal, generateShiftedPOSCARs, singleDisp)

    implicit none

    ! Input variables:
    integer, intent(in) :: disp2AtomInd(2)
      !! Index of atoms to check displacement
      !! between if calcMaxDisp

    real(kind=dp), intent(in) :: freqThresh
      !! Threshold for frequency to determine if the mode 
      !! should be skipped or not
    real(kind=dp), intent(in) :: shift
      !! Magnitude of shift along phonon eigenvectors
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold to count modes above
    real(kind=dp), intent(in) :: temperature

    character(len=300), intent(inout) :: basePOSCARFName
      !! File name for intial POSCAR to calculate shift from
    character(len=300), intent(in) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(in) :: dqFName
      !! File name for generalized-coordinate norms
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table to get allowed transitions
    character(len=300), intent(in) :: phononFName, phononPrimeFName
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
    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(in) :: dqEigvecsFinal
      !! If final phonon eigenvectors should be used for
      !! Delta q, if applicable
    logical, intent(in) :: generateShiftedPOSCARs
      !! If shifted POSCARs should be generated
    logical, intent(in) :: singleDisp
      !! If there is just a single displacement to consider

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    ! Must-haves for reading phonons and calculating nj:
    abortExecution = checkFileInitialization('phononFName', phononFName)
    abortExecution = checkDoubleInitialization('freqThres', freqThresh, 0.0_dp, 1.0_dp) .or. abortExecution
    abortExecution = checkDoubleInitialization('temperature', temperature, 0.0_dp, 1500.0_dp) .or. abortExecution


    ! Need for calculating Sj:
    write(*,'("calcSj = ''",L1,"''")') calcSj
    if(calcSj) then

      abortExecution = checkDoubleInitialization('SjThresh', SjThresh, 0.0_dp, 1.0_dp) .or. abortExecution

      write(*,'("singleDisp = ''",L1,"''")') singleDisp
      if(singleDisp) then

        abortExecution = checkFileInitialization('initPOSCARFName', initPOSCARFName) .or. abortExecution
        abortExecution = checkFileInitialization('finalPOSCARFName', finalPOSCARFName) .or. abortExecution

      else

        write(*,'("CONTCARsBaseDir = ''",a,"''")') trim(CONTCARsBaseDir)
        abortExecution = checkFileInitialization('energyTableDir', trim(energyTableDir)//'/energyTable.1') .or. abortExecution

      endif


      write(*,'("diffOmega = ''",L1,"''")') diffOmega
      if(diffOmega) then
        abortExecution = checkFileInitialization('phononPrimeFName', phononPrimeFName)
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
      if(diffOmega) write(*,'("dqEigvecsFinal = ''",L1,"''")') dqEigvecsFinal

      abortExecution = checkDoubleInitialization('shift', shift, 0.0_dp, 10.0_dp) .or. abortExecution
      
      if(trim(basePOSCARFName) == '') then
        write(*,'("No input detected for basePOSCARFName! Defaulting to input value for initPOSCARFName.")')
        basePOSCARFName = initPOSCARFName
      endif

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
  subroutine readPhonons(freqThresh, phononFName, nAtoms, nModes, coordFromPhon, eigenvector, mass, omega)

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
      !! Coordinates from phonon file
    real(kind=dp), allocatable, intent(out) :: eigenvector(:,:,:)
      !! Eigenvectors for each atom for each mode
    real(kind=dp), allocatable, intent(out) :: mass(:)
      !! Masses of atoms
    real(kind=dp), allocatable, intent(out) :: omega(:)
      !! Frequency for each mode

    ! Local variables:
    integer :: j, ia, ix, jx
      !! Loop index
    integer :: nUsed
      !! Number of phonon modes used based on frequency threshold

    real(kind=dp), allocatable :: centerOfMassCoords(:,:)
      !! Coordinates of the center of mass
    real(kind=dp) :: omega_j
      !! Frequency for a single mode
    real(kind=dp) :: qPos(3)
      !! Phonon q position
    real(kind=dp) :: realLattVec(3,3)
      !! Real space lattice vectors

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

      ! Ignore line with "lattice:"
      read(57,*)

      do ix = 1,3
        read(57,'(A)') line 
        read(line(4:len(trim(line))),*) (realLattVec(jx,ix), jx=1,3)
      enddo

    endif


    call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
    allocate(coordFromPhon(3,nAtoms))
    allocate(centerOfMassCoords(3,nAtoms))
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

      call standardizeCoordinates(nAtoms, mass, realLattVec, coordFromPhon, centerOfMassCoords)

    endif

    call MPI_BCAST(coordFromPhon, size(coordFromPhon), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(mass, size(mass), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(centerOfMassCoords, size(centerOfMassCoords), MPI_DOUBLE_PRECISION, root, worldComm, ierr)


    nModes = 3*nAtoms - 3

    allocate(eigenvector(3,nAtoms,nModes))
    allocate(omega(nModes))


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
        read(line(15:len(trim(line))),*) omega_j
        if(abs(omega_j) > freqThresh) then
          nUsed = nUsed + 1
          omega(nUsed) = omega_j

          if(omega_j < 0) write(*,'("WARNING: Negative frequency detected in mode ", i7,&
            "! Re-relax the structure using the POSCAR shifted along this mode''s eigenvector.")') j
        else
          write(*,'("Skipping mode ", i7, " with freqency ", ES24.15E3,".")') j, omega_j
        endif


        read(57,'(A)') ! Ignore eigenvector section header

        do ia = 1, nAtoms

          read(57,'(A)') line ! Ignore atom number
  
          do ix = 1, 3

            read(57,'(A)') line
            if(abs(omega_j) > freqThresh)  read(line(10:len(trim(line))),*) eigenvector(ix,ia,nUsed)
              !! Only store the eigenvectors if the mode isn't skipped
              !! because those are the translational modes 

          enddo
        enddo
      enddo

      if(nUsed /= nModes) call exitError('readPhonons', 'Did not skip 3 modes based on threshold set. Check threshold.', 1)

      omega(:) = omega(:)*2*pi

      close(57)

    endif

    call MPI_BCAST(eigenvector, size(eigenvector), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omega, size(omega), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine readPhonons

!----------------------------------------------------------------------------
  subroutine lineUpModes(nAtoms, nModes, eigenvector, eigenvectorPrime, omegaPrime)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms in intial system
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each initial-state mode

    ! Output variables:
    real(kind=dp), intent(inout) :: eigenvectorPrime(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each final-state mode
      !! to be resorted
    real(kind=dp), intent(inout) :: omegaPrime(nModes)
      !! Frequency for each mode in the order originally passed
      !! in

    ! Local variables:
    integer :: jPrime
      !! Loop indices
    integer, allocatable :: maxMode(:)
      !! Mode index to match initial and final

    real(kind=dp), allocatable :: eigDotEigGlobal(:,:)
      !! Dot product of all eigenvectors
    real(kind=dp), allocatable :: eigenvectorPrime_in(:,:,:)
      !! Eigenvectors for each atom for each final-state mode
      !! in the order originally passed in
    real(kind=dp), allocatable :: omegaPrime_in(:)
      !! Frequency for each mode in the order originally passed
      !! in


    allocate(eigenvectorPrime_in(3,nAtoms,nModes), omegaPrime_in(nModes), maxMode(nModes))

    call getAllDotProds(nAtoms, nModes, eigenvector, eigenvectorPrime, eigDotEigGlobal)

    if(ionode) then
      call hungarianOptimalMatching(nModes, eigDotEigGlobal, .false., maxMode)

      ! Output the dot products to a file for visual confirmation
      open(17,file="optimalPairs.out")

      write(17,'("# Initial-mode index, Final-mode index, Dot product")')


      ! Store the original order
      eigenvectorPrime_in = eigenvectorPrime
      eigenvectorPrime = 0.0_dp
      omegaPrime_in = omegaPrime
      omegaPrime = 0.0_dp


      do jPrime = 1, nModes ! Loop over final
        write(17,'(2i10,ES13.4E3)') maxMode(jPrime), jPrime, eigDotEigGlobal(maxMode(jPrime),jPrime) 
        eigenvectorPrime(:,:,maxMode(jPrime)) = eigenvectorPrime_in(:,:,jPrime)
        omegaPrime(maxMode(jPrime)) = omegaPrime_in(jPrime)
      enddo

    endif

    call MPI_BARRIER(worldComm,ierr)

    call MPI_BCAST(eigenvectorPrime, size(eigenvectorPrime), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omegaPrime, size(omegaPrime), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine lineUpModes

!----------------------------------------------------------------------------
  subroutine getAllDotProds(nAtoms, nModes, eigenvector, eigenvectorPrime, eigDotEigGlobal)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms in intial system
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: eigenvector(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each initial-state mode
    real(kind=dp), intent(in) :: eigenvectorPrime(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each final-state mode

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: eigDotEigGlobal(:,:)
      !! Dot product of all eigenvectors

    ! Local variables:
    integer, allocatable :: displacement(:)
      !! Offset from beginning of array for
      !! scattering/gathering coefficients 
    integer :: iproc, jPrimeLocal, jPrimeGlobal, j, ia
      !! Loop index
    integer, allocatable :: sendCount(:)
      !! Number of items to send/recieve to/from each process

    real(kind=dp) :: eig(3), eigPrime(3)
      !! Local storage of eigenvectors
    real(kind=dp), allocatable :: eigDotEigLocal(:,:)
      !! Dot product of local eigenvectors

    
    allocate(eigDotEigLocal(nModes,nModesLocal))

    if(ionode) then
      allocate(eigDotEigGlobal(nModes,nModes))
    else
      allocate(eigDotEigGlobal(1,1))
    endif

    ! Create array with the number of dot products calculated by each process
    allocate(sendCount(nProcs))
    sendCount = 0
    sendCount(myid+1) = nModesLocal*nModes
    call mpiSumIntV(sendCount, worldComm)

    ! Create array with the starting point for each process
    allocate(displacement(nProcs))
    if(ionode) then
      displacement(1) = 0
      do iproc = 2, nProcs
        displacement(iproc) = displacement(iproc-1) + sendCount(iproc-1)
      enddo
    endif


    do jPrimeLocal = 1, nModesLocal ! Loop over final phonons
      jPrimeGlobal = jPrimeLocal+iModeStart-1
      do j = 1, nModes
        eigDotEigLocal(j,jPrimeLocal) = 0.0_dp
        do ia = 1, nAtoms

          eig = eigenvector(:,ia,j)
          eigPrime = eigenvectorPrime(:,ia,jPrimeGlobal)

          eigDotEigLocal(j,jPrimeLocal) = eigDotEigLocal(j,jPrimeLocal) + dot_product(eig,eigPrime)

        enddo
      enddo
    enddo

    call MPI_GATHERV(eigDotEigLocal, nModes*nModesLocal, MPI_DOUBLE_PRECISION, &
                     eigDotEigGlobal, sendCount, displacement, MPI_DOUBLE_PRECISION, 0, worldComm, ierr)

    deallocate(eigDotEigLocal)

    return

  end subroutine getAllDotProds

!----------------------------------------------------------------------------
  subroutine calcAndWriteNj(nModes, omega, temperature)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: omega(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: temperature

    ! Local variables:
    integer :: j
      !! Loop index

    real(kind=dp) :: beta
      !! 1/kb*T
    real(kind=dp) :: nj(nModes)
      !! \(n_j\) occupation number

    
    if(ionode) then

      beta = 1.0_dp/(kB_atomic*temperature)

      nj(:) = 1.0_dp/(exp(hbar_atomic*omega(:)*THzToHartree*beta) - 1.0_dp)

      open(17,file='njThermal.out')

      write(17,'("# Mode index, nj from thermal Bose-Einstein")')

      write(17,'("# Temperature (K): ")')
      write(17,'(f7.1)') temperature

      write(17,'("# Mode index, nj")')

      do j = 1, nModes
        write(17,'(1i7, 1ES24.15E3)') j, nj(j)
      enddo

      close(17)

    endif

    return

  end subroutine

!----------------------------------------------------------------------------
  subroutine calculateSj(nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, omega, omegaPrime, SjThresh, diffOmega, singleDisp, &
          CONTCARsBaseDir, energyTableDir, initPOSCARFName, finalPOSCARFName)

    implicit none

    ! Input variables:
    integer, intent(inout) :: nAtoms
      !! Number of atoms in intial system
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: coordFromPhon(3,nAtoms)
      !! Coordinates from phonon file
    real(kind=dp), intent(in) :: dqEigenvectors(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode to be
      !! used to calculate Delta q_j and displacements
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(in) :: SjThresh
      !! Threshold to count modes above

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
    logical, intent(in) :: singleDisp
      !! If there is just a single displacement to consider

    character(len=300), intent(in) :: CONTCARsBaseDir
      !! Base dir for sets of relaxed files if not captured
    character(len=300), intent(in) :: energyTableDir
      !! Path to energy table
    character(len=300), intent(inout) :: initPOSCARFName, finalPOSCARFName
      !! File name for POSCAR for relaxed initial and final charge states

    ! Local variables:
    integer :: iE
      !! Loop index
    integer, allocatable :: iki(:), ikf(:), ibi(:), ibf(:)
      !! Allowed state indices
    integer :: nTransitions
      !! Total number of transitions 
    integer :: modeIndex(nModes)
      !! Track mode indices after sorting

    real(kind=dp), allocatable :: dE(:,:)
      !! All energy differences from energy table
    real(kind=dp) :: Sj(nModes), SjPrime(nModes)
      !! Sj and Sj' Huang-Rhys factors using omega
      !! and omega'

    logical :: fileExists
      !! If a file exists

    character(len=300) :: SjFName
      !! File name for the Sj output


    if(ionode) then
      open(43,file='Sj.analysis.out')
      write(43,'("# SjThresh = ",ES11.3E2)') SjThresh
    endif

    if(singleDisp) then
  
      SjFName = 'Sj.out'
      call getAndWriteSingleSj(nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, omega, omegaPrime, diffOmega, &
                initPOSCARFName, finalPOSCARFName, SjFName, modeIndex, Sj, SjPrime)

      if(ionode) then
        if(.not. diffOmega) then
          write(43,'("# Max Sj (mode index, maxval), Num. Sj >= SjThresh")')
          write(43,'(i10,ES24.15E3,i10)') modeIndex(maxloc(Sj)), maxval(Sj), count(Sj >= SjThresh)
        else
          write(43,'("# Max Sj (mode index, maxval), Num. Sj >= SjThresh, Max Sj'' (maxval), Num. Sj'' >= SjThresh")')
          write(43,'(i10,ES24.15E3,i10,ES24.15E3,i10)') modeIndex(maxloc(Sj)), maxval(Sj), count(Sj >= SjThresh), &
                                                         maxval(SjPrime), count(SjPrime >= SjThresh) 
        endif
      endif

    else

      call readScatterEnergyTable(1, energyTableDir, ibi, ibf, iki, ikf, nTransitions, dE)
        ! Assume that the transitions allowed is the same for both
        ! spin channels

      if(ionode) then
        if(.not. diffOmega) then
          write(43,'("# iki, ibi, ikf, ibf, Max Sj (mode index, maxval), Num. Sj >= SjThresh")')
        else
          write(43,'("# iki, ibi, ikf, ibf, Max Sj (mode index, maxval), Num. Sj >= SjThresh, Max Sj'' (maxval), Num. Sj'' >= SjThresh")')
        endif
      endif

      do iE = 1, nTransitions

        ! Set initial and final file names based on state indices read
        ! from energy table. Make sure both files exist.
        if(ionode) then
          initPOSCARFName = trim(CONTCARsBaseDir)//'/k'//trim(int2str(iki(iE)))//'_b'//trim(int2str(ibi(iE)))//'/CONTCAR'
          inquire(file=trim(initPOSCARFName), exist=fileExists)
          if(.not. fileExists) call exitError('calculateSj', 'File does not exist!! '//trim(initPOSCARFName), 1)

          finalPOSCARFName = trim(CONTCARsBaseDir)//'/k'//trim(int2str(ikf(iE)))//'_b'//trim(int2str(ibf(iE)))//'/CONTCAR'
          inquire(file=trim(finalPOSCARFName), exist=fileExists)
          if(.not. fileExists) call exitError('calculateSj', 'File does not exist!! '//trim(finalPOSCARFName), 1)
        endif


        ! Set Sj output file name
        SjFName = 'Sj.k'//trim(int2str(iki(iE)))//'_b'//trim(int2str(ibi(iE)))//'.k'&
                        //trim(int2str(ikf(iE)))//'_b'//trim(int2str(ibf(iE)))//'.out'

        call getAndWriteSingleSj(nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, omega, omegaPrime, diffOmega, &
                  initPOSCARFName, finalPOSCARFName, SjFName, modeIndex, Sj, SjPrime)

        if(ionode) then
          if(.not. diffOmega) then
            write(43,'(5i10,ES24.15E3,i10)') iki(iE), ibi(iE), ikf(iE), ibf(iE), &
                                             modeIndex(maxloc(Sj)), maxval(Sj), count(Sj >= SjThresh)
          else
            write(43,'(5i10,ES24.15E3,i10,ES24.15E3,i10)') iki(iE), ibi(iE), ikf(iE), ibf(iE), &
                                                            modeIndex(maxloc(Sj)), maxval(Sj), count(Sj >= SjThresh), &
                                                            maxval(SjPrime), count(SjPrime >= SjThresh) 
          endif
        endif

      enddo

      deallocate(iki)
      deallocate(ikf)
      deallocate(ibi)
      deallocate(ibf)
      deallocate(dE)

    endif

    if(ionode) close(43)

    return

  end subroutine calculateSj

!----------------------------------------------------------------------------
  subroutine getAndWriteSingleSj(nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, omega, omegaPrime, diffOmega, &
            initPOSCARFName, finalPOSCARFName, SjFName, modeIndex, Sj, SjPrime)

    use generalComputations, only: direct2cart
    use cell, only: cartDispProjOnPhononEigsNorm

    implicit none

    ! Input variables:
    integer, intent(inout) :: nAtoms
      !! Number of atoms in intial system
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: coordFromPhon(3,nAtoms)
      !! Coordinates from phonon file
    real(kind=dp), intent(in) :: dqEigenvectors(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode to be
      !! used to calculate Delta q_j and displacements
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms
    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(in) :: initPOSCARFName, finalPOSCARFName
      !! File names for CONTCARs for different states
    character(len=300), intent(in) :: SjFName
      !! File name for the Sj output

    ! Output variables:
    integer, intent(out) :: modeIndex(nModes)
      !! Track mode indices after sorting

    real(kind=dp), intent(out) :: Sj(nModes), SjPrime(nModes)
      !! Sj and Sj' Huang-Rhys factors using omega
      !! and omega'

    ! Local variables:
    integer :: j
      !! Loop index

    real(kind=dp), allocatable :: atomPositionsDirInit(:,:), atomPositionsDirFinal(:,:)
      !! Atom positions in difference relaxed positions
    real(kind=dp) :: centerOfMassCoords(3,nAtoms)
      !! Coordinates of the center of mass
    real(kind=dp) :: displacement(3,nAtoms)
      !! Displacement 
    real(kind=dp) :: volumeInit, volumeFinal
      !! Volume of supercell
    real(kind=dp) :: projNorm(nModes)
      !! Generalized norms after displacement
    real(kind=dp) :: realLattVec(3,3)
      !! Real space lattice vectors

    
    call getCoordsAndDisp(nAtoms, coordFromPhon, mass, initPOSCARFName, atomPositionsDirInit, centerOfMassCoords, &
            displacement, realLattVec, volumeInit)
      ! Don't need displacement here. We only check compatibility with
      ! coordinates from phonons. 

    call getCoordsAndDisp(nAtoms, atomPositionsDirInit, mass, finalPOSCARFName, atomPositionsDirFinal, centerOfMassCoords, &
            displacement, realLattVec, volumeFinal)
      ! It is assumed for now that the lattice vectors are the same

    displacement = direct2cart(nAtoms, displacement, realLattVec)
    
    if(ionode) then
      if(abs(volumeInit - volumeFinal) > 1e-8) call exitError('calculateAndWriteSingleSj', 'volumes don''t match', 1)
    endif
  
    projNorm = 0.0_dp
    do j = iModeStart, iModeEnd

      projNorm(j) = cartDispProjOnPhononEigsNorm(nAtoms, displacement, dqEigenvectors(:,:,j), mass, realLattVec)

    enddo

    projNorm = projNorm*angToBohr*sqrt(daltonToElecM)
    call calcAndWriteSj(nModes, omega, omegaPrime, projNorm, diffOmega, SjFName, modeIndex, Sj, SjPrime)


    deallocate(atomPositionsDirInit)
    deallocate(atomPositionsDirFinal)

    return

  end subroutine getAndWriteSingleSj

!----------------------------------------------------------------------------
  subroutine getCoordsAndDisp(nAtoms, atomPositionsDirStart, mass, POSCARFName, atomPositionsDirEnd, centerOfMassCoords, &
          displacement, realLattVec, volume)

    use cell, only: readPOSCAR

    implicit none

    ! Input variables:
    integer, intent(inout) :: nAtoms
      !! Number of atoms in intial system

    real(kind=dp), intent(in) :: atomPositionsDirStart(3,nAtoms)
      !! Atom positions in direct coordinates to serve as the
      !! starting point for calculating the displacement (end-start)
    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms

    character(len=300), intent(in) :: POSCARFName
      !! File name for POSCAR

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: atomPositionsDirEnd(:,:)
      !! Atom positions in direct coordinates to serve as the
      !! ending point for calculating the displacement (end-start)
    real(kind=dp), intent(out) :: centerOfMassCoords(3,nAtoms)
      !! Coordinates of the center of mass
    real(kind=dp), intent(out) :: displacement(3,nAtoms)
      !! Displacement vector; not used here, only output 
      !! from check compatibility
    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: volume
      !! Volume of supercell

    ! Local variables:

    
    if(ionode) then
      call readPOSCAR(POSCARFName, nAtoms, atomPositionsDirEnd, realLattVec, volume)
      call standardizeCoordinates(nAtoms, mass, realLattVec, atomPositionsDirEnd, centerOfMassCoords)
    endif

    call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)

    if(.not. ionode) allocate(atomPositionsDirEnd(3,nAtoms))
    call MPI_BCAST(atomPositionsDirEnd, size(atomPositionsDirEnd), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(volume, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(centerOfMassCoords, size(centerOfMassCoords), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call getRelaxDispAndCheckCompatibility(nAtoms, atomPositionsDirEnd, atomPositionsDirStart, displacement)

    return

  end subroutine getCoordsAndDisp

!----------------------------------------------------------------------------
  subroutine standardizeCoordinates(nAtoms, mass, realLattVec, atomPositionsDir, centerOfMassCoords)

    use generalComputations, only: direct2cart, cart2direct

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms in system

    real(kind=dp), intent(in) :: mass(nAtoms)
      !! Mass of atoms
    real(kind=dp), intent(in) :: realLattVec(3,3)
      !! Real space lattice vectors

    ! Output variables:
    real(kind=dp), intent(inout) :: atomPositionsDir(3,nAtoms)
      !! Atom positions in direct coordinates
    real(kind=dp), intent(out) :: centerOfMassCoords(3,nAtoms)
      !! Coordinates of the center of mass

    ! Local variables: 
    real(kind=dp) :: atomPositionsCart(3,nAtoms)
      !! Atom positions in Cartesian coordinates


    ! Standardize to make all coordinates positive
    where(atomPositionsDir(:,:) < 0)
      atomPositionsDir(:,:) = 1.0_dp + atomPositionsDir(:,:)
    end where


    ! Convert to Cartesian coordinates to calculate center of mass
    atomPositionsCart = direct2cart(nAtoms, atomPositionsDir, realLattVec)

    ! Subtract center of mass
    centerOfMassCoords = centerOfMass(nAtoms, atomPositionsCart, mass)
    atomPositionsCart = atomPositionsCart - centerOfMassCoords

    ! Convert back to direct for easy standardization
    atomPositionsDir = cart2direct(nAtoms, atomPositionsCart, realLattVec)

    return

  end subroutine standardizeCoordinates

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
  subroutine getRelaxDispAndCheckCompatibility(nAtoms, atomPositionsDirEnd, atomPositionsDirStart, displacement)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: atomPositionsDirStart(3,nAtoms), atomPositionsDirEnd(3,nAtoms)
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
    real(kind=dp) :: posStart, posEnd
      !! Local storage of single coordinate

    logical :: abortExecution
      !! Whether or not to abort execution


    if(ionode) then

      abortExecution = .false.

      do ia= 1, nAtoms

        do ix = 1, 3

          posStart = atomPositionsDirStart(ix,ia)
          posEnd = atomPositionsDirEnd(ix,ia)

          dispInner = posEnd - posStart

          if(posEnd > posStart) then
            dispOuter = -posStart + posEnd - 1
          else
            dispOuter = 1 - posStart + posEnd
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
  subroutine calcAndWriteSj(nModes, omega, omegaPrime, projNorm, diffOmega, SjFName, modeIndex, Sj, SjPrime)

    use miscUtilities, only: hpsort_eps

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    real(kind=dp), intent(in) :: omega(nModes), omegaPrime(nModes)
      !! Frequency for each mode
    real(kind=dp), intent(inout) :: projNorm(nModes)
      !! Generalized norms after displacement

    logical, intent(in) :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different

    character(len=300), intent(in) :: SjFName
      !! File name for the Sj output

    ! Output variables:
    integer, intent(out) :: modeIndex(nModes)
      !! Track mode indices after sorting

    real(kind=dp), intent(out) :: Sj(nModes), SjPrime(nModes)
      !! Sj and Sj' Huang-Rhys factors using omega
      !! and omega'

    ! Local variables:
    integer :: j, jSort
      !! Loop index


    call mpiSumDoubleV(projNorm, worldComm)
      !! * Get the generalized-displacement norms
      !!   from all processes

    if(ionode) then

      Sj = 0.0_dp
      SjPrime = 0.0_dp

      do j = 1, nModes

        modeIndex(j) = j

        Sj(j) = projNorm(j)**2*omega(j)*THzToHartree/(2.0_dp*hbar_atomic)

        if(diffOmega) SjPrime(j) = projNorm(j)**2*omegaPrime(j)*THzToHartree/(2.0_dp*hbar_atomic)
          ! The way the algebra works, they are calculated using the same
          ! Delta q_j (projNorm)

      enddo

      call hpsort_eps(nModes, Sj, modeIndex, 1e-14_dp)
        ! Sort in ascending order

      open(60, file=trim(SjFName))

      write(60,'("# Number of modes")')

      write(60,'(1i7)') nModes

      write(60,'("# Tabulated for different initial and final frequencies?")')
      write(60,'(L5)') diffOmega


      ! Currently always sort by initial Sj
      if(diffOmega) then
        write(60,'("# Mode index, Sj (highest to lowest), omega_j, Sj'', omega_j''")')
      else
        write(60,'("# Mode index, Sj, (highest to lowest), omega_j")')
      endif


      do j = 1, nModes

        jSort = modeIndex(nModes-(j-1))

        if(diffOmega) then
          write(60,'(1i7, 4ES24.15E3)') jSort, Sj(nModes-(j-1)), omega(jSort), &
                                               SjPrime(jSort), omegaPrime(jSort)
        else
          write(60,'(1i7, 2ES24.15E3)') jSort, Sj(jSort), omega(jSort)
        endif

      enddo

      close(60)

    endif


    return

  end subroutine calcAndWriteSj

!----------------------------------------------------------------------------
  subroutine calculateShiftAndDq(disp2AtomInd, nAtoms, nModes, coordFromPhon, dqEigenvectors, mass, shift, calcDq, calcMaxDisp, &
        generateShiftedPOSCARs, basePOSCARFName, dqFName, prefix)

    use generalComputations, only: direct2cart
    use cell, only: cartDispProjOnPhononEigsNorm, readPOSCAR, writePOSCARNewPos

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
      !! Coordinates from phonon file
    real(kind=dp), intent(in) :: dqEigenvectors(3,nAtoms,nModes)
      !! Eigenvectors for each atom for each mode to be
      !! used to calculate Delta q_j and displacements
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
    real(kind=dp) :: centerOfMassCoords(3,nAtoms)
      !! Coordinates of the center of mass
    real(kind=dp) :: displacement(3,nAtoms)
      !! Atom displacements in angstrom
    real(kind=dp) :: volume
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


    call getCoordsAndDisp(nAtoms, coordFromPhon, mass, basePOSCARFName, atomPositionsDirBase, centerOfMassCoords, &
          displacement, realLattVec, volume)


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

      displacement = getShiftDisplacement(nAtoms, dqEigenvectors(:,:,j), realLattVec, mass, shift)

      if(calcMaxDisp) then
        relDisp = displacement(:,disp2AtomInd(1)) - displacement(:,disp2AtomInd(2))

        relDispMag(j) = sqrt(dot_product(relDisp,relDisp))
      endif

      if(calcDq) projNorm(j) = cartDispProjOnPhononEigsNorm(nAtoms, displacement, dqEigenvectors(:,:,j), mass, realLattVec)

      if(generateShiftedPOSCARs) then

        shiftedPositions = direct2cart(nAtoms, atomPositionsDirBase, realLattVec) + displacement + centerOfMassCoords

        shiftedPOSCARFName = trim(prefix)//"_"//trim(int2strLeadZero(j,suffixLength))

        call writePOSCARNewPos(nAtoms, shiftedPositions, basePOSCARFName, shiftedPOSCARFName, trim(memoLine)//' j = '//int2str(j), .true.)

      endif

    enddo


    deallocate(atomPositionsDirBase)
    deallocate(shiftedPositions)

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
  function getShiftDisplacement(nAtoms, dqEigenvectors, realLattVec, mass, shift) result(displacement)

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms

    real(kind=dp), intent(in) :: dqEigenvectors(3,nAtoms)
      !! Eigenvectors for each atom for each mode to be
      !! used to calculate Delta q_j and displacements
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

      eig = dqEigenvectors(:,ia)/sqrt(mass(ia))

      displacement(:,ia) = eig

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

!----------------------------------------------------------------------------
  subroutine readSjOneFreq(SjInput, nModes, omega, Sj)

    implicit none

    ! Input variables
    character(len=300), intent(in) :: SjInput
      !! Path to Sj.out file

    ! Output variables:
    integer, intent(out) :: nModes
      !! Number of phonon modes

    real(kind=dp),allocatable, intent(out) :: omega(:)
      !! Frequency for each mode
    real(kind=dp), allocatable, intent(out) :: Sj(:)
      !! Huang-Rhys factor for each mode

    ! Local variables:
    integer :: j, jSort
      !! Loop index

    real(kind=dp) :: Sj_, omega_
      !! Input variables

    logical :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
  
  
    if(ionode) then

      open(12,file=trim(SjInput))

      read(12,*)
      read(12,*) nModes

      read(12,*)
      read(12,*) diffOmega
      read(12,*)

      if(diffOmega) call exitError('readSjOneFreq', 'called one-frequency Sj, but this file tabulated for two frequencies!', 1)

    endif

    call MPI_BCAST(nModes, 1, MPI_INTEGER, root, worldComm, ierr)


    allocate(Sj(nModes))
    allocate(omega(nModes))

    
    if(ionode) then

      do j = 1, nModes
        read(12,*) jSort, Sj_, omega_! freq read from Sj.out is f(in Thz)*2pi
        Sj(jSort) = Sj_
        omega(jSort) = omega_
      end do

      omega(:) = omega(:)*THzToHartree
        ! Convert to Hz*2pi

      close(12)

    endif

    call MPI_BCAST(omega, size(omega), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(Sj, size(Sj), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return 

  end subroutine readSjOneFreq

!----------------------------------------------------------------------------
  subroutine readSjTwoFreq(SjInput, nModes, omega, omegaPrime, Sj, SjPrime)

    implicit none

    ! Input variables
    character(len=300), intent(in) :: SjInput
      !! Path to Sj.out file

    ! Output variables:
    integer, intent(out) :: nModes
      !! Number of phonon modes

    real(kind=dp),allocatable, intent(out) :: omega(:), omegaPrime(:)
      !! Frequency for each mode
    real(kind=dp), allocatable, intent(out) :: Sj(:), SjPrime(:)
      !! Huang-Rhys factor for each mode

    ! Local variables:
    integer :: j, jSort
      !! Loop index

    real(kind=dp) :: rDum
      !! Dummy variable to ignore input
    real(kind=dp) :: Sj_, SjPrime_, omega_, omegaPrime_
      !! Input variables

    logical :: diffOmega
      !! If initial- and final-state frequencies 
      !! should be treated as different
  
  
    if(ionode) then

      open(12,file=trim(SjInput))

      read(12,*)
      read(12,*) nModes

      read(12,*)
      read(12,*) diffOmega
      read(12,*)

      if(.not. diffOmega) call exitError('readSjTwoFreq', 'called two-frequency Sj, but this file tabulated for a single frequency!', 1)

    endif

    call MPI_BCAST(nModes, 1, MPI_INTEGER, root, worldComm, ierr)


    allocate(Sj(nModes), SjPrime(nModes))
    allocate(omega(nModes), omegaPrime(nModes))

    
    if(ionode) then

      do j = 1, nModes
        read(12,*) jSort, Sj_, omega_, SjPrime_, omegaPrime_! freq read from Sj.out is f(in Thz)*2pi
        Sj(jSort) = Sj_
        omega(jSort) = omega_
        SjPrime(jSort) = SjPrime_
        omegaPrime(jSort) = omegaPrime_
      end do

      omega(:) = omega(:)*THzToHartree
      omegaPrime(:) = omegaPrime(:)*THzToHartree
        ! Convert to Hz*2pi

      close(12)

    endif

    call MPI_BCAST(omega, size(omega), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(Sj, size(Sj), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omegaPrime, size(omegaPrime), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(SjPrime, size(SjPrime), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return 

  end subroutine readSjTwoFreq

!----------------------------------------------------------------------------
  subroutine readNj(nModes, njInput, nj)

    implicit none

    ! Input variables:
    integer, intent(in) :: nModes
      !! Number of modes

    character(len=300), intent(in) :: njInput
      !! Path to nj file

    ! Output variables:
    real(kind=dp), intent(out) :: nj(nModes)
      !! \(n_j\) occupation number

    ! Local variables:
    integer :: iDum
      !! Dummy integer to ignore input
    integer :: j
      !! Loop index

    
    if(ionode) then

      open(17,file=trim(njInput))

      read(17,*)
      read(17,*)
      read(17,*)
      read(17,*)

      do j = 1, nModes
        read(17,'(1i7, 1ES24.15E3)') iDum, nj(j)
      enddo

      close(17)

    endif

    call MPI_BCAST(nj, size(nj), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    return

  end subroutine

end module PhononPPMod
