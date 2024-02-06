module wfcExportVASPMod
  
  use constants, only: dp, angToBohr, eVToRy, RyToHartree, eVToHartree, pi, twopi
  use base, only: nKPoints, nSpins, ispSelect, loopSpins
  use errorsAndMPI
  use cell
  use mpi

  implicit none

  ! Parameters:
  integer, parameter :: mainOutFileUnit = 50
    !! Main output file unit
  integer, parameter :: potcarUnit = 71
    !! POTCAR unit for I/O
  integer, parameter :: nonlPseudoGridSize = 100
    !! Size of non-local pseudopotential grid
  integer, parameter :: wavecarUnit = 72
    !! WAVECAR unit for I/O
  integer, parameter :: maxNumDispKPerCoord = 6
    !! Max number of displaced k-points per coordinate
    !! for group velocity calculations

  real(kind = dp), parameter :: twoPiSquared = (2.0_dp*pi)**2
    !! This is used in place of \(2\pi/a\) which assumes that \(a=1\)


  ! Global variables not passed as arguments:
  integer, allocatable :: iGkEnd_pool(:)
    ! Ending index for G+k vectors on
    ! single process in a given pool
  integer, allocatable :: iGkStart_pool(:)
    ! Starting index for G+k vectors on
    ! single process in a given pool

  real(kind=dp) :: t1, t2, t0
    !! Timers


  ! Variables that should be passed as arguments:
  real(kind=dp) :: eFermi
    !! Fermi energy
  real(kind=dp) :: eTot
    !! Total energy
  real(kind=dp), allocatable :: gVecInCart(:,:)
    !! G-vectors in Cartesian coordinates
  real(kind=dp), allocatable :: bandOccupation(:,:,:)
    !! Occupation of band
  real(kind=dp) :: wfcVecCut
    !! Energy cutoff converted to vector cutoff
  real(kind=dp), allocatable :: kPosition(:,:)
    !! Position of k-points in reciprocal space
  real(kind=dp), allocatable :: kWeight(:)
    !! Weight of k-points
  real(kind=dp), allocatable :: patternArr(:)
    !! Displacement pattern for groups of
    !! k-points for group velocity calculations

  complex*16, allocatable :: eigenE(:,:,:)
    !! Band eigenvalues
  
  integer, allocatable :: gIndexLocalToGlobal(:)
    !! Converts local index `ig` to global index
  integer, allocatable :: gKIndexLocalToGlobal(:,:)
    !! Local to global indices for \(G+k\) vectors 
    !! ordered by magnitude at a given k-point
  integer, allocatable :: gKSort(:,:)
    !! Indices to recover sorted order on reduced
    !! \(G+k\) grid
  integer, allocatable :: gVecMillerIndicesGlobalOrig(:,:)
    !! Integer coefficients for G-vectors on all processors (original order)
  integer, allocatable :: gVecMillerIndicesGlobalSort(:,:)
    !! Integer coefficients for G-vectors on all processors (sorted)
  integer, allocatable :: igkSort2OrigLocal(:,:)
    !! Indices of \(G+k\) vectors in just this pool
    !! and for local PWs in the original order
  integer, allocatable :: iMill(:)
    !! Indices of miller indices after sorting
  integer, allocatable :: iType(:)
    !! Atom type index
  integer :: fftGridSize(3)
    !! Max number of points on the FFT grid in each direction
  integer :: nBands
    !! Total number of bands
  integer :: nDispkPerCoord
    !! Number of displaced k-points per coordinate
  integer, allocatable :: nGkLessECutGlobal(:)
    !! Global number of \(G+k\) vectors with magnitude
    !! less than `wfcVecCut` for each k-point
  integer, allocatable :: nGkVecsLocal(:)
    !! Local number of G+k-vectors on this processor
  integer :: nGVecsGlobal
    !! Global number of G-vectors
  integer :: nGVecsLocal
    !! Local number of G-vectors on this processor
  integer, allocatable :: nPWs1kGlobal(:)
    !! Input number of plane waves for a single k-point for all processors
  integer :: maxGIndexGlobal
    !! Maximum G-vector index among all \(G+k\)
    !! and processors
  integer :: maxGkVecsLocal
    !! Max number of G+k vectors across all k-points
    !! in this pool
  integer :: reclenWav
    !! Number of records in WAVECAR file

  logical :: energiesOnly
    !! If only energy-related files should be exported
  logical :: gammaOnly
    !! If the gamma only VASP code is used
  logical :: groupForGroupVelocity
    !! If there are groups of k-points for group
    !! velocity calculation
  
  character(len=256) :: exportDir
    !! Directory to be used for export
  character(len=256) :: mainOutputFile
    !! Main output file
  character(len=300) :: pattern
    !! Character input for displacement pattern
  character(len=256) :: VASPDir
    !! Directory with VASP files

  type potcar
    integer :: angMom(16) = 0
      !! Angular momentum of projectors
    integer :: iRAugMax
      !! Max index of augmentation sphere
    integer :: lmmax
      !! Total number of nlm channels
    integer :: nChannels
      !! Number of l channels;
      !! also number of projectors
    integer :: nmax
      !! Number of radial grid points

    real(kind=dp), allocatable :: dRadGrid(:)
      !! Derivative of radial grid
    real(kind=dp) :: maxGkNonlPs
      !! Max \(|G+k|\) for non-local potential
    real(kind=dp) :: psRMax
      !! Max r for non-local contribution
    real(kind=dp), allocatable :: radGrid(:)
      !! Radial grid points
    real(kind=dp) :: rAugMax
      !! Maximum radius of augmentation sphere
    real(kind=dp) :: recipProj(16,0:nonlPseudoGridSize)
      !! Reciprocal-space projectors
    real(kind=dp), allocatable :: wae(:,:)
      !! AE wavefunction
    real(kind=dp), allocatable :: wps(:,:)
      !! PS wavefunction

    character(len=2) :: element
  end type potcar

  type (potcar), allocatable :: pot(:)

  namelist /inputParams/ VASPDir, exportDir, gammaOnly, energiesOnly, groupForGroupVelocity, nDispkPerCoord, pattern, ispSelect


  contains

!----------------------------------------------------------------------------
  subroutine readInputParams(ispSelect, nDispkPerCoord, patternArr, energiesOnly, gammaOnly, groupForGroupVelocity, loopSpins, &
        exportDir, pattern, VASPDir)

    implicit none

    ! Output variables:
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: nDispkPerCoord
      !! Number of displaced k-points per coordinate
      
    real(kind=dp), allocatable, intent(out) :: patternArr(:)
      !! Displacement pattern for groups of
      !! k-points for group velocity calculations

    logical, intent(out) :: energiesOnly
      !! If only energy-related files should be exported
    logical, intent(out) :: gammaOnly
      !! If the gamma only VASP code is used
    logical, intent(out) :: groupForGroupVelocity
      !! If there are groups of k-points for group
      !! velocity calculation
    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    character(len=256), intent(out) :: exportDir
      !! Directory to be used for export
    character(len=300), intent(out) :: pattern
      ! Character input for displacement pattern
    character(len=256), intent(out) :: VASPDir
      !! Directory with VASP files

    ! Local variables:
    integer :: idk
      !! Loop index


    if(ionode) then
    
      call initialize(ispSelect, nDispkPerCoord, energiesOnly, gammaOnly, groupForGroupVelocity, exportDir, pattern, VASPDir)

      read(5, inputParams, iostat=ierr)
    
      if(ierr /= 0) call exitError('readInputParams', 'reading inputParams namelist', abs(ierr))

      call checkInitialization(ispSelect, nDispkPerCoord, energiesOnly, gammaOnly, groupForGroupVelocity, exportDir, &
            pattern, VASPDir, loopSpins)

    endif

    call MPI_BCAST(exportDir, len(exportDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(VASPDir, len(VASPDir), MPI_CHARACTER, root, worldComm, ierr)
    call MPI_BCAST(energiesOnly, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(gammaOnly, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(loopSpins, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(groupForGroupVelocity, 1, MPI_LOGICAL, root, worldComm, ierr)
    call MPI_BCAST(ispSelect, 1, MPI_INTEGER, root, worldComm, ierr)


    if(groupForGroupVelocity) then

      call MPI_BCAST(nDispkPerCoord, 1, MPI_INTEGER, root, worldComm, ierr)

      allocate(patternArr(nDispkPerCoord))

      if(ionode) read(pattern,*) (patternArr(idk), idk=1,nDispkPerCoord)

      call MPI_BCAST(patternArr, size(patternArr), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    endif

    return

  end subroutine readInputParams

!----------------------------------------------------------------------------
  subroutine initialize(ispSelect, nDispkPerCoord, energiesOnly, gammaOnly, groupForGroupVelocity, exportDir, pattern, VASPDir)
    
    implicit none

    ! Output variables:
    integer, intent(out) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(out) :: nDispkPerCoord
      !! Number of displaced k-points per coordinate

    logical, intent(out) :: energiesOnly
      !! If only energy-related files should be exported
    logical, intent(out) :: gammaOnly
      !! If the gamma only VASP code is used
    logical, intent(out) :: groupForGroupVelocity
      !! If there are groups of k-points for group
      !! velocity calculation

    character(len=256), intent(out) :: exportDir
      !! Directory to be used for export
    character(len=300), intent(out) :: pattern
      !! Character input for displacement pattern
    character(len=256), intent(out) :: VASPDir
      !! Directory with VASP files


    nDispkPerCoord = -1
    ispSelect = -1

    energiesOnly = .false.
    gammaOnly = .false.
    groupForGroupVelocity = .false.

    exportDir = './export'
    pattern = ''
    VASPDir = './'

    return 

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine checkInitialization(ispSelect, nDispkPerCoord, energiesOnly, gammaOnly, groupForGroupVelocity, exportDir, pattern, &
        VASPDir, loopSpins)
    
    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nDispkPerCoord
      !! Number of displaced k-points per coordinate
      
    logical, intent(in) :: energiesOnly
      !! If only energy-related files should be exported
    logical, intent(in) :: gammaOnly
      !! If the gamma only VASP code is used
    logical, intent(in) :: groupForGroupVelocity
      !! If there are groups of k-points for group
      !! velocity calculation

    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export
    character(len=300), intent(in) :: pattern
      !! Character input for displacement pattern
    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    ! Output variables:
    logical, intent(out) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    ! Local variables:
    logical :: abortExecution
      !! Whether or not to abort the execution


    abortExecution = .false.

    if(gammaOnly .and. .not. energiesOnly) then
      write(*,'("ERROR: gamma-only version only currently implemented for exporting energies only")')
      abortExecution = .true.
    endif

    if(groupForGroupVelocity) then
      abortExecution = checkIntInitialization('nDispkPerCoord', nDispkPerCoord, 1, maxNumDispkPerCoord) .or. abortExecution
      abortExecution = checkStringInitialization('pattern', pattern) .or. abortExecution
    endif

    if(ispSelect < 1 .or. ispSelect > 2) then
      write(*,*) "No valid choice for spin channel selection given. Looping over spin."
      loopSpins = .true.
    else
      write(*,'("Only exporting spin channel ", i2)') ispSelect
      loopSpins = .false.
    endif

    abortExecution = checkDirInitialization('VASPDir', VASPDir, 'OUTCAR') .or. abortExecution
  
    call execute_command_line('mkdir -p '//trim(exportDir))
    
    mainOutputFile = trim(exportDir)//"/input"
    
    write(*,*) "Opening file "//trim(mainOutputFile)
    open(mainOutFileUnit, file=trim(mainOutputFile))

    if(abortExecution) then
      write(*, '(" Program stops!")')
      stop
    endif

    return 

  end subroutine checkInitialization

!----------------------------------------------------------------------------
  subroutine readWAVECAR(ispSelect, loopSpins, VASPDir, realLattVec, recipLattVec, bandOccupation, omega, wfcVecCut, &
        kPosition, fftGridSize, nBands, nKPoints, nPWs1kGlobal, nSpins, reclenWav, eigenE)
    !! Read cell and wavefunction data from the WAVECAR file
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user

    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files
    
    ! Output variables:
    real(kind=dp), intent(out) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), allocatable, intent(out) :: bandOccupation(:,:,:)
      !! Occupation of band
    real(kind=dp), intent(out) :: omega
      !! Volume of unit cell
    real(kind=dp), intent(out) :: wfcVecCut
      !! Energy cutoff converted to vector cutoff
    real(kind=dp), allocatable, intent(out) :: kPosition(:,:)
      !! Position of k-points in reciprocal space

    integer, intent(out) :: fftGridSize(3)
      !! Max number of points on the FFT grid in each direction
    integer, intent(out) :: nBands
      !! Total number of bands
    integer, intent(out) :: nKPoints
      !! Total number of k-points
    integer, allocatable, intent(out) :: nPWs1kGlobal(:)
      !! Input number of plane waves for a single k-point 
      !! for all processors
    integer, intent(out) :: nSpins
      !! Number of spins
    integer, intent(out) :: reclenWav
      !! Number of records in WAVECAR file

    complex*16, allocatable, intent(out) :: eigenE(:,:,:)
      !! Band eigenvalues


    ! Local variables:
    real(kind=dp) :: c = 0.26246582250210965422
      !! \(2m/\hbar^2\) converted from J\(^{-1}\)m\(^{-2}\)
      !! to eV\(^{-1}\)A\(^{-2}\)
    real(kind=dp) :: reclenWav_real, nspin_real, prec_real, nkstot_real 
      !! Real version of integers for reading from file
    real(kind=dp) :: nbnd_real
      !! Real version of integers for reading from file
    real(kind=dp) :: wfcECut
      !! Plane wave energy cutoff in Ry

    integer :: j
      !! Index used for reading lattice vectors
    integer :: prec
      !! Precision of plane wave coefficients

    character(len=256) :: fileName
      !! Full WAVECAR file name including path

    
    if(ionode) then

      fileName = trim(VASPDir)//'/WAVECAR'

      reclenWav = 24
        ! Set a starting value for the number of records

      open(unit=wavecarUnit, file=fileName, access='direct', recl=reclenWav, iostat=ierr, status='old')
      if (ierr .ne. 0) write(*,*) 'open error - iostat =', ierr
        !! * If root node, open the `WAVECAR` file

      read(unit=wavecarUnit,rec=1) reclenWav_real, nspin_real, prec_real
        !! @note Must read in as real first then convert to integer @endnote

      close(unit=wavecarUnit)

      reclenWav = nint(reclenWav_real)
      nSpins = nint(nspin_real)
      prec = nint(prec_real)
        ! Convert input variables to integers


      ! If not looping spins (i.e., a single spin channel was selected),
      ! make sure that the selected spin channel is not larger than
      ! the number of spin channels available in the system
      if(.not. loopSpins) then

        if(checkIntInitialization('ispSelect', ispSelect, 1, nSpins)) &
          call exitError('readWAVECAR', 'selected spin larger than number of spin channels available', 1)
          ! Put this inside the other if statement and not as an `.and.`
          ! condition because fortran does not short-circuit within
          ! conditionals, and the call to `checkIntInitialization` will
          ! output an error message even if `loopSpins = .true.`
      endif


      if(prec .eq. 45210) call exitError('readWAVECAR', 'WAVECAR_double requires complex*16', 1)

      open(unit=wavecarUnit, file=fileName, access='direct', recl=reclenWav, iostat=ierr, status='old')
      if (ierr .ne. 0) write(*,*) 'open error - iostat =', ierr
        !! * Reopen WAVECAR with correct number of records

      read(unit=wavecarUnit,rec=2) nkstot_real, nbnd_real, wfcECut,(realLattVec(j,1),j=1,3),&
          (realLattVec(j,2),j=1,3), (realLattVec(j,3),j=1,3)
        !! * Read total number of k-points, plane wave cutoff energy, and real
        !!   space lattice vectors
      !read(unit=wavecarUnit,rec=2) nkstot_real, nbnd_real, wfcECut,((realLattVec(i,j),j=1,3),i=1,3)
        !! @todo Test this more compact form @endtodo

      close(wavecarUnit)

      wfcVecCut = sqrt(wfcECut*c)/angToBohr
        !! * Calculate vector cutoff from energy cutoff

      realLattVec = realLattVec*angToBohr

      nKPoints = nint(nkstot_real)
      nBands = nint(nbnd_real)
        ! Convert input variables to integers

      call calculateOmega(realLattVec, omega)
        !! * Calculate the cell volume as \(a_1\cdot a_2\times a_3\)

      call getReciprocalVectors(realLattVec, omega, recipLattVec)
        !! * Calculate the reciprocal lattice vectors from the real-space
        !!   lattice vectors and the cell volume

      call estimateMaxNumPlanewaves(recipLattVec, wfcECut*eVToRy, fftGridSize)
        !! * Get the maximum number of plane waves

      !> * Write out total number of k-points, number of bands, 
      !>   the energy cutoff, the real-space-lattice vectors,
      !>   the cell volume, and the reciprocal lattice vectors
      write(*,*) 'no. k points =', nKPoints
      write(*,*) 'no. bands =', nBands
      write(*,*) 'max. energy (eV) =', sngl(wfcECut)
        !! @note 
        !!  The energy cutoff is currently output to stdout
        !!  in eV to compare with output from WaveTrans.
        !! @endnote
      write(*,*) 'real space lattice vectors:'
      write(*,*) 'a1 =', (sngl(realLattVec(j,1)),j=1,3)
      write(*,*) 'a2 =', (sngl(realLattVec(j,2)),j=1,3)
      write(*,*) 'a3 =', (sngl(realLattVec(j,3)),j=1,3)
      write(*,*) 
      write(*,*) 'volume unit cell =', sngl(omega)
      write(*,*) 
      write(*,*) 'reciprocal lattice vectors:'
      write(*,*) 'b1 =', (sngl(recipLattVec(j,1)),j=1,3)
      write(*,*) 'b2 =', (sngl(recipLattVec(j,2)),j=1,3)
      write(*,*) 'b3 =', (sngl(recipLattVec(j,3)),j=1,3)
      write(*,*) 
        !! @note
        !!  I made an intentional choice to stick with the unscaled lattice
        !!  vectors until I see if it will be convenient to scale them down.
        !!  QE uses the `alat` and `tpiba` scaling quite a bit though, so I
        !!  will have to be careful with the scaling/units.
        !! @endnote

      write(mainOutFileUnit, '("# Cell volume (a.u.)^3. Format: ''(ES24.15E3)''")')
      write(mainOutFileUnit, '(ES24.15E3)' ) omega
      flush(mainOutFileUnit)

    endif

    call MPI_BCAST(reclenWav, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nSpins, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nKPoints, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nBands, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(wfcVecCut, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(omega, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(fftGridSize, 3, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(realLattVec, size(realLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(recipLattVec, size(recipLattVec), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    call preliminaryWAVECARScan(ispSelect, nBands, nKPoints, nSpins, reclenWav, loopSpins, bandOccupation, kPosition, nPWs1kGlobal, eigenE)
      !! * For each spin and k-point, read the number of
      !!   \(G+k\) vectors below the energy cutoff, the
      !!   position of the k-point in reciprocal space, 
      !!   and the eigenvalue and occupation for each band

    return
  end subroutine readWAVECAR

!----------------------------------------------------------------------------
  subroutine estimateMaxNumPlanewaves(recipLattVec, wfcECut, fftGridSize)
    !! Get the maximum number of plane waves. I'm not sure how 
    !! this is done completely. It seems to be just basic vector
    !! stuff, but I haven't been able to make sense of it.
    !! 
    !! @todo Figure out how `estimateMaxNumPlanewaves` works @endtodo
    !!
    !! <h2>Walkthrough</h2>
    !!

    use generalComputations, only: vcross

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: wfcECut
      !! Plane wave energy cutoff in Ry

    ! Output variables:
    integer, intent(out) :: fftGridSize(3)
      !! Max number of points on the FFT grid in each direction

    ! Local variables:
    real(kind=dp) :: b1mag, b2mag, b3mag
      !! Reciprocal vector magnitudes
    real(kind=dp) :: c = 0.26246582250210965422
      !! \(2m/\hbar^2\) converted from J\(^{-1}\)m\(^{-2}\)
      !! to eV\(^{-1}\)A\(^{-2}\)
    real(kind=dp) :: phi12, phi13, phi23
      !! Angle between vectors
    real(kind=dp) :: sinphi123
      !! \(\sin\phi_{123}\)
    real(kind=dp) :: vmag
      !! Magnitude of temporary vector
    real(kind=dp) :: vtmp(3)
      !! Temporary vector for calculating angles

    integer :: nb1max, nb2max, nb3max
      !! Max magnitude of G-vectors in each direction
    integer :: nb1maxA, nb2maxA, nb3maxA
      !! First estimate of max
    integer :: nb1maxB, nb2maxB, nb3maxB
      !! Second estimate of max
    integer :: nb1maxC, nb2maxC, nb3maxC
      !! Third estimate of max


    b1mag = sqrt(sum(recipLattVec(:,1)**2))
    b2mag = sqrt(sum(recipLattVec(:,2)**2))
    b3mag = sqrt(sum(recipLattVec(:,3)**2))

    write(*,*) 'reciprocal lattice vector magnitudes:'
    write(*,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)
      !! * Calculate and output reciprocal vector magnitudes


    phi12 = acos(dot_product(recipLattVec(:,1),recipLattVec(:,2))/(b1mag*b2mag))
      !! * Calculate angle between \(b_1\) and \(b_2\)

    call vcross(recipLattVec(:,1), recipLattVec(:,2), vtmp)
    vmag = sqrt(dot_product(vtmp,vtmp))
    sinphi123 = dot_product(recipLattVec(:,3),vtmp(:))/(vmag*b3mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxA = (dsqrt(wfcECut/eVToRy*c)/(b1mag*abs(sin(phi12)))) + 1
    nb2maxA = (dsqrt(wfcECut/eVToRy*c)/(b2mag*abs(sin(phi12)))) + 1
    nb3maxA = (dsqrt(wfcECut/eVToRy*c)/(b3mag*abs(sinphi123))) + 1
      !! * Get first set of max values


    phi13 = acos(dot_product(recipLattVec(:,1),recipLattVec(:,3))/(b1mag*b3mag))
      !! * Calculate angle between \(b_1\) and \(b_3\)

    call vcross(recipLattVec(:,1), recipLattVec(:,3), vtmp)
    vmag = sqrt(dot_product(vtmp,vtmp))
    sinphi123 = dot_product(recipLattVec(:,2),vtmp(:))/(vmag*b2mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxB = (dsqrt(wfcECut/eVToRy*c)/(b1mag*abs(sin(phi13)))) + 1
    nb2maxB = (dsqrt(wfcECut/eVToRy*c)/(b2mag*abs(sinphi123))) + 1
    nb3maxB = (dsqrt(wfcECut/eVToRy*c)/(b3mag*abs(sin(phi13)))) + 1
      !! * Get first set of max values


    phi23 = acos(dot_product(recipLattVec(:,2),recipLattVec(:,3))/(b2mag*b3mag))
      !! * Calculate angle between \(b_2\) and \(b_3\)

    call vcross(recipLattVec(:,2), recipLattVec(:,3), vtmp)
    vmag = sqrt(dot_product(vtmp,vtmp))
    sinphi123 = dot_product(recipLattVec(:,1),vtmp(:))/(vmag*b1mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxC = (dsqrt(wfcECut/eVToRy*c)/(b1mag*abs(sinphi123))) + 1
    nb2maxC = (dsqrt(wfcECut/eVToRy*c)/(b2mag*abs(sin(phi23)))) + 1
    nb3maxC = (dsqrt(wfcECut/eVToRy*c)/(b3mag*abs(sin(phi23)))) + 1
      !! * Get first set of max values


    nb1max = max0(nb1maxA,nb1maxB,nb1maxC)
    nb2max = max0(nb2maxA,nb2maxB,nb2maxC)
    nb3max = max0(nb3maxA,nb3maxB,nb3maxC)

    fftGridSize(1) = 2*nb1max + 1
    fftGridSize(2) = 2*nb2max + 1
    fftGridSize(3) = 2*nb3max + 1

    write(*,*) 'max. no. G values; 1,2,3 =', nb1max, nb2max, nb3max
    write(*,*) ' '

    return
  end subroutine estimateMaxNumPlanewaves

!----------------------------------------------------------------------------
  subroutine preliminaryWAVECARScan(ispSelect, nBands, nKPoints, nSpins, reclenWav, loopSpins, bandOccupation, kPosition, &
          nPWs1kGlobal, eigenE)
    !! For each spin and k-point, read the number of
    !! \(G+k\) vectors below the energy cutoff, the
    !! position of the k-point in reciprocal space, 
    !! and the eigenvalue and occupation for each band
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nSpins
      !! Number of spins
    integer, intent(in) :: reclenWav
      !! Number of records in the WAVECAR file

    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: bandOccupation(:,:,:)
      !! Occupation of band
    real(kind=dp), allocatable, intent(out) :: kPosition(:,:)
      !! Position of k-points in reciprocal space

    integer, allocatable, intent(out) :: nPWs1kGlobal(:)
      !! Input number of plane waves for a single k-point 
      !! for all processors

    complex*16, allocatable, intent(out) :: eigenE(:,:,:)
      !! Band eigenvalues


    ! Local variables:
    real(kind=dp) :: nPWs1kGlobal_real
      !! Real version of integers for reading from file

    integer :: irec, isp, ik, i, iband
      !! Loop indices

    character(len=256) :: fileName
      !! Full WAVECAR file name including path


    allocate(bandOccupation(nSpins, nBands, nKPoints))
    allocate(kPosition(3,nKPoints))
    allocate(nPWs1kGlobal(nKPoints))
    allocate(eigenE(nSpins,nKPoints,nBands))

    bandOccupation = -1.0_dp
    eigenE = 0.0_dp
    
    fileName = trim(VASPDir)//'/WAVECAR'

    if(ionode) then
      open(unit=wavecarUnit, file=fileName, access='direct', recl=reclenWav, iostat=ierr, status='old')
      if (ierr .ne. 0) write(*,*) 'open error - iostat =', ierr

      irec=2

      do isp = 1, nSpins
        !! * For each spin:
        !!    * Go through each k-point
        !!       * Read in the number of \(G+k\) plane wave
        !!         vectors below the energy cutoff
        !!       * Read the position of the k-point in 
        !!         reciprocal space
        !!       * Read in the eigenvalue and occupation for
        !!         each band

        do ik = 1, nKPoints
        
          irec = irec + 1


          if(loopSpins .or. isp == ispSelect) then
            ! After adding this, I realized that the spin-independent
            ! values like number of plane waves and k position are
            ! being overwritten. This seems to be working okay, so
            ! VASP must just repeat the same information for each spin 
            ! channel


            ! Read in the number of \(G+k\) plane wave vectors below the energy
            ! cutoff, the position of the k-point in reciprocal space, and
            ! the eigenvalue and occupation for each band.
            read(unit=wavecarUnit,rec=irec) nPWs1kGlobal_real, (kPosition(i,ik),i=1,3), &
                     (eigenE(isp,ik,iband), bandOccupation(isp, iband, ik), iband=1,nBands)


            nPWs1kGlobal(ik) = nint(nPWs1kGlobal_real)
              !! @note
              !!  `nPWs1kGlobal(ik)` corresponds to `WDES%NPLWKP_TOT(K)` in VASP (see 
              !!  subroutine `OUTWAV` in `fileio.F`). In the `GEN_INDEX` subroutine in
              !!  `wave.F`, `WDES%NGVECTOR(NK) = WDES%NPLWKP(NK)/WDES%NRSPINORS`, where
              !!  `NRSPINORS=1` for our case. `WDES%NPLWKP(NK)` and `WDES%NPLWKP_TOT(K)`
              !!  are set the same way and neither `NGVECTOR(NK)` nor `NPLWKP_TOT(K)` are
              !!  changed anywhere else, so I am treating them as equivalent. I am not
              !!  sure why there are two separate variables defined.
              !! @endnote

          endif

          irec = irec + nBands
            ! Skip the records for the plane-wave coefficients for now.
            ! Those are read in the `readAndWriteWavefunction` subroutine.

        enddo
      enddo

      close(wavecarUnit)

      eigenE(:,:,:) = eigenE(:,:,:)*eVToHartree

    endif

    call MPI_BCAST(kPosition, size(kPosition), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(bandOccupation, size(bandOccupation), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(eigenE, size(eigenE), MPI_COMPLEX, root, worldComm, ierr)
    call MPI_BCAST(nPWs1kGlobal, size(nPWs1kGlobal), MPI_INTEGER, root, worldComm, ierr)

    return

  end subroutine preliminaryWAVECARScan

!----------------------------------------------------------------------------
  subroutine read_vasprun_xml(nKPoints, VASPDir, eFermi, eTot, kWeight, iType, nAtoms, nAtomsEachType, nAtomTypes)
    !! Read the k-point weights and cell info from the `vasprun.xml` file
    !!
    !! <h2>Walkthrough</h2>
    !!

    use miscUtilities, only: getFirstLineWithKeyword, ignoreNextNLinesFromFile

    implicit none

    ! Input variables:
    integer, intent(in) :: nKPoints
      !! Total number of k-points

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files


    ! Output variables:
    real(kind=dp), intent(out) :: eFermi
      !! Fermi energy
    real(kind=dp), intent(out) :: eTot
      !! Total energy
    real(kind=dp), allocatable, intent(out) :: kWeight(:)
      !! K-point weights

    integer, allocatable, intent(out) :: iType(:)
      !! Atom type index
    integer, intent(out) :: nAtoms
      !! Number of atoms
    integer, allocatable, intent(out) :: nAtomsEachType(:)
      !! Number of atoms of each type
    integer, intent(out) :: nAtomTypes
      !! Number of types of atoms


    ! Local variables:
    integer :: ik, ia
      !! Loop indices

    character(len=256) :: cDum
      !! Dummy variable to ignore input
    character(len=256) :: fileName
      !! `vasprun.xml` with path
    character(len=300) :: line
      !! Line read from file

    logical :: fileExists
      !! If the `vasprun.xml` file exists
    logical :: orbitalMag
      !! If can safely ignore `VKPT_SHIFT`
    logical :: spinSpiral
      !! If spin spirals considered (LSPIRAL)
    logical :: useRealProj
      !! If real-space projectors are used (LREAL)

    allocate(kWeight(nKPoints))

    if (ionode) then

      fileName = trim(VASPDir)//'/vasprun.xml'

      inquire(file = fileName, exist = fileExists)

      if (.not. fileExists) call exitError('read_vasprun_xml', 'Required file vasprun.xml does not exist', 1)

      open(57, file=fileName)
        !! * If root node, open `vasprun.xml`


      line = getFirstLineWithKeyword(57,'weights')
        !! * Ignore everything until you get to a
        !!   line with `'weights'`, indicating the
        !!   tag surrounding the k-point weights

      do ik = 1, nKPoints
        !! * Read in the weight for each k-point

        read(57,*) cDum, kWeight(ik), cDum

      enddo


      line = getFirstLineWithKeyword(57,'LREAL')
        !! * Ignore everything until you get to a
        !!   line with `'LREAL'`, indicating the
        !!   tag that determines if real-space 
        !!   projectors are used

      read(line,'(a35,L4,a4)') cDum, useRealProj, cDum

      if(useRealProj) call exitError('read_vasprun_xml', &
        '*** error - expected LREAL = F but got T', 1)


      line = getFirstLineWithKeyword(57,'LSPIRAL')
        !! * Ignore everything until you get to a
        !!   line with `'LSPIRAL'`, indicating the
        !!   tag that determines if spin spirals are
        !!   included

      read(line,'(a37,L4,a4)') cDum, spinSpiral, cDum

      if(spinSpiral) call exitError('read_vasprun_xml', &
        '*** error - expected LSPIRAL = F but got T', 1)


      line = getFirstLineWithKeyword(57,'ORBITALMAG')
        !! * Ignore everything until you get to a
        !!   line with `'ORBITALMAG'`, indicating the
        !!   tag that determines if can ignore `VKPT_SHIFT`

      read(line,'(a39,L4,a4)') cDum, orbitalMag, cDum

      if(orbitalMag) call exitError('read_vasprun_xml', &
        '*** error - expected ORBITALMAG = F but got T', 1)


      line = getFirstLineWithKeyword(57,'atominfo')
        !! * Ignore everything until you get to a
        !!   line with `'atominfo'`, indicating the
        !!   tag surrounding the cell info

      read(57,*) cDum, nAtoms, cDum
      read(57,*) cDum, nAtomTypes, cDum

      call ignoreNextNLinesFromFile(57, 5)

      allocate(iType(nAtoms), nAtomsEachType(nAtomTypes))

      nAtomsEachType = 0
      do ia = 1, nAtoms
        !! * Read in the atom type index for each atom
        !!   and calculate the number of atoms of each
        !!   type

        read(57,'(a21,i3,a9)') cDum, iType(ia), cDum

        nAtomsEachType(iType(ia)) = nAtomsEachType(iType(ia)) + 1

      enddo


      line = getFirstLineWithKeyword(57,'alphaZ')
        !! * Read until first occurence of `'alphaZ'`. 
        !!   This is the first step in the electronic
        !!   minimization. 

      line = getFirstLineWithKeyword(57,'alphaZ')
        !! * Read until next (last) occurence of `'alphaZ'`. 
        !!   This is the last step in the electronic
        !!   minimization. 
        !! @note
        !!    The `Export` code assumes that there are only
        !!    two occurences of `alphaZ` in the `vasprun.xml`
        !!    file indicating the beginning and end of the
        !!    electronic minimization. This may not be the
        !!    case for structural relaxations, but this code
        !!    should always be run on the SCF calculation that
        !!    follows a potential relaxation because that is
        !!    required for VASP to set the correct energy 
        !!    density.
        !! @endnote

      line = getFirstLineWithKeyword(57,'e_wo_entrp')
        !! * Ignore everything until you get to a
        !!   line with `'e_wo_entrp'`, indicating the
        !!   tag with the total energy

      read(line,'(a25,f16.8,a5)') cDum, eTot, cDum
      eTot = eTot*eVToRy


      line = getFirstLineWithKeyword(57,'efermi')
        !! * Ignore everything until you get to a
        !!   line with `'efermi'`, indicating the
        !!   tag with the Fermi energy

      read(line,*) cDum, cDum, eFermi, cDum
      eFermi = eFermi*eVToRy

      close(57)

    endif

    call MPI_BCAST(kWeight, size(kWeight), MPI_DOUBLE_PRECISION, root, worldComm, ierr)
    call MPI_BCAST(nAtoms, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nAtomTypes, 1, MPI_INTEGER, root, worldComm, ierr)

    if (.not. ionode) then
      allocate(iType(nAtoms))
      allocate(nAtomsEachType(nAtomTypes))
    endif

    call MPI_BCAST(iType, size(iType), MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(nAtomsEachType, size(nAtomsEachType), MPI_INTEGER, root, worldComm, ierr)

    return
  end subroutine read_vasprun_xml

!----------------------------------------------------------------------------
  subroutine calculateGvecs(fftGridSize, recipLattVec, gVecInCart, gIndexLocalToGlobal, gVecMillerIndicesGlobalOrig, gVecMillerIndicesGlobalSort, &
      iMill, nGVecsGlobal, nGVecsLocal)
    !! Calculate Miller indices and G-vectors and split
    !! over processors
    !!
    !! <h2>Walkthrough</h2>
    !!

    use generalComputations, only: direct2cart
    use miscUtilities, only: hpsort_eps

    implicit none

    ! Input variables:
    integer, intent(in) :: fftGridSize(3)
      !! Max number of points on the FFT grid in each direction

    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: gVecInCart(:,:)
      !! G-vectors in Cartesian coordinates

    integer, allocatable, intent(out) :: gIndexLocalToGlobal(:)
      !! Converts local index `ig` to global index
    integer, allocatable :: gVecMillerIndicesGlobalOrig(:,:)
      !! Integer coefficients for G-vectors on all processors (original order)
    integer, allocatable, intent(out) :: gVecMillerIndicesGlobalSort(:,:)
      !! Integer coefficients for G-vectors on all processors (sorted)
    integer, allocatable, intent(out) :: iMill(:)
      !! Indices of miller indices after sorting
    integer, intent(out) :: nGVecsGlobal
      !! Global number of G-vectors
    integer, intent(out) :: nGVecsLocal
      !! Local number of G-vectors on this processor


    ! Local variables:
    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp), allocatable :: millSum(:)
      !! Sum of integer coefficients for G-vectors
    real(kind=dp), allocatable :: realMillLocal(:,:)
      !! Real version of integer coefficients for G-vectors

    integer :: igx, igy, igz, ig
      !! Loop indices
    integer, allocatable :: gVecMillerIndicesLocal(:,:)
      !! Integer coefficients for G-vectors
    integer :: millX, millY, millZ
      !! Miller indices for each direction; in order
      !! 0,1,...,((fftGridSize(:)-1)/2),-((fftGridSize(:)-1)/2-1),...,-1
    integer :: npmax
      !! Max number of plane waves


    npmax = fftGridSize(1)*fftGridSize(2)*fftGridSize(3) 
    allocate(gVecMillerIndicesGlobalSort(3,npmax))
    allocate(millSum(npmax))

    if(ionode) then

      allocate(gVecMillerIndicesGlobalOrig(3,npmax))

      nGVecsGlobal = 0
      gVecMillerIndicesGlobalOrig = 0

      !> * Generate Miller indices for every possible G-vector
      !>   regardless of the \(|G+k|\) cutoff
      do igz = 1, fftGridSize(3)

        millZ = igz - 1

        if (igz - 1 .gt. (fftGridSize(3)-1)/2) millZ = igz - 1 - fftGridSize(3)

        do igy = 1, fftGridSize(2)

          millY = igy - 1

          if (igy - 1 .gt. (fftGridSize(2)-1)/2) millY = igy - 1 - fftGridSize(2)

          do igx = 1, fftGridSize(1)

            millX = igx - 1

            if (igx - 1 .gt. (fftGridSize(1)-1)/2) millX = igx - 1 - fftGridSize(1)

            nGVecsGlobal = nGVecsGlobal + 1

            gVecMillerIndicesGlobalOrig(1,nGVecsGlobal) = millX
            gVecMillerIndicesGlobalOrig(2,nGVecsGlobal) = millY
            gVecMillerIndicesGlobalOrig(3,nGVecsGlobal) = millZ
              !! * Calculate Miller indices

            millSum(nGVecsGlobal) = sqrt(real(millX**2 + millY**2 + millZ**2))
              !! * Calculate the sum of the Miller indices
              !!   for sorting

          enddo
        enddo
      enddo

      if (nGVecsGlobal .ne. npmax) call exitError('calculateGvecs', & 
        '*** error - computed no. of G-vectors != estimated number of plane waves', 1)
        !! * Check that number of G-vectors are the same as the number of plane waves

      allocate(iMill(nGVecsGlobal))

      do ig = 1, nGVecsGlobal
        !! * Initialize the index array that will track elements
        !!   after sorting

        iMill(ig) = ig

      enddo

      call hpsort_eps(nGVecsGlobal, millSum, iMill, eps8)
        !! * Order indices `iMill` by the G-vector length `millSum`

    endif

    deallocate(millSum)

    if(ionode) then

      do ig = 1, nGVecsGlobal
        !! * Rearrange the miller indices to match order of `millSum`

        gVecMillerIndicesGlobalSort(:,ig) = gVecMillerIndicesGlobalOrig(:,iMill(ig))

      enddo

    endif

    call MPI_BCAST(nGVecsGlobal, 1, MPI_INTEGER, root, worldComm, ierr)
    call MPI_BCAST(gVecMillerIndicesGlobalSort, size(gVecMillerIndicesGlobalSort), MPI_INTEGER, root, worldComm, ierr)


    call distributeGvecsOverProcessors(nGVecsGlobal, gVecMillerIndicesGlobalSort, gIndexLocalToGlobal, gVecMillerIndicesLocal, nGVecsLocal)
      !! * Split up the G-vectors and Miller indices over processors 

    allocate(gVecInCart(3,nGVecsLocal))
    allocate(realMillLocal(3,nGVecsLocal))

    realMillLocal = real(gVecMillerIndicesLocal)
    deallocate(gVecMillerIndicesLocal)

    gVecInCart = direct2cart(nGVecsLocal, realMillLocal, recipLattVec)

    deallocate(realMillLocal)

    return
  end subroutine calculateGvecs

!----------------------------------------------------------------------------
  subroutine distributeGvecsOverProcessors(nGVecsGlobal, gVecMillerIndicesGlobalSort, gIndexLocalToGlobal, gVecMillerIndicesLocal, nGVecsLocal)
    !! Figure out how many G-vectors there should be per processor.
    !! G-vectors are split up in a round robin fashion over processors
    !! in a single k-point pool.
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors
    !integer, intent(in) :: nProcPerPool
      ! Number of processes per pool
    integer, intent(in) :: gVecMillerIndicesGlobalSort(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors

    
    ! Output variables:
    integer, allocatable, intent(out) :: gIndexLocalToGlobal(:)
      !! Converts local index `ig` to global index
    integer, allocatable, intent(out) :: gVecMillerIndicesLocal(:,:)
      !! Integer coefficients for G-vectors
    integer, intent(out) :: nGVecsLocal
      !! Local number of G-vectors on this processor


    ! Local variables:
    integer :: ig_l, ig_g
      !! Loop indices
    integer :: ngr
      !! Number of G-vectors left over after evenly divided across processors


    if( nGVecsGlobal > 0 ) then
      nGVecsLocal = nGVecsGlobal/nProcPerPool
        !!  * Calculate number of G-vectors per processor

      ngr = nGVecsGlobal - nGVecsLocal*nProcPerPool 
        !! * Calculate the remainder

      if( indexInPool < ngr ) nGVecsLocal = nGVecsLocal + 1
        !! * Assign the remainder to the first `ngr` processors

      !> * Generate an array to map a local index
      !>   (`ig` passed to `gIndexLocalToGlobal`) to a global
      !>   index (the value stored at `gIndexLocalToGlobal(ig)`)
      !>   and get local miller indices
      allocate(gIndexLocalToGlobal(nGVecsLocal))
      allocate(gVecMillerIndicesLocal(3,nGVecsLocal))

      ig_l = 0
      do ig_g = 1, nGVecsGlobal

        if(indexInPool == mod(ig_g-1,nProcPerPool)) then
        
          ig_l = ig_l + 1
          gIndexLocalToGlobal(ig_l) = ig_g
          gVecMillerIndicesLocal(:,ig_l) = gVecMillerIndicesGlobalSort(:,ig_g)

        endif

      enddo

      if (ig_l /= nGVecsLocal) call exitError('distributeGvecsOverProcessors', 'unexpected number of G-vecs for this processor', 1)

    endif

    return
  end subroutine distributeGvecsOverProcessors

!----------------------------------------------------------------------------
  subroutine reconstructFFTGrid(nGVecsLocal, gIndexLocalToGlobal, nKPoints, nPWs1kGlobal, kPosition, gVecInCart, recipLattVec, &
      wfcVecCut, igkSort2OrigLocal, maxGIndexGlobal, maxGkVecsLocal, nGkLessECutGlobal, nGkVecsLocal)
    !! Determine which G-vectors result in \(G+k\)
    !! below the energy cutoff for each k-point and
    !! sort the indices based on \(|G+k|^2\)
    !!
    !! <h2>Walkthrough</h2>
    !!

    use miscUtilities, only: hpsort_eps

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsLocal
      !! Number of G-vectors on this processor
    !integer, intent(in) :: nGVecsGlobal
      ! Global number of G-vectors
    integer, intent(in) :: gIndexLocalToGlobal(nGVecsLocal)
      ! Converts local index `ig` to global index
    !integer, intent(in) :: iMill(nGVecsGlobal)
      ! Indices of miller indices after sorting
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: nPWs1kGlobal(nKPoints)
      !! Input number of plane waves for a single k-point

    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space
    real(kind=dp), intent(in) :: gVecInCart(3,nGVecsLocal)
      !! G-vectors in Cartesian coordinates
    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: wfcVecCut
      !! Energy cutoff converted to vector cutoff


    ! Output variables:
    integer, allocatable, intent(out) :: igkSort2OrigLocal(:,:)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order
    integer, intent(out) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(out) :: maxGkVecsLocal
      !! Max number of G+k vectors across all k-points
      !! in this pool
    integer, intent(out) :: nGkLessECutGlobal(nKPoints)
      !! Global number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` for each k-point
    integer, allocatable, intent(out) :: nGkVecsLocal(:)
      !! Local number of G+k-vectors on this processor


    ! Local variables:
    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp) :: gkMod(nkPerPool,nGVecsLocal)
      !! \(|G+k|^2\);
      !! only stored if less than `wfcVecCut`
    real(kind=dp) :: q
      !! \(|q|^2\) where \(q = G+k\)
    real(kind=dp), allocatable :: realiMillGk(:)
      !! Indices of miller indices after sorting
    real(kind=dp) :: xkCart(3)
      !! Cartesian coordinates for given k-point

    integer :: ik, igk, ig
      !! Loop indices
    integer, allocatable :: igk2igGlobal(:,:)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors
    integer, allocatable :: gKIndexLocalToGlobal(:,:)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point
    integer :: igk2igLocal(nkPerPool,nGVecsLocal)
      !! Map from local \(G+k\) index to global
      !! \(G\); indexed up to `nGVecsLocal` which
      !! is greater than `maxNumPWsPool` and
      !! stored for each k-point
    integer, allocatable :: igk2igLocal_ik(:)
      !! Map from local \(G+k\) index to global
      !! \(G\) index for single k-point
    integer, allocatable :: igkSort2OrigGlobal(:,:)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors in the original order
    integer :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` among all k-points
    integer :: maxNumPWsPool
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! pool
    integer :: ngk_tmp
      !! Temporary variable to hold `nGkLessECutLocal`
      !! value so that don't have to keep accessing
      !! array
    integer :: nGkLessECutLocal(nkPerPool)
      !! Number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` for each
      !! k-point, on this processor
    integer :: maxGIndexLocal
      !! Maximum G-vector index among all \(G+k\)
      !! for just this processor
    integer :: maxNumPWsLocal
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    
    maxNumPWsLocal = 0
    nGkLessECutLocal(:) = 0
    igk2igLocal(:,:) = 0

    do ik = 1, nkPerPool
      !! * For each \(G+k\) combination, calculate the 
      !!   magnitude and, if it is less than the energy
      !!   cutoff, store the G index and magnitude and 
      !!   increment the number of \(G+k\) vectors at
      !!   this k-point. Also, keep track of the maximum 
      !!   number of \(G+k\) vectors among all k-points
      !!
      !! @note
      !!  All of the above calculations are local to a single
      !!  processor.
      !! @endnote

      xkCart = matmul(recipLattVec,kPosition(:,ik+ikStart_pool-1))

      ngk_tmp = 0

      do ig = 1, nGVecsLocal

        q = sqrt(sum((xkCart(:) + gVecInCart(:,ig))**2))
          ! Calculate \(|G+k|\)

        if (q <= eps8) q = 0.d0

        if (q <= wfcVecCut) then

          ngk_tmp = ngk_tmp + 1
            ! If \(|G+k| \leq \) `wfcVecCut` increment the count for
            ! this k-point

          gkMod(ik,ngk_tmp) = q
            ! Store the modulus for sorting

          igk2igLocal(ik,ngk_tmp) = ig
            ! Store the index for this G-vector

        !else

          !if (sqrt(sum(gVecInCart(:, ig)**2)) .gt. &
            !sqrt(sum(kPosition(:,ik+ikStart_pool-1)**2) + sqrt(wfcVecCut))) goto 100
            ! if |G| > |k| + sqrt(Ecut)  stop search
            !! @todo Figure out if there is valid exit check for `ig` loop @endtodo

        endif
      enddo

      if (ngk_tmp == 0) call exitError('reconstructFFTGrid', 'no G+k vectors on this processor', 1) 

100   maxNumPWsLocal = max(maxNumPWsLocal, ngk_tmp)
        ! Track the maximum number of \(G+k\)
        ! vectors among all k-points

      nGkLessECutLocal(ik) = ngk_tmp
        ! Store the total number of \(G+k\)
        ! vectors for this k-point

    enddo

    !nGkLessECutGlobal = 0
      ! This is initialized in the main program, but I wanted
      ! to include this here for clarity on how the sum across
      ! processes works. 
    nGkLessECutGlobal(ikStart_pool:ikEnd_pool) = nGkLessECutLocal(1:nkPerPool)
    CALL mpiSumIntV(nGkLessECutGlobal, worldComm)
      !! * Calculate the global number of \(G+k\) 
      !!   vectors for each k-point
      
    if (ionode) then

      do ik = 1, nKPoints

        if (nGkLessECutGlobal(ik) .ne. nPWs1kGlobal(ik)) call exitError('reconstructFFTGrid', &
          'computed no. of G-vectors != input no. of plane waves', 1)
          !! * Make sure that number of G-vectors isn't higher than the calculated maximum

      enddo
    endif

    if (maxNumPWsLocal <= 0) call exitError('reconstructFFTGrid', &
                'No plane waves found: running on too many processors?', 1)
      !! * Make sure that each processor gets some \(G+k\) vectors. If not,
      !!   should rerun with fewer processors.

    call MPI_ALLREDUCE(maxNumPWsLocal, maxNumPWsPool, 1, MPI_INTEGER, MPI_MAX, intraPoolComm, ierr)
    if(ierr /= 0) call exitError('reconstructFFTGrid', 'error in mpi_allreduce 1', ierr)
      !! * When using pools, set `maxNumPWsPool` to the maximum value of `maxNumPWsLocal` 
      !!   in the pool 


    allocate(gKIndexLocalToGlobal(maxNumPWsPool,nkPerPool))
    allocate(igk2igLocal_ik(maxNumPWsPool))

    gKIndexLocalToGlobal = 0
    igk2igLocal_ik = 0


    do ik = 1, nkPerPool
      !! * Reorder the indices of the G-vectors so that
      !!   they are sorted by \(|G+k|^2\) for each k-point

      ngk_tmp = nGkLessECutLocal(ik)

      igk2igLocal_ik(1:ngk_tmp) = igk2igLocal(ik,1:ngk_tmp)

      call hpsort_eps(ngk_tmp, gkMod(ik,:), igk2igLocal_ik, eps8)
        ! Order vector `igk` by \(|G+k|\) (`gkMod`)

      do igk = 1, ngk_tmp
        
        gKIndexLocalToGlobal(igk,ik) = gIndexLocalToGlobal(igk2igLocal_ik(igk))
        
      enddo
     
      gKIndexLocalToGlobal(ngk_tmp+1:maxNumPWsPool, ik) = 0

    enddo


    deallocate(igk2igLocal_ik)


    maxGIndexLocal = maxval(gKIndexLocalToGlobal(:,:))
    call MPI_ALLREDUCE(maxGIndexLocal, maxGIndexGlobal, 1, MPI_INTEGER, MPI_MAX, worldComm, ierr)
    if(ierr /= 0) call exitError('reconstructFFTGrid', 'error in mpi_allreduce 2', ierr)
      !! * Calculate the maximum G-vector index 
      !!   among all \(G+k\) and processors

    maxNumPWsGlobal = maxval(nGkLessECutGlobal(1:nKPoints))
      !! * Calculate the maximum number of G-vectors 
      !!   among all k-points

    allocate(igk2igGlobal(maxNumPWsGlobal, nKPoints))

  
    igk2igGlobal(:,:) = 0
    do ik = 1, nKPoints

      call getGlobalGkIndices(nKPoints, maxNumPWsPool, gKIndexLocalToGlobal, ik, maxGIndexGlobal, maxNumPWsGlobal, nGkLessECutGlobal, &
          nGkLessECutLocal, igk2igGlobal)
        !! * For each k-point, gather all of the \(G+k\) indices
        !!   among all processors in a single global array
    
    enddo

    deallocate(gKIndexLocalToGlobal)

    allocate(igkSort2OrigGlobal(maxNumPWsGlobal, nKPoints))

    igkSort2OrigGlobal = igk2igGlobal

    if(ionode) then

      allocate(realiMillGk(maxNumPWsGlobal))

      do ik = 1, nKPoints

        realiMillGk = 0._dp
        ngk_tmp = nGkLessECutGlobal(ik)

        do ig = 1, ngk_tmp

          realiMillGk(ig) = real(iMill(igk2igGlobal(ig,ik)))
            !! * Get only the original indices that correspond
            !!   to G vectors s.t. \(|G+k|\) is less than the
            !!   cutoff

        enddo

        call hpsort_eps(ngk_tmp, realiMillGk(1:ngk_tmp), igkSort2OrigGlobal(1:ngk_tmp,ik), eps8)
          !! * Order the \(G+k\) indices by the original indices `realiMillGk`. 
          !!   This will allow us to recover only specific G-vectors in the 
          !!   original ordering. Have to cast to `real` because that is what the 
          !!   sorting algorithm expects. Would be better to have a different
          !!   interface for different types, but we don't really need that here
          !!   and it shouldn't affect the results. 

      enddo

      deallocate(iMill)
      deallocate(realiMillGk)

    endif

    deallocate(igk2igGlobal)

    call MPI_BCAST(igkSort2OrigGlobal, size(igkSort2OrigGlobal), MPI_INTEGER, root, worldComm, ierr)

    call distributeGkVecsInPool(maxNumPWsGlobal, igkSort2OrigGlobal, nKPoints, nGkLessECutGlobal, igkSort2OrigLocal, maxGkVecsLocal, nGkVecsLocal)
      !! * Distribute the G+k vectors evenly across the processes in a band group within each pool

    deallocate(igkSort2OrigGlobal)

    return
  end subroutine reconstructFFTGrid

!----------------------------------------------------------------------------
  subroutine getGlobalGkIndices(nKPoints, maxNumPWsPool, gKIndexLocalToGlobal, ik, maxGIndexGlobal, maxNumPWsGlobal, nGkLessECutGlobal, &
      nGkLessECutLocal, igk2igGlobal)
    !! Gather the \(G+k\) vector indices in single, global 
    !! array
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: maxNumPWsPool
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    integer, intent(in) :: gKIndexLocalToGlobal(maxNumPWsPool, nkPerPool)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point;
      !! the first index goes up to `maxNumPWsPool`,
      !! but only valid values are up to `nGkLessECutLocal`
    integer, intent(in) :: ik
      !! Index of current k-point
    integer, intent(in) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(in) :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` among all k-points
    integer, intent(in) :: nGkLessECutGlobal(nKPoints)
      !! Global number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` for each k-point
    integer, intent(in) :: nGkLessECutLocal(nkPerPool)
      !! Number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` for each
      !! k-point, on this processor


    ! Output variables:
    integer, intent(out) :: igk2igGlobal(maxNumPWsGlobal, nKPoints)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors


    ! Local variables:
    integer :: ig
    integer, allocatable :: itmp1(:)
      !! Global \(G+k\) indices for single
      !! k-point with zeros for G-vector indices
      !! where \(G+k\) was greater than the cutoff
    integer :: ngg 
      !! Counter for \(G+k\) vectors for given
      !! k-point; should equal `nGkLessECutGlobal`

    
    allocate(itmp1(maxGIndexGlobal), stat=ierr)
    if (ierr/= 0) call exitError('getGlobalGkIndices','allocating itmp1', abs(ierr))

    itmp1 = 0
    if(ik >= ikStart_pool .and. ik <= ikEnd_pool) then

      do ig = 1, nGkLessECutLocal(ik-ikStart_pool+1)

        itmp1(gKIndexLocalToGlobal(ig, ik-ikStart_pool+1)) = gKIndexLocalToGlobal(ig, ik-ikStart_pool+1)
          !! * For each k-point and \(G+k\) vector for this processor,
          !!   store the local to global indices (`gKIndexLocalToGlobal`) in an
          !!   array that will later be combined globally
          !!
          !! @note
          !!  This will leave zeros in spots where the \(G+k\) 
          !!  combination for this k-point was greater than the energy 
          !!  cutoff.
          !! @endnote

      enddo
    endif

    call mpiSumIntV(itmp1, worldComm)

    ngg = 0
    do  ig = 1, maxGIndexGlobal

      if(itmp1(ig) == ig) then
        !! * Go through and find all of the non-zero
        !!   indices in the now-global `itmp1` array,
        !!   and store them in a new array that won't
        !!   have the extra zeros

        ngg = ngg + 1

        igk2igGlobal(ngg, ik) = ig

      endif
    enddo


    if(ionode .and. ngg /= nGkLessECutGlobal(ik)) call exitError('getGlobalGkIndices', 'Unexpected number of G+k vectors', 1)
      !! * Make sure that the total number of non-zero
      !!   indices matches the global number of \(G+k\)
      !!   vectors for this k-point
    
    deallocate( itmp1 )

    return
  end subroutine getGlobalGkIndices

!----------------------------------------------------------------------------
  subroutine distributeGkVecsInPool(maxNumPWsGlobal, igkSort2OrigGlobal, nKPoints, nGkLessECutGlobal, igkSort2OrigLocal, maxGkVecsLocal, nGkVecsLocal)
    !! Distribute the G+k vectors across band groups in each pool by 
    !! splitting up the `igkSort2OrigGlobal` array
    !! into local arrays
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: maxNumPWsGlobal
      !! Max number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` among all k-points
    integer, intent(in) :: igkSort2OrigGlobal(maxNumPWsGlobal, nKPoints)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors in the original order
    !integer, intent(in) :: ikEnd_pool
      ! Ending index for k-points in single pool 
    !integer, intent(in) :: ikStart_pool
      ! Starting index for k-points in single pool 
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: nGkLessECutGlobal(nKPoints)
      !! Global number of G+k-vectors
    !integer, intent(in) :: nProcPerPool
      ! Number of processes per pool

    
    ! Output variables:
    !integer, allocatable, intent(out) :: iGkEnd_pool(:)
      ! Ending index for G+k vectors on
      ! single process in a given pool
    !integer, allocatable, intent(out) :: iGkStart_pool(:)
      ! Starting index for G+k vectors on
      ! single process in a given pool
    integer, allocatable, intent(out) :: igkSort2OrigLocal(:,:)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order
    integer, intent(out) :: maxGkVecsLocal
      !! Max number of G+k vectors across all k-points
      !! in this pool
    integer, allocatable, intent(out) :: nGkVecsLocal(:)
      !! Local number of G+k-vectors on this processor


    ! Local variables:
    integer :: ik
      !! Loop indices
    integer :: ngkr
      !! Number of G+k vectors left over after evenly 
      !! divided across processors in pool


    allocate(nGkVecsLocal(nkPerPool), iGkStart_pool(nkPerPool), iGkEnd_pool(nkPerPool))

    do ik = 1, nkPerPool
      nGkVecsLocal(ik) = nGkLessECutGlobal(ik+ikStart_pool-1)/nProcPerBgrp
        !!  * Calculate the number of G+k vectors per processors
        !!    in a band group

      ngkr = nGkLessECutGlobal(ik+ikStart_pool-1) - nGkVecsLocal(ik)*nProcPerBgrp
        !! * Calculate the remainder

      if( indexInBgrp < ngkr ) nGkVecsLocal(ik) = nGkVecsLocal(ik) + 1
        !! * Assign the remainder to the first `ngr` processors

      !>  * Calculate the index of the first G+k vector for this process
      iGkStart_pool(ik) = nGkVecsLocal(ik) * indexInBgrp + 1
      if( indexInBgrp >= ngkr ) iGkStart_pool(ik) = iGkStart_pool(ik) + ngkr

      iGkEnd_pool(ik) = iGkStart_pool(ik) + nGkVecsLocal(ik) - 1
        !!  * Calculate the index of the last G+k vector in a band group in this pool

    enddo

    maxGkVecsLocal = maxval(nGkVecsLocal)
      !! * Get the max number of G+k vectors across
      !!   all k-points in this pool

    allocate(igkSort2OrigLocal(maxGkVecsLocal, nkPerPool))

    do ik = 1, nkPerPool

      igkSort2OrigLocal(1:nGkVecsLocal(ik),ik) = igkSort2OrigGlobal(iGkStart_pool(ik):iGkEnd_pool(ik),ik+ikStart_pool-1)
        !! * Split up the PWs `igkSort2OrigGlobal` across processors and 
        !!   store the G-vector indices locally

    enddo
      
    return
  end subroutine distributeGkVecsInPool

!----------------------------------------------------------------------------
  subroutine readPOTCAR(nAtomTypes, VASPDir, pot)
    !! Read PAW pseudopotential information from POTCAR
    !! file
    !!
    !! <h2>Walkthrough</h2>
    !!

    use miscUtilities, only: getFirstLineWithKeyword

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    
    ! Output variables:
    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR


    ! Local variables:
    real(kind=dp) :: dummyD(1000)
      !! Dummy variable to ignore input
    real(kind=dp), allocatable :: dummyDA1(:), dummyDA2(:,:)
      !! Allocatable dummy variable to ignore input
    real(kind=dp) :: H
      !! Factor for generating derivative of 
      !! radial grid

    integer :: angMom
      !! Angular momentum of projectors
    integer :: iT, ip, ir, i
      !! Loop indices
    integer :: nProj
      !! Number of projectors with given angular momentum

    character(len=1) :: charSwitch
      !! Switch to determine what section reading
    character(len=256) :: dummyC
      !! Dummy character to ignore input
    character(len=256) :: fileName
      !! Full WAVECAR file name including path
    character(len=300) :: line
      !! Line read from file

    logical :: found
      !! If min index has been found


    if(ionode) then
      fileName = trim(VASPDir)//'/POTCAR'

      open(unit=potcarUnit, file=fileName, iostat=ierr, status='old')
      if (ierr .ne. 0) write(*,*) 'open error - iostat =', ierr
        !! * If root node, open the `POTCAR` file

      do iT = 1, nAtomTypes

        pot(iT)%nChannels = 0
        pot(iT)%lmmax = 0

        read(potcarUnit,*) dummyC, pot(iT)%element, dummyC
          !! * Read in the header
        read(potcarUnit,*)
          !! * Ignore the valence line
        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read in character switch to determine if there is a 
          !!   PSCRT section (switch not used)
          !! @note
          !!  Some of the switches do not actually seem to be used
          !!  as a switch because the following code does not include
          !!  any logic to process the switch. If the POTCAR files 
          !!  ever have a different form than assumed here, the logic
          !!  will need to be updated.
          !! @endnote

        line = getFirstLineWithKeyword(potcarUnit,'END')
          !! * Ignore all lines until you get to the `END` of
          !!   the PSCRT section

        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch (switch not used)
        read(potcarUnit,*)
          !! * Ignore the max G for local potential
        read(potcarUnit,*) (dummyD(i), i=1,1000)
          !! * Ignore the local pseudopotential in reciprocal
          !!   space
        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch

        if (charSwitch == 'g') then
          !! * Ignore gradient correction

          read(potcarUnit,*)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        if (charSwitch == 'c') then
          !! * Ignore core charge density

          read(potcarUnit,*) (dummyD(i), i=1,1000)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        if (charSwitch == 'k') then
          !! * Ignore partial kinetic energy density

          read(potcarUnit,*) (dummyD(i), i=1,1000)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        if (charSwitch == 'K') then
          !! * Ignore kinetic energy density

          read(potcarUnit,*) (dummyD(i), i=1,1000)
          read(potcarUnit,'(1X,A1)') charSwitch
          
        endif

        read(potcarUnit,*) (dummyD(i), i=1,1000)
          !! * Ignore the atomic pseudo charge density

        read(potcarUnit,*) pot(iT)%maxGkNonlPs, dummyC
          !! * Read the max \(|G+k|\) for non-local potential 
          !!   and ignore unused boolean (`LDUM` in VASP)

        pot(iT)%maxGkNonlPs = pot(iT)%maxGkNonlPs/angToBohr

        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch

        allocate(dummyDA1(nonlPseudoGridSize))

        do while (charSwitch /= 'D' .and. charSwitch /= 'A' .and. charSwitch /= 'P' &
          .and. charSwitch /= 'E')
            !! * Until you have read in all of the momentum channels
            !!   (i.e. you get to a character switch that is not `'N'`)
            !!     * Read in the angular momentum, the number of 
            !!       projectors at this angular momentum, and the max
            !!       r for the non-local contribution
            !!     * Increment the number of nlm channels
            !!     * Ignore non-local strength multipliers
            !!     * Read in the reciprocal-space projectors and set
            !!       boundary
            !!     * Increment the number of l channels
            !!     * Read the next character switch

          read(potcarUnit,*) angMom, nProj, pot(iT)%psRMax
            ! Read in angular momentum, the number of projectors
            ! at this angular momentum, and the max r for the 
            ! non-local contribution

          pot(iT)%lmmax = pot(iT)%lmmax + (2*angMom+1)*nProj
            ! Increment the number of nlm channels

          allocate(dummyDA2(nProj,nProj))

          read(potcarUnit,*) dummyDA2(:,:)
            ! Ignore non-local strength multipliers

          do ip = 1, nProj
            ! Read in the reciprocal-space and real-space
            ! projectors

            pot(iT)%angMom(pot(iT)%nChannels+ip) = angMom

            read(potcarUnit,*) 
            read(potcarUnit,*) (pot(iT)%recipProj(pot(iT)%nChannels+ip,i), i=1,nonlPseudoGridSize)
              ! Read in reciprocal-space projector
              ! I believe these units are Ang^(3/2). When multiplied by `1/sqrt(omega)`,
              ! the projectors are then unitless. 

            ! Not really sure what the purpose of this is. Seems to be setting the grid boundary,
            ! but I'm not sure on the logic.
            if(mod(angMom,2) == 0) then
              pot(iT)%recipProj(pot(iT)%nChannels+ip, 0) = pot(iT)%recipProj(pot(iT)%nChannels+ip, 2) 
            else
              pot(iT)%recipProj(pot(iT)%nChannels+ip, 0) = -pot(iT)%recipProj(pot(iT)%nChannels+ip, 2) 
            endif

            read(potcarUnit,*) 
            read(potcarUnit,*) (dummyDA1(i), i=1,nonlPseudoGridSize)
              ! Ignore real-space projector

          enddo

          pot(iT)%nChannels = pot(iT)%nChannels + nProj
            ! Increment the number of l channels

          deallocate(dummyDA2)

          read(potcarUnit,'(1X,A1)') charSwitch
            ! Read character switch

        enddo

        deallocate(dummyDA1)

        if (charSwitch /= 'P') then
          !! * Ignore depletion charges

          read(potcarUnit,*)
          read(potcarUnit,*)
          
        else

          read(potcarUnit,*) pot(iT)%nmax, pot(iT)%rAugMax  
            !! * Read the number of mesh grid points and
            !!   the maximum radius in the augmentation sphere

          pot(iT)%rAugMax = pot(iT)%rAugMax*angToBohr 
            ! Convert units 

          read(potcarUnit,*)
            !! * Ignore format specifier
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch (not used)

          allocate(pot(iT)%radGrid(pot(iT)%nmax))
          allocate(pot(iT)%dRadGrid(pot(iT)%nmax))
          allocate(pot(iT)%wps(pot(iT)%nChannels,pot(iT)%nmax))
          allocate(pot(iT)%wae(pot(iT)%nChannels,pot(iT)%nmax))
          allocate(dummyDA2(pot(iT)%nChannels, pot(iT)%nChannels))

          read(potcarUnit,*) dummyDA2(:,:)
            !! * Ignore augmentation charges
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch

          if (charSwitch == 't') then
            !! * Ignore total charge in each channel 

            read(potcarUnit,*) dummyDA2(:,:)
            read(potcarUnit,*) 

          endif

          read(potcarUnit,*) dummyDA2
            !! * Ignore initial occupancies in atom

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'g') call exitError('readPOTCAR', 'expected grid section', 1)

          read(potcarUnit,*) (pot(iT)%radGrid(i), i=1,pot(iT)%nmax)

          pot(iT)%radGrid(:) = pot(iT)%radGrid(:)*angToBohr

          H = log(pot(iT)%radGrid(pot(iT)%nmax)/pot(iT)%radGrid(1))/(pot(iT)%nmax - 1)
            !! * Calculate \(H\) which is used to generate the derivative of the grid
            !!
            !! @note
            !!  The grid in VASP is defined as \(R_i = R_0e^{H(i-1)}\), so we define the
            !!  derivative as \(dR_i = R_0He^{H(i-1)}\)
            !! @endnote
          
          found = .false.
          do ir = 1, pot(iT)%nmax
            !! * Calculate the max index of the augmentation sphere and
            !!   the derivative of the radial grid

            if (.not. found .and. pot(iT)%radGrid(ir) > pot(iT)%rAugMax) then
              pot(iT)%iRAugMax = ir - 1
              found = .true.
            endif

            pot(iT)%dRadGrid(ir) = pot(iT)%radGrid(1)*H*exp(H*(ir-1))

          enddo

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'a') call exitError('readPOTCAR', 'expected aepotential section', 1)

          allocate(dummyDA1(pot(iT)%nmax))

          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore AE potential

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'c') call exitError('readPOTCAR', 'expected core charge-density section', 1)

          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore the frozen core charge
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch

          if (charSwitch == 'k') then
            !! * Ignore kinetic energy density

            read(potcarUnit,*) dummyDA1(:)
            read(potcarUnit,'(1X,A1)') charSwitch

          endif

          if (charSwitch == 'm') then
            !! * Ignore pseudo-ized kinetic energy density

            read(potcarUnit,*) dummyDA1(:)
            read(potcarUnit,'(1X,A1)') charSwitch

          endif

          if (charSwitch == 'l') then
            !! * Ignore local pseudopotential core

            read(potcarUnit,*) dummyDA1(:)
            read(potcarUnit,'(1X,A1)') charSwitch

          endif

          if (charSwitch /= 'p') call exitError('readPOTCAR', 'expected pspotential section', 1)
          
          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore PS potential

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'c') call exitError('readPOTCAR', 'expected core charge-density section', 1)
          
          read(potcarUnit,*) dummyDA1(:)
            !! * Ignore core charge density

          do ip = 1, pot(iT)%nChannels
            !! * Read the AE and PS partial waves for each projector
            
            read(potcarUnit,'(1X,A1)') charSwitch
            if (charSwitch /= 'p') call exitError('readPOTCAR', 'expected pseudowavefunction section', 1)
            read(potcarUnit,*) (pot(iT)%wps(ip,i), i=1,pot(iT)%nmax)

            read(potcarUnit,'(1X,A1)') charSwitch
            if (charSwitch /= 'a') call exitError('readPOTCAR', 'expected aewavefunction section', 1)
            read(potcarUnit,*) (pot(iT)%wae(ip,i), i=1,pot(iT)%nmax)

          enddo

          pot(iT)%wps(:,:) = pot(iT)%wps(:,:)/sqrt(angToBohr)
          pot(iT)%wae(:,:) = pot(iT)%wae(:,:)/sqrt(angToBohr)
            !! @note
            !!  Based on the fact that this does not have an x/y/z
            !!  dimension and that these values get multiplied by 
            !!  `radGrid` and `dRadGrid`, which we treat as being 
            !!  one dimensional, I think these are one dimensional 
            !!  and should just have `1/sqrt(angToBohr)`. That has 
            !!  worked well in the past too.
            !! @endnote

          deallocate(dummyDA1)
          deallocate(dummyDA2)

        endif

        line = getFirstLineWithKeyword(potcarUnit,'End of Dataset')
          !! * Ignore all lines until you get to the `End of Dataset`

      enddo

    endif    

    do iT = 1, nAtomTypes

      call MPI_BCAST(pot(iT)%nChannels, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(pot(iT)%lmmax, 1, MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(pot(iT)%maxGkNonlPs, 1, MPI_DOUBLE_PRECISION, root, worldComm, ierr)
      call MPI_BCAST(pot(iT)%angmom, size(pot(iT)%angmom), MPI_INTEGER, root, worldComm, ierr)
      call MPI_BCAST(pot(iT)%recipProj, size(pot(iT)%recipProj), MPI_DOUBLE_PRECISION, root, worldComm, ierr)

    enddo

    return
  end subroutine readPOTCAR

!----------------------------------------------------------------------------
  subroutine projAndWav(ispSelect, maxGkVecsLocal, nAtoms, nAtomTypes, nBands, nGkVecsLocal, nGVecsGlobal, nKPoints, nSpins, &
      gVecMillerIndicesGlobalSort, igkSort2OrigLocal, nPWs1kGlobal, reclenWav, atomPositionsDir, kPosition, omega, recipLattVec, &
      exportDir, VASPDir, gammaOnly, loopSpins, pot)

    implicit none

    ! Input variables: 
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: maxGkVecsLocal
      !! Max number of G+k vectors across all k-points
      !! in this pool
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nGkVecsLocal(nkPerPool)
      !! Local number of G+k-vectors on this processor
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nSpins
      !! Number of spins
    integer, intent(in) :: gVecMillerIndicesGlobalSort(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(in) :: igkSort2OrigLocal(maxGkVecsLocal,nkPerPool)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order
    integer, intent(in) :: nPWs1kGlobal(nKPoints)
      !! Input number of plane waves for a single k-point
    integer, intent(in) :: reclenWav
      !! Number of records in the WAVECAR file

    real(kind=dp), intent(in) :: atomPositionsDir(3,nAtoms)
      !! Atom positions
    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space
    real(kind=dp), intent(in) :: omega
      !! Volume of unit cell
    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors

    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export
    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    logical, intent(in) :: gammaOnly
      !! If the gamma only VASP code is used
    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR

    ! Local variables:
    integer, allocatable :: igkSort2OrigLocal_ik(:)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order for a
      !! given k-point
    integer :: nGkVecsLocal_ik
      !! Number of G+k vectors locally for 
      !! a given k-point
    integer :: nProj
      !! Number of projectors across all atom types
    integer :: nPWs1k
      !! Input number of plane waves for the given k-point
    integer :: irec
      !! Record number in WAVECAR file;
      !! needed for shared access
    integer :: ikLocal, ikGlobal, isp, isk
      !! Loop indices

    real(kind=dp), allocatable :: realProjWoPhase(:,:,:)
      !! Real projectors without phase
    real(kind=dp) :: t1, t2
      !! Timers

    complex*8, allocatable :: coeffLocal(:,:)
      !! Plane wave coefficients
    complex(kind=dp) :: compFact(64,nAtomTypes)
      !! Complex "phase" factor
    complex(kind=dp), allocatable :: phaseExp(:,:)
      !! Complex phase exponential

    character(len=256) :: fileName
      !! Full WAVECAR file name including path

    
    if(indexInBgrp == 0) then
      !! Have the root node in each band group open the WAVECAR file

      fileName = trim(VASPDir)//'/WAVECAR'

      open(unit=wavecarUnit, file=fileName, access='direct', recl=reclenWav, iostat=ierr, status='old', SHARED)
      if (ierr .ne. 0) write(*,*) 'open error - iostat =', ierr

    endif


    if(ionode) write(*,'("Max number of k-points in pools = ", i4)') nkPerPool


    do ikLocal = 1, nkPerPool
      nGkVecsLocal_ik = nGkVecsLocal(ikLocal)

      allocate(phaseExp(nGkVecsLocal_ik, nAtoms))
      allocate(realProjWoPhase(nGkVecsLocal_ik, 64, nAtomTypes))
      allocate(coeffLocal(nGkVecsLocal_ik, ibStart_bgrp:ibEnd_bgrp))
      allocate(igkSort2OrigLocal_ik(nGkVecsLocal_ik))

      igkSort2OrigLocal_ik = igkSort2OrigLocal(1:nGkVecsLocal_ik,ikLocal)

      ikGlobal = ikLocal+ikStart_pool-1
        !! Get the global `ik` index from the local one

      nPWs1k = nPWs1kGlobal(ikGlobal)

      if(ionode) &
        write(*, '("   k-point ",i4," in pool: [ ] Phase  [ ] Real(projector)  [ ] Write projectors")') ikLocal
      call cpu_time(t1)

      !> Calculate the projectors and phase only once for each k-point
      !> because they are not dependent on spin. Calculate in one band 
      !> group and broadcast to the others because doesn't depend on 
      !> band index.
      if(myBgrpId == 0) call calculatePhase(nAtoms, nGkVecsLocal_ik, nGVecsGlobal, gVecMillerIndicesGlobalSort, igkSort2OrigLocal_ik, &
            atomPositionsDir, phaseExp)

      call MPI_BCAST(phaseExp, size(phaseExp), MPI_DOUBLE_COMPLEX, 0, interBgrpComm, ierr)


      call cpu_time(t2)
      if(ionode) &
        write(*, '("   k-point ",i4," in pool: [X] Phase  [ ] Real(projector)  [ ] Write projectors (",f7.2," secs)")') &
              ikLocal, t2-t1
      call cpu_time(t1)


      if(myBgrpId == 0) call calculateRealProjWoPhase(ikLocal, nAtomTypes, nGkVecsLocal_ik, nKPoints, gVecMillerIndicesGlobalSort, &
            igkSort2OrigLocal_ik, kPosition, omega, recipLattVec, gammaOnly, pot, realProjWoPhase, compFact)

      call MPI_BCAST(realProjWoPhase, size(realProjWoPhase), MPI_DOUBLE_PRECISION, 0, interBgrpComm, ierr)
      call MPI_BCAST(compFact, size(compFact), MPI_DOUBLE_COMPLEX, 0, interBgrpComm, ierr)


      call cpu_time(t2)
      if(ionode) &
        write(*, '("   k-point ",i4," in pool: [X] Phase  [X] Real(projector)  [ ] Write projectors (",f7.2," secs)")') &
              ikLocal, t2-t1
      call cpu_time(t1)


      if(myBgrpId == 0) call writeProjectors(ikLocal, nAtoms, iType, nAtomTypes, nAtomsEachType, nGkVecsLocal_ik, nPWs1k, realProjWoPhase, &
                compFact, phaseExp, exportDir, pot, nProj)

      call MPI_BCAST(nProj, 1, MPI_INTEGER, root, intraPoolComm, ierr)


      call cpu_time(t2)
      if(ionode) &
        write(*, '("   k-point ",i4," in pool: [X] Phase  [X] Real(projector)  [X] Write projectors (",f7.2," secs)")') &
              ikLocal, t2-t1
      call cpu_time(t1)


      do isp = 1, nSpins

        if(loopSpins .or. isp == ispSelect) then
          isk = ikGlobal + (isp - 1)*nKPoints
            !! Define index to combine k-point and spin

          irec = 2 + isk + (isk - 1)*nBands
            ! Have all processes increment the record number so
            ! they know where they are supposed to access the WAVECAR
            ! once/if they are the I/O node

          if(ionode) &
            write(*, '("      k-point ",i4," in pool, spin ",i1,": [ ] Wavefunctions  [ ] Projections")') ikLocal, isp
          call cpu_time(t1)

          call readAndWriteWavefunction(ikLocal, isp, nGkVecsLocal_ik, nPWs1k, exportDir, irec, coeffLocal)


          call cpu_time(t2)
          if(ionode) &
            write(*, '("      k-point ",i4," in pool, spin ",i1,": [X] Wavefunctions  [ ] Projections (",f7.2," secs)")') &
                  ikLocal, isp, t2-t1
          call cpu_time(t1)


          call getAndWriteProjections(ikGlobal, isp, nAtoms, nAtomTypes, nAtomsEachType, nGkVecsLocal_ik, nProj, realProjWoPhase, compFact, & 
                    phaseExp, coeffLocal, exportDir, pot)


          call cpu_time(t2)
          if(ionode) &
            write(*, '("      k-point ",i4," in pool, spin ",i1,": [X] Wavefunctions  [X] Projections (",f7.2," secs)")') &
                  ikLocal, isp, t2-t1
          call cpu_time(t1)


        endif

      enddo

      deallocate(phaseExp, realProjWoPhase, coeffLocal, igkSort2OrigLocal_ik)

      if(indexInPool == 0) write(*,'("k-point ", i4," complete!!")') ikGlobal
    enddo

    if(indexInBgrp == 0) close(wavecarUnit)

    return
  end subroutine projAndWav

!----------------------------------------------------------------------------
  subroutine calculatePhase(nAtoms, nGkVecsLocal_ik, nGVecsGlobal, gVecMillerIndicesGlobalSort, igkSort2OrigLocal_ik, &
                atomPositionsDir, phaseExp)

    implicit none

    ! Input variables: 
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: gVecMillerIndicesGlobalSort(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(in) :: igkSort2OrigLocal_ik(nGkVecsLocal_ik)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order for a 
      !! given k-point

    real(kind=dp), intent(in) :: atomPositionsDir(3,nAtoms)
      !! Atom positions

    ! Output variables:
    complex(kind=dp), intent(out) :: phaseExp(nGkVecsLocal_ik,nAtoms)

    ! Local variables:
    integer :: ia, ipw
      !! Loop indices

    real(kind=dp) :: atomPosDir(3)
      !! Direct coordinates for current atom

    complex(kind=dp) :: expArg
      !! Argument for phase exponential
    complex(kind=dp) :: itwopi = (0._dp, 1._dp)*twopi
      !! Complex phase exponential
    
    do ia = 1, nAtoms

      atomPosDir = atomPositionsDir(:,ia)
        !! Store positions locally so don't have to access 
        !! array every loop over plane waves

      do ipw = 1, nGkVecsLocal_ik

        expArg = itwopi*dot_product(atomPosDir, gVecMillerIndicesGlobalSort(:,igkSort2OrigLocal_ik(ipw)))
          !! \(2\pi i (\mathbf{G} \cdot \mathbf{r})\)

        phaseExp(ipw, ia) = exp(expArg)

      enddo
    enddo

    return
  end subroutine calculatePhase

!----------------------------------------------------------------------------
  subroutine calculateRealProjWoPhase(ik, nAtomTypes, nGkVecsLocal_ik, nKPoints, gVecMillerIndicesGlobalSort, igkSort2OrigLocal_ik, &
      kPosition, omega, recipLattVec, gammaOnly, pot, realProjWoPhase, compFact)
    implicit none

    ! Input variables:
    integer, intent(in) :: ik
      !! Current k-point 
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: gVecMillerIndicesGlobalSort(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(in) :: igkSort2OrigLocal_ik(nGkVecsLocal_ik)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order for a
      !! given k-point

    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space
    real(kind=dp), intent(in) :: omega
      !! Volume of unit cell
    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors

    logical, intent(in) :: gammaOnly
      !! If the gamma only VASP code is used

    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR

    ! Output variables:
    real(kind=dp), intent(out) :: realProjWoPhase(nGkVecsLocal_ik,64,nAtomTypes)
      !! Real projectors without phase

    complex(kind=dp), intent(out) :: compFact(64,nAtomTypes)
      !! Complex "phase" factor

    ! Local variables:
    integer :: angMom
      !! Angular momentum of projectors
    integer :: ilm
      !! Index to track going over l and m
    integer :: ilmBase
      !! Starting index over l and m based 
      !! on the current angular momentum
    integer :: ilmMax
      !! Max index over l and m based 
      !! on the current angular momentum
    integer :: imMax
      !! Max index of magnetic quantum number;
      !! loop from 0 to `imMax=2*angMom` because
      !! \(m_l\) can go from \(-l, \dots, l \)
    integer :: YDimL
      !! L dimension of spherical harmonics;
      !! max l quantum number across all
      !! pseudopotentials
    integer :: YDimLM
      !! Total number of lm combinations
    integer :: iT, ip, im, ipw
      !! Loop index
      
    real(kind=dp) :: gkMod(nGkVecsLocal_ik)
      !! \(|G+k|^2\)
    real(kind=dp) :: gkUnit(3,nGkVecsLocal_ik)
      !! \( (G+k)/|G+k| \)
    real(kind=dp) :: multFact(nGkVecsLocal_ik)
      !! Multiplicative factor for the pseudopotential;
      !! only used in the Gamma-only version
    real(kind=dp), allocatable :: pseudoV(:)
      !! Pseudopotential
    real(kind=dp), allocatable :: Ylm(:,:)
      !! Spherical harmonics


    call generateGridTable(nGkVecsLocal_ik, nKPoints, gVecMillerIndicesGlobalSort, igkSort2OrigLocal_ik, ik, kPosition, &
          recipLattVec, gammaOnly, gkMod, gkUnit, multFact)

    YDimL = maxL(nAtomTypes, pot)
      !! Get the L dimension for the spherical harmonics by
      !! finding the max l quantum number across all pseudopotentials

    YDimLM = (YDimL + 1)**2
      !! Calculate the total number of lm pairs

    allocate(Ylm(nGkVecsLocal_ik, YDimLM))

    call getYlm(nGkVecsLocal_ik, YDimL, YDimLM, gkUnit, Ylm)

    do iT = 1, nAtomTypes
      ilm = 1

      do ip = 1, pot(iT)%nChannels

        call getPseudoV(ip, nGkVecsLocal_ik, gkMod, multFact, omega, pot(iT), pseudoV)

        angMom = pot(iT)%angMom(ip)
        imMax = 2*angMom
        ilmBase = angMom**2 + 1

        !> Store complex phase factor
        ilmMax = ilm + imMax
        if(angMom == 0) then
          compFact(ilm:ilmMax,iT) = 1._dp
        else if(angMom == 1) then
          compFact(ilm:ilmMax,iT) = (0._dp, 1._dp)
        else if(angMom == 2) then
          compFact(ilm:ilmMax,iT) = -1._dp
        else if(angMom == 3) then
          compFact(ilm:ilmMax,iT) = (0._dp, -1._dp)
        endif

        do im = 0, imMax
          
          do ipw = 1, nGkVecsLocal_ik

            realProjWoPhase(ipw,ilm+im,iT) = pseudoV(ipw)*Ylm(ipw,ilmBase+im)
              !! @note
              !!  This code does not work with spin spirals! For that to work, would need 
              !!  an additional index at the end of the array for `ISPINOR`.
              !! @endnote
              !!
              !! @todo Add test to kill job if `NONL_S%LSPIRAL = .TRUE.` @endtodo
              !!
              !! @note
              !!  `realProjWoPhase` corresponds to `QPROJ`, but it only stores as much as needed 
              !!  for our application.
              !! @endnote
              !!
              !! @note
              !!  `QPROJ` is accessed as `QPROJ(ipw,ilm,iT,ik,1)`, where `ipw` is over the number
              !!  of plane waves at a specfic k-point, `ilm` goes from 1 to `WDES%LMMAX(iT)` and
              !!  `iT` is the atom-type index.
              !! @endnote
              !!
              !! @note
              !!  At the end of the subroutine `STRENL` in `nonl.F` that calculates the forces,
              !!  `SPHER` is called with `IZERO=1` along with the comment "relalculate the 
              !!  projection operators (the array was used as a workspace)." `SPHER` is what is
              !!  used to calculate `QPROJ`.
              !!
              !!  Based on this comment, I am going to assume `IZERO = 1`, which means that
              !!  `realProjWoPhase` is calculated directly rather than being added to what was
              !!  previously stored in the array, as is done in the `SPHER` subroutine.
              !! @endnote
          enddo
        enddo
        
        deallocate(pseudoV)

        ilm = ilm + imMax + 1

      enddo

      if(ilm - 1 /= pot(iT)%lmmax) call exitError('calculateRealProjWoPhase', 'LMMAX is wrong', 1)

    enddo

    deallocate(Ylm)

    return
  end subroutine calculateRealProjWoPhase

!----------------------------------------------------------------------------
  subroutine generateGridTable(nGkVecsLocal_ik, nKPoints, gVecMillerIndicesGlobalSort, igkSort2OrigLocal_ik, ik, kPosition, &
        recipLattVec, gammaOnly, gkMod, gkUnit, multFact)
    implicit none

    ! Input variables:
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    !integer, intent(in) :: nkPerPool
      ! Number of k-points in each pool
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: gVecMillerIndicesGlobalSort(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(in) :: igkSort2OrigLocal_ik(nGkVecsLocal_ik)
      !! Indices of \(G+k\) vectors in just this pool
      !! and for local PWs in the original order for a
      !! given k-point
    integer, intent(in) :: ik
      !! Current k-point 

    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space
    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors

    logical, intent(in) :: gammaOnly
      !! If the gamma only VASP code is used

    ! Output variables:
    real(kind=dp), intent(out) :: gkMod(nGkVecsLocal_ik)
      !! \(|G+k|^2\)
    real(kind=dp), intent(out) :: gkUnit(3,nGkVecsLocal_ik)
      !! \( (G+k)/|G+k| \)
    real(kind=dp), intent(out) :: multFact(nGkVecsLocal_ik)
      !! Multiplicative factor for the pseudopotential;
      ! only used in the Gamma-only version

    ! Local variables:
    integer :: gVec(3)
      !! Local storage of this G-vector
    integer :: ipw
      !! Loop indices

    real(kind=dp) :: gkCart(3)
      !! \(G+k\) in Cartesian coordinates for only
      !! vectors that satisfy the cutoff
    real(kind=dp) :: gkDir(3)
      !! \(G+k\) in direct coordinates for only
      !! vectors that satisfy the cutoff


    multFact(:) = 1._dp
      !! Initialize the multiplicative factor to 1

    do ipw = 1, nGkVecsLocal_ik

      !N1 = MOD(WDES%IGX(ipw,ik) + fftGridSize(1), fftGridSize(1)) + 1
      !N2 = MOD(WDES%IGY(ipw,ik) + fftGridSize(2), fftGridSize(2)) + 1
      !N3 = MOD(WDES%IGZ(ipw,ik) + fftGridSize(3), fftGridSize(3)) + 1

      !G1 = (GRID%LPCTX(N1) + kPosition(1,ik))
      !G2 = (GRID%LPCTY(N2) + kPosition(2,ik))
      !G3 = (GRID%LPCTZ(N3) + kPosition(3,ik))
        ! @note
        !  The original code from VASP is left commented out above. 
        !  `GRID%LPCT*` corresponds to our `gVecMillerIndicesGlobalOrig`
        !  variable that holds the unsorted G-vectors in Miller indices. 
        !  `WDES%IGX/IGY/IGZ` holds only the G-vectors s.t. \(|G+k| <\)
        !  cutoff. Their values do not seem to be sorted, so I use 
        !  `igkSort2OrigGlobal` to recreate the unsorted order from
        !  the sorted array.
        ! @endnote


      gVec(:) = gVecMillerIndicesGlobalSort(:,igkSort2OrigLocal_ik(ipw))
      gkDir(:) = gVec(:) + kPosition(:,ik+ikStart_pool-1)
        ! I belive this recreates the purpose of the original lines 
        ! above by getting \(G+k\) in direct coordinates for only
        ! the \(G+k\) combinations that satisfy the cutoff

      !IF (ASSOCIATED(NONL_S%VKPT_SHIFT)) THEN
        ! @note 
        !  `NONL_S%VKPT_SHIFT` is only set in `us.F::SETDIJ_AVEC_`.
        !  That subroutine is only called in `nmr.F::SETDIJ_AVEC`, but
        !  that call is only reached if `ORBITALMAG = .TRUE.`. This value
        !  can be set in the `INCAR` file, but there is no wiki entry,
        !  so it looks like a legacy option. Not sure how it relates
        !  to other magnetic switches like `MAGMOM`. However, `ORBITALMAG`
        !  is written to the `vasprun.xml` file, so will just test
        !  to make sure that `ORBITALMAG = .FALSE.` and throw an error 
        !  if not.
        ! @endnote

      if(gammaOnly .and. (gVec(1) /= 0 .or. gVec(2) /= 0 .or. gVec(3) /= 0)) multFact(ipw) = sqrt(2._dp)

      gkCart = matmul(recipLattVec, gkDir)
        ! VASP has a factor of `twopi` here, but I removed
        ! it because the vectors in `reconstructFFTGrid` 
        ! do not have that factor, but they result in the
        ! number of \(G+k\) vectors for each k-point matching
        ! the value input from the WAVECAR.
        !! @note
        !!  There was originally a subtraction within the parentheses of `QX`/`QY`/`QZ`
        !!  representing the spin spiral propagation vector. It was removed here because 
        !!  we do not consider spin spirals.
        !! @endnote


      gkMod(ipw) = max(sqrt(dot_product(gkCart,gkCart)), 1e-10_dp)
        !! * Get magnitude of G+k vector 

      gkUnit(:,ipw)  = gkCart(:)/gkMod(ipw)
        !! * Calculate unit vector in direction of \(G+k\)

      !IF (PRESENT(DK)) THEN
        !! @note
        !!  At the end of the subroutine `STRENL` in `nonl.F` that calculates the forces,
        !!  `SPHER` is called without `DK`  along with the comment "relalculate the 
        !!  projection operators (the array was used as a workspace)." `SPHER` is what is
        !!  used to calculate the real projectors without the complex phase.
        !!
        !!  Based on this comment, I am going to assume that `DK` isn't present, which means
        !!  that this section in the original `SPHER` subroutine is skipped. 
        !! @endnote
    enddo
  
    return
  end subroutine generateGridTable

!----------------------------------------------------------------------------
  function maxL(nAtomTypes, pot)
    !! Get the maximum L quantum number across all
    !! pseudopotentials

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms

    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR

    ! Output variables:
    integer :: maxL
      !! The maximum L quantum number across all 
      !! pseudopotentials

    ! Local variables
    integer :: maxLTmp
      !! Max L in all channels of single atom type
    integer :: iT, ip
      !! Loop indices

    maxL = 0

    do iT = 1, nAtomTypes
      maxLTmp = 0

      do ip = 1, pot(iT)%nChannels
            
        maxLTmp = max(pot(iT)%angMom(ip), maxLTmp)

      enddo

      maxL = max(maxL, maxLTmp)
      
    enddo

  end function maxL

!----------------------------------------------------------------------------
  subroutine getYlm(nGkVecsLocal_ik, YDimL, YDimLM, gkUnit, Ylm)
    implicit none

    ! Input variables:
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    integer, intent(in) :: YDimL
      !! L dimension of spherical harmonics;
      !! max l quantum number across all
      !! pseudopotentials
    integer, intent(in) :: YDimLM
      !! Total number of lm combinations

    real(kind=dp), intent(in) :: gkUnit(3,nGkVecsLocal_ik)
      !! \( (G+k)/|G+k| \)

    ! Output variables:
    real(kind=dp), intent(out) :: Ylm(nGkVecsLocal_ik,YDimLM)
      !! Spherical harmonics

    ! Local variables:
    real(kind=dp) :: multFact
      !! Factor that is multiplied in front
      !! of all spherical harmonics
    real(kind=dp) :: multFactTmp
      !! Multiplication factor for a specific
      !! calculation


    if(YDimL < 0) return
      !! Return if there is no angular momentum dimension.
      !! This shouldn't happen, but a check just in case

    Ylm(:,:) = 0._dp
      !! Initialize all spherical harmonics to zero

    multFact = 1/(2._dp*sqrt(pi))
      !! Set factor that is in front of all spherical
      !! harmonics

    Ylm(:,1) = multFact
      !! Directly calculate L=0 case

    if(YDimL < 1) return
      !! Return if the max L quantum number is 0


    !> Directly calculate L=1 case
    multFactTmp = multFact*sqrt(3._dp)

    Ylm(:,2)  = multFactTmp*gkUnit(2,:)
    Ylm(:,3)  = multFactTmp*gkUnit(3,:)
    Ylm(:,4)  = multFactTmp*gkUnit(1,:)

    if(YDimL < 2) return
      !! Return if the max L quantum number is 1


    !> Directly calculate L=2 case
    multFactTmp = multFact*sqrt(15._dp)

    Ylm(:,5)= multFactTmp*gkUnit(1,:)*gkUnit(2,:)
    Ylm(:,6)= multFactTmp*gkUnit(2,:)*gkUnit(3,:)
    Ylm(:,7)= (multFact*sqrt(5._dp)/2._dp)*(3._dp*gkUnit(3,:)*gkUnit(3,:) - 1)
    Ylm(:,8)= multFactTmp*gkUnit(1,:)*gkUnit(3,:)
    Ylm(:,9)= (multFactTmp/2._dp)*(gkUnit(1,:)*gkUnit(1,:) - gkUnit(2,:)*gkUnit(2,:))

    if(YDimL < 3) return
      !! Return if the max L quantum number is 2


    call exitError('getYlm', &
        '*** error - expected YDimL < 3', 1)
      !> @note
      !>  This code only considers up to d electrons! The spherical
      !>  harmonics are much more complicated to calculate past that 
      !>  point, but we have no use for it right now, so I am just 
      !>  going to skip it. 
      !> @endnote


    return
  end subroutine getYlm

!----------------------------------------------------------------------------
  subroutine getPseudoV(ip, nGkVecsLocal_ik, gkMod, multFact, omega, pot, pseudoV)
    implicit none

    ! Input variables:
    integer, intent(in) :: ip
      !! Channel index
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point

    real(kind=dp), intent(in) :: gkMod(nGkVecsLocal_ik)
      !! \(|G+k|^2\)
    real(kind=dp), intent(in) :: multFact(nGkVecsLocal_ik)
      !! Multiplicative factor for the pseudopotential;
      !! only used in the Gamma-only version
    real(kind=dp), intent(in) :: omega
      !! Volume of unit cell

    type (potcar) :: pot
      !! Holds all information needed from POTCAR
      !! for the specific atom type considered

    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: pseudoV(:)
      !! Pseudopotential

    ! Local variables:
    integer :: iPsGr
      !! Index on pseudopotential grid
    integer :: ipw
      !! Loop index

    real(kind=dp) :: a_ipw, b_ipw, c_ipw, d_ipw
      !! Cubic spline coefficients for recreating
      !! pseudopotential
    real(kind=dp) :: GkLenToPseudoGrid
      !! Factor to scale from \(G+k\) length scale
      !! to non-linear grid of size `nonlPseudoGridSize`
    real(kind=dp) :: divSqrtOmega
      !! 1/sqrt(omega) for multiplying pseudopotential
    real(kind=dp) :: pseudoGridLoc
      !! Location of \(G+k\) vector on pseudopotential 
      !! grid, scaled by the \(|G+k|\)
    real(kind=dp) :: rem
      !! Decimal part of `pseudoGridLoc`, used for recreating
      !! pseudopotential from cubic spline interpolation
    real(kind=dp) :: rp1, rp2, rp3, rp4
      !! Compressed recipocal-space projectors as read from
      !! the POTCAR file


    allocate(pseudoV(nGkVecsLocal_ik))

    divSqrtOmega = 1/sqrt(omega/angToBohr**3)

    GkLenToPseudoGrid = nonlPseudoGridSize/pot%maxGkNonlPs
      !! * Define a scale factor for the argument based on the
      !!   length of the G-vector. Convert from continous G-vector
      !!   length scale to discrete scale of size `nonlPseudoGridSize`.

    do ipw = 1, nGkVecsLocal_ik

      pseudoGridLoc = gkMod(ipw)*GkLenToPseudoGrid + 1
        !! * Get a location of this \(G+k\) vector scaled to
        !!   the size of the non-linear pseudopotential grid.
        !!   This value is real.

      iPsGr = int(pseudoGridLoc)
        !! * Get the integer part of the location, which will be 
        !!   used to index the reciprocal-space projectors as read
        !!   in from the POTCAR file

      rem = mod(pseudoGridLoc,1.0_dp)
        !! * Get the remainder, which is used for recreating the
        !!   pseudopotential from a cubic spline interpolation

      pseudoV(ipw) = 0._dp

      !IF (ASSOCIATED(P(NT)%PSPNL_SPLINE)) THEN
        !! @note
        !!  Default `NLSPLINE = .FALSE.`, which is recommended except for specific
        !!  applications not relevant to our purposes, so this section in the original
        !!  `SPHER` subroutine is skipped.
        !! @endnote

      rp1 = pot%recipProj(ip, iPsGr-1)
      rp2 = pot%recipProj(ip, iPsGr)
      rp3 = pot%recipProj(ip, iPsGr+1)
      rp4 = pot%recipProj(ip, iPsGr+2)

      a_ipw = rp2
      b_ipw = (6*rp3 - 2*rp1 - 3*rp2 - rp4)/6._dp
      c_ipw = (rp1 + rp3 - 2*rp2)/2._dp
      d_ipw = (rp4 - rp1 + 3*(rp2 - rp3))/6._dp
        !! * Decode the spline coefficients from the compressed reciprocal
        !!   projectors read from the POTCAR file

      pseudoV(ipw) = (a_ipw + rem*(b_ipw + rem*(c_ipw + rem*d_ipw)))*divSqrtOmega*multFact(ipw)
        !! * Recreate full pseudopotential from cubic spline coefficients:
        !!   \( \text{pseudo} = a_i + dx\cdot b_i + dx^2\cdot c_i + dx^3\cdot d_i \)
        !!   where the \(i\) index is the plane-wave index, and \(dx\) is the decimal
        !!   part of the pseudopotential-grid location

      !IF (VPS(IND) /= 0._dp .AND. PRESENT(DK)) THEN
        !! @note
        !!  At the end of the subroutine `STRENL` in `nonl.F` that calculates the forces,
        !!  `SPHER` is called without `DK`  along with the comment "relalculate the 
        !!  projection operators (the array was used as a workspace)." `SPHER` is what is
        !!  used to calculate the real projectors without the complex phase.
        !!
        !!  Based on this comment, I am going to assume that `DK` isn't present, which means
        !!  that this section in the original `SPHER` subroutine is skipped. 
        !! @endnote
    enddo

    return
  end subroutine getPseudoV

!----------------------------------------------------------------------------
  subroutine writeProjectors(ik, nAtoms, iType, nAtomTypes, nAtomsEachType, nGkVecsLocal_ik, nPWs1k, realProjWoPhase, compFact, &
      phaseExp, exportDir, pot, nProj)

    use miscUtilities, only: int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: ik
      !! Current k-point
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms
    integer, intent(in) :: nAtomsEachType(nAtomTypes)
      !! Number of atoms of each type
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    integer, intent(in) :: nPWs1k
      !! Input number of plane waves for the given k-point

    real(kind=dp), intent(in) :: realProjWoPhase(nGkVecsLocal_ik,64,nAtomTypes)
      !! Real projectors without phase

    complex(kind=dp), intent(in) :: compFact(64,nAtomTypes)
      !! Complex "phase" factor
    complex(kind=dp), intent(in) :: phaseExp(nGkVecsLocal_ik,nAtoms)
      !! Exponential phase factor

    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export

    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR

    ! Output variables:
    integer, intent(out) :: nProj
      !! Number of projectors across all atom types

    ! Local variables:
    integer :: reclen
      !! Record length for projectors file
    integer :: iT, ia, ilm, igkLocal, igkGlobal, ikGlobal, ipr
      !! Loop indices

    real(kind=dp) :: nProj_real, nPWs1k_real
      !! Real versions of output variables for better 
      !! compatibility of binary file

    complex(kind=dp), allocatable :: beta(:)
      !! Projectors for single plane wave


    !> Set up base file name
    ikGlobal = ik+ikStart_pool-1


    nProj = 0
    do iT = 1, nAtomTypes
      !! Calculate the total number of projectors across all
      !! atom types

      nProj = nProj + pot(iT)%lmmax*nAtomsEachType(iT)

    enddo


    allocate(beta(nProj))

    inquire(iolength=reclen) beta(:)
      !! Get the record length needed to write a double complex
      !! array of size `nProj`


    open(63, file=trim(exportDir)//'/projectors.'//trim(int2str(ikGlobal)), access='direct', form='unformatted', recl=reclen)
      !! Open output file with direct access


    if(indexInBgrp == 0) then
      !! Have root process write header values

      nProj_real = nProj
      nPWs1k_real = nPWs1k
      write(63,rec=1) nProj_real, nPWs1k_real
        !! Write out the number of projectors and number of
        !! \(G+k\) vectors at this k-point below the energy
        !! cutoff

    endif


    do igkLocal = 1, nGkVecsLocal_ik

      ipr = 0
      do ia = 1, nAtoms

        iT = iType(ia)
        do ilm = 1, pot(iT)%lmmax

          ipr = ipr + 1
          beta(ipr) = conjg(realProjWoPhase(igkLocal,ilm,iT)*phaseExp(igkLocal,ia)*compFact(ilm,iT))
            !! @note
            !!    The projectors are stored as \(\langle\beta|\), so need to take the complex conjugate
            !!    to output \(|\beta\rangle.
            !! @endnote
            !! @note
            !!    The projectors should have units inverse to those of the coefficients. That was
            !!    previously listed as (a.u.)^(-3/2), but the `TME` code seems to expect both the
            !!    projectors and the wave function coefficients to be unitless, so there should be
            !!    no unit conversion here.
            !! @endnote
            !! @note
            !!    `NONL_S%LSPIRAL = .FALSE.`, so spin spirals are not calculated, which makes
            !!    `NONL_S%QPROJ` spin-independent. This is why there is no spin index on `realProjWoPhase`.
            !! @endnote
        enddo
      enddo


      igkGlobal = igkLocal+iGkStart_pool(ik)-1

      write(63,rec=igkGlobal+1) (beta(ipr), ipr=1,nProj)

    enddo

    close(63)

    return
  end subroutine writeProjectors

!----------------------------------------------------------------------------
  subroutine readAndWriteWavefunction(ik, isp, nGkVecsLocal_ik, nPWs1k, exportDir, irec, coeffLocal)
    !! For each spin and k-point, read and write the plane
    !! wave coefficients for each band
    !!
    !! <h2>Walkthrough</h2>
    !!

    use miscUtilities, only: int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: ik
      !! Current k-point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    integer, intent(in) :: nPWs1k
      !! Input number of plane waves for the given k-point
      
    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export

    ! Output variables:
    integer, intent(inout) :: irec

    complex*8, intent(out) :: coeffLocal(nGkVecsLocal_ik, ibStart_bgrp:ibEnd_bgrp)
      !! Plane wave coefficients

    ! Local variables:
    integer :: reclen
      !! Record length for projectors file
    integer, allocatable :: sendCount(:)
      !! Number of items to send/recieve to/from each process
    integer, allocatable :: displacement(:)
      !! Offset from beginning of array for
      !! scattering/gathering coefficients 
    integer :: ib, ipw, ikGlobal
      !! Loop indices

    complex*8, allocatable :: coeff(:)
      !! Plane wave coefficients

    character(len=300) :: fNameExport
      !! Output file name


    ikGlobal = ik+ikStart_pool-1


    if(indexInBgrp == 0) then

      allocate(coeff(nPWs1k))

      fNameExport = trim(exportDir)//"/wfc."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)) 

      inquire(iolength=reclen) coeff(:)
        !! Get the record length needed to write a double complex
        !! array of length nPWs1k

      open(63, file=trim(fNameExport), access='direct', form='unformatted', recl=reclen)
        !! Open output file with direct access

    else

      allocate(coeff(1))

    endif


    allocate(sendCount(nProcPerBgrp))
    sendCount = 0
    sendCount(indexInBgrp+1) = nGkVecsLocal_ik
    call mpiSumIntV(sendCount, intraBgrpComm)
      !! * Put the number of G+k vectors on each process
      !!   in a single array per band group

    allocate(displacement(nProcPerBgrp))
    displacement = 0
    displacement(indexInBgrp+1) = iGkStart_pool(ik)-1
    call mpiSumIntV(displacement, intraBgrpComm)
      !! * Put the displacement from the beginning of the array
      !!   for each process in a single array per band group


    irec = irec + ibStart_bgrp - 1

    do ib = ibStart_bgrp, ibEnd_bgrp

      irec = irec + 1

      if(indexInBgrp == 0) then

        read(unit=wavecarUnit,rec=irec) (coeff(ipw), ipw=1,nPWs1k)
          ! Read in the plane wave coefficients for each band

        write(63,rec=ib) (coeff(ipw), ipw=1,nPWs1k)
          !! @note
          !!  I was trying to convert these coefficients based
          !!  on the units previously listed in the `wfc.ik` file, 
          !!  but I don't think those are accurate. Based on the 
          !!  `TME` code, it seems like these coefficients are 
          !!  actually treated as unitless, so there should be no 
          !!  unit conversion here.
          !! @endnote

      endif

      call MPI_SCATTERV(coeff(:), sendCount, displacement, MPI_COMPLEX, coeffLocal(1:nGkVecsLocal_ik,ib), nGkVecsLocal_ik, &
          MPI_COMPLEX, 0, intraBgrpComm, ierr)
      !! * For each band, scatter the coefficients across all 
      !!   of the processes in the pool

    enddo


    deallocate(coeff)
    deallocate(sendCount)
    deallocate(displacement)

    return
  end subroutine readAndWriteWavefunction

!----------------------------------------------------------------------------
  subroutine getAndWriteProjections(ik, isp, nAtoms, nAtomTypes, nAtomsEachType, nGkVecsLocal_ik, nProj, realProjWoPhase, compFact, &
          phaseExp, coeffLocal, exportDir, pot)

    use miscUtilities, only: int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: ik
      !! Current k-point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms
    integer, intent(in) :: nAtomsEachType(nAtomTypes)
      !! Number of atoms of each type
    integer, intent(in) :: nGkVecsLocal_ik
      !! Local number of G+k-vectors on this processor
      !! for a given k-point
    integer, intent(in) :: nProj
      !! Number of projectors across all atom types

    real(kind=dp), intent(in) :: realProjWoPhase(nGkVecsLocal_ik,64,nAtomTypes)
      !! Real projectors without phase

    complex(kind=dp), intent(in) :: compFact(64,nAtomTypes)
      !! Complex "phase" factor
    complex(kind=dp), intent(in) :: phaseExp(nGkVecsLocal_ik,nAtoms)

    complex*8, intent(in) :: coeffLocal(nGkVecsLocal_ik, ibStart_bgrp:ibEnd_bgrp)
      !! Plane wave coefficients
      
    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export

    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR

    ! Local variables:
    integer :: ib, iT, ia, iaBase, ilm, ipr
      !! Loop indices
    integer :: reclen
      !! Record length for projectors file

    character(len=300) :: fNameExport
      !! Output file name

    complex(kind=dp) :: projection(nProj)
      !! Projection for current atom/band/lm channel


    if(indexInBgrp == 0) then
      !! Have the root node within the band group handle I/O

      fNameExport = trim(exportDir)//"/projections."//trim(int2str(isp))//"."//trim(int2str(ik)) 

      inquire(iolength=reclen) projection(:)

      open(83, file=trim(fNameExport), access='direct', form='unformatted', recl=reclen)
        !! Open output file with direct access

    endif

    do ib = ibStart_bgrp, ibEnd_bgrp
      iaBase = 1
      ipr = 0
      
      do iT = 1, nAtomTypes
        do ia = iaBase, nAtomsEachType(iT)+iaBase-1
          do ilm = 1, pot(iT)%lmmax

            ipr = ipr + 1

            projection(ipr) = compFact(ilm,iT)*sum(realProjWoPhase(:,ilm,iT)*phaseExp(:,ia)*coeffLocal(:,ib))
              ! Calculate projection (sum over plane waves)
              ! Don't need to worry about sorting because projection
              ! has sum over plane waves.

          enddo
        enddo

        iaBase = iaBase + nAtomsEachType(iT)

      enddo

      call mpiSumComplexV(projection, intraBgrpComm)

      if(indexInBgrp == 0) write(83,rec=ib) projection(:)
    enddo

    if(indexInBgrp == 0) close(83)

    return
  end subroutine getAndWriteProjections

!----------------------------------------------------------------------------
  subroutine writeKInfo(nKPoints, nGkLessECutGlobal, nSpins, kWeight, kPosition)
    !! Calculate the highest occupied band for each k-point
    !! and write out k-point information
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nGkLessECutGlobal(nKPoints)
      !! Global number of \(G+k\) vectors with magnitude
      !! less than `wfcVecCut` for each k-point
    integer, intent(in) :: nSpins
      !! Number of spins

    real(kind=dp), intent(in) :: kWeight(nKPoints)
      !! K-point weights
    real(kind=dp), intent(in) :: kPosition(3,nKPoints)
      !! Position of k-points in reciprocal space

    ! Local variables:
    integer :: ik
      !! Loop index


    if(ionode) then

      write(mainOutFileUnit, '("# Number of spins. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nSpins

      write(mainOutFileUnit, '("# Number of K-points. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nKPoints

      write(mainOutFileUnit, '("# ik, nGkLessECutGlobal(ik), wk(ik), xk(1:3,ik). Format: ''(2i10,4ES24.15E3)''")')
      flush(mainOutFileUnit)

      if(energiesOnly) then
        do ik = 1, nKPoints
    
          write(mainOutFileUnit, '(2i10,4ES24.15E3)') ik, -1, kWeight(ik), kPosition(1:3,ik)
            ! Skip writing out grid info

        enddo
      else
        do ik = 1, nKPoints
    
          write(mainOutFileUnit, '(2i10,4ES24.15E3)') ik, nGkLessECutGlobal(ik), kWeight(ik), kPosition(1:3,ik)

        enddo
      endif

    endif

    return
  end subroutine writeKInfo

!----------------------------------------------------------------------------
  subroutine writeGridInfo(nGVecsGlobal, gVecMillerIndicesGlobalOrig, maxGIndexGlobal, exportDir)
    !! Write out grid boundaries and miller indices
    !! for just \(G+k\) combinations below cutoff energy
    !! in one file and all miller indices in another 
    !! file
    !!
    !! <h2>Walkthrough</h2>
    !!

    use miscUtilities, only: int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: nGVecsGlobal
      !! Global number of G-vectors

    integer, intent(in) :: gVecMillerIndicesGlobalOrig(3,nGVecsGlobal)
      !! Integer coefficients for G-vectors on all processors (original order)
    integer, intent(in) :: maxGIndexGlobal
      !! Maximum G-vector index among all \(G+k\)
      !! and processors

    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export

    ! Local variables:
    integer :: ig
      !! Loop index


    if (ionode) then
    
      !> * Write the global number of G-vectors, the maximum
      !>   G-vector index, and the max/min miller indices
      write(mainOutFileUnit, '("# Number of G-vectors. Format: ''(i10)''")')
      if(energiesOnly) then
        write(mainOutFileUnit,'(i10)') -1
      else
        write(mainOutFileUnit, '(i10)') nGVecsGlobal
      endif
    
      write(mainOutFileUnit, '("# Number of PW-vectors. Format: ''(i10)''")')
      if(energiesOnly) then
        write(mainOutFileUnit,'(i10)') -1
      else
        write(mainOutFileUnit, '(i10)') maxGIndexGlobal
      endif
    
      write(mainOutFileUnit, '("# Number of min - max values of fft grid in x, y and z axis. Format: ''(6i10)''")')
      if(energiesOnly) then
        write(mainOutFileUnit,'(6i10)') -1, -1, -1, -1, -1, -1
      else
        write(mainOutFileUnit, '(6i10)') minval(gVecMillerIndicesGlobalOrig(1,1:nGVecsGlobal)), maxval(gVecMillerIndicesGlobalOrig(1,1:nGVecsGlobal)), &
                          minval(gVecMillerIndicesGlobalOrig(2,1:nGVecsGlobal)), maxval(gVecMillerIndicesGlobalOrig(2,1:nGVecsGlobal)), &
                          minval(gVecMillerIndicesGlobalOrig(3,1:nGVecsGlobal)), maxval(gVecMillerIndicesGlobalOrig(3,1:nGVecsGlobal))
      endif

    endif
    
    if(.not. energiesOnly .and. ionode) then

      !> * Output all miller indices in `mgrid` file
      open(72, file=trim(exportDir)//"/mgrid")
      write(72, '("# Full G-vectors grid")')
      write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
    
      do ig = 1, nGVecsGlobal
        write(72, '(4i10)') ig, gVecMillerIndicesGlobalOrig(1:3,ig)
      enddo
    
      close(72)

    endif

    return
  end subroutine writeGridInfo


!----------------------------------------------------------------------------
  subroutine writeCellInfo(iType, nAtoms, nBands, nAtomTypes, realLattVec, recipLattVec, atomPositionsDir)
    !! Write out the real- and reciprocal-space lattice vectors, 
    !! the number of atoms, the number of types of atoms, the
    !! final atom positions, number of bands, and number of spins

    use generalComputations, only: direct2cart

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtoms
      !! Number of atoms
    integer, intent(in) :: iType(nAtoms)
      !! Atom type index
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms

    real(kind=dp), intent(in) :: realLattVec(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: recipLattVec(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: atomPositionsDir(3,nAtoms)
      !! Atom positions

    ! Local variables:
    real(kind=dp) :: atomPositionsCart(3,nAtoms)
      !! Position of given atom in cartesian coordinates

    integer :: ia
      !! Loop indices


    if (ionode) then
    
      write(mainOutFileUnit, '("# Cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
      write(mainOutFileUnit, '("# a1 ",3ES24.15E3)') realLattVec(:,1)
      write(mainOutFileUnit, '("# a2 ",3ES24.15E3)') realLattVec(:,2)
      write(mainOutFileUnit, '("# a3 ",3ES24.15E3)') realLattVec(:,3)
    
      write(mainOutFileUnit, '("# Reciprocal cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
      write(mainOutFileUnit, '("# b1 ",3ES24.15E3)') recipLattVec(:,1)
      write(mainOutFileUnit, '("# b2 ",3ES24.15E3)') recipLattVec(:,2)
      write(mainOutFileUnit, '("# b3 ",3ES24.15E3)') recipLattVec(:,3)
    
      write(mainOutFileUnit, '("# Number of Atoms. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nAtoms
    
      write(mainOutFileUnit, '("# Number of Types. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nAtomTypes
    
      write(mainOutFileUnit, '("# Atoms type, position(1:3) (a.u.). Format: ''(i10,3ES24.15E3)''")')

      atomPositionsCart = direct2cart(nAtoms, atomPositionsDir, realLattVec)

      do ia = 1, nAtoms

        write(mainOutFileUnit,'(i10,3ES24.15E3)') iType(ia), atomPositionsCart(:,ia)

      enddo
    
      write(mainOutFileUnit, '("# Number of Bands. Format: ''(i10)''")')
      write(mainOutFileUnit, '(i10)') nBands

    endif

    return
  end subroutine writeCellInfo

!----------------------------------------------------------------------------
  subroutine writePseudoInfo(nAtomTypes, nAtomsEachType, pot)
    !! For each atom type, write out the element name,
    !! number of atoms of this type, projector info,
    !! radial grid info, and partial waves

    implicit none

    ! Input variables:
    integer, intent(in) :: nAtomTypes
      !! Number of types of atoms

    integer, intent(in) :: nAtomsEachType(nAtomTypes)
      !! Number of atoms of each type

    type (potcar) :: pot(nAtomTypes)
      !! Holds all information needed from POTCAR


    ! Output variables:


    ! Local variables:
    integer :: iT, ip, ir
      !! Loop index

  
    if (ionode) then

      do iT = 1, nAtomTypes
        
        write(mainOutFileUnit, '("# Element")')
        write(mainOutFileUnit, *) trim(pot(iT)%element)
        write(mainOutFileUnit, '("# Number of Atoms of this type. Format: ''(i10)''")')
        write(mainOutFileUnit, '(i10)') nAtomsEachType(iT)
        write(mainOutFileUnit, '("# Number of projectors. Format: ''(i10)''")')
        write(mainOutFileUnit, '(i10)') pot(iT)%nChannels
        
        write(mainOutFileUnit, '("# Angular momentum, index of the projectors. Format: ''(2i10)''")')
        do ip = 1, pot(iT)%nChannels

          write(mainOutFileUnit, '(2i10)') pot(iT)%angMom(ip), ip

        enddo
        
        write(mainOutFileUnit, '("# Number of channels. Format: ''(i10)''")')
        write(mainOutFileUnit, '(i10)') pot(iT)%lmmax
        
        write(mainOutFileUnit, '("# Number of radial mesh points. Format: ''(2i10)''")')
        write(mainOutFileUnit, '(2i10)') pot(iT)%nmax, pot(iT)%iRAugMax
          ! Number of points in the radial mesh, number of points inside the aug sphere
        
        write(mainOutFileUnit, '("# Radial grid, Integratable grid. Format: ''(2ES24.15E3)''")')
        do ir = 1, pot(iT)%nmax
          write(mainOutFileUnit, '(2ES24.15E3)') pot(iT)%radGrid(ir), pot(iT)%dRadGrid(ir) 
            ! Radial grid, derivative of radial grid
        enddo
        
        write(mainOutFileUnit, '("# AE, PS radial wfc for each beta function. Format: ''(2ES24.15E3)''")')
        do ip = 1, pot(iT)%nChannels
          do ir = 1, pot(iT)%nmax
            write(mainOutFileUnit, '(2ES24.15E3)') pot(iT)%wae(ip,ir), pot(iT)%wps(ip,ir)
          enddo
        enddo
      
      enddo
    
    endif

    return
  end subroutine writePseudoInfo

!----------------------------------------------------------------------------
  subroutine writeEigenvalues(ispSelect, nBands, nKPoints, nSpins, bandOccupation, eFermi, eTot, eigenE, loopSpins)
    !! Write Fermi energy and eigenvalues and occupations for each band

    use miscUtilities, only: int2str

    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nSpins
      !! Number of spins
      
    real(kind=dp), intent(in) :: bandOccupation(nSpins, nBands,nKPoints)
      !! Occupation of band
    real(kind=dp), intent(in) :: eFermi
      !! Fermi energy
    real(kind=dp), intent(in) :: eTot
      !! Total energy

    complex*16, intent(in) :: eigenE(nSpins,nKPoints,nBands)
      !! Band eigenvalues

    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    ! Local variables:
    integer :: ik, ib, isp
      !! Loop indices


    if(ionode) then
    
      write(mainOutFileUnit, '("# Fermi Energy (Hartree). Format: ''(ES24.15E3)''")')
      write(mainOutFileUnit, '(ES24.15E3)') eFermi*ryToHartree
      write(mainOutFileUnit, '("# Total Energy (Hartree). Format: ''(ES24.15E3)''")')
      write(mainOutFileUnit, '(ES24.15E3)') eTot*ryToHartree
      flush(mainOutFileUnit)
    
      do isp = 1, nSpins
        if(loopSpins .or. isp == ispSelect) then
          do ik = 1, nKPoints

            open(72, file=trim(exportDir)//"/eigenvalues."//trim(int2str(isp))//"."//trim(int2str(ik)))
      
            write(72, '("# Spin, k-point index: ",2i10, " Format: ''(a23, 2i10)''")') isp, ik
            write(72, '("# Eigenvalues (Hartree), band occupation number. Format: ''(2ES24.15E3)''")')
      
            do ib = 1, nBands
  
              write(72, '(i10, 2ES24.15E3)') ib, real(eigenE(isp,ik,ib)), bandOccupation(isp,ib,ik)

            enddo
      
            close(72)
          enddo
        endif
      enddo

    endif

    return
  end subroutine writeEigenvalues

!----------------------------------------------------------------------------
  subroutine writeGroupedEigenvalues(ispSelect, nBands, nDispkPerCoord, nKPoints, nSpins, bandOccupation, patternArr, eigenE, &
        loopSpins, pattern)

    use miscUtilities, only: int2str, hpsort_eps

    implicit none

    ! Input variables:
    integer, intent(in) :: ispSelect
      !! Selection of a single spin channel if input
      !! by the user
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nDispkPerCoord
      !! Number of displaced k-points per coordinate
    integer, intent(in) :: nKPoints
      !! Total number of k-points
    integer, intent(in) :: nSpins
      !! Number of spins
      
    real(kind=dp), intent(in) :: bandOccupation(nSpins, nBands,nKPoints)
      !! Occupation of band
    real(kind=dp), intent(in) :: patternArr(nDispkPerCoord)
      !! Displacement pattern for groups of
      !! k-points for group velocity calculations

    complex*16, intent(in) :: eigenE(nSpins,nKPoints,nBands)
      !! Band eigenvalues

    logical, intent(in) :: loopSpins
      !! Whether to loop over available spin channels;
      !! otherwise, use selected spin channel

    character(len=300), intent(in) :: pattern
      ! Character input for displacement pattern

    ! Local variables
    integer :: ikGroup, ib, isp, ikx, ix, ik
      !! Loop indices
    integer :: iPattSort(nDispkPerCoord+1)
      !! Sorted index order for patternArr from 
      !! negative to positive
    integer :: nKGroups
      !! Number of k-point groups if grouping for
      !! group velocity calculation
    integer :: nkPerGroup
      !! Number of k-points per group

    real(kind=dp) :: eigsForOutput(nBands,nDispkPerCoord+1,3)
      !! Sorted real eigenvalues in correct units
      !! to be output
    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp) :: patternArrSort(nDispkPerCoord+1)
      !! Sorted displacement pattern for groups of
      !! k-points for group velocity calculations

    character(len=300) :: formatString
      !! String to dynamically determine the output
      !! format based on `nkPerGroup`


    if(ionode) then

      nkPerGroup = nDispkPerCoord*3 + 1

      if(mod(nKPoints,nkPerGroup) /= 0) &
        call exitError('writeEigenvalues', 'groups of '//trim(int2str(nkPerGroup))//' k-points expected', 1)

      nKGroups = nKPoints/nkPerGroup

      formatString = "("//trim(int2str(nkPerGroup+1))//"ES19.10E3)"
        ! Eigenvalues for all k-points in group plus the occupation

      do ikx = 1, nDispkPerCoord+1
        iPattSort(ikx) = ikx
      enddo
      patternArrSort(1) = 0.0_dp
      patternArrSort(2:nDispkPerCoord+1) = patternArr(:)
      call hpsort_eps(nDispkPerCoord+1, patternArrSort, iPattSort, eps8)
        ! Get sorted index order for patternArr from negative to positive
    
      do isp = 1, nSpins
        if(loopSpins .or. isp == ispSelect) then
          do ikGroup = 1, nKGroups
            do ix = 1, 3
              do ikx = 1, nDispkPerCoord+1

                if(iPattSort(ikx) == 1) then
                  ! This is the zero point. We add it to the
                  ! sorting array in case the points aren't 
                  ! symmetric about zero displacement.

                  ik = (ikGroup-1)*nkPerGroup + & ! Skip all k-points from previous group
                       1                          ! Skip to base k-point

                else

                  ik = (ikGroup-1)*nkPerGroup + & ! Skip all k-points from previous group
                       1 + &                      ! Skip base k-point
                       (ix-1)*nDispkPerCoord + &  ! Skip previous displaced coordinates
                       iPattSort(ikx)-1           ! Skip to sorted value from displacement pattern minus zero-displacement point
                endif

                eigsForOutput(:,ikx,ix) = real(eigenE(isp,ik,:))*ryToHartree

              enddo

              if(ix == 1) then
                open(72, file=trim(exportDir)//"/groupedEigenvaluesX."//trim(int2str(isp))//"."//trim(int2str(ikGroup)))
              else if(ix == 2) then
                open(72, file=trim(exportDir)//"/groupedEigenvaluesY."//trim(int2str(isp))//"."//trim(int2str(ikGroup)))
              else if(ix == 3) then
                open(72, file=trim(exportDir)//"/groupedEigenvaluesZ."//trim(int2str(isp))//"."//trim(int2str(ikGroup)))
              endif
        
              write(72, '("# Spin : ",i10, " Format: ''(a9, i10)''")') isp
              write(72,'("# Input displacement pattern:")')
              write(72,'(a)') trim(pattern)

              if(ix == 1) then
                write(72, '("# Eigenvalues (Hartree) - to + for x, band occupation number Format: ''",a"''")') trim(formatString)
              else if(ix == 2) then
                write(72, '("# Eigenvalues (Hartree) - to + for y, band occupation number Format: ''",a"''")') trim(formatString)
              else if(ix == 3) then
                write(72, '("# Eigenvalues (Hartree) - to + for z, band occupation number Format: ''",a"''")') trim(formatString)
              endif
      
              do ib = 1, nBands

                write(72,formatString) &
                  (eigsForOutput(ib,ikx,ix), ikx=1,nDispkPerCoord+1), & ! k-points for this coordinate from - to +
                  bandOccupation(isp,ib,(ikGroup-1)*nkPerGroup+1)       ! band occupation from base k-point           

              enddo
      
              close(72)

            enddo
          enddo
        endif
     enddo
    endif
    
    return
  end subroutine writeGroupedEigenvalues

!----------------------------------------------------------------------------
  subroutine subroutineTemplate()
    implicit none


    return
  end subroutine subroutineTemplate

end module wfcExportVASPMod
