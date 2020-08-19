module wfcExportVASPMod

  USE wrappers,      ONLY : f_mkdir_safe

  !USE pwcom
  USE constants, ONLY : e2, rytoev, tpi, fpi
  USE cell_base, ONLY : celldm, ibrav
  USE ener, ONLY : ef
  USE wvfct, ONLY : et
  USE lsda_mod, ONLY : isk

  USE io_files,  ONLY : prefix, outdir, tmp_dir
  USE iotk_module
  use mpi
  USE mp,        ONLY: mp_max, mp_get, mp_bcast, mp_rank
  USE mp_wave, ONLY : mergewf

  implicit none

  ! Parameters:
  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision
  integer, parameter :: root = 0
    !! ID of the root node
  integer, parameter :: root_pool = 0
    !! Index of the root process within each pool
  integer, parameter :: mainout = 50
    !! Main output file unit
  integer, parameter :: potcarUnit = 86
    !! POTCAR unit for I/O
  integer, parameter :: stdout = 6
    !! Standard output unit
  integer, parameter :: wavecarUnit = 86
    !! WAVECAR unit for I/O

  real(kind = dp), parameter :: eVToRy = 0.073498618_dp
    !! Conversion factor from eV Rydberg
  real(kind = dp), parameter :: pi = 3.141592653589793_dp
    !! \(\pi\)
  real(kind = dp), parameter :: ryToHartree = 0.5_dp
    !! Conversion factor from Rydberg to Hartree
  real(kind = dp), parameter :: twoPiSquared = (2.0_dp*pi)**2
    !! This is used in place of \(2\pi/a\) which assumes that \(a=1\)


  ! Global variables not passed as arguments:
  integer :: ierr
    !! Error returned by MPI
  integer :: ikEnd
    !! Ending index for k-points in single pool 
  integer :: ikStart
    !! Starting index for k-points in single pool 
  integer :: ios
    !! Error for input/output
  integer :: indexInPool
    !! Process index within pool
  integer :: inter_pool_comm_local = 0
    !! Inter-pool communicator
    !! @todo Change back to `inter_pool_comm` once extracted from QE #end @endtodo
  integer :: intra_pool_comm_local = 0
    !! Intra-pool communicator
    !! @todo Change back to `intra_pool_comm` once extracted from QE #end @endtodo
  integer :: myid
    !! ID of this process
  integer :: myPoolId
    !! Pool index for this process
  integer :: nk_Pool
    !! Number of k-points in each pool
  integer :: npool_local = 1
    !! Number of pools for k-point parallelization
    !! @todo Change back to `npool` once extracted from QE #end @endtodo
  integer :: nproc_local
    !! Number of processes
    !! @todo Change back to `nproc` once extracted from QE #end @endtodo
  integer :: nproc_pool_local
    !! Number of processes per pool
    !! @todo Change back to `nproc_pool` once extracted from QE #end @endtodo
  integer :: world_comm_local
    !! World communicator
    !! @todo Change back to `world_comm` once extracted from QE #end @endtodo

  logical :: ionode_local
    !! If this node is the root node
    !! @todo Change back to `ionode` once extracted from QE #end @endtodo


  ! Variables that should be passed as arguments:
  real(kind=dp) :: at_local(3,3)
    !! Real space lattice vectors
    !! @todo Change back to `at` once extracted from QE #end @endtodo
  real(kind=dp) :: bg_local(3,3)
    !! Reciprocal lattice vectors
    !! @todo Change back to `bg` once extracted from QE #end @endtodo
  real(kind=dp) :: ecutwfc_local
    !! Plane wave energy cutoff in Ry
    !! @todo Change back to `ecutwfc` once extracted from QE #end @endtodo
  real(kind=dp), allocatable :: gCart_local(:,:)
    !! G-vectors in Cartesian coordinates
  real(kind=dp), allocatable :: occ(:,:)
    !! Occupation of band
  real(kind=dp) :: omega_local
    !! Volume of unit cell
    !! @todo Change back to `omega` once extracted from QE #end @endtodo
  real(kind=dp), allocatable :: tau(:,:)
    !! Atom positions
  real(kind=dp) :: tStart
    !! Start time
  real(kind=dp) :: vcut_local
    !! Energy cutoff converted to vector cutoff
  real(kind=dp), allocatable :: xk_local(:,:)
    !! Position of k-points in reciprocal space
  real(kind=dp), allocatable :: wk_local(:)
    !! Weight of k-points
  
  integer, allocatable :: ig_l2g(:)
    !! Converts local index `ig` to global index
  integer, allocatable :: igk_l2g(:,:)
    !! Local to global indices for \(G+k\) vectors 
    !! ordered by magnitude at a given k-point
  integer, allocatable :: igk_large(:,:)
    !! Index map from \(G\) to \(G+k\);
    !! indexed up to `ngm_local` which
    !! is greater than `npwx_local` and
    !! stored for each k-point
  integer, allocatable :: igwk(:,:)
    !! Indices of \(G+k\) vectors for each k-point
    !! and all processors
  integer, allocatable :: ityp(:)
    !! Atom type index
  integer, allocatable :: mill_g(:,:)
    !! Integer coefficients for G-vectors on all processors
  integer :: nb1max, nb2max, nb3max
    !! Not sure what this is??
  integer :: nat
    !! Number of atoms
  integer :: nbnd_local
    !! Total number of bands
    !! @todo Change back to `nbnd` once extracted from QE #end @endtodo
  integer, allocatable :: ngk_g(:)
    !! Global number of \(G+k\) vectors with energy
    !! less than `ecutwfc_local` for each k-point
  integer, allocatable :: ngk_local(:)
    !! Number of \(G+k\) vectors with energy
    !! less than `ecutwfc_local` for each
    !! k-point, on this processor
    !! @todo Change back to `ngk` once extracted from QE #end @endtodo
  integer :: ngk_max
    !! Maximum number of \(G+k\) combinations
  integer :: ngm_g_local
    !! Global number of G-vectors
    !! @todo Change back to `ngm_g` once extracted from QE #end @endtodo
  integer :: ngm_local
    !! Local number of G-vectors on this processor
    !! @todo Change back to `ngm` once extracted from QE #end @endtodo
  integer :: nkstot_local
    !! Total number of k-points
    !! @todo Change back to `nkstot` once extracted from QE #end @endtodo
  integer, allocatable :: nnTyp(:)
    !! Number of atoms of each type
  integer, allocatable :: nplane_g(:)
    !! Input number of plane waves for a single k-point for all processors
  integer :: npw_g
    !! Maximum G-vector index among all \(G+k\)
    !! and processors
  integer :: npwx_g
    !! Max number of \(G+k\) vectors with energy
    !! less than `ecutwfc_local` among all k-points
  integer :: npwx_local
    !! Maximum number of \(G+k\) vectors
    !! across all k-points for just this
    !! processor
    !! @todo Change back to `npwx` once extracted from QE #end @endtodo
  integer :: nsp
    !! Number of types of atoms
  integer :: nspin_local
    !! Number of spins
    !! @todo Change back to `nspin` once extracted from QE #end @endtodo
  
  character(len=256) :: exportDir
    !! Directory to be used for export
  character(len=256) :: mainOutputFile
    !! Main output file
  character(len=256) :: QEDir
    !! Directory with QE files
  character(len=256) :: VASPDir
    !! Directory with VASP files

  type pseudo
    integer :: angmom(16) = 0
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
    real(kind=dp), allocatable :: radGrid(:)
      !! Radial grid points
    real(kind=dp) :: rAugMax
      !! Maximum radius of augmentation sphere
    real(kind=dp) :: recipProj(16,100)
      !! Reciprocal-space projectors
    real(kind=dp) :: realProj(16,100)
      !! Real-space projectors
    real(kind=dp), allocatable :: wae(:,:)
      !! AE wavefunction
    real(kind=dp), allocatable :: wps(:,:)
      !! PS wavefunction

    character(len=2) :: element
  end type pseudo

  type (pseudo), allocatable :: ps(:)

  namelist /inputParams/ prefix, QEDir, VASPDir, exportDir


  contains

!----------------------------------------------------------------------------
  subroutine mpiInitialization()
    !! Generate MPI processes and communicators 
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Output variables:
    !logical, intent(out) :: ionode_local
      ! If this node is the root node
    !integer, intent(out) :: inter_pool_comm_local = 0
      ! Inter-pool communicator
    !integer, intent(out) :: intra_pool_comm_local = 0
      ! Intra-pool communicator
    !integer, intent(out) :: indexInPool
      ! Process index within pool
    !integer, intent(out) :: myid
      ! ID of this process
    !integer, intent(out) :: myPoolId
      ! Pool index for this process
    !integer, intent(out) :: npool_local
      ! Number of pools for k-point parallelization
    !integer, intent(out) :: nproc_local
      ! Number of processes
    !integer, intent(out) :: nproc_pool_local
      ! Number of processes per pool
    !integer, intent(out) :: world_comm_local
      ! World communicator


    call MPI_Init(ierr)
    if (ierr /= 0) call mpiExitError( 8001 )

    world_comm_local = MPI_COMM_WORLD

    call MPI_COMM_RANK(world_comm_local, myid, ierr)
    if (ierr /= 0) call mpiExitError( 8002 )
      !! * Determine the rank or ID of the calling process
    call MPI_COMM_SIZE(world_comm_local, nproc_local, ierr)
    if (ierr /= 0) call mpiExitError( 8003 )
      !! * Determine the size of the MPI pool (i.e., the number of processes)

    ionode_local = (myid == root)
      ! Set a boolean for if this is the root process

    call getCommandLineArguments()
      !! * Get the number of pools from the command line

    call setUpImages()
      ! This sets up variables only used by QE, so it will be removed
      ! once extracted from QE. Only one image will be used.

    call setUpPools()
      !! * Split up processors between pools and generate MPI
      !!   communicators for pools

    call setUpBands()
      ! This sets up variables only used by QE, so it will be removed
      ! once extracted from QE. 

    call setUpDiag()
      ! This sets up variables only used by QE, so it will be removed
      ! once extracted from QE. 

    call setGlobalVariables()
      ! This sets up variables only used by QE, so it will be removed
      ! once extracted from QE.


    return
  end subroutine mpiInitialization

!----------------------------------------------------------------------------
  subroutine getCommandLineArguments()
    !! Get the command line arguments. This currently
    !! only processes the number of pools
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Output variables:
    !integer, intent(out) :: npool_local
      ! Number of pools for k-point parallelization


    ! Local variables:
    integer :: narg = 0
      !! Arguments processed
    integer :: nargs
      !! Total number of command line arguments
    integer :: npool_ = 1
      !! Number of k point pools for parallelization

    character(len=256) :: arg = ' '
      !! Command line argument
    character(len=256) :: command_line = ' '
      !! Command line arguments that were not processed


    nargs = command_argument_count()
      !! * Get the number of arguments input at command line

    call MPI_BCAST(nargs, 1, MPI_INTEGER, root, world_comm_local, ierr)

    if(ionode_local) then

      do while (narg <= nargs)
        call get_command_argument(narg, arg)
          !! * Get the flag
          !! @note
          !!  This program only currently processes the number of pools,
          !!  represented by `-nk`/`-npool`/`-npools`. All other flags 
          !!  will be ignored.
          !! @endnote

        narg = narg + 1

        !> * Process the flag and store the following value
        select case (trim(arg))
          case('-nk', '-npool', '-npools') 
            call get_command_argument(narg, arg)
            read(arg, *) npool_
            narg = narg + 1
          case default
            command_line = trim(command_line) // ' ' // trim(arg)
        end select
      enddo

      write(stdout,*) 'Unprocessed command line arguments: ' // trim(command_line)
    endif

    call MPI_BCAST(npool_, 1, MPI_INTEGER, root, world_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8005)

    npool_local = npool_

    return
  end subroutine getCommandLineArguments

!----------------------------------------------------------------------------
  subroutine setUpImages()
    !! @todo Remove this once extracted from QE #end @endtodo

    use mp_images, only : nproc_image, me_image, intra_image_comm

    implicit none

    intra_image_comm = world_comm_local

    nproc_image = nproc_local
      !! * Calculate how many processes there are per image

    me_image = myid
      !! * Get the index of the process within the image

    return
  end subroutine setUpImages

!----------------------------------------------------------------------------
  subroutine setUpPools()
    !! Split up processors between pools and generate MPI
    !! communicators for pools
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    !integer, intent(in) :: myid
      ! ID of this process
    !integer, intent(in) :: npool_local
      ! Number of pools for k-point parallelization
    !integer, intent(in) :: nproc_local
      ! Number of processes


    ! Output variables:
    !integer, intent(out) :: inter_pool_comm_local = 0
      ! Inter-pool communicator
    !integer, intent(out) :: intra_pool_comm_local = 0
      ! Intra-pool communicator
    !integer, intent(out) :: indexInPool
      ! Process index within pool
    !integer, intent(out) :: myPoolId
      ! Pool index for this process
    !integer, intent(out) :: nproc_pool_local
      ! Number of processes per pool


    if(npool_local < 1 .or. npool_local > nproc_local) call exitError('mpiInitialization', &
      'invalid number of pools, out of range', 1)
      !! * Verify that the number of pools is between 1 and the number of processes

    if(mod(nproc_local, npool_local) /= 0) call exitError('mpiInitialization', &
      'invalid number of pools, mod(nproc_local,npool_local) /=0 ', 1)
      !! * Verify that the number of processes is evenly divisible by the number of pools

    nproc_pool_local = nproc_local / npool_local
      !! * Calculate how many processes there are per pool

    myPoolId = myid / nproc_pool_local
      !! * Get the pool index for this process

    indexInPool = mod(myid, nproc_pool_local)
      !! * Get the index of the process within the pool

    call MPI_BARRIER(world_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8007)

    call MPI_COMM_SPLIT(world_comm_local, myPoolId, myid, intra_pool_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8008)
      !! * Create intra-pool communicator

    call MPI_BARRIER(world_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8009)

    call MPI_COMM_SPLIT(world_comm_local, indexInPool, myid, inter_pool_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8010)
      !! * Create inter-pool communicator

    return
  end subroutine setUpPools

!----------------------------------------------------------------------------
  subroutine setUpBands()
    !! @todo Remove this once extracted from QE #end @endtodo

    use mp_bands, only : nproc_bgrp, me_bgrp, intra_bgrp_comm, inter_bgrp_comm, &
        my_bgrp_id

    implicit none

    nproc_bgrp = nproc_local

    me_bgrp = myid

    call MPI_BARRIER(intra_pool_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8010)

    call MPI_COMM_SPLIT(intra_pool_comm_local, my_bgrp_id, myid, intra_bgrp_comm, ierr)
    if(ierr /= 0) call mpiExitError(8011)

    call MPI_BARRIER(intra_pool_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8012)

    call MPI_COMM_SPLIT(intra_pool_comm_local, me_bgrp, myid, inter_bgrp_comm, ierr)
    if(ierr /= 0) call mpiExitError(8013)

    return
  end subroutine setUpBands

!----------------------------------------------------------------------------
  subroutine setUpDiag()
    !! @todo Remove this once extracted from QE #end @endtodo

    use mp_diag, only : np_ortho, me_ortho, nproc_ortho, leg_ortho, &
        ortho_comm, ortho_parent_comm, me_ortho1, ortho_comm_id, &
        ortho_row_comm, ortho_col_comm

    implicit none

    integer :: nproc_all, color, key

    call MPI_COMM_SIZE(intra_pool_comm_local, nproc_all, ierr)
    if (ierr /= 0) call mpiExitError( 8014 )

    np_ortho = 1 

    if( nproc_all >= 4*nproc_ortho ) then
       !  Here we choose a processor every 4, in order not to stress memory BW
       !  on multi core procs, for which further performance enhancements are
       !  possible using OpenMP BLAS inside regter/cegter/rdiaghg/cdiaghg
       !  (to be implemented)

       color = 0
       if( indexInPool < 4*nproc_ortho .AND. MOD( indexInPool, 4 ) == 0 ) color = 1
       
       leg_ortho = 4
       
    else if( nproc_all >= 2*nproc_ortho ) then
       !  here we choose a processor every 2, in order not to stress memory BW
       
       color = 0
       if( indexInPool < 2*nproc_ortho .AND. MOD( indexInPool, 2 ) == 0 ) color = 1
       
       leg_ortho = 2
       
    else
       !  here we choose the first processors
       
       color = 0
       if( indexInPool < nproc_ortho ) color = 1
       !
       leg_ortho = 1
       !
    end if

    key = indexInPool
    
    call MPI_COMM_SPLIT(intra_pool_comm_local, color, key, ortho_comm, ierr)
    if(ierr /= 0) call mpiExitError(8015)
      !  initialize the communicator for the new group by splitting the input communicator
    
    ortho_parent_comm = intra_pool_comm_local
      ! and remember where it comes from
    
    me_ortho1 = mp_rank( ortho_comm )
      !  Computes coordinates of the processors, in row maior order
    
    IF( indexInPool == 0 .AND. me_ortho1 /= 0 ) &
         CALL exitError( " init_ortho_group ", " wrong root task in ortho group ", ierr )
    !
    if( color == 1 ) then
      ortho_comm_id = 1
      CALL GRID2D_COORDS( 'R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
      CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
      if( ierr /= me_ortho1 ) &
          CALL exitError( " init_ortho_group ", " wrong task coordinates in ortho group ", ierr )
      if( me_ortho1*leg_ortho /= indexInPool ) &
          CALL exitError( " init_ortho_group ", " wrong rank assignment in ortho group ", ierr )

      call MPI_COMM_SPLIT(ortho_comm, me_ortho(2), me_ortho(1), ortho_col_comm, ierr)
      if(ierr /= 0) call mpiExitError(8017)
      call MPI_COMM_SPLIT(ortho_comm, me_ortho(1), me_ortho(2), ortho_row_comm, ierr)
      if(ierr /= 0) call mpiExitError(8018)

    else
       ortho_comm_id = 0
       me_ortho(1) = me_ortho1
       me_ortho(2) = me_ortho1
    endif

    return
  end subroutine setUpDiag

!----------------------------------------------------------------------------
  subroutine setGlobalVariables()
    !! @todo Remove this once extracted from QE #end @endtodo

    use mp_world, only : world_comm, nproc, mpime
    use mp_pools, only : npool, nproc_pool, me_pool, my_pool_id, intra_pool_comm, inter_pool_comm
    use io_global, only : ionode, meta_ionode, meta_ionode_id

    implicit none

    ionode = ionode_local
    meta_ionode = (myid == root)
    meta_ionode_id = root
    world_comm = world_comm_local
    nproc = nproc_local
    mpime = myid
    npool = npool_local
    nproc_pool = nproc_pool_local
    me_pool = indexInPool
    my_pool_id = myPoolId
    inter_pool_comm = inter_pool_comm_local
    intra_pool_comm = intra_pool_comm_local

    return
  end subroutine setGlobalVariables

!----------------------------------------------------------------------------
  subroutine initialize(exportDir, QEDir, VASPDir)
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>
    !!

    use io_files, only : nd_nmbr
      !! @todo Remove this once extracted from QE #end @endtodo
    
    implicit none

    ! Input variables:
    !integer, intent(in) :: npool_local
      ! Number of pools for k-point parallelization
    !integer, intent(in) :: nproc_local
      ! Number of processes


    ! Output variables:
    character(len=256), intent(out) :: exportDir
      !! Directory to be used for export
    character(len=256), intent(out) :: QEDir
      !! Directory with QE files
    character(len=256), intent(out) :: VASPDir
      !! Directory with VASP files


    ! Local variables:
    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time

    character(len=6), external :: int_to_char
      !! @todo Remove this once extracted from QE #end @endtodo


    prefix = ''
      !! @todo Create local version of `prefix` #end @endtodo
    QEDir = './'
    VASPDir = './'
    exportDir = './Export'

#ifdef __INTEL_COMPILER
    call remove_stack_limit()
      !! * Removed the stack limit because Intel compiler allocates a lot of stack space
      !!   which leads to seg faults and crash. This always works unlike `ulimit -s unlimited`
#endif

    call cpu_time(tStart)

    call date_and_time(cdate, ctime)

#ifdef __MPI
    nd_nmbr = trim(int_to_char(myid+1))
      !! @todo Create local version of `nd_nmbr` #end @endtodo
#else
    nd_nmbr = ' '
#endif

    if(ionode_local) then

      write(stdout, '(/5X,"VASP wavefunction export program starts on ",A9," at ",A9)') &
             cdate, ctime

#ifdef __MPI
      write(stdout, '(/5X,"Parallel version (MPI), running on ",I5," processors")') nproc_local

      if(npool_local > 1) write(stdout, '(5X,"K-points division:     npool_local     = ",I7)') npool_local
#else
      write(stdout, '(/5X,"Serial version")')
#endif

    else

      open(unit = stdout, file='/dev/null', status='unknown')
        ! Make the stdout unit point to null for non-root processors
        ! to avoid tons of duplicate output

    endif

  end subroutine initialize

!----------------------------------------------------------------------------
  subroutine mpiSumIntV(msg, comm)
    implicit none

    ! Input/output variables:
    integer, intent(in) :: comm
      !! MPI communicator
    integer, intent(inout) :: msg(:)
      !! Message to be sent


#if defined(__MPI)
    ! Local variables:
    integer, parameter :: maxb = 100000
      !! Max buffer size

    integer :: ib
      !! Loop index
    integer :: buff(maxb)
      !! Buffer
    integer :: msglen
      !! Length of message to be sent
    integer :: nbuf
      !! Number of buffers

    msglen = size(msg)

    nbuf = msglen/maxb
      !! * Get the number of buffers of size `maxb` needed
  
    do ib = 1, nbuf
      !! * Send message in buffers of size `maxb`
     
        call MPI_ALLREDUCE(msg(1+(ib-1)*maxb), buff, maxb, MPI_INTEGER, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumIntV', 'error in mpi_allreduce 1', ierr)

        msg((1+(ib-1)*maxb):(ib*maxb)) = buff(1:maxb)

    enddo

    if((msglen - nbuf*maxb) > 0 ) then
      !! * Send any data left of size less than `maxb`

        call MPI_ALLREDUCE(msg(1+nbuf*maxb), buff, (msglen-nbuf*maxb), MPI_INTEGER, MPI_SUM, comm, ierr)
        if(ierr /= 0) call exitError('mpiSumIntV', 'error in mpi_allreduce 2', ierr)

        msg((1+nbuf*maxb):msglen) = buff(1:(msglen-nbuf*maxb))
    endif

#endif

    return
  end subroutine mpiSumIntV

!----------------------------------------------------------------------------
  subroutine mpiExitError(code)
    !! Exit on error with MPI communication

    implicit none
    
    integer, intent(in) :: code

    write( stdout, '( "*** MPI error ***")' )
    write( stdout, '( "*** error code: ",I5, " ***")' ) code

    call MPI_ABORT(world_comm_local,code,ierr)
    
    stop

    return
  end subroutine mpiExitError

!----------------------------------------------------------------------------
  subroutine exitError(calledFrom, message, ierror)
    !! Output error message and abort if ierr > 0
    !!
    !! Can ensure that error will cause abort by
    !! passing abs(ierror)
    !!
    !! <h2>Walkthrough</h2>
    !!
    
    implicit none

    integer, intent(in) :: ierror
      !! Error

    character(len=*), intent(in) :: calledFrom
      !! Place where this subroutine was called from
    character(len=*), intent(in) :: message
      !! Error message

    integer :: id
      !! ID of this process
    integer :: mpierr
      !! Error output from MPI

    character(len=6) :: cerr
      !! String version of error


    if ( ierr <= 0 ) return
      !! * Do nothing if the error is less than or equal to zero

    write( cerr, fmt = '(I6)' ) ierror
      !! * Write ierr to a string
    write(unit=*, fmt = '(/,1X,78("%"))' )
      !! * Output a dividing line
    write(unit=*, fmt = '(5X,"Error in ",A," (",A,"):")' ) trim(calledFrom), trim(adjustl(cerr))
      !! * Output where the error occurred and the error
    write(unit=*, fmt = '(5X,A)' ) TRIM(message)
      !! * Output the error message
    write(unit=*, fmt = '(1X,78("%"),/)' )
      !! * Output a dividing line

    write( *, '("     stopping ...")' )
  
    call flush( stdout )

#if defined (__MPI)
  
    id = 0
  
    !> * For MPI, get the id of this process and abort
    call MPI_COMM_RANK( world_comm_local, id, mpierr )
    call MPI_ABORT( world_comm_local, mpierr, ierr )
    call MPI_FINALIZE( mpierr )

#endif

    stop 2

    return

  end subroutine exitError

!----------------------------------------------------------------------------
  subroutine readWAVECAR(VASPDir, at_local, bg_local, ecutwfc_local, occ, omega_local, vcut_local, &
        xk_local, nb1max, nb2max, nb3max, nbnd_local, ngk_max, nkstot_local, nplane_g, nspin_local)
    !! Read cell and wavefunction data from the WAVECAR file
    !!
    !! <h2>Walkthrough</h2>
    !!

    use cell_base, only : at, bg, omega, alat, tpiba
    use wvfct, only : ecutwfc, nbnd
    use klist, only : nkstot
    use lsda_mod, only : nspin
      !! @todo Remove this once extracted from QE #end @endtodo

    implicit none

    ! Input variables:
    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    
    ! Output variables:
    real(kind=dp), intent(out) :: at_local(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: bg_local(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(out) :: ecutwfc_local
      !! Plane wave energy cutoff in Ry
    real(kind=dp), allocatable, intent(out) :: occ(:,:)
      !! Occupation of band
    real(kind=dp), intent(out) :: omega_local
      !! Volume of unit cell
    real(kind=dp), intent(out) :: vcut_local
      !! Energy cutoff converted to vector cutoff
    real(kind=dp), allocatable, intent(out) :: xk_local(:,:)
      !! Position of k-points in reciprocal space

    integer, intent(out) :: nb1max, nb2max, nb3max
      !! Not sure what this is??
    integer, intent(out) :: nbnd_local
      !! Total number of bands
    integer, intent(out) :: ngk_max
      !! Maximum number of \(G+k\) combinations
    integer, intent(out) :: nkstot_local
      !! Total number of k-points
    integer, allocatable, intent(out) :: nplane_g(:)
      !! Input number of plane waves for a single k-point 
      !! for all processors
    integer, intent(out) :: nspin_local
      !! Number of spins


    ! Local variables:
    real(kind=dp) :: c = 0.26246582250210965422
      !! \(2m/\hbar^2\) converted from J\(^{-1}\)m\(^{-2}\)
      !! to eV\(^{-1}\)A\(^{-2}\)
    real(kind=dp) :: nRecords_real, nspin_real, prec_real, nkstot_real 
      !! Real version of integers for reading from file
    real(kind=dp) :: nbnd_real
      !! Real version of integers for reading from file

    integer :: j
      !! Index used for reading lattice vectors
    integer :: nRecords
      !! Number of records in WAVECAR file
    integer :: prec
      !! Precision of plane wave coefficients

    character(len=256) :: fileName
      !! Full WAVECAR file name including path


    if(ionode_local) then
      fileName = trim(VASPDir)//'/WAVECAR'

      nRecords = 24
        ! Set a starting value for the number of records

      open(unit=wavecarUnit, file=fileName, access='direct', recl=nRecords, iostat=ierr, status='old')
      if (ierr .ne. 0) write(stdout,*) 'open error - iostat =', ierr
        !! * If root node, open the `WAVECAR` file

      read(unit=wavecarUnit,rec=1) nRecords_real, nspin_real, prec_real
        !! @note Must read in as real first then convert to integer @endnote

      close(unit=wavecarUnit)

      nRecords = nint(nRecords_real)
      nspin_local = nint(nspin_real)
      prec = nint(prec_real)
        ! Convert input variables to integers

      !if(prec .eq. 45210) call exitError('readWAVECAR', 'WAVECAR_double requires complex*16', 1)

      open(unit=wavecarUnit, file=fileName, access='direct', recl=nRecords, iostat=ierr, status='old')
      if (ierr .ne. 0) write(stdout,*) 'open error - iostat =', ierr

      read(unit=wavecarUnit,rec=2) nkstot_real, nbnd_real, ecutwfc_local,(at_local(j,1),j=1,3),&
          (at_local(j,2),j=1,3), (at_local(j,3),j=1,3)
        !! * Read total number of k-points, plane wave cutoff energy, and real
        !!   space lattice vectors

      ecutwfc_local = ecutwfc_local*eVToRy
        !! * Convert energy from VASP from eV to Rydberg to match QE/TME expectation

      vcut_local = sqrt(ecutwfc_local/evToRy*c)
        !! * Calculate vector cutoff from energy cutoff

      nkstot_local = nint(nkstot_real)
      nbnd_local = nint(nbnd_real)
        ! Convert input variables to integers

      call calculateOmega(at_local, omega_local)
        !! * Calculate the cell volume as \(a_1\cdot a_2\times a_3\)

      call getReciprocalVectors(at_local, omega_local, bg_local)
        !! * Calculate the reciprocal lattice vectors from the real-space
        !!   lattice vectors and the cell volume

      call estimateMaxNumPlanewaves(bg_local, ecutwfc_local, nb1max, nb2max, nb3max, ngk_max)
        !! * Get the maximum number of plane waves

      !> * Write out total number of k-points, number of bands, 
      !>   the energy cutoff, the real-space-lattice vectors,
      !>   the cell volume, and the reciprocal lattice vectors
      write(stdout,*) 'no. k points =', nkstot_local
      write(stdout,*) 'no. bands =', nbnd_local
      write(stdout,*) 'max. energy (eV) =', sngl(ecutwfc_local/eVToRy)
        !! @note 
        !!  The energy cutoff is currently output to the `stdout` file
        !!  in eV to compare with output from WaveTrans.
        !! @endnote
      write(stdout,*) 'real space lattice vectors:'
      write(stdout,*) 'a1 =', (sngl(at_local(j,1)),j=1,3)
      write(stdout,*) 'a2 =', (sngl(at_local(j,2)),j=1,3)
      write(stdout,*) 'a3 =', (sngl(at_local(j,3)),j=1,3)
      write(stdout,*) 
      write(stdout,*) 'volume unit cell =', sngl(omega_local)
      write(stdout,*) 
      write(stdout,*) 'reciprocal lattice vectors:'
      write(stdout,*) 'b1 =', (sngl(bg_local(j,1)),j=1,3)
      write(stdout,*) 'b2 =', (sngl(bg_local(j,2)),j=1,3)
      write(stdout,*) 'b3 =', (sngl(bg_local(j,3)),j=1,3)
      write(stdout,*) 
        !! @note
        !!  I made an intentional choice to stick with the unscaled lattice
        !!  vectors until I see if it will be convenient to scale them down.
        !!  QE uses the `alat` and `tpiba` scaling quite a bit though, so I
        !!  will have to be careful with the scaling/units.
        !! @endnote

      write(mainout, '("# Cell volume (a.u.)^3. Format: ''(ES24.15E3)''")')
      write(mainout, '(ES24.15E3)' ) omega_local
      flush(mainout)

    endif

    call MPI_BCAST(nb1max, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nb2max, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nb3max, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(ngk_max, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nspin_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nkstot_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nbnd_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(ecutwfc_local, 1, MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(vcut_local, 1, MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(omega_local, 1, MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(at_local, size(at_local), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(bg_local, size(bg_local), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)

    if (ionode_local) then
      write(stdout,*) 
      write(stdout,*) '******'
      write(stdout,*) 'Starting to read wavefunction'
    endif

    call readWavefunction(nbnd_local, ngk_max, nkstot_local, nspin_local, occ, xk_local, nplane_g)
      !! * Get the position of each k-point in reciprocal space 
      !!   and the number of \(G+k) vectors below the cutoff 
      !!   energy for each k-point

    if(ionode_local) close(wavecarUnit)

    if (ionode_local) then
      write(stdout,*) 'Finished reading wavefunction'
      write(stdout,*) '******'
      write(stdout,*) 
    endif

    at = at_local/alat
    bg = bg_local/tpiba
    omega = omega_local
    ecutwfc = ecutwfc_local
    nbnd = nbnd_local
    nkstot = nkstot_local
    nspin = nspin_local
      !! @todo Remove QE variable assignment once extracted from QE #end @endtodo

    return
  end subroutine readWAVECAR

!----------------------------------------------------------------------------
  subroutine calculateOmega(at_local, omega_local)
    !! Calculate the cell volume as \(a_1\cdot a_2\times a_3\)

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: at_local(3,3)
      !! Real space lattice vectors


    ! Output variables:
    real(kind=dp), intent(out) :: omega_local
      !! Volume of unit cell


    ! Local variables:
    real(kind=dp) :: vtmp(3)
      !! \(a_2\times a_3\)


    call vcross(at_local(:,2), at_local(:,3), vtmp)

    omega_local = sum(at_local(:,1)*vtmp(:))

    return
  end subroutine calculateOmega

!----------------------------------------------------------------------------
  subroutine getReciprocalVectors(at_local, omega_local, bg_local)
    !! Calculate the reciprocal lattice vectors from the real-space
    !! lattice vectors and the cell volume

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: at_local(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: omega_local
      !! Volume of unit cell


    ! Output variables:
    real(kind=dp), intent(out) :: bg_local(3,3)
      !! Reciprocal lattice vectors


    ! Local variables:
    integer :: i
      !! Loop index
    

    call vcross(2.0d0*pi*at_local(:,2)/omega_local, at_local(:,3), bg_local(:,1))
      ! \(b_1 = 2\pi/\Omega a_2\times a_3\)
    call vcross(2.0d0*pi*at_local(:,3)/omega_local, at_local(:,1), bg_local(:,2))
      ! \(b_2 = 2\pi/\Omega a_3\times a_1\)
    call vcross(2.0d0*pi*at_local(:,1)/omega_local, at_local(:,2), bg_local(:,3))
      ! \(b_3 = 2\pi/\Omega a_1\times a_2\)


    return
  end subroutine getReciprocalVectors

!----------------------------------------------------------------------------
  subroutine vcross(vec1, vec2, crossProd)
    !! Calculate the cross product `crossProd` of 
    !! two vectors `vec1` and `vec2`

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: vec1(3), vec2(3)
      !! Input vectors


    ! Output variables:
    real(kind=dp), intent(out) :: crossProd(3)
      !! Cross product of input vectors


    crossProd(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    crossProd(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    crossProd(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

    return
  end subroutine vcross

!----------------------------------------------------------------------------
  subroutine estimateMaxNumPlanewaves(bg_local, ecutwfc_local, nb1max, nb2max, nb3max, ngk_max)
    !! Get the maximum number of plane waves. I'm not sure how 
    !! this is done completely. It seems to be just basic vector
    !! stuff, but I haven't been able to make sense of it.
    !! 
    !! @todo Figure out how `estimateMaxNumPlanewaves` works @endtodo
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: bg_local(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: ecutwfc_local
      !! Plane wave energy cutoff in Ry


    ! Output variables:
    integer, intent(out) :: nb1max, nb2max, nb3max
      !! Not sure what this is??
    integer, intent(out) :: ngk_max
      !! Maximum number of \(G+k\) combinations


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

    integer :: nb1maxA, nb2maxA, nb3maxA
      !! Not sure what this is??
    integer :: nb1maxB, nb2maxB, nb3maxB
      !! Not sure what this is??
    integer :: nb1maxC, nb2maxC, nb3maxC
      !! Not sure what this is??
    integer :: npmaxA, npmaxB, npmaxC
      !! Not sure what this is??


    b1mag = sqrt(sum(bg_local(:,1)**2))
    b2mag = sqrt(sum(bg_local(:,2)**2))
    b3mag = sqrt(sum(bg_local(:,3)**2))

    write(stdout,*) 'reciprocal lattice vector magnitudes:'
    write(stdout,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)
      !! * Calculate and output reciprocal vector magnitudes


    phi12 = acos(sum(bg_local(:,1)*bg_local(:,2))/(b1mag*b2mag))
      !! * Calculate angle between \(b_1\) and \(b_2\)

    call vcross(bg_local(:,1), bg_local(:,2), vtmp)
    vmag = sqrt(sum(vtmp(:)**2))
    sinphi123 = sum(bg_local(:,3)*vtmp(:))/(vmag*b3mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxA = (dsqrt(ecutwfc_local/eVToRy*c)/(b1mag*abs(sin(phi12)))) + 1
    nb2maxA = (dsqrt(ecutwfc_local/eVToRy*c)/(b2mag*abs(sin(phi12)))) + 1
    nb3maxA = (dsqrt(ecutwfc_local/eVToRy*c)/(b3mag*abs(sinphi123))) + 1
    npmaxA = nint(4.0*pi*nb1maxA*nb2maxA*nb3maxA/3.0)
      !! * Get first set of max values


    phi13 = acos(sum(bg_local(:,1)*bg_local(:,3))/(b1mag*b3mag))
      !! * Calculate angle between \(b_1\) and \(b_3\)

    call vcross(bg_local(:,1), bg_local(:,3), vtmp)
    vmag = sqrt(sum(vtmp(:)**2))
    sinphi123 = sum(bg_local(:,2)*vtmp(:))/(vmag*b2mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxB = (dsqrt(ecutwfc_local/eVToRy*c)/(b1mag*abs(sin(phi13)))) + 1
    nb2maxB = (dsqrt(ecutwfc_local/eVToRy*c)/(b2mag*abs(sinphi123))) + 1
    nb3maxB = (dsqrt(ecutwfc_local/eVToRy*c)/(b3mag*abs(sin(phi13)))) + 1
    npmaxB = nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      !! * Get first set of max values


    phi23 = acos(sum(bg_local(:,2)*bg_local(:,3))/(b2mag*b3mag))
      !! * Calculate angle between \(b_2\) and \(b_3\)

    call vcross(bg_local(:,2), bg_local(:,3), vtmp)
    vmag = sqrt(sum(vtmp(:)**2))
    sinphi123 = sum(bg_local(:,1)*vtmp(:))/(vmag*b1mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxC = (dsqrt(ecutwfc_local/eVToRy*c)/(b1mag*abs(sinphi123))) + 1
    nb2maxC = (dsqrt(ecutwfc_local/eVToRy*c)/(b2mag*abs(sin(phi23)))) + 1
    nb3maxC = (dsqrt(ecutwfc_local/eVToRy*c)/(b3mag*abs(sin(phi23)))) + 1
    npmaxC = nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)
      !! * Get first set of max values


    nb1max = max0(nb1maxA,nb1maxB,nb1maxC)
    nb2max = max0(nb2maxA,nb2maxB,nb2maxC)
    nb3max = max0(nb3maxA,nb3maxB,nb3maxC)
    ngk_max = min0(npmaxA,npmaxB,npmaxC)

    write(stdout,*) 'max. no. G values; 1,2,3 =', nb1max, nb2max, nb3max
    write(stdout,*) ' '

    write(stdout,*) 'estimated max. no. plane waves =', ngk_max

    return
  end subroutine estimateMaxNumPlanewaves

!----------------------------------------------------------------------------
  subroutine readWavefunction(nbnd_local, ngk_max, nkstot_local, nspin_local, occ, xk_local, nplane_g)
    !! For each spin and k-point, read the number of
    !! \(G+k\) vectors below the energy cutoff, the
    !! position of the k-point in reciprocal space, 
    !! and the eigenvalue, occupation, and plane
    !! wave coefficients for each band
    !!
    !! <h2>Walkthrough</h2>
    !!

    use klist, only : xk
      !! @todo Remove this once extracted from QE #end @endtodo

    implicit none

    ! Input variables:
    integer, intent(in) :: nbnd_local
      !! Total number of bands
    integer, intent(in) :: ngk_max
      !! Maximum number of \(G+k\) combinations
    integer, intent(in) :: nkstot_local
      !! Total number of k-points
    integer, intent(in) :: nspin_local
      !! Number of spins


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: occ(:,:)
      !! Occupation of band
    real(kind=dp), allocatable, intent(out) :: xk_local(:,:)
      !! Position of k-points in reciprocal space

    integer, allocatable, intent(out) :: nplane_g(:)
      !! Input number of plane waves for a single k-point 
      !! for all processors


    ! Local variables:
    real(kind=dp) :: nplane_g_real
      !! Real version of integers for reading from file

    complex*16, allocatable :: cener(:)
      !! Band eigenvalues
    complex*8, allocatable :: coeff(:,:)
      !! Plane wave coefficients

    integer :: irec, isp, ik, i, iband, iplane
      !! Loop indices


    allocate(occ(nbnd_local, nkstot_local))
    allocate(xk_local(3,nkstot_local))
    allocate(nplane_g(nkstot_local))

    if(ionode_local) then

      allocate(cener(nbnd_local))
      allocate(coeff(ngk_max,nbnd_local))
    
      irec=2

      do isp = 1, nspin_local
        !! * For each spin:
        !!    * Go through each k-point
        !!       * Read in the number of \(G+k\) plane wave
        !!         vectors below the energy cutoff
        !!       * Read the position of the k-point in 
        !!         reciprocal space
        !!       * Read in the eigenvalue and occupation for
        !!         each band
        !!       * Read in the plane wave coefficients for
        !!         each band

       write(stdout,*) '  Reading spin ', isp

       do ik = 1, nkstot_local

        write(stdout,*) '    Reading k-point ', ik

        irec = irec + 1
       
        read(unit=wavecarUnit,rec=irec) nplane_g_real, (xk_local(i,ik),i=1,3), &
               (cener(iband), occ(iband, ik), iband=1,nbnd_local)
          ! Read in the number of \(G+k\) plane wave vectors below the energy
          ! cutoff, the position of the k-point in reciprocal space, and
          ! the eigenvalue and occupation for each band

        nplane_g(ik) = nint(nplane_g_real)

        do iband = 1, nbnd_local

          irec = irec + 1

          read(unit=wavecarUnit,rec=irec) (coeff(iplane,iband), iplane=1,nplane_g(ik))
            ! Read in the plane wave coefficients for each band

          write(45+ik,*) cener(iband)
            !! @todo 
            !!  Figure out how `cener` and `eigF`/`eigI` relate to `et`
            !! @endtodo

        enddo

        write(45+ik,*) "--------------------------------------------------------"

      enddo
    enddo

      deallocate(cener)
      deallocate(coeff)
        !! @note 
        !!  The band eigenvalues and the plane wave coefficients are 
        !!  not currently used anywhere.
        !! @endnote

    endif

    call MPI_BCAST(xk_local, size(xk_local), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(occ, size(occ), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(nplane_g, size(nplane_g), MPI_INTEGER, root, world_comm_local, ierr)

    xk = xk_local
      !! @todo Remove this once extracted from QE #end @endtodo

    return
  end subroutine readWavefunction

!----------------------------------------------------------------------------
  subroutine distributeKpointsInPools(nkstot_local)
    !! Figure out how many k-points there should be per pool
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nkstot_local
      !! Total number of k-points


    ! Output variables:
    !integer, intent(out) :: ikEnd
      ! Ending index for k-points in single pool 
    !integer, intent(out) :: ikStart
      ! Starting index for k-points in single pool 
    !integer, intent(out) :: nk_Pool
      ! Number of k-points in each pool


    ! Local variables:
    integer :: nkr
      !! Number of k-points left over after evenly divided across pools


    if( nkstot_local > 0 ) then

      IF( ( nproc_pool_local > nproc_local ) .or. ( mod( nproc_local, nproc_pool_local ) /= 0 ) ) &
        CALL exitError( 'distributeKpointsInPools','nproc_pool_local', 1 )

      nk_Pool = nkstot_local / npool_local
        !!  * Calculate k-points per pool

      nkr = nkstot_local - nk_Pool * npool_local 
        !! * Calculate the remainder `nkr`

      IF( myPoolId < nkr ) nk_Pool = nk_Pool + 1
        !! * Assign the remainder to the first `nkr` pools

      !>  * Calculate the index of the first k-point in this pool
      ikStart = nk_Pool * myPoolId + 1
      IF( myPoolId >= nkr ) ikStart = ikStart + nkr

      ikEnd = ikStart + nk_Pool - 1
        !!  * Calculate the index of the last k-point in this pool

    endif

    return
  end subroutine distributeKpointsInPools

!----------------------------------------------------------------------------
  subroutine calculateGvecs(nb1max, nb2max, nb3max, bg_local, gCart_local, ig_l2g, mill_g, &
      ngm_g_local, ngm_local)
    !! Calculate Miller indices and G-vectors and split
    !! over processors
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nb1max, nb2max, nb3max
      !! Not sure what this is??

    real(kind=dp), intent(in) :: bg_local(3,3)
      !! Reciprocal lattice vectors


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: gCart_local(:,:)
      !! G-vectors in Cartesian coordinates

    integer, allocatable, intent(out) :: ig_l2g(:)
      !! Converts local index `ig` to global index
    integer, allocatable, intent(out) :: mill_g(:,:)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(out) :: ngm_g_local
      !! Global number of G-vectors
    integer, intent(out) :: ngm_local
      !! Local number of G-vectors on this processor


    ! Local variables:
    integer :: ig1, ig2, ig3, ig1p, ig2p, ig3p, ig, ix
      !! Loop indices
    integer :: igEnd
      !! Ending index for G-vectors across processors 
    integer :: igStart
      !! Starting index for G-vectors across processors 
    integer, allocatable :: mill_local(:,:)
      !! Integer coefficients for G-vectors
    integer :: npmax
      !! Max number of plane waves


    npmax = (2*nb1max+1)*(2*nb2max+1)*(2*nb3max+1) 
    allocate(mill_g(3,npmax))

    if(ionode_local) then
      write(stdout,*)
      write(stdout,*) "***************"
      write(stdout,*) "Calculating miller indices"

      ngm_g_local = 0
      mill_g = 0

      do ig3 = 0, 2*nb3max

        ig3p = ig3

        if (ig3 .gt. nb3max) ig3p = ig3 - 2*nb3max - 1

        !write(stdout,*) " Outer miller index: ", ig3p

        do ig2 = 0, 2*nb2max

          ig2p = ig2

          if (ig2 .gt. nb2max) ig2p = ig2 - 2*nb2max - 1

          do ig1 = 0, 2*nb1max

            ig1p = ig1

            if (ig1 .gt. nb1max) ig1p = ig1 - 2*nb1max - 1

            ngm_g_local = ngm_g_local + 1

            mill_g(1,ngm_g_local) = ig1p
            mill_g(2,ngm_g_local) = ig2p
            mill_g(3,ngm_g_local) = ig3p
              !! * Calculate Miller indices

          enddo
        enddo
      enddo

      if (ngm_g_local .ne. npmax) call exitError('calculateGvecs', & 
        '*** error - computed no. of G-vectors != estimated number of plane waves', 1)
        !! * Check that number of G-vectors are the same as the number of plane waves

      write(stdout,*) "Done calculating miller indices"
      write(stdout,*) "***************"
      write(stdout,*)
    endif

    call MPI_BCAST(ngm_g_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(mill_g, size(mill_g), MPI_INTEGER, root, world_comm_local, ierr)

    if (ionode_local) then
      write(stdout,*)
      write(stdout,*) "***************"
      write(stdout,*) "Distributing G-vecs over processors"
    endif

    call distributeGvecsOverProcessors(ngm_g_local, mill_g, ig_l2g, igEnd, igStart, &
        mill_local, ngm_local)
      !! * Split up the G-vectors and Miller indices over processors 

    if (ionode_local) write(stdout,*) "Calculating G-vectors"

    allocate(gCart_local(3,ngm_local))

    do ig = 1, ngm_local

      do ix = 1, 3
        !! * Calculate \(G = m_1b_1 + m_2b_2 + m_3b_3\)

        gCart_local(ix,ig) = sum(mill_local(:,ig)*bg_local(ix,:))

      enddo
      
    enddo

    if (ionode_local) then
      write(stdout,*) "***************"
      write(stdout,*)
    endif

    deallocate(mill_local)

    return
  end subroutine calculateGvecs

!----------------------------------------------------------------------------
  subroutine distributeGvecsOverProcessors(ngm_g_local, mill_g, ig_l2g, igEnd, igStart, &
      mill_local, ngm_local)
    !! Figure out how many G-vectors there should be per processor
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: ngm_g_local
      !! Global number of G-vectors
      
    integer, intent(in) :: mill_g(3,ngm_g_local)
      !! Integer coefficients for G-vectors on all processors

    
    ! Output variables:
    integer, allocatable, intent(out) :: ig_l2g(:)
      ! Converts local index `ig` to global index
    integer, intent(out) :: igEnd
      !! Ending index for G-vectors across processors 
    integer, intent(out) :: igStart
      !! Starting index for G-vectors across processors 
    integer, allocatable, intent(out) :: mill_local(:,:)
      !! Integer coefficients for G-vectors
    integer, intent(out) :: ngm_local
      !! Local number of G-vectors on this processor


    ! Local variables:
    integer :: ig
      !! Loop index
    integer :: ngr
      !! Number of G-vectors left over after evenly divided across processors


#if defined (__MPI)

    if( ngm_g_local > 0 ) then
      ngm_local = ngm_g_local/nproc_local
        !!  * Calculate number of G-vectors per processor

      ngr = ngm_g_local - ngm_local*nproc_local 
        !! * Calculate the remainder

      IF( myid < ngr ) ngm_local = ngm_local + 1
        !! * Assign the remainder to the first `ngr` processors

      !>  * Calculate the index of the first k point in this pool
      igStart = ngm_local * myid + 1
      IF( myid >= ngr ) igStart = igStart + ngr

      igEnd = igStart + ngm_local - 1
        !!  * Calculate the index of the last k point in this pool


      !> * Generate an array to map a local index
      !>   (`ig` passed to `ig_l2g`) to a global
      !>   index (the value stored at `ig_l2g(ig)`)
      !>   and get local miller indices
      allocate(ig_l2g(ngm_local))
      allocate(mill_local(3,ngm_local))

      do ig = 1, ngm_local

        ig_l2g(ig) = igStart + ig - 1 
        mill_local(:,ig) = mill_g(:,ig_l2g(ig))

      enddo

    endif

#else

    ngm_local = ngm_g_local

#endif

    return
  end subroutine distributeGvecsOverProcessors

!----------------------------------------------------------------------------
  subroutine reconstructFFTGrid(ngm_local, ig_l2g, ngk_max, nkstot_local, nplane_g, bg_local, gCart_local, &
      vcut_local, xk_local, igk_l2g, igk_large, ngk_local, ngk_g, npw_g, npwx_g, npwx_local)
    !! Determine which G-vectors result in \(G+k\)
    !! below the energy cutoff for each k-point and
    !! sort the indices based on \(|G+k|^2\)
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: ngm_local
      !! Number of G-vectors on this processor

    integer, intent(in) :: ig_l2g(ngm_local)
      ! Converts local index `ig` to global index
    integer, intent(in) :: ngk_max
      !! Maximum number of \(G+k\) combinations
    integer, intent(in) :: nkstot_local
      !! Total number of k-points
    integer, intent(in) :: nplane_g(nkstot_local)
      !! Input number of plane waves for a single k-point

    real(kind=dp), intent(in) :: bg_local(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: gCart_local(3,ngm_local)
      !! G-vectors in Cartesian coordinates
    real(kind=dp), intent(in) :: vcut_local
      !! Energy cutoff converted to vector cutoff
    real(kind=dp), intent(in) :: xk_local(3,nkstot_local)
      !! Position of k-points in reciprocal space


    ! Output variables:
    integer, allocatable, intent(out) :: igk_l2g(:,:)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point
    integer, allocatable, intent(out) :: igk_large(:,:)
      !! Index map from \(G\) to \(G+k\);
      !! indexed up to `ngm_local` which
      !! is greater than `npwx_local` and
      !! stored for each k-point
    integer, allocatable, intent(out) :: ngk_local(:)
      !! Number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each
      !! k-point, on this processor
    integer, allocatable, intent(out) :: ngk_g(:)
      !! Global number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each k-point
    integer, intent(out) :: npw_g
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(out) :: npwx_g
      !! Max number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` among all k-points
    integer, intent(out) :: npwx_local
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor


    ! Local variables:
    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp) :: gkMod(nk_Pool,ngm_local)
      !! \(|G+k|^2\);
      !! only stored if less than `vcut_local`
    real(kind=dp) :: q
      !! \(|q|^2\) where \(q = G+k\)
    real(kind=dp) :: xkCart(3)
      !! Cartesian coordinates for given k-point

    integer :: ik, ig, ix
      !! Loop indices
    integer, allocatable :: igk(:)
      !! Index map from \(G\) to \(G+k\)
      !! indexed up to `npwx_local`
    integer :: ngk_tmp
      !! Temporary variable to hold `ngk_local`
      !! value so that don't have to keep accessing
      !! array

    allocate(ngk_local(nk_Pool))
    allocate(igk_large(nk_Pool,ngm_local))
    
    npwx_local = 0
    ngk_local(:) = 0
    igk_large(:,:) = 0

    if (ionode_local) then
      write(stdout,*)
      write(stdout,*) "***************"
      write(stdout,*) "Determining G+k combinations less than energy cutoff"
    endif


    do ik = 1, nk_Pool
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

      if (ionode_local) write(stdout,*) "Processing k-point ", ik

      do ix = 1, 3
        xkCart(ix) = sum(xk_local(:,ik+ikStart-1)*bg_local(ix,:))
      enddo

      ngk_tmp = 0

      do ig = 1, ngm_local

        q = sqrt(sum((xkCart(:) + gCart_local(:,ig))**2))
          ! Calculate \(|G+k|\)

        if (q <= eps8) q = 0.d0
   
        if (q <= vcut_local) then

          ngk_tmp = ngk_tmp + 1
            ! If \(|G+k| \leq \) `vcut` increment the count for
            ! this k-point

          gkMod(ik,ngk_tmp) = q
            ! Store the modulus for sorting

          igk_large(ik,ngk_tmp) = ig
            ! Store the index for this G-vector

        !else

          !if (sqrt(sum(gCart_local(:, ig)**2)) .gt. &
            !sqrt(sum(xk_local(:,ik+ikStart-1)**2) + sqrt(vcut_local))) goto 100
            ! if |G| > |k| + sqrt(Ecut)  stop search
            !! @todo Figure out if there is valid exit check for `ig` loop @endtodo

        endif
      enddo

100   npwx_local = max(npwx_local, ngk_tmp)
        ! Track the maximum number of \(G+k\)
        ! vectors among all k-points

      ngk_local(ik) = ngk_tmp
        ! Store the total number of \(G+k\)
        ! vectors for this k-point

    enddo

    allocate(ngk_g(nkstot_local))
    ngk_g = 0
    ngk_g(ikStart:ikEnd) = ngk_local(1:nk_Pool)
    CALL mpiSumIntV(ngk_g, world_comm_local)
      !! * Calculate the global number of \(G+k\) 
      !!   vectors for each k-point
      
    if (ionode_local) then

      do ik = 1, nkstot_local

        if (ngk_g(ik) .ne. nplane_g(ik)) call exitError('reconstructFFTGrid', &
          'computed no. of G-vectors != input no. of plane waves', 1)
          !! * Make sure that number of G-vectors isn't higher than the calculated maximum

        if (ngk_g(ik) .gt. ngk_max) call exitError('reconstructFFTGrid', &
          'G-vector count exceeds estimate', 1)
          !! * Make sure that number of G-vectors isn't higher than the calculated maximum

      enddo
    endif

    if (ionode_local) then
      write(stdout,*) "Done determining G+k combinations less than energy cutoff"
      write(stdout,*) "***************"
      write(stdout,*)
    endif

    if (npwx_local <= 0) call exitError('reconstructFFTGrid', &
                'No plane waves found: running on too many processors?', 1)
      !! * Make sure that each processor gets some \(G+k\) vectors. If not,
      !!   should rerun with fewer processors.

    CALL mp_max(npwx_local, inter_pool_comm_local)
      !! * When using pools, set `npwx_local` to the maximum value across pools
      !! @todo Change call to QE `mp_max` to actual call to `MPI_ALLREDUCE` here @endtodo


    allocate(igk_l2g(npwx_local,nk_Pool))
    allocate(igk(npwx_local))

    igk_l2g = 0
    igk = 0

    if (ionode_local) then
      write(stdout,*)
      write(stdout,*) "***************"
      write(stdout,*) "Sorting G+k combinations by magnitude"
    endif

    do ik = 1, nk_Pool
      !! * Reorder the indices of the G-vectors so that
      !!   they are sorted by \(|G+k|^2\) for each k-point

      ngk_tmp = ngk_local(ik)

      igk(1:ngk_tmp) = igk_large(ik,1:ngk_tmp)

      call hpsort_eps(ngk_tmp, gkMod(ik,:), igk, eps8)
        ! Order vector `gkMod` keeping initial position in `igk`

      do ig = 1, ngk_tmp
        
        igk_l2g(ig,ik) = ig_l2g(igk(ig))
        
      enddo
     
      igk_l2g(ngk_tmp+1:npwx_local, ik) = 0

    enddo

    if (ionode_local) then
      write(stdout,*) "Done sorting G+k combinations by magnitude"
      write(stdout,*) "***************"
      write(stdout,*)
    endif

    deallocate(igk)

    npw_g = maxval(igk_l2g(:,:))
    call mp_max( npw_g, world_comm_local )
      !! * Calculate the maximum G-vector index 
      !!   among all \(G+k\) and processors
      !! @todo Change call to QE `mp_max` to actual call to `MPI_ALLREDUCE` here @endtodo


    npwx_g = maxval(ngk_g(1:nkstot_local))
      !! * Calculate the maximum number of G-vectors 
      !!   among all k-points

    return
  end subroutine reconstructFFTGrid

!----------------------------------------------------------------------------
  subroutine hpsort_eps(n, ra, ind, eps)
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! considering two elements equal if their difference is less than `eps`
    !!
    !! n is input, ra is replaced on output by its sorted rearrangement.
    !! Create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (ra).
    !! In case of equal values in the data array (ra) the values in the
    !! index array (ind) are used to order the entries.
    !!
    !! if on input ind(1)  = 0 then indices are initialized in the routine,
    !! if on input ind(1) != 0 then indices are assumed to have been
    !!                initialized before entering the routine and these
    !!                indices are carried around during the sorting process
    !!
    !!
    !! From QE code, adapted from Numerical Recipes pg. 329 (new edition)
    !!

    implicit none

    ! Input/Output variables:
    real(kind=dp), intent(in) :: eps
    integer, intent(in) :: n

    integer, intent(inout) :: ind(:)
    real(kind=dp), intent(inout) :: ra (:)


    ! Local variables
    integer :: i, ir, j, l, iind
    real(kind=dp) :: rra

    ! initialize index array
    if (ind (1) .eq.0) then
      do i = 1, n
        ind (i) = i
      enddo
    endif
    ! nothing to order
    if (n.lt.2) return
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1

    ir = n

    sorting: do

      ! still in hiring phase
      if ( l .gt. 1 ) then
        l    = l - 1
        rra  = ra (l)
        iind = ind (l)
        ! in retirement-promotion phase.
      else
        ! clear a space at the end of the array
        rra  = ra (ir)
        !
        iind = ind (ir)
        ! retire the top of the heap into it
        ra (ir) = ra (1)
        !
        ind (ir) = ind (1)
        ! decrease the size of the corporation
        ir = ir - 1
        ! done with the last promotion
        if ( ir .eq. 1 ) then
          ! the least competent worker at all !
          ra (1)  = rra
          !
          ind (1) = iind
          exit sorting
        endif
      endif
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      do while ( j .le. ir )
        if ( j .lt. ir ) then
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then
            if (ra(j).lt.ra(j+1)) j = j + 1
          else
            ! this means ra(j) == ra(j+1) within tolerance
            if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
        endif
        ! demote rra
        if ( abs(rra - ra(j)).ge.eps ) then
          if (rra.lt.ra(j)) then
            ra (i) = ra (j)
            ind (i) = ind (j)
            i = j
            j = j + j
          else
            ! set j to terminate do-while loop
            j = ir + 1
          end if
        else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind (j) ) then
            ra (i) = ra (j)
            ind (i) = ind (j)
            i = j
            j = j + j
          else
            ! set j to terminate do-while loop
            j = ir + 1
          endif
        end if
      enddo
      ra (i) = rra
      ind (i) = iind

    end do sorting

    return 
  end subroutine hpsort_eps

!----------------------------------------------------------------------------
  subroutine read_vasprun_xml(at_local, nkstot_local, VASPDir, wk_local, ityp, nat, nsp)
    !! Read the k-point weights and cell info from the `vasprun.xml` file
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    real(kind=dp), intent(in) :: at_local(3,3)
      !! Real space lattice vectors

    integer, intent(in) :: nkstot_local
      !! Total number of k-points

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files


    ! Output variables:
    real(kind=dp), allocatable, intent(out) :: wk_local(:)
      !! K-point weights

    integer, allocatable, intent(out) :: ityp(:)
      !! Atom type index
    integer, intent(out) :: nat
      !! Number of atoms
    integer, intent(out) :: nsp
      !! Number of types of atoms


    ! Local variables:
    real(kind=dp) :: dir(3)
      !! Direct coordinates read from file

    integer :: ik, ia, i, ix
      !! Loop indices

    character(len=256) :: cDum
      !! Dummy variable to ignore input
    character(len=256) :: fileName
      !! `vasprun.xml` with path
    character(len=256) :: line
      !! Line read from file

    logical :: fileExists
      !! If the `vasprun.xml` file exists
    logical :: found
      !! If the required tag was found

    allocate(wk_local(nkstot_local))

    if (ionode_local) then

      fileName = trim(VASPDir)//'/vasprun.xml'

      inquire(file = fileName, exist = fileExists)

      if (.not. fileExists) call exitError('getKPointWeights', 'Required file vasprun.xml does not exist', 1)

      open(57, file=fileName)
        !! * If root node, open `vasprun.xml`


      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'weights'`, indicating the
        !!   tag surrounding the k-point weights
        
        read(57, '(A)') line

        if (index(line,'weights') /= 0) found = .true.
        
      enddo

      do ik = 1, nkstot_local
        !! * Read in the weight for each k-point

        read(57,*) cDum, wk_local(ik), cDum

      enddo


      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'atominfo'`, indicating the
        !!   tag surrounding the cell info
        
        read(57, '(A)') line

        if (index(line,'atominfo') /= 0) found = .true.
        
      enddo

      read(57,*) cDum, nat, cDum
      read(57,*) cDum, nsp, cDum
      read(57,*) 
      read(57,*) 
      read(57,*) 
      read(57,*) 
      read(57,*) 

      allocate(ityp(nat))

      do ia = 1, nat
        !! * Read in the atom type index for each atom

        read(57,'(a21,i3,a9)') cDum, ityp(ia), cDum

      enddo


      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'finalpos'`, indicating the
        !!   tag surrounding the final cell parameters
        !!   and positions
        
        read(57, '(A)') line

        if (index(line,'finalpos') /= 0) found = .true.
        
      enddo

      found = .false.
      do while (.not. found)
        !! * Ignore everything until you get to a
        !!   line with `'positions'`, indicating the
        !!   tag surrounding the final positions
        
        read(57, '(A)') line

        if (index(line,'positions') /= 0) found = .true.
        
      enddo

      allocate(tau(3,nat))

      do ia = 1, nat
        !! * Read in the final position for each atom

        read(57,*) cDum, (dir(i),i=1,3), cDum
          !! @note
          !!  I assume that the coordinates are always direct
          !!  in the `vasprun.xml` file and that the scaling
          !!  factor is already included as I cannot find it 
          !!  listed anywhere in that file. Extensive testing
          !!  needs to be done to confirm this assumption.
          !! @endnote
          

        do ix = 1, 3
          tau(ix,ia) = sum(dir(:)*at_local(ix,:))
            !! @todo Test logic of direct to cartesian coordinates with scaling factor @endtodo
        enddo

      enddo

    endif

    call MPI_BCAST(wk_local, size(wk_local), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(nat, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nsp, 1, MPI_INTEGER, root, world_comm_local, ierr)

    if (.not. ionode_local) allocate(ityp(nat))
    if (.not. ionode_local) allocate(tau(3,nat))

    call MPI_BCAST(ityp, size(ityp), MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(tau, size(tau), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)

    return
  end subroutine read_vasprun_xml

!----------------------------------------------------------------------------
  subroutine readPOTCAR(nsp, VASPDir, ps)

    implicit none

    ! Input variables:
    integer, intent(in) :: nsp
      !! Number of types of atoms

    character(len=256), intent(in) :: VASPDir
      !! Directory with VASP files

    
    ! Output variables:
    type (pseudo) :: ps(nsp)


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
    integer :: ityp, i, j, ip, ir
      !! Loop indices
    integer :: nProj
      !! Number of projectors with given angular momentum

    character(len=1) :: charSwitch
      !! Switch to determine what section reading
    character(len=256) :: dummyC
      !! Dummy character to ignore input
    character(len=256) :: fileName
      !! Full WAVECAR file name including path

    logical :: found
      !! If the required tag was found


    if(ionode_local) then
      fileName = trim(VASPDir)//'/POTCAR'

      open(unit=potcarUnit, file=fileName, iostat=ierr, status='old')
      if (ierr .ne. 0) write(stdout,*) 'open error - iostat =', ierr
        !! * If root node, open the `POTCAR` file

      do ityp = 1, nsp

        ps(ityp)%nChannels = 0
        ps(ityp)%lmmax = 0

        read(potcarUnit,*) dummyC, ps(ityp)%element, dummyC
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

        found = .false.
        do while (.not. found)
          !! * Ignore all lines until you get to the `END` of
          !!   the PSCRT section
        
          read(potcarUnit, '(A)') dummyC

          if (dummyC(1:3) == 'END') found = .true.
        
        enddo

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

        read(potcarUnit,*) 
          !! * Ignore the max G for non-local potential and 
          !!   unused boolean (`LDUM` in VASP)
        read(potcarUnit,'(1X,A1)') charSwitch
          !! * Read character switch

        do while (charSwitch /= 'D' .and. charSwitch /= 'A' .and. charSwitch /= 'P' &
          .and. charSwitch /= 'E')
            !! * Until you have read in all of the momentum channels
            !!   (i.e. you get to a character switch that is not `'N'`)
            !!     * Read in the angular momentum and the number of 
            !!       projectors at this angular momentum
            !!     * Increment the number of nlm channels
            !!     * Ignore non-local strength multipliers
            !!     * Read in the reciprocal-space and real-space
            !!       projectors
            !!     * Increment the number of l channels
            !!     * Read the next character switch

          read(potcarUnit,*) angMom, nProj, dummyC
            ! Read in angular momentum and the number of projectors
            ! at this angular momentum

          ps(ityp)%lmmax = ps(ityp)%lmmax + (2*angMom+1)*nProj
            ! Increment the number of nlm channels

          allocate(dummyDA2(nProj,nProj))

          read(potcarUnit,*) dummyDA2(:,:)
            ! Ignore non-local strength multipliers

          do ip = 1, nProj
            ! Read in the reciprocal-space and real-space
            ! projectors

            ps(ityp)%angmom(ps(ityp)%nChannels+ip) = angmom

            read(potcarUnit,*) 
            read(potcarUnit,*) (ps(ityp)%recipProj(ps(ityp)%nChannels+ip,i), i=1,100)
            read(potcarUnit,*) 
            read(potcarUnit,*) (ps(ityp)%realProj(ps(ityp)%nChannels+ip,i), i=1,100)

          enddo

          ps(ityp)%nChannels = ps(ityp)%nChannels + nProj
            ! Increment the number of l channels

          deallocate(dummyDA2)

          read(potcarUnit,'(1X,A1)') charSwitch
            ! Read character switch

        enddo

        if (charSwitch /= 'P') then
          !! * Ignore depletion charges

          read(potcarUnit,*)
          read(potcarUnit,*)
          
        else

          read(potcarUnit,*) ps(ityp)%nmax, ps(ityp)%rAugMax  
            !! * Read the number of mesh grid points and
            !!   the maximum radius in the augmentation sphere
          read(potcarUnit,*)
            !! * Ignore format specifier
          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch (not used)

          allocate(ps(ityp)%radGrid(ps(ityp)%nmax))
          allocate(ps(ityp)%wps(ps(ityp)%nChannels,ps(ityp)%nmax))
          allocate(ps(ityp)%wae(ps(ityp)%nChannels,ps(ityp)%nmax))
          allocate(dummyDA2(ps(ityp)%nChannels, ps(ityp)%nChannels))

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

          read(potcarUnit,*) (ps(ityp)%radGrid(i), i=1,ps(ityp)%nmax)

          H = log(ps(ityp)%radGrid(ps(ityp)%nmax)/ps(ityp)%radGrid(1))/(ps(ityp)%nmax - 1)
            !! * Calculate \(H\) which is used to generate the derivative of the grid
            !! @note
            !!  The grid in VASP is defined as \(R_i = R_0e^{H(i-1)}\), so we define the
            !!  derivative as \(dR_i = R_0He^{H(i-1)}\)
            !! @endnote
          
          do ir = 1, ps(ityp)%nmax
            !! * Calculate the max index of the augmentation sphere and
            !!   the derivative of the radial grid

            if (ps(ityp)%radGrid(ir) > ps(ityp)%rAugMax) ps(ityp)%iRAugMax = ir - 1
            ps(ityp)%dRadGrid(ir) = ps(ityp)%radGrid(1)*H*exp(H*(ir-1))

          enddo

          read(potcarUnit,'(1X,A1)') charSwitch
            !! * Read character switch
          if (charSwitch /= 'a') call exitError('readPOTCAR', 'expected aepotential section', 1)

          allocate(dummyDA1(ps(ityp)%nmax))

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

          do ip = 1, ps(ityp)%nChannels
            
            read(potcarUnit,'(1X,A1)') charSwitch
            if (charSwitch /= 'p') call exitError('readPOTCAR', 'expected pseudowavefunction section', 1)
            read(potcarUnit,*) (ps(ityp)%wps(ip,i), i=1,ps(ityp)%nmax)

            read(potcarUnit,'(1X,A1)') charSwitch
            if (charSwitch /= 'a') call exitError('readPOTCAR', 'expected aewavefunction section', 1)
            read(potcarUnit,*) (ps(ityp)%wae(ip,i), i=1,ps(ityp)%nmax)

          enddo

          deallocate(dummyDA1)
          deallocate(dummyDA2)

        endif

        found = .false.
        do while (.not. found)
          !! * Ignore all lines until you get to the `END` of
          !!   the PSCRT section
        
          read(potcarUnit, '(A)') dummyC

          if (index(dummyC,'End of Dataset') /= 0) found = .true.
        
        enddo

      enddo

    endif    

    return
  end subroutine readPOTCAR

!----------------------------------------------------------------------------
  subroutine writeKInfo(nkstot_local, npwx_local, igk_l2g, nbnd_local, ngk_g, ngk_local, &
      npw_g, npwx_g, occ, wk_local, xk_local, igwk)
    !! Calculate the highest occupied band for each k-point,
    !! gather the \(G+k\) vector indices in single, global 
    !! array, and write out k-point information
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nkstot_local
      !! Total number of k-points
    integer, intent(in) :: npwx_local
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    integer, intent(in) :: igk_l2g(npwx_local, nk_Pool)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point;
      !! the first index goes up to `npwx_local`,
      !! but only valid values are up to `ngk_local`
    integer, intent(in) :: nbnd_local
      !! Total number of bands
    integer, intent(in) :: ngk_g(nkstot_local)
      !! Global number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each k-point
    integer, intent(in) :: ngk_local(nk_Pool)
      !! Number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each
      !! k-point, on this processor
    integer, intent(in) :: npw_g
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(in) :: npwx_g
      !! Max number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` among all k-points

    real(kind=dp), intent(in) :: occ(nbnd_local, nkstot_local)
      !! Occupation of band
    real(kind=dp), intent(in) :: wk_local(nkstot_local)
      !! K-point weights
    real(kind=dp), intent(in) :: xk_local(3,nkstot_local)
      !! Position of k-points in reciprocal space


    ! Output variables:
    integer, allocatable, intent(out) :: igwk(:,:)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors


    ! Local variables:
    integer, allocatable :: groundState(:)
      !! Holds the highest occupied band
      !! for each k-point
    integer :: ik, ig
      !! Loop indices


    if(ionode_local) then

      write(stdout,*)
      write(stdout,*) "***************"
      write(stdout,*) "Getting ground state bands"
    
      write(mainout, '("# Number of K-points. Format: ''(i10)''")')
      write(mainout, '(i10)') nkstot_local
      write(mainout, '("# ik, groundState, ngk_g(ik), wk(ik), xk(1:3,ik). Format: ''(3i10,4ES24.15E3)''")')
      flush(mainout)
    
      allocate(groundState(nkstot_local))

      call getGroundState(nbnd_local, nkstot_local, occ, groundState)
        !! * For each k-point, find the index of the 
        !!   highest occupied band
        !!
        !! @note
        !!  Although `groundState` is written out in `Export`,
        !!  it is not currently used by the `TME` program.
        !! @endtodo

      write(stdout,*) "Done getting ground state bands"
      write(stdout,*) "***************"
      write(stdout,*)

    endif

    if(ionode_local) then

      write(stdout,*)
      write(stdout,*) "***************"
      write(stdout,*) "Getting global G+k indices"

    endif
  
    allocate(igwk(npwx_g, nkstot_local))
  
    igwk(:,:) = 0
    do ik = 1, nkstot_local

      if (ionode_local) write(stdout,*) "Processing k-point ", ik

      call getGlobalGkIndices(nkstot_local, npwx_local, igk_l2g, ik, ngk_g, ngk_local, npw_g, &
          npwx_g, igwk)
        !! * For each k-point, gather all of the \(G+k\) indices
        !!   among all processors in a single global array
    
      if (ionode_local) write(mainout, '(3i10,4ES24.15E3)') ik, groundState(ik), ngk_g(ik), wk_local(ik), xk_local(1:3,ik)
      if (ionode_local) flush(mainout)
    
    enddo

    if(ionode_local) then

      write(stdout,*) "Done getting global G+k indices"
      write(stdout,*) "***************"
      write(stdout,*)
      flush(stdout)

    endif

    if (ionode_local) deallocate(groundState)

    return
  end subroutine writeKInfo

!----------------------------------------------------------------------------
  subroutine getGroundState(nbnd_local, nkstot_local, occ, groundState)
    !! * For each k-point, find the index of the 
    !!   highest occupied band

    implicit none

    ! Input variables:
    integer, intent(in) :: nbnd_local
      !! Total number of bands
    integer, intent(in) :: nkstot_local
      !! Total number of k-points

    real(kind=dp), intent(in) :: occ(nbnd_local, nkstot_local)
      !! Occupation of band

    
    ! Output variables:
    integer, intent(out) :: groundState(nkstot_local)
      !! Holds the highest occupied band
      !! for each k-point


    ! Local variables:
    integer :: ik, ibnd
      !! Loop indices


    groundState(:) = 0
    do ik = 1, nkstot_local

      do ibnd = 1, nbnd_local

        if (occ(ibnd,ik) < 0.5_dp) then
          !! @todo Figure out if boundary for "occupied" should be 0.5 or less @endtodo
        !if (et(ibnd,ik) > ef) then

          groundState(ik) = ibnd - 1
          goto 10

        endif
      enddo

10    continue

    enddo

    return
  end subroutine getGroundState

!----------------------------------------------------------------------------
  subroutine getGlobalGkIndices(nkstot_local, npwx_local, igk_l2g, ik, ngk_g, ngk_local, npw_g, &
      npwx_g, igwk)
    !! Gather the \(G+k\) vector indices in single, global 
    !! array
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: nkstot_local
      !! Total number of k-points
    integer, intent(in) :: npwx_local
      !! Maximum number of \(G+k\) vectors
      !! across all k-points for just this 
      !! processor

    integer, intent(in) :: igk_l2g(npwx_local, nk_Pool)
      !! Local to global indices for \(G+k\) vectors 
      !! ordered by magnitude at a given k-point;
      !! the first index goes up to `npwx_local`,
      !! but only valid values are up to `ngk_local`
    integer, intent(in) :: ik
      !! Index of current k-point
    integer, intent(in) :: ngk_g(nkstot_local)
      !! Global number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each k-point
    integer, intent(in) :: ngk_local(nk_Pool)
      !! Number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each
      !! k-point, on this processor
    integer, intent(in) :: npw_g
      !! Maximum G-vector index among all \(G+k\)
      !! and processors
    integer, intent(in) :: npwx_g
      !! Max number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` among all k-points


    ! Output variables:
    integer, intent(out) :: igwk(npwx_g, nkstot_local)
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
      !! k-point; should equal `ngk_g`

    
    allocate(itmp1(npw_g), stat=ierr)
    if (ierr/= 0) call exitError('getGlobalGkIndices','allocating itmp1', abs(ierr))

    itmp1 = 0
    if(ik >= ikStart .and. ik <= ikEnd) then

      do ig = 1, ngk_local(ik-ikStart+1)

        itmp1(igk_l2g(ig, ik-ikStart+1)) = igk_l2g(ig, ik-ikStart+1)
          !! * For each k-point and \(G+k\) vector for this processor,
          !!   store the local to global indices (`igk_l2g`) in an
          !!   array that will later be combined globally
          !!
          !! @note
          !!  This will leave zeros in spots where the \(G+k\) 
          !!  combination for this k-point was greater than the energy 
          !!  cutoff.
          !! @endnote

      enddo
    endif

    call mpiSumIntV(itmp1, world_comm_local)

    ngg = 0
    do  ig = 1, npw_g

      if(itmp1(ig) == ig) then
        !! * Go through and find all of the non-zero
        !!   indices in the now-global `itmp1` array,
        !!   and store them in a new array that won't
        !!   have the extra zeros

        ngg = ngg + 1

        igwk(ngg, ik) = ig

      endif
    enddo


    if(ionode_local .and. ngg /= ngk_g(ik)) call exitError('writeKInfo', 'Unexpected number of G+k vectors', 1)
      !! * Make sure that the total number of non-zero
      !!   indices matches the global number of \(G+k\)
      !!   vectors for this k-point
    
    deallocate( itmp1 )

    return
  end subroutine getGlobalGkIndices

!----------------------------------------------------------------------------
  subroutine writeGridInfo(ngm_g_local, nkstot_local, npwx_g, igwk, mill_g, ngk_g, npw_g, exportDir)
    !! Write out grid boundaries and miller indices
    !! for just \(G+k\) combinations below cutoff energy
    !! in one file and all miller indices in another 
    !! file
    !!
    !! <h2>Walkthrough</h2>
    !!

    implicit none

    ! Input variables:
    integer, intent(in) :: ngm_g_local
      !! Global number of G-vectors
    integer, intent(in) :: nkstot_local
      !! Total number of k-points
    integer, intent(in) :: npwx_g
      !! Max number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` among all k-points

    integer, intent(in) :: igwk(npwx_g, nkstot_local)
      !! Indices of \(G+k\) vectors for each k-point
      !! and all processors
    integer, intent(in) :: mill_g(3,ngm_g_local)
      !! Integer coefficients for G-vectors on all processors
    integer, intent(in) :: ngk_g(nkstot_local)
      !! Global number of \(G+k\) vectors with energy
      !! less than `ecutwfc_local` for each k-point
    integer, intent(in) :: npw_g
      !! Maximum G-vector index among all \(G+k\)
      !! and processors

    character(len=256), intent(in) :: exportDir
      !! Directory to be used for export


    ! Output variables:


    ! Local variables:
    integer :: ik, ig, igk
      !! Loop indices


    if (ionode_local) then
    
      write(mainout, '("# Number of G-vectors. Format: ''(i10)''")')
      write(mainout, '(i10)') ngm_g_local
    
      write(mainout, '("# Number of PW-vectors. Format: ''(i10)''")')
      write(mainout, '(i10)') npw_g
    
      write(mainout, '("# Number of min - max values of fft grid in x, y and z axis. Format: ''(6i10)''")')
      write(mainout, '(6i10)') minval(mill_g(1,1:ngm_g_local)), maxval(mill_g(1,1:ngm_g_local)), &
                          minval(mill_g(2,1:ngm_g_local)), maxval(mill_g(2,1:ngm_g_local)), &
                          minval(mill_g(3,1:ngm_g_local)), maxval(mill_g(3,1:ngm_g_local))
      flush(mainout)
    
      do ik = 1, nkstot_local
        !! * For each k-point, write out the miller indices
        !!   resulting in \(G+k\) vectors less than the energy
        !!   cutoff in a `grid.ik` file
      
        open(72, file=trim(exportDir)//"/grid"//iotk_index(ik))
          !! @todo Move `iotk_index` here #thisbranch @endtodo
        write(72, '("# Wave function G-vectors grid")')
        write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
      
        do igk = 1, ngk_g(ik)
          write(72, '(4i10)') igwk(igk,ik), mill_g(1:3,igwk(igk,ik))
          flush(72)
        enddo
      
        close(72)
      
      enddo

      !> * Output all miller indices
      open(72, file=trim(exportDir)//"/mgrid")
      write(72, '("# Full G-vectors grid")')
      write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
    
      do ig = 1, ngm_g_local
        write(72, '(4i10)') ig, mill_g(1:3,ig)
        flush(72)
      enddo
    
      close(72)

    endif

    return
  end subroutine writeGridInfo


!----------------------------------------------------------------------------
  subroutine writeCellInfo(ityp, nat, nbnd_local, nsp, nspin_local, at_local, bg_local, tau, nnTyp)
    implicit none

    ! Input variables:
    integer, intent(in) :: ityp(nat)
      !! Atom type index
    integer, intent(in) :: nat
      !! Number of atoms
    integer, intent(in) :: nbnd_local
      !! Total number of bands
    integer, intent(in) :: nsp
      !! Number of types of atoms
    integer, intent(in) :: nspin_local
      !! Number of spins

    real(kind=dp), intent(in) :: at_local(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(in) :: bg_local(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: tau(3,nat)
      !! Atom positions


    ! Output variables:
    integer, allocatable, intent(out) :: nnTyp(:)
      !! Number of atoms of each type


    ! Local variables:
    integer :: i
      !! Loop index


    if (ionode_local) then
    
      write(mainout, '("# Cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
      write(mainout, '("# a1 ",3ES24.15E3)') at_local(:,1)
      write(mainout, '("# a2 ",3ES24.15E3)') at_local(:,2)
      write(mainout, '("# a3 ",3ES24.15E3)') at_local(:,3)
    
      write(mainout, '("# Reciprocal cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
      write(mainout, '("# b1 ",3ES24.15E3)') bg_local(:,1)
      write(mainout, '("# b2 ",3ES24.15E3)') bg_local(:,2)
      write(mainout, '("# b3 ",3ES24.15E3)') bg_local(:,3)
    
      write(mainout, '("# Number of Atoms. Format: ''(i10)''")')
      write(mainout, '(i10)') nat
    
      write(mainout, '("# Number of Types. Format: ''(i10)''")')
      write(mainout, '(i10)') nsp
    
      write(mainout, '("# Atoms type, position(1:3) (a.u.). Format: ''(i10,3ES24.15E3)''")')
      do i = 1, nat
        write(mainout,'(i10,3ES24.15E3)') ityp(i), tau(:,i)
      enddo
    
      write(mainout, '("# Number of Bands. Format: ''(i10)''")')
      write(mainout, '(i10)') nbnd_local

      write(mainout, '("# Spin. Format: ''(i10)''")')
      write(mainout, '(i10)') nspin_local
    
      allocate( nnTyp(nsp) )
      nnTyp = 0
      do i = 1, nat
        nnTyp(ityp(i)) = nnTyp(ityp(i)) + 1
      enddo

    endif

    return
  end subroutine writeCellInfo

!----------------------------------------------------------------------------
  subroutine writePseudoInfo(nsp, nnTyp, ps)

    implicit none

    ! Input variables:
    integer, intent(in) :: nsp
      !! Number of types of atoms

    integer, intent(in) :: nnTyp(nsp)
      !! Number of atoms of each type

    type (pseudo) :: ps(nsp)


    ! Output variables:


    ! Local variables:
    integer :: ityp, ip, ir
      !! Loop index

  
    if (ionode_local) then

      do ityp = 1, nsp
        
        write(mainout, '("# Element")')
        write(mainout, *) trim(ps(ityp)%element)
        write(mainout, '("# Number of Atoms of this type. Format: ''(i10)''")')
        write(mainout, '(i10)') nnTyp(ityp)
        write(mainout, '("# Number of projectors. Format: ''(i10)''")')
        write(mainout, '(i10)') ps(ityp)%nChannels
        
        write(mainout, '("# Angular momentum, index of the projectors. Format: ''(2i10)''")')
        do ip = 1, ps(ityp)%nChannels

          write(mainout, '(2i10)') ps(ityp)%angmom(ip), ip

        enddo
        
        write(mainout, '("# Number of channels. Format: ''(i10)''")')
        write(mainout, '(i10)') ps(ityp)%lmmax
        
        write(mainout, '("# Number of radial mesh points. Format: ''(2i10)''")')
        write(mainout, '(2i10)') ps(ityp)%nmax, ps(ityp)%iRAugMax
          ! Number of points in the radial mesh, number of points inside the aug sphere
        
        write(mainout, '("# Radial grid, Integratable grid. Format: ''(2ES24.15E3)''")')
        do ir = 1, ps(ityp)%nmax
          write(mainout, '(2ES24.15E3)') ps(ityp)%radGrid(ir), ps(ityp)%dRadGrid(ir) 
            ! Radial grid, derivative of radial grid
        enddo
        
        write(mainout, '("# AE, PS radial wfc for each beta function. Format: ''(2ES24.15E3)''")')
        do ip = 1, ps(ityp)%nChannels
          do ir = 1, ps(ityp)%nmax
            write(mainout, '(2ES24.15E3)') ps(ityp)%wae(ip,ir), ps(ityp)%wps(ip,ir)
          enddo
        enddo
      
      enddo
    
    endif

    return
  end subroutine writePseudoInfo

!----------------------------------------------------------------------------
  subroutine subroutineTemplate()
    implicit none


    return
  end subroutine subroutineTemplate

!----------------------------------------------------------------------------
! ..  This subroutine write wavefunctions to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
    !CALL mp_max( npw_g, world_comm_local )

!
    SUBROUTINE write_restart_wfc(iuni, exportDir, &
      ik, nk, ispin, nspin_local, scal, wf0, t0, wfm, tm, ngw, gamma_only, nbnd_local, igl, ngwl )
!
!
      USE mp,             ONLY : mp_sum
      IMPLICIT NONE
!
      INTEGER, INTENT(in) :: iuni
      character(len = 256), intent(in) :: exportDir
      INTEGER, INTENT(in) :: ik, nk, ispin, nspin_local
      COMPLEX(DP), INTENT(in) :: wf0(:,:)
      COMPLEX(DP), INTENT(in) :: wfm(:,:)
      INTEGER, INTENT(in) :: ngw   !
      LOGICAL, INTENT(in) :: gamma_only
      INTEGER, INTENT(in) :: nbnd_local
      INTEGER, INTENT(in) :: ngwl
      INTEGER, INTENT(in) :: igl(:)
      REAL(DP), INTENT(in) :: scal
      LOGICAL, INTENT(in) :: t0, tm

      INTEGER :: i, j, ierr, idum = 0
      INTEGER :: nkt, ikt, igwx, ig
      INTEGER :: ipmask( nproc_local ), ipsour
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER, ALLOCATABLE :: igltot(:)

      CHARACTER(len=20) :: section_name = 'wfc'

      LOGICAL :: twrite = .true.

      INTEGER :: ierr_iotk
      CHARACTER(len=iotk_attlenx) :: attr

!
! ... Subroutine Body
!

        ! set working variables for k point index (ikt) and k points number (nkt)
        ikt = ik
        nkt = nk

        ipmask = 0
        ipsour = root

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool_local > 1 ) THEN
          IF( ( ikt >= ikStart ) .and. ( ikt <= ikEnd ) ) THEN
            IF( indexInPool == root_pool ) ipmask( myid + 1 ) = 1
          ENDIF
          CALL mp_sum( ipmask, world_comm_local )
          DO i = 1, nproc_local
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
          ENDDO
        ENDIF

        igwx = 0
        ierr = 0
        IF( ( ikt >= ikStart ) .and. ( ikt <= ikEnd ) ) THEN
          IF( ngwl > size( igl ) ) THEN
            ierr = 1
          ELSE
            igwx = maxval( igl(1:ngwl) )
          ENDIF
        ENDIF

        ! get the maximum index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm_local )

        ! now notify all procs if an error has been found
        !
        CALL mp_max( ierr, world_comm_local )

        IF( ierr > 0 ) &
          CALL exitError(' write_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipsour /= root ) THEN
          CALL mp_get( igwx, igwx, myid, root, ipsour, 1, world_comm_local )
        ENDIF

        ALLOCATE( wtmp( max(igwx,1) ) )
        wtmp = cmplx(0.0_dp, 0.0_dp, kind=dp)

        DO j = 1, nbnd_local
          IF( t0 ) THEN
            IF( npool_local > 1 ) THEN
              IF( ( ikt >= ikStart ) .and. ( ikt <= ikEnd ) ) THEN
                CALL mergewf(wf0(:,j), wtmp, ngwl, igl, indexInPool, &
                             nproc_pool_local, root_pool, intra_pool_comm_local)
              ENDIF
              IF( ipsour /= root ) THEN
                CALL mp_get( wtmp, wtmp, myid, root, ipsour, j, world_comm_local )
              ENDIF
            ELSE
              CALL mergewf(wf0(:,j), wtmp, ngwl, igl, myid, nproc_local, &
                           root, world_comm_local )
            ENDIF

            IF( ionode_local ) THEN
              do ig = 1, igwx
                write(iuni, '(2ES24.15E3)') wtmp(ig)
              enddo
              !
!              do j = 1, nbnd_local
!                do i = 1, igwx ! ngk_g(ik)
!                  write(74,'(2ES24.15E3)') wf0(i,j) ! wf0 is the local array for evc(i,j)
!                enddo
!              enddo
              !
            ENDIF
          ELSE
          ENDIF
        ENDDO

!        DO j = 1, nbnd_local
!          IF( tm ) THEN
!            IF( npool_local > 1 ) THEN
!              IF( ( ikt >= ikStart ) .and. ( ikt <= ikEnd ) ) THEN
!                CALL mergewf(wfm(:,j), wtmp, ngwl, igl, indexInPool, &
!                             nproc_pool_local, root_pool, intra_pool_comm_local)
!              ENDIF
!              IF( ipsour /= root ) THEN
!                CALL mp_get( wtmp, wtmp, myid, root, ipsour, j, world_comm_local )
!              ENDIF
!            ELSE
!              CALL mergewf(wfm(:,j), wtmp, ngwl, igl, myid, nproc_local, root, world_comm_local )
!            ENDIF
!            IF( ionode_local ) THEN
!              CALL iotk_write_dat(iuni,"Wfcm"//iotk_index(j),wtmp(1:igwx))
!            ENDIF
!          ELSE
!          ENDIF
!        ENDDO
        IF(ionode_local) then
          close(iuni)
          !CALL iotk_write_end  (iuni,"Kpoint"//iotk_index(ik))
        endif
      
        DEALLOCATE( wtmp )

      RETURN
    END SUBROUTINE

  SUBROUTINE write_export (mainOutputFile, exportDir)
    !-----------------------------------------------------------------------
    !
    USE iotk_module

    use gvect, only : g, ngm, ngm_g
    use klist, only : nks, xk, ngk, wk
    use cell_base, only : tpiba2, alat

    USE kinds,          ONLY : DP
    USE start_k,        ONLY : nk1, nk2, nk3, k1, k2, k3
    USE control_flags,  ONLY : gamma_only
    USE global_version, ONLY : version_number
    USE becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type

    USE uspp,          ONLY : nkb, vkb
    USE wavefunctions_module,  ONLY : evc
    USE io_files,       ONLY : outdir, prefix, iunwfc, nwordwfc
    USE io_files,       ONLY : psfile
    USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
    USE mp,             ONLY : mp_sum, mp_max
  
    USE upf_module,     ONLY : read_upf
  
    USE pseudo_types, ONLY : pseudo_upf
    USE radial_grids, ONLY : radial_grid_type
    
    USE wvfct,         ONLY : wg, npw, g2kin
  
    USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
    USE paw_onecenter,        ONLY : PAW_potential
    USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
    USE uspp_param,           ONLY : nh, nhm ! used for PAW
    USE uspp,                 ONLY : qq_so, dvan_so, qq, dvan
    USE scf,                  ONLY : rho

    IMPLICIT NONE
  
    CHARACTER(5), PARAMETER :: fmt_name="QEXPT"
    CHARACTER(5), PARAMETER :: fmt_version="1.1.0"

    CHARACTER(256), INTENT(in) :: mainOutputFile, exportDir

    INTEGER :: i, j, k, ig, ik, ibnd, na, ngg,ig_, ierr
    real(DP) :: xyz(3), tmp(3)
    INTEGER :: im, ink, inb, ms
    INTEGER :: ispin, local_pw
    INTEGER, ALLOCATABLE :: l2g_new( : )
  
    character(len = 300) :: text
  

    real(DP) :: wfc_scal
    LOGICAL :: twf0, twfm, file_exists
    CHARACTER(iotk_attlenx) :: attr
    TYPE(pseudo_upf) :: upf       ! the pseudo data
    TYPE(radial_grid_type) :: grid

    integer, allocatable :: nnTyp(:)

    integer, allocatable :: igk(:)
      !! Index map from \(G\) to \(G+k\)
      !! indexed up to `npwx_local`


      !------------------------------------------------------------------------------------

#ifdef __MPI
  CALL poolrecover (et, nbnd_local, nkstot_local, nk_Pool)
#endif

      !------------------------------------------------------------------------------------

    WRITE(stdout,*) "Writing Eigenvalues"

    IF( ionode_local ) THEN
    
      write(mainout, '("# Fermi Energy (Hartree). Format: ''(ES24.15E3)''")')
      write(mainout, '(ES24.15E3)') ef*ryToHartree
      flush(mainout)
    
      DO ik = 1, nkstot_local
      
        ispin = isk( ik )
      
        open(72, file=trim(exportDir)//"/eigenvalues"//iotk_index(ik))
      
        write(72, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
        write(72, '("# Eigenvalues (Hartree), band occupation number. Format: ''(2ES24.15E3)''")')
      
        do ibnd = 1, nbnd_local
          if ( wk(ik) == 0.D0 ) then
              write(72, '(2ES24.15E3)') et(ibnd,ik)*ryToHartree, wg(ibnd,ik)
           else
            write(72, '(2ES24.15E3)') et(ibnd,ik)*ryToHartree, wg(ibnd,ik)/wk(ik)
          endif
        enddo
      
        close(72)
      
      ENDDO
    
    endif

      !------------------------------------------------------------------------------------
  
    if ( ionode_local ) WRITE(stdout,*) "Writing Wavefunctions"
  
    wfc_scal = 1.0d0
    twf0 = .true.
    twfm = .false.
  
    IF ( nkb > 0 ) THEN
    
      CALL init_us_1
      CALL init_at_1
    
      CALL allocate_bec_type (nkb,nbnd_local, becp)

      allocate(igk(npwx_local))
      
      igk = 0
    
      DO ik = 1, nkstot_local
      
        local_pw = 0
        IF ( (ik >= ikStart) .and. (ik <= ikEnd) ) THEN
          CALL gk_sort (xk_local(1, ik+ikStart-1), ngm_local, g, vcut_local, npw, igk, g2kin)
            !! @note
            !!  The call to `gk_sort` here returns values for `igk` and `g2kin` 
            !! @endnote
            !! @note 
            !!  I am assuming that this call to `gk_sort` will not actually be needed. I think 
            !!  I should be able to use values already calculated instead of calling `gk_sort`.
            !!  If I'm wrong and this turns out to be needed, `g` should be `gCart_local` here 
            !! @endnote
            !! @note 
            !!  I changed the way the `vcut_local` was calculated, so this will not give the
            !!  expected result. If this is ultimately needed, will need to fix what is passed
            !!  in
            !! @endnote
          CALL davcio (evc, nwordwfc, iunwfc, (ik-ikStart+1), - 1)

          igk(1:ngk_local(ik-ikStart+1)) = igk_large(ik-ikStart+1,1:ngk_local(ik-ikStart+1))

          CALL init_us_2(npw, igk, xk_local(1,ik), vkb)
          local_pw = ngk(ik-ikStart+1)

          IF ( gamma_only ) THEN
            CALL calbec ( ngk_g(ik), vkb, evc, becp )
            WRITE(0,*) 'Gamma only PW_EXPORT not yet tested'
          ELSE
            CALL calbec ( npw, vkb, evc, becp )
            if ( ionode_local ) then

              WRITE(stdout,*) "Writing projectors of kpt", ik

              file_exists = .false.
              inquire(file =trim(exportDir)//"/projections"//iotk_index(ik), exist = file_exists)
              if ( .not. file_exists ) then
                open(72, file=trim(exportDir)//"/projections"//iotk_index(ik))
                write(72, '("# Complex projections <beta|psi>. Format: ''(2ES24.15E3)''")')
                do j = 1,  becp%nbnd ! number of bands
                  do i = 1, nkb      ! number of projections
                    write(72,'(2ES24.15E3)') becp%k(i,j)
                  enddo
                enddo
              
                close(72)
              
              endif
            endif
          ENDIF
        ENDIF

        ALLOCATE(l2g_new(local_pw))

        l2g_new = 0
        DO ig = 1, local_pw
          ngg = igk_l2g(ig,ik-ikStart+1)
          DO ig_ = 1, ngk_g(ik)
            IF(ngg == igwk(ig_,ik)) THEN
              l2g_new(ig) = ig_
              exit
            ENDIF
          ENDDO
        ENDDO
        
        ispin = isk( ik )
        
        if ( ionode_local ) then

          file_exists = .false.
          inquire(file =trim(exportDir)//"/wfc"//iotk_index(ik), exist = file_exists)
          if ( .not. file_exists ) then
            open (72, file=trim(exportDir)//"/wfc"//iotk_index(ik))
            write(72, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
            write(72, '("# Complex : wavefunction coefficients (a.u.)^(-3/2). Format: ''(2ES24.15E3)''")')
            
            open(73, file=trim(exportDir)//"/projectors"//iotk_index(ik))
            write(73, '("# Complex projectors |beta>. Format: ''(2ES24.15E3)''")')
            write(73,'(2i10)') nkb, ngk_g(ik)
!            WRITE(stdout,*) "Writing Wavefunctions of kpt", ik
!            open(74, file=trim(exportDir)//"/evc"//iotk_index(ik))
!            write(74, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
!            write(74, '("# Complex : wavefunction coefficients (a.u.)^(-3/2). Format: ''(2ES24.15E3)''")')
          endif
        endif
        
        call MPI_BCAST(file_exists, 1, MPI_LOGICAL, root, world_comm_local, ierr)
        
        if ( .not. file_exists ) then
          CALL write_restart_wfc(72, exportDir, ik, nkstot_local, ispin, nspin_local, &
                                 wfc_scal, evc, twf0, evc, twfm, npw_g, gamma_only, nbnd_local, &
                                 l2g_new(:),local_pw )
          CALL write_restart_wfc(73, exportDir, ik, nkstot_local, ispin, nspin_local, &
                                 wfc_scal, vkb, twf0, evc, twfm, npw_g, gamma_only, nkb, &
                                 l2g_new(:), local_pw )
        endif
      
        if ( .not. file_exists .and. ionode_local ) then
          close(72)
          close(73)
!          close(74)
        endif
      
        DEALLOCATE(l2g_new)
      ENDDO

      deallocate(igk)
    
      CALL deallocate_bec_type ( becp )
    
    ENDIF

      !------------------------------------------------------------------------------------

    deallocate(xk_local)
    deallocate(igk_large)
    deallocate(igk_l2g)
    DEALLOCATE( igwk )
    deallocate(ngk_g)
END SUBROUTINE write_export

end module wfcExportVASPMod
