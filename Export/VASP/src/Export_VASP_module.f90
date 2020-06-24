module wfcExportVASPMod

  USE wrappers,      ONLY : f_mkdir_safe

  !USE pwcom
  USE constants, ONLY : e2, rytoev, pi, tpi, fpi
  USE cell_base, ONLY : celldm, ibrav
  USE klist, ONLY : wk
  USE ener, ONLY : ef
  USE wvfct, ONLY : igk, et
  USE lsda_mod, ONLY : isk

  USE io_files,  ONLY : prefix, outdir, tmp_dir
  USE ions_base, ONLY : ntype => nsp
  USE iotk_module
  use mpi
  USE mp,        ONLY: mp_sum, mp_max, mp_get, mp_bcast, mp_rank
  USE mp_wave, ONLY : mergewf

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision
  integer, parameter :: root = 0
    !! ID of the root node
  integer, parameter :: root_pool = 0
    !! Index of the root process within each pool
  integer, parameter :: mainout = 50
    !! Main output file unit
  integer, parameter :: stdout = 6
    !! Standard output unit
  integer, parameter :: wavecarUnit = 86
    !! WAVECAR unit for I/O

  real(kind = dp), parameter :: eVToRy = 0.073498618_dp
    !! Conversion factor from eV Rydberg
  real(kind = dp), parameter :: ryToHartree = 0.5_dp
    !! Conversion factor from Rydberg to Hartree

  real(kind=dp) :: at_local(3,3)
    !! Real space lattice vectors
    !! @todo Change back to `at` once extracted from QE #end @endtodo
  real(kind=dp) :: bg_local(3,3)
    !! Reciprocal lattice vectors
    !! @todo Change back to `bg` once extracted from QE #end @endtodo
  real(kind=dp) :: ecutwfc_local
    !! Plane wave energy cutoff
    !! @todo Change back to `ecutwfc` once extracted from QE #end @endtodo
  real(kind=dp) :: omega_local
    !! Volume of unit cell
    !! @todo Change back to `omega` once extracted from QE #end @endtodo
  real(kind=dp) :: tStart
    !! Start time
  
  INTEGER :: ik, i
    !! @todo Move these variables to be local #thisbranch @endtodo
  integer :: ierr
    !! Error returned by MPI
  integer :: ikEnd
    !! Ending index for kpoints in single pool 
  integer :: ikStart
    !! Starting index for kpoints in single pool 
  integer :: ios
    !! Error for input/output
  integer :: indexInPool
    !! Process index within pool
  integer :: inter_pool_comm_local = 0
    !! Inter pool communicator
    !! @todo Change back to `inter_pool_comm` once extracted from QE #end @endtodo
  integer :: intra_pool_comm_local = 0
    !! Intra pool communicator
    !! @todo Change back to `intra_pool_comm` once extracted from QE #end @endtodo
  integer :: myPoolId
    !! Pool index for this process
  integer :: nbnd_local
    !! Total number of bands
    !! @todo Change back to `nbnd` once extracted from QE #end @endtodo
  integer :: nkstot_local
    !! Total number of kpoints
    !! @todo Change back to `nkstot` once extracted from QE #end @endtodo
  integer :: npool_local = 1
    !! Number of pools for kpoint parallelization
    !! @todo Change back to `npool` once extracted from QE #end @endtodo
  integer :: nproc_local
    !! Number of processes
    !! @todo Change back to `nproc` once extracted from QE #end @endtodo
  integer :: nproc_pool_local
    !! Number of processes per pool
    !! @todo Change back to `nproc_pool` once extracted from QE #end @endtodo
  integer :: npw_g
    !! ??Not sure what this is
  integer :: npwx_g
    !! ??Not sure what this is
  integer :: nspin_local
    !! Number of spins
    !! @todo Change back to `nspin` once extracted from QE #end @endtodo
  integer :: myid
    !! ID of this process
  integer :: world_comm_local
    !! World communicator
    !! @todo Change back to `world_comm` once extracted from QE #end @endtodo
  
  character(len=256) :: exportDir
    !! Directory to be used for export
  character(len=256) :: mainOutputFile
    !! Main output file
  character(len=256) :: QEDir
    !! Directory with QE files
  character(len=256) :: VASPDir
    !! Directory with VASP files

  logical :: ionode_local
    !! If this node is the root node
    !! @todo Change back to `ionode` once extracted from QE #end @endtodo

  real(kind=dp), allocatable :: xk_local(:,:)

  integer, allocatable :: igall(:,:)
    !! Integer coefficients for G-vectors
  integer, allocatable :: igk_l2g( :, : )
    !! ??Not sure what this is
  integer, allocatable :: itmp_g( :, : )
    !! ??Not sure what this is
  integer, allocatable :: ngk_g( : )
    !! ??Not sure what this is

  NAMELIST /inputParams/ prefix, QEDir, VASPDir, exportDir


  contains

!----------------------------------------------------------------------------
  subroutine mpiInitialization()
    !! Initialize MPI processes and split into pools
    !!
    !! <h2>Walkthrough</h2>

    implicit none

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

    call getCommandLineArguments()

    call setUpImages()

    call setUpPools()

    call setUpBands()

    call setUpDiag()

    call setGlobalVariables()

    return
  end subroutine mpiInitialization

!----------------------------------------------------------------------------
  subroutine getCommandLineArguments()

    implicit none

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
          !! * Get the flag (currently only processes number of pools)

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
  subroutine setUpPools()

    implicit none

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
      !! * Create intra pool communicator

    call MPI_BARRIER(world_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8009)

    call MPI_COMM_SPLIT(world_comm_local, indexInPool, myid, inter_pool_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8010)
      !! * Create inter pool communicator

    return
  end subroutine setUpPools

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
      !! * Create intra pool communicator

    call MPI_BARRIER(intra_pool_comm_local, ierr)
    if(ierr /= 0) call mpiExitError(8012)

    call MPI_COMM_SPLIT(intra_pool_comm_local, me_bgrp, myid, inter_bgrp_comm, ierr)
    if(ierr /= 0) call mpiExitError(8013)
      !! * Create inter pool communicator

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
  subroutine initialize()
    !! Set the default values for input variables, open output files,
    !! and start timer
    !!
    !! <h2>Walkthrough</h2>

    use io_files, only : nd_nmbr
      !! @todo Remove this once extracted from QE #end @endtodo
    
    implicit none

    character(len=8) :: cdate
      !! String for date
    character(len=10) :: ctime
      !! String for time

    character(len=6), external :: int_to_char
      !! @todo Remove this once extracted from QE #end @endtodo

    prefix = ''
    QEDir = './'
    VASPDir = './'
    exportDir = './Export'

#ifdef __INTEL_COMPILER
    call remove_stack_limit ( )
      !! * Removed the stack limit because Intel compiler allocates a lot of stack space
      !!   which leads to seg faults and crash. This always works unlike `ulimit -s unlimited`
#endif

    call cpu_time(tStart)

    call date_and_time( cdate, ctime )

#ifdef __MPI
    nd_nmbr = trim(int_to_char(myid+1))
#else
    nd_nmbr = ' '
#endif

    if(ionode_local) then

      write( stdout, '(/5X,"VASP wavefunction export program starts on ",A9," at ",A9)' ) &
             cdate, ctime

#ifdef __MPI
      write( stdout, '(/5X,"Parallel version (MPI), running on ",I5," processors")' ) nproc_local

      if(npool_local > 1) write( stdout, '(5X,"K-points division:     npool_local     = ",I7)' ) npool_local
#else
      write( stdout, '(/5X,"Serial version")' )
#endif

    else

      open( unit = stdout, file='/dev/null', status='unknown' )
        ! Make the stdout unit point to null for non-root processors
        ! to avoid tons of duplicate output

    endif

  end subroutine initialize

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
  subroutine readInputFiles()
    implicit none


    return
  end subroutine readInputFiles

!----------------------------------------------------------------------------
  subroutine readWAVECAR()
    !! Read data from the WAVECAR file
    !!
    !! @todo Update this to have input/output variables #thisbranch @endtodo
    !! <h2>Walkthrough</h2>

    use cell_base, only : at, bg, omega, alat, tpiba
    use wvfct, only : ecutwfc, nbnd
    use klist, only : nkstot
    use lsda_mod, only : nspin
      !! @todo Remove this once extracted from QE #end @endtodo

    implicit none

    real(kind=dp) :: nRecords_real, nspin_real, prec_real, nkstot_real 
      !! Real version of integers for reading from file
    real(kind=dp) :: nbnd_real
      !! Real version of integers for reading from file

    integer :: j
      !! Index used for reading lattice vectors
    integer :: prec
      !! Precision of plane wave coefficients
    integer :: npmax
      !! Maximum number of plane waves
    integer :: nRecords
      !! Number of records in WAVECAR file

    character(len=256) :: fileName
      !! Full WAVECAR file name including path

    if(ionode_local) then
      fileName = trim(VASPDir)//'/WAVECAR'

      nRecords = 24
        ! Set a starting value for the number of records

      open(unit=wavecarUnit, file=fileName, access='direct', recl=nRecords, iostat=ierr, status='old')
      if (ierr .ne. 0) write(stdout,*) 'open error - iostat =',ierr

      read(unit=wavecarUnit,rec=1) nRecords_real, nspin_real, prec_real
        !! @note Must read in as real first then convert to integer @endnote

      close(unit=wavecarUnit)

      nRecords = nint(nRecords_real)
      nspin_local = nint(nspin_real)
      prec = nint(prec_real)
        ! Convert input variables to integers

      !if(prec .eq. 45210) call exitError('readWAVECAR', 'WAVECAR_double requires complex*16', 1)

      open(unit=wavecarUnit, file=fileName, access='direct', recl=nRecords, iostat=ierr, status='old')
      if (ierr .ne. 0) write(stdout,*) 'open error - iostat =',ierr

      read(unit=wavecarUnit,rec=2) nkstot_real, nbnd_real, ecutwfc_local,(at_local(j,1),j=1,3),&
          (at_local(j,2),j=1,3), (at_local(j,3),j=1,3)
        !! * Read total number of kpoints, plane wave cutoff energy, and real
        !!   space lattice vectors

      ecutwfc_local = ecutwfc_local*eVToRy
        !! * Convert energy from VASP to Rydberg to match QE expectation
        !! @todo Remove this once extracted from QE #end @endtodo

      nkstot_local = nint(nkstot_real)
      nbnd_local = nint(nbnd_real)
        ! Convert input variables to integers

      write(stdout,*) 'no. k points =', nkstot_local
      write(stdout,*) 'no. bands =', nbnd_local
      write(stdout,*) 'max. energy =', sngl(ecutwfc_local/eVToRy)
      write(stdout,*) 'real space lattice vectors:'
      write(stdout,*) 'a1 =', (sngl(at_local(j,1)),j=1,3)
      write(stdout,*) 'a2 =', (sngl(at_local(j,2)),j=1,3)
      write(stdout,*) 'a3 =', (sngl(at_local(j,3)),j=1,3)
      write(stdout,*) 

      call calculateOmega(at_local, omega_local)
        !! * Calculate unit cell volume

      write(stdout,*) 'volume unit cell =', sngl(omega)
      write(stdout,*) 

      call getReciprocalVectors(at_local, omega_local, bg_local)

      write(stdout,*) 'reciprocal lattice vectors:'
      write(stdout,*) 'b1 =', (sngl(bg_local(j,1)),j=1,3)
      write(stdout,*) 'b2 =', (sngl(bg_local(j,2)),j=1,3)
      write(stdout,*) 'b3 =', (sngl(bg_local(j,3)),j=1,3)
      write(stdout,*) 


      call estimateMaxNumPlanewaves(bg_local, npmax)

      call mainDoLoop(npmax)
        !! @todo Rename this when figure out what loop does #thisbranch @endtodo

      close(wavecarUnit)

    endif

    call MPI_BCAST(nspin_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nkstot_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(nbnd_local, 1, MPI_INTEGER, root, world_comm_local, ierr)
    call MPI_BCAST(ecutwfc_local, 1, MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(omega_local, 1, MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(at_local, size(at_local), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)
    call MPI_BCAST(bg_local, size(bg_local), MPI_DOUBLE_PRECISION, root, world_comm_local, ierr)

    at = at_local/alat
    bg = bg_local/tpiba
    omega = omega_local
    ecutwfc = ecutwfc_local
    nbnd = nbnd_local
    nkstot = nkstot_local
    nspin = nspin_local
      !! @todo Remove this once extracted from QE #end @endtodo

    return
  end subroutine readWAVECAR

!----------------------------------------------------------------------------
  subroutine calculateOmega(at_local, omega_local)
    implicit none

    real(kind=dp), intent(in) :: at_local(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: omega_local
      !! Volume of unit cell
    real(kind=dp) :: vtmp(3)
      !! \(a_2\times a_3\)

    call vcross(at_local(:,2), at_local(:,3), vtmp)

    omega_local = at_local(1,1)*vtmp(1) + at_local(2,1)*vtmp(2) + at_local(3,1)*vtmp(3)

    return
  end subroutine calculateOmega

!----------------------------------------------------------------------------
  subroutine getReciprocalVectors(at_local, omega_local, bg_local)
    implicit none

    real(kind=dp), intent(in) :: at_local(3,3)
      !! Real space lattice vectors
    real(kind=dp), intent(out) :: bg_local(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp), intent(in) :: omega_local
      !! Volume of unit cell

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
    implicit none

    real(kind=dp) :: vec1(3)
    real(kind=dp) :: vec2(3)
    real(kind=dp) :: crossProd(3)
    
    crossProd(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    crossProd(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    crossProd(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

    return
  end subroutine vcross

!----------------------------------------------------------------------------
  subroutine estimateMaxNumPlanewaves(bg_local, npmax)
    implicit none

    real(kind=dp) :: b1mag, b2mag, b3mag
      !! Reciprocal vector magnitudes
    real(kind=dp), intent(in) :: bg_local(3,3)
      !! Reciprocal lattice vectors
    real(kind=dp) :: c = 0.26246582250210965422
    real(kind=dp) :: phi12, phi13, phi23
      !! Angle between vectors
    real(kind=dp) :: sinphi123
      !! \(\sin\phi_{123}\)
    real(kind=dp) :: vmag
      !! Magnitude of temporary vector
    real(kind=dp) :: vtmp(3)
      !! Temporary vector for calculating angles

    integer :: nb1maxA, nb2maxA, nb3maxA
    integer :: nb1maxB, nb2maxB, nb3maxB
    integer :: nb1maxC, nb2maxC, nb3maxC
    integer :: npmaxA, npmaxB, npmaxC
    integer :: nb1max, nb2max, nb3max
    integer, intent(out) :: npmax
      !! Maximum number of plane waves

    b1mag = sqrt(bg_local(1,1)**2 + bg_local(2,1)**2 + bg_local(3,1)**2)
    b2mag = sqrt(bg_local(1,2)**2 + bg_local(2,2)**2 + bg_local(3,2)**2)
    b3mag = sqrt(bg_local(1,3)**2 + bg_local(2,3)**2 + bg_local(3,3)**2)
      !! * Calculate reciprocal vector magnitudes

    write(stdout,*) 'reciprocal lattice vector magnitudes:'
    write(stdout,*) sngl(b1mag),sngl(b2mag),sngl(b3mag)

    phi12 = acos((bg_local(1,1)*bg_local(1,2) + bg_local(2,1)*bg_local(2,2) + &
        bg_local(3,1)*bg_local(3,2))/(b1mag*b2mag))
      !! * Calculate angle between \(b_1\) and \(b_2\)

    call vcross(bg_local(:,1), bg_local(:,2), vtmp)
    vmag = sqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
    sinphi123 = (bg_local(1,3)*vtmp(1) + bg_local(2,3)*vtmp(2) + &
        bg_local(3,3)*vtmp(3))/(vmag*b3mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxA = (dsqrt(ecutwfc_local/eVToRy*c)/(b1mag*abs(sin(phi12)))) + 1
    nb2maxA = (dsqrt(ecutwfc_local/eVToRy*c)/(b2mag*abs(sin(phi12)))) + 1
    nb3maxA = (dsqrt(ecutwfc_local/eVToRy*c)/(b3mag*abs(sinphi123))) + 1
    npmaxA = nint(4.0*pi*nb1maxA*nb2maxA*nb3maxA/3.0)
      !! * Get first set of max values


    phi13 = acos((bg_local(1,1)*bg_local(1,3) + bg_local(2,1)*bg_local(2,3) + &
        bg_local(3,1)*bg_local(3,3))/(b1mag*b3mag))
      !! * Calculate angle between \(b_1\) and \(b_3\)

    call vcross(bg_local(:,1), bg_local(:,3), vtmp)
    vmag = sqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
    sinphi123 = (bg_local(1,2)*vtmp(1) + bg_local(2,2)*vtmp(2) + &
        bg_local(3,2)*vtmp(3))/(vmag*b3mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxB = (dsqrt(ecutwfc_local/eVToRy*c)/(b1mag*abs(sin(phi13)))) + 1
    nb2maxB = (dsqrt(ecutwfc_local/eVToRy*c)/(b2mag*abs(sinphi123))) + 1
    nb3maxB = (dsqrt(ecutwfc_local/eVToRy*c)/(b3mag*abs(sin(phi13)))) + 1
    npmaxB = nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      !! * Get first set of max values


    phi23 = acos((bg_local(1,2)*bg_local(1,3) + bg_local(2,2)*bg_local(2,3) + &
        bg_local(3,2)*bg_local(3,3))/(b2mag*b3mag))
      !! * Calculate angle between \(b_2\) and \(b_3\)

    call vcross(bg_local(:,2), bg_local(:,3), vtmp)
    vmag = sqrt(vtmp(1)**2 + vtmp(2)**2 + vtmp(3)**2)
    sinphi123 = (bg_local(1,1)*vtmp(1) + bg_local(2,1)*vtmp(2) + &
        bg_local(3,1)*vtmp(3))/(vmag*b1mag)
      !! * Get \(\sin\phi_{123}\)

    nb1maxC = (dsqrt(ecutwfc_local/eVToRy*c)/(b1mag*abs(sinphi123))) + 1
    nb2maxC = (dsqrt(ecutwfc_local/eVToRy*c)/(b2mag*abs(sin(phi23)))) + 1
    nb3maxC = (dsqrt(ecutwfc_local/eVToRy*c)/(b3mag*abs(sin(phi23)))) + 1
    npmaxC = nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)
      !! * Get first set of max values


    nb1max = max0(nb1maxA,nb1maxB,nb1maxC)
    nb2max = max0(nb2maxA,nb2maxB,nb2maxC)
    nb3max = max0(nb3maxA,nb3maxB,nb3maxC)
    npmax = min0(npmaxA,npmaxB,npmaxC)

    write(stdout,*) 'max. no. G values; 1,2,3 =', nb1max, nb2max, nb3max
    write(stdout,*) ' '

    write(stdout,*) 'estimated max. no. plane waves =', npmax

    return
  end subroutine estimateMaxNumPlanewaves

!----------------------------------------------------------------------------
  subroutine mainDoLoop(npmax)

    use klist, only : xk
      !! @todo Remove this once extracted from QE #end @endtodo

    implicit none

    real(kind=dp), allocatable :: cener(:)
      !! Band eigenvalues
    real(kind=dp), allocatable :: coeff(:,:)
      !! Plane wave coefficients
    real(kind=dp), allocatable :: occ(:)
      !! Occupation of band
    real(kind=dp) :: nPlane_real
      !! Real version of integers for reading from file

    integer :: irec, isp, ik, i, ig1, ig2, ig3, iband, iplane
      !! Loop indices
    integer :: nPlane
      !! Input number of plane waves
    integer, intent(in) :: npmax
      !! Maximum number of plane waves

    irec=2

    allocate(occ(nbnd_local))
    allocate(cener(nbnd_local))
    allocate(xk_local(3,nkstot_local))
    allocate(igall(3,npmax))
    allocate(coeff(npmax,nbnd_local))
      !! @todo Make sure that these variables are also deallocated #thisbranch @endtodo

    do isp = 1, nspin_local

       write(stdout,*) ' '
       write(stdout,*) '******'
       write(stdout,*) 'Reading spin ', isp

       do ik = 1, nkstot_local

          irec = irec + 1
       
          read(unit=wavecarUnit,rec=irec) nPlane_real, (xk_local(i,ik),i=1,3), &
               (cener(iband), occ(iband), iband=1,nbnd_local)

          nplane = nint(nPlane_real)
            !! @todo Figure out the difference between this and npw/npwx #thisbranch @endtodo
            !! @note 
            !!  `nplane` is read within the k-point loop, so it could be the 
            !!  number of plane waves per k-point, `npws`/`ngk_g`. May need to 
            !!  add a sum within the loop to get the total number of plane waves.
            !! @endnote

          write(stdout,*) 'Number of plane waves at k-point ', ik, ' (VASP): ', nplane

!    !!$*   Calculate plane waves
!          ncnt=0
!          do ig3=0,2*nb3max
!             ig3p=ig3
!             if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
!             do ig2=0,2*nb2max
!                ig2p=ig2
!                if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
!                do ig1=0,2*nb1max
!                   ig1p=ig1
!                   if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
!                   do j=1,3
!                      sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
!                           (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
!                   enddo
!                   gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
!                   etot=gtot**2/c
!                   if (etot.lt.ecut) then
!                      ncnt=ncnt+1
!                      igall(1,ncnt)=ig1p
!                      igall(2,ncnt)=ig2p
!                      igall(3,ncnt)=ig3p
!                   end if
!                enddo
!             enddo
!          enddo
!
!          if (ncnt.ne.nplane) then
!             write(stdout,*) '*** error - computed no. != input no.'
!             stop
!          endif
!          if (ncnt.gt.npmax) then
!             write(stdout,*) '*** error - plane wave count exceeds estimate'
!             stop
!          endif
!
          do iband = 1, nbnd_local

            irec = irec + 1

            read(unit=wavecarUnit,rec=irec) (coeff(iplane,iband), iplane=1,nplane)

            write(46+ik,*) cener(iband)
              !! @todo 
              !!  Figure out how this and `eigF`/`eigI` relates to `et` #thisbranch @endtodo

          enddo
       enddo
    enddo

    xk = xk_local
      !! @todo Remove this once extracted from QE #end @endtodo

    return
  end subroutine mainDoLoop

!----------------------------------------------------------------------------
  subroutine distributeKpointsInPools()
    !! Figure out how many kpoints there should be per pool
    !!
    !! <h2>Walkthrough</h2>

    implicit none

    integer :: nk_Pool
      !! Number of kpoints in each pool
    integer :: nkr
      !! Number of kpoints left over after evenly divided across pools


    if( nkstot_local > 0 ) then

      IF( ( nproc_pool_local > nproc_local ) .or. ( mod( nproc_local, nproc_pool_local ) /= 0 ) ) &
        CALL exitError( 'distributeKpointsInPools','nproc_pool_local', 1 )

      nk_Pool = nkstot_local / npool_local
        !!  * Calculate k points per pool

      nkr = nkstot_local - nk_Pool * npool_local 
        !! * Calculate the remainder

      IF( myPoolId < nkr ) nk_Pool = nk_Pool + 1
        !! * Assign the remainder to the first `nkr` pools

      !>  * Calculate the index of the first k point in this pool
      ikStart = nk_Pool * myPoolId + 1
      IF( myPoolId >= nkr ) ikStart = ikStart + nkr

      ikEnd = ikStart + nk_Pool - 1
        !!  * Calculate the index of the last k point in this pool

    endif

    return
  end subroutine distributeKpointsInPools

!----------------------------------------------------------------------------
  subroutine reconstructMainGrid()
    use gvect, only : g, ngm, ngm_g, ig_l2g, mill
    use wvfct, only : npwx, npw, g2kin
    use klist, only : nks, xk, ngk
    use cell_base, only : tpiba2

    implicit none

    integer :: ig, ik
      !! Loop indices

    real(kind=dp), allocatable :: rtmp_g( :, : )
      !! ??Not sure what this is

    integer, allocatable :: kisort(:)
      !! ??Not sure what this is

    write(stdout,*) 'Max. no. plane waves (QE) =', npwx

    ! find out the global number of G vectors: ngm_g
    ngm_g = ngm
      !! @todo Figure out how to get this value from VASP files #thisbranch @endtodo

    call MPI_ALLREDUCE(ngm, ngm_g, 1, MPI_INTEGER, MPI_SUM, intra_pool_comm_local, ierr)
    if( ierr /= 0 ) call exitError( 'reconstructMainGrid', 'error in mpi_allreduce 1', ierr)

    if( ionode_local ) then 
    
      write(stdout,*) "Reconstructing the main grid"
    
    endif

    ! collect all G vectors across processors within the pools
    ! and compute their modules
  
    ALLOCATE( itmp_g( 3, ngm_g ) )
    ALLOCATE( rtmp_g( 3, ngm_g ) )

    itmp_g = 0
    DO  ig = 1, ngm
      itmp_g( 1, ig_l2g( ig ) ) = mill(1,ig )
      itmp_g( 2, ig_l2g( ig ) ) = mill(2,ig )
      itmp_g( 3, ig_l2g( ig ) ) = mill(3,ig )
    ENDDO
  
    CALL mp_sum( itmp_g , intra_pool_comm_local )
  
    ! here we are in crystal units
    rtmp_g(1:3,1:ngm_g) = REAL( itmp_g(1:3,1:ngm_g) )
  
    ! go to cartesian units (tpiba)
    CALL cryst_to_cart( ngm_g, rtmp_g, bg_local , 1 )
  
    DEALLOCATE( rtmp_g )

    ! build the G+k array indexes
    ALLOCATE ( igk_l2g ( npwx, nks ) )
    ALLOCATE ( kisort( npwx ) )
    DO ik = 1, nks
      kisort = 0
      npw = npwx
      CALL gk_sort (xk (1, ik+ikStart-1), ngm, g, ecutwfc_local / tpiba2, npw, kisort(1), g2kin)

      ! mapping between local and global G vector index, for this kpoint
     
      DO ig = 1, npw
        
        igk_l2g(ig,ik) = ig_l2g( kisort(ig) )
        
      ENDDO
     
      igk_l2g( npw+1 : npwx, ik ) = 0
     
      ngk (ik) = npw
      write(stdout,*) 'Number of plane waves at k-point ', ik, ' (QE): ', npw

    ENDDO
    DEALLOCATE (kisort)

    ! compute the global number of G+k vectors for each k point
    ALLOCATE( ngk_g( nkstot_local ) )
    ngk_g = 0
    ngk_g( ikStart:ikEnd ) = ngk( 1:nks )
    CALL mp_sum( ngk_g, world_comm_local )

    ! compute the Maximum G vector index among all G+k and processors
    npw_g = maxval( igk_l2g(:,:) )
    CALL mp_max( npw_g, world_comm_local )

    ! compute the Maximum number of G vector among all k points
    npwx_g = maxval( ngk_g( 1:nkstot_local ) )

    return 
  end subroutine reconstructMainGrid

!----------------------------------------------------------------------------
! ..  This subroutine write wavefunctions to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_wfc(iuni, exportDir, &
      ik, nk, ispin, nspin_local, scal, wf0, t0, wfm, tm, ngw, gamma_only, nbnd_local, igl, ngwl )
!
!
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
    use klist, only : nks, xk, ngk
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
    INTEGER, ALLOCATABLE :: itmp1( : )
    INTEGER, ALLOCATABLE :: igwk( :, : )
    INTEGER, ALLOCATABLE :: l2g_new( : )
  
    character(len = 300) :: text
  

    real(DP) :: wfc_scal
    LOGICAL :: twf0, twfm, file_exists
    CHARACTER(iotk_attlenx) :: attr
    TYPE(pseudo_upf) :: upf       ! the pseudo data
    TYPE(radial_grid_type) :: grid

    integer, allocatable :: nnTyp(:), groundState(:)

    IF( ionode_local ) THEN
    

      write(mainout, '("# Cell volume (a.u.)^3. Format: ''(ES24.15E3)''")')
      write(mainout, '(ES24.15E3)' ) omega_local
    
      write(mainout, '("# Number of K-points. Format: ''(i10)''")')
      write(mainout, '(i10)') nkstot_local
    
      write(mainout, '("# ik, groundState, ngk_g(ik), wk(ik), xk(1:3,ik). Format: ''(3i10,4ES24.15E3)''")')
    
      allocate ( groundState(nkstot_local) )

      groundState(:) = 0
      DO ik=1,nkstot_local
        do ibnd = 1, nbnd_local
          if ( wg(ibnd,ik)/wk(ik) < 0.5_dp ) then
          !if (et(ibnd,ik) > ef) then
            groundState(ik) = ibnd - 1
            goto 10
          endif
        enddo
10      continue
      enddo
    
    endif
  
    ALLOCATE( igwk( npwx_g, nkstot_local ) )
  
    DO ik = 1, nkstot_local
      igwk(:,ik) = 0
    
      ALLOCATE( itmp1( npw_g ), STAT= ierr )
      IF ( ierr/=0 ) CALL exitError('pw_export','allocating itmp1', abs(ierr) )
      itmp1 = 0
    
      IF( ik >= ikStart .and. ik <= ikEnd ) THEN
        DO  ig = 1, ngk( ik-ikStart+1 )
          itmp1( igk_l2g( ig, ik-ikStart+1 ) ) = igk_l2g( ig, ik-ikStart+1 )
        ENDDO
      ENDIF
    
      CALL mp_sum( itmp1, world_comm_local )
    
      ngg = 0
      DO  ig = 1, npw_g
        IF( itmp1( ig ) == ig ) THEN
          ngg = ngg + 1
          igwk( ngg , ik) = ig
        ENDIF
      ENDDO
      IF( ngg /= ngk_g( ik ) ) THEN
        if ( ionode_local ) WRITE(mainout, *) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
      ENDIF
    
      DEALLOCATE( itmp1 )
    
      if ( ionode_local ) write(mainout, '(3i10,4ES24.15E3)') ik, groundState(ik), ngk_g(ik), wk(ik), xk(1:3,ik)
    
    ENDDO
  
    if ( ionode_local ) then
    
      write(mainout, '("# Number of G-vectors. Format: ''(i10)''")')
      write(mainout, '(i10)') ngm_g
    
      write(mainout, '("# Number of PW-vectors. Format: ''(i10)''")')
      write(mainout, '(i10)') npw_g
    
      write(mainout, '("# Number of min - max values of fft grid in x, y and z axis. Format: ''(6i10)''")')
      write(mainout, '(6i10)') minval(itmp_g(1,1:ngm_g)), maxval(itmp_g(1,1:ngm_g)), &
                          minval(itmp_g(2,1:ngm_g)), maxval(itmp_g(2,1:ngm_g)), &
                          minval(itmp_g(3,1:ngm_g)), maxval(itmp_g(3,1:ngm_g))
    
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
      DO i = 1, nat
        xyz = tau(:,i)
        write(mainout,'(i10,3ES24.15E3)') ityp(i), tau(:,i)*alat
      ENDDO
    
      write(mainout, '("# Number of Bands. Format: ''(i10)''")')
      write(mainout, '(i10)') nbnd_local
    
      DO ik = 1, nkstot_local
      
        open(72, file=trim(exportDir)//"/grid"//iotk_index(ik))
        write(72, '("# Wave function G-vectors grid")')
        write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
      
        do ink = 1, ngk_g(ik)
          write(72, '(4i10)') igwk(ink,ik), itmp_g(1:3,igwk(ink,ik))
        enddo
      
        close(72)
      
      ENDDO
    
      open(72, file=trim(exportDir)//"/mgrid")
      write(72, '("# Full G-vectors grid")')
      write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
    
      do ink = 1, ngm_g
        write(72, '(4i10)') ink, itmp_g(1:3,ink)
      enddo
    
      close(72)

      write(mainout, '("# Spin. Format: ''(i10)''")')
      write(mainout, '(i10)') nspin_local
    
      allocate( nnTyp(nsp) )
      nnTyp = 0
      do i = 1, nat
        nnTyp(ityp(i)) = nnTyp(ityp(i)) + 1
      enddo

      DO i = 1, nsp
      
        call read_upf(upf, grid, ierr, 71, trim(outdir)//'/'//trim(prefix)//'.save/'//trim(psfile(i)))
      
        if (  upf%typ == 'PAW' ) then
        
          write(stdout, *) ' PAW type pseudopotential found !'
        
          write(mainout, '("# Element")')
          write(mainout, *) trim(atm(i))
          write(mainout, '("# Number of Atoms of this type. Format: ''(i10)''")')
          write(mainout, '(i10)') nnTyp(i)
          write(mainout, '("# Number of projectors. Format: ''(i10)''")')
          write(mainout, '(i10)') upf%nbeta              ! number of projectors
        
          write(mainout, '("# Angular momentum, index of the projectors. Format: ''(2i10)''")')
          ms = 0
          do inb = 1, upf%nbeta
            write(mainout, '(2i10)') upf%lll(inb), inb
            ms = ms + 2*upf%lll(inb) + 1
          enddo
        
          write(mainout, '("# Number of channels. Format: ''(i10)''")')
          write(mainout, '(i10)') ms
        
          write(mainout, '("# Number of radial mesh points. Format: ''(2i10)''")')
          write(mainout, '(2i10)') upf%mesh, upf%kkbeta ! number of points in the radial mesh, number of point inside the aug sphere
        
          write(mainout, '("# Radial grid, Integratable grid. Format: ''(2ES24.15E3)''")')
          do im = 1, upf%mesh
            write(mainout, '(2ES24.15E3)') upf%r(im), upf%rab(im) ! r(mesh) radial grid, rab(mesh) dr(x)/dx (x=linear grid)
          enddo
        
          write(mainout, '("# AE, PS radial wfc for each beta function. Format: ''(2ES24.15E3)''")')
          if ( upf%has_wfc ) then   ! if true, UPF contain AE and PS wfc for each beta
            do inb = 1, upf%nbeta
              do im = 1, upf%mesh
                write(mainout, '(2ES24.15E3)') upf%aewfc(im, inb), upf%pswfc(im, inb)
                                          ! wfc(mesh,nbeta) AE wfc, wfc(mesh,nbeta) PS wfc
              enddo
            enddo
          else
            write(mainout, *) 'UPF does not contain AE and PS wfcs!!'
            stop
          endif
        
        endif
      
      enddo
    
    ENDIF

#ifdef __MPI
  CALL poolrecover (et, nbnd_local, nkstot_local, nks)
#endif


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
  
    if ( ionode_local ) WRITE(stdout,*) "Writing Wavefunctions"
  
    wfc_scal = 1.0d0
    twf0 = .true.
    twfm = .false.
  
    IF ( nkb > 0 ) THEN
    
      CALL init_us_1
      CALL init_at_1
    
      CALL allocate_bec_type (nkb,nbnd_local, becp)
    
      DO ik = 1, nkstot_local
      
        local_pw = 0
        IF ( (ik >= ikStart) .and. (ik <= ikEnd) ) THEN
          CALL gk_sort (xk (1, ik+ikStart-1), ngm, g, ecutwfc_local / tpiba2, npw, igk, g2kin)
          CALL davcio (evc, nwordwfc, iunwfc, (ik-ikStart+1), - 1)

          CALL init_us_2(npw, igk, xk(1, ik), vkb)
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
    
      CALL deallocate_bec_type ( becp )
    
    ENDIF

    DEALLOCATE( igk_l2g )
    DEALLOCATE( igwk )
    DEALLOCATE ( ngk_g )
END SUBROUTINE write_export

end module wfcExportVASPMod
