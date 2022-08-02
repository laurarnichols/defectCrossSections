!
! Copyright (C) 2003-2009 Andrea Ferretti and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE io_base_export
!=----------------------------------------------------------------------------=!

! do i = 1, nk             !                                                   !
!   WAVEFUNCTIONS( i )     !  write_restart_wfc         read_restart_wfc       !
! end do                   !                                                   !

  USE io_global,  ONLY : stdout
  USE kinds
  USE parameters, ONLY: nsx

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: file_version = 202
  INTEGER :: restart_module_verbosity = 0

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!
!
!=----------------------------------------------------------------------------=!

! ..  This subroutine write wavefunctions to the disk
! .. Where:
! iuni    = Restart file I/O fortran unit
!
    SUBROUTINE write_restart_wfc(iuni, exportDir, &
      ik, nk, kunit, ispin, nspin, scal, wf0, t0, wfm, tm, ngw, gamma_only, nbnd, igl, ngwl )
!
      USE mp_wave
      USE mp, ONLY: mp_sum, mp_get, mp_bcast, mp_max
      USE mp_global, ONLY: mpime, nproc, root, me_pool, my_pool_id, &
        nproc_pool, intra_pool_comm, root_pool, world_comm
      USE io_global, ONLY: ionode, ionode_id
      USE iotk_module
!
      IMPLICIT NONE
!
      INTEGER, INTENT(in) :: iuni
      character(len = 256), intent(in) :: exportDir
      INTEGER, INTENT(in) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP), INTENT(in) :: wf0(:,:)
      COMPLEX(DP), INTENT(in) :: wfm(:,:)
      INTEGER, INTENT(in) :: ngw   !
      LOGICAL, INTENT(in) :: gamma_only
      INTEGER, INTENT(in) :: nbnd
      INTEGER, INTENT(in) :: ngwl
      INTEGER, INTENT(in) :: igl(:)
      REAL(DP), INTENT(in) :: scal
      LOGICAL, INTENT(in) :: t0, tm

      INTEGER :: i, j, ierr, idum = 0
      INTEGER :: nkl, nkr, nkbl, iks, ike, nkt, ikt, igwx, ig
      INTEGER :: npool, ipmask( nproc ), ipsour
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

        !  find out the number of pools
        npool = nproc / nproc_pool

        !  find out number of k points blocks
        nkbl = nkt / kunit

        !  k points per pool
        nkl = kunit * ( nkbl / npool )

        !  find out the reminder
        nkr = ( nkt - nkl * npool ) / kunit

        !  Assign the reminder to the first nkr pools
        IF( my_pool_id < nkr ) nkl = nkl + kunit

        !  find out the index of the first k point in this pool
        iks = nkl * my_pool_id + 1
        IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

        !  find out the index of the last k point in this pool
        ike = iks + nkl - 1

        ipmask = 0
        ipsour = ionode_id

        !  find out the index of the processor which collect the data in the pool of ik
        IF( npool > 1 ) THEN
          IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
            IF( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
          ENDIF
          CALL mp_sum( ipmask )
          DO i = 1, nproc
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
          ENDDO
        ENDIF

        igwx = 0
        ierr = 0
        IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
          IF( ngwl > size( igl ) ) THEN
            ierr = 1
          ELSE
            igwx = maxval( igl(1:ngwl) )
          ENDIF
        ENDIF

        ! get the maximum index within the pool
        !
        CALL mp_max( igwx, intra_pool_comm )

        ! now notify all procs if an error has been found
        !
        CALL mp_max( ierr )

        IF( ierr > 0 ) &
          CALL errore(' write_restart_wfc ',' wrong size ngl ', ierr )

        IF( ipsour /= ionode_id ) THEN
          CALL mp_get( igwx, igwx, mpime, ionode_id, ipsour, 1 )
        ENDIF

        ALLOCATE( wtmp( max(igwx,1) ) )
        wtmp = cmplx(0.0_dp, 0.0_dp, kind=dp)

        DO j = 1, nbnd
          IF( t0 ) THEN
            IF( npool > 1 ) THEN
              IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
                CALL mergewf(wf0(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm)
              ENDIF
              IF( ipsour /= ionode_id ) THEN
                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j, world_comm )
              ENDIF
            ELSE
              CALL mergewf(wf0(:,j), wtmp, ngwl, igl, mpime, nproc, &
                           ionode_id, world_comm )
            ENDIF

            IF( ionode ) THEN
              do ig = 1, igwx
                write(iuni, '(2ES24.15E3)') wtmp(ig)
              enddo
              !
!              do j = 1, nbnd
!                do i = 1, igwx ! ngk_g(ik)
!                  write(74,'(2ES24.15E3)') wf0(i,j) ! wf0 is the local array for evc(i,j)
!                enddo
!              enddo
              !
            ENDIF
          ELSE
          ENDIF
        ENDDO

!        DO j = 1, nbnd
!          IF( tm ) THEN
!            IF( npool > 1 ) THEN
!              IF( ( ikt >= iks ) .and. ( ikt <= ike ) ) THEN
!                CALL mergewf(wfm(:,j), wtmp, ngwl, igl, me_pool, &
!                             nproc_pool, root_pool, intra_pool_comm)
!              ENDIF
!              IF( ipsour /= ionode_id ) THEN
!                CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j, world_comm )
!              ENDIF
!            ELSE
!              CALL mergewf(wfm(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id, world_comm )
!            ENDIF
!            IF( ionode ) THEN
!              CALL iotk_write_dat(iuni,"Wfcm"//iotk_index(j),wtmp(1:igwx))
!            ENDIF
!          ELSE
!          ENDIF
!        ENDDO
        IF(ionode) then
          close(iuni)
          !CALL iotk_write_end  (iuni,"Kpoint"//iotk_index(ik))
        endif
      
        DEALLOCATE( wtmp )

      RETURN
    END SUBROUTINE

  END MODULE
!=----------------------------------------------------------------------------=!



!-----------------------------------------------------------------------
PROGRAM pw_export_for_TME
  !-----------------------------------------------------------------------
  !
  ! writes PWSCF data for postprocessing purposes in XML format using IOTK lib
  ! Wave-functions are collected and written using IO_BASE module.
  !
  ! input:  namelist "&inputpp", with variables
  !   prefix       prefix of input files saved by program pwscf
  !   outdir       temporary directory where files resides
  !   exportDir    output directory. A directory
  !                "exportDir" is created and
  !                output files are put there. All the data
  !                are accessible through the ""exportDir"/input" file.
  !
  USE wrappers,      ONLY : f_mkdir
  USE pwcom
  USE io_global,     ONLY : ionode, ionode_id
  USE io_files,      ONLY : prefix, tmp_dir, outdir
  USE ions_base, ONLY : ntype => nsp
  USE iotk_module
  USE mp_global,     ONLY : mp_startup, kunit
  USE mp,            ONLY : mp_bcast
  USE environment,   ONLY : environment_start
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ik, i, kunittmp, ios
  !
  real(kind = dp), parameter :: ryToHartree = 0.5_dp
  !
  CHARACTER(len=256) :: pp_file, exportDir
  LOGICAL :: writeWFC
  !
  NAMELIST /inputpp/ prefix, outdir, exportDir, writeWFC
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW_EXPORT' )
  !
  !   set default values for variables in namelist
  !
  prefix = ''
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  exportDir = './Export'
  !
  writeWFC  = .true.        ! gdb : by default the wavefunctions are needed, 
                            !       this gives the user the ability not to write the wavefunctions
  !
  !    Reading input file
  !
  IF ( ionode ) THEN
    !
    CALL input_from_file ( )
    !
    READ(5, inputpp, IOSTAT=ios)
    !
    IF (ios /= 0) CALL errore ('pw_export', 'reading inputpp namelist', abs(ios) )
    !
    ios = f_mkdir( trim(exportDir) )
    !
    pp_file = trim(exportDir)//"/input"
    !
    !
  ENDIF
  !
  ! ... Broadcasting variables
  !
  tmp_dir = trimcheck( outdir )
  CALL mp_bcast( outdir, ionode_id )
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  !

#if defined __MPI
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  CALL write_export (pp_file, exportDir, kunittmp )

  CALL stop_pp
 
  !
CONTAINS

!
!-----------------------------------------------------------------------
SUBROUTINE write_export (pp_file, exportDir, kunit )
  !-----------------------------------------------------------------------
  !
  USE iotk_module


  USE kinds,          ONLY : DP
  USE pwcom
  USE start_k,        ONLY : nk1, nk2, nk3, k1, k2, k3
  USE control_flags,  ONLY : gamma_only
  USE global_version, ONLY : version_number
  USE becmod,         ONLY : bec_type, becp, calbec, &
                             allocate_bec_type, deallocate_bec_type
  USE uspp,          ONLY : nkb, vkb
  USE wavefunctions_module,  ONLY : evc
  USE io_files,       ONLY : outdir, prefix, iunwfc, nwordwfc
  USE io_files,       ONLY : psfile
  USE io_base_export, ONLY : write_restart_wfc
  USE io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : atm, nat, ityp, tau, nsp
  USE mp_global,      ONLY : nproc, nproc_pool, mpime
  USE mp_global,      ONLY : my_pool_id, intra_pool_comm, inter_pool_comm
  USE mp,             ONLY : mp_sum, mp_max
  !
  USE upf_module,     ONLY : read_upf
  !
  USE pseudo_types, ONLY : pseudo_upf
  USE radial_grids, ONLY : radial_grid_type
  !
  USE wvfct,         ONLY : wg
  !
  USE paw_variables,        ONLY : okpaw, ddd_paw, total_core_energy, only_paw
  USE paw_onecenter,        ONLY : PAW_potential
  USE paw_symmetry,         ONLY : PAW_symmetrize_ddd
  USE uspp_param,           ONLY : nh, nhm ! used for PAW
  USE uspp,                 ONLY : qq_so, dvan_so, qq, dvan
  USE scf,                  ONLY : rho
  !
  IMPLICIT NONE
  !
  CHARACTER(5), PARAMETER :: fmt_name="QEXPT"
  CHARACTER(5), PARAMETER :: fmt_version="1.1.0"

  INTEGER, INTENT(in) :: kunit
  CHARACTER(256), INTENT(in) :: pp_file, exportDir

  INTEGER :: i, j, k, ig, ik, ibnd, na, ngg,ig_, ierr
  INTEGER, ALLOCATABLE :: kisort(:)
  real(DP) :: xyz(3), tmp(3)
  INTEGER :: npool, nkbl, nkl, nkr, npwx_g, im, ink, inb, ms
  INTEGER :: ike, iks, npw_g, ispin, local_pw
  INTEGER, ALLOCATABLE :: ngk_g( : )
  INTEGER, ALLOCATABLE :: itmp_g( :, : )
  real(DP),ALLOCATABLE :: rtmp_g( :, : )
  real(DP),ALLOCATABLE :: rtmp_gg( : )
  INTEGER, ALLOCATABLE :: itmp1( : )
  INTEGER, ALLOCATABLE :: igwk( :, : )
  INTEGER, ALLOCATABLE :: l2g_new( : )
  INTEGER, ALLOCATABLE :: igk_l2g( :, : )
  !


  !
  character(len = 300) :: text
  !

  real(DP) :: wfc_scal
  LOGICAL :: twf0, twfm, file_exists
  CHARACTER(iotk_attlenx) :: attr
  TYPE(pseudo_upf) :: upf       ! the pseudo data
  TYPE(radial_grid_type) :: grid

  integer, allocatable :: nTyp(:), groundState(:)

  IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .or. ( mod( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' write_export ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .or. ( mod( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' write_export ',' nproc_pool ', 1 )

     !  find out the number of pools
     npool = nproc / nproc_pool

     !  find out number of k points blocks
     nkbl = nkstot / kunit

     !  k points per pool
     nkl = kunit * ( nkbl / npool )

     !  find out the reminder
     nkr = ( nkstot - nkl * npool ) / kunit

     !  Assign the reminder to the first nkr pools
     IF( my_pool_id < nkr ) nkl = nkl + kunit

     !  find out the index of the first k point in this pool
     iks = nkl * my_pool_id + 1
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

     !  find out the index of the last k point in this pool
     ike = iks + nkl - 1

  ENDIF

  ! find out the global number of G vectors: ngm_g
  ngm_g = ngm
  CALL mp_sum( ngm_g , intra_pool_comm )


  !  Open file PP_FILE

  IF( ionode ) THEN
    !
    WRITE(stdout,*) "Opening file "//trim(pp_file)
    !
    open(50, file=trim(pp_file))
    !
    WRITE(stdout,*) "Reconstructing the main grid"
    !
  endif

  ! collect all G vectors across processors within the pools
  ! and compute their modules
  !
  ALLOCATE( itmp_g( 3, ngm_g ) )
  ALLOCATE( rtmp_g( 3, ngm_g ) )
  ALLOCATE( rtmp_gg( ngm_g ) )

  itmp_g = 0
  DO  ig = 1, ngm
    itmp_g( 1, ig_l2g( ig ) ) = mill(1,ig )
    itmp_g( 2, ig_l2g( ig ) ) = mill(2,ig )
    itmp_g( 3, ig_l2g( ig ) ) = mill(3,ig )
  ENDDO
  !
  CALL mp_sum( itmp_g , intra_pool_comm )
  !
  ! here we are in crystal units
  rtmp_g(1:3,1:ngm_g) = REAL( itmp_g(1:3,1:ngm_g) )
  !
  ! go to cartesian units (tpiba)
  CALL cryst_to_cart( ngm_g, rtmp_g, bg , 1 )
  !
  ! compute squared moduli
  DO  ig = 1, ngm_g
     rtmp_gg(ig) = rtmp_g(1,ig)**2 + rtmp_g(2,ig)**2 + rtmp_g(3,ig)**2
  ENDDO
  DEALLOCATE( rtmp_g )

  ! build the G+k array indexes
  ALLOCATE ( igk_l2g ( npwx, nks ) )
  ALLOCATE ( kisort( npwx ) )
  DO ik = 1, nks
     kisort = 0
     npw = npwx
     CALL gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, kisort(1), g2kin)
     !
     ! mapping between local and global G vector index, for this kpoint
     !
     DO ig = 1, npw
        !
        igk_l2g(ig,ik) = ig_l2g( kisort(ig) )
        !
     ENDDO
     !
     igk_l2g( npw+1 : npwx, ik ) = 0
     !
     ngk (ik) = npw
  ENDDO
  DEALLOCATE (kisort)

  ! compute the global number of G+k vectors for each k point
  ALLOCATE( ngk_g( nkstot ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g )

  ! compute the Maximum G vector index among all G+k and processors
  npw_g = maxval( igk_l2g(:,:) )
  CALL mp_max( npw_g )

  ! compute the Maximum number of G vector among all k points
  npwx_g = maxval( ngk_g( 1:nkstot ) )

  IF( ionode ) THEN
    !

    write(50, '("# Cell volume (a.u.)^3. Format: ''(ES24.15E3)''")')
    write(50, '(ES24.15E3)' ) omega
    !
    write(50, '("# Number of K-points. Format: ''(i10)''")')
    write(50, '(i10)') nkstot 
    !
    write(50, '("# ik, groundState, ngk_g(ik), wk(ik), xk(1:3,ik). Format: ''(3i10,4ES24.15E3)''")')
    !
    allocate ( groundState(nkstot) )
    !
    groundState = 0
    DO ik=1,nkstot
      do ibnd = 1, nbnd
        if ( wg(ibnd,ik)/wk(ik) < 0.5_dp ) then
        !if (et(ibnd,ik) > ef) then
          groundState(ik) = ibnd - 1
          goto 10
        endif
      enddo
 10   continue
    enddo
    !
  endif
  !
  ALLOCATE( igwk( npwx_g, nkstot ) )
  !
  DO ik = 1, nkstot
    igwk(:,ik) = 0
    !
    ALLOCATE( itmp1( npw_g ), STAT= ierr )
    IF ( ierr/=0 ) CALL errore('pw_export','allocating itmp1', abs(ierr) )
    itmp1 = 0
    !
    IF( ik >= iks .and. ik <= ike ) THEN
      DO  ig = 1, ngk( ik-iks+1 )
        itmp1( igk_l2g( ig, ik-iks+1 ) ) = igk_l2g( ig, ik-iks+1 )
      ENDDO
    ENDIF
    !
    CALL mp_sum( itmp1 )
    !
    ngg = 0
    DO  ig = 1, npw_g
      IF( itmp1( ig ) == ig ) THEN
        ngg = ngg + 1
        igwk( ngg , ik) = ig
      ENDIF
    ENDDO
    IF( ngg /= ngk_g( ik ) ) THEN
      if ( ionode ) WRITE(50, *) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    ENDIF
    !
    DEALLOCATE( itmp1 )
    !
    if ( ionode ) write(50, '(3i10,4ES24.15E3)') ik, groundState(ik), ngk_g(ik), wk(ik), xk(1:3,ik)
    !
  ENDDO
  !
  if ( ionode ) then
    !  
    write(50, '("# Number of G-vectors. Format: ''(i10)''")')
    write(50, '(i10)') ngm_g
    !
    write(50, '("# Number of PW-vectors. Format: ''(i10)''")')
    write(50, '(i10)') npw_g
    !
    write(50, '("# Number of min - max values of fft grid in x, y and z axis. Format: ''(6i10)''")')
    write(50, '(6i10)') minval(itmp_g(1,1:ngm_g)), maxval(itmp_g(1,1:ngm_g)), &
                        minval(itmp_g(2,1:ngm_g)), maxval(itmp_g(2,1:ngm_g)), &
                        minval(itmp_g(3,1:ngm_g)), maxval(itmp_g(3,1:ngm_g))
    !
    write(50, '("# Cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
    write(50, '("# a1 ",3ES24.15E3)') at(:,1)*alat
    write(50, '("# a2 ",3ES24.15E3)') at(:,2)*alat
    write(50, '("# a3 ",3ES24.15E3)') at(:,3)*alat
    !
    write(50, '("# Reciprocal cell (a.u.). Format: ''(a5, 3ES24.15E3)''")')
    write(50, '("# b1 ",3ES24.15E3)') bg(:,1)*tpiba
    write(50, '("# b2 ",3ES24.15E3)') bg(:,2)*tpiba
    write(50, '("# b3 ",3ES24.15E3)') bg(:,3)*tpiba
    !
    write(50, '("# Number of Atoms. Format: ''(i10)''")')
    write(50, '(i10)') nat
    !
    write(50, '("# Number of Types. Format: ''(i10)''")')
    write(50, '(i10)') nsp
    !
    write(50, '("# Atoms type, position(1:3) (a.u.). Format: ''(i10,3ES24.15E3)''")')
    DO i = 1, nat
      xyz = tau(:,i)
      write(50,'(i10,3ES24.15E3)') ityp(i), tau(:,i)*alat
    ENDDO
    !
    write(50, '("# Number of Bands. Format: ''(i10)''")')
    write(50, '(i10)') nbnd
    !
    DO ik = 1, nkstot
      !
      open(72, file=trim(exportDir)//"/grid"//iotk_index(ik))
      write(72, '("# Wave function G-vectors grid")')
      write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
      !
      do ink = 1, ngk_g(ik)
        write(72, '(4i10)') igwk(ink,ik), itmp_g(1:3,igwk(ink,ik))
      enddo
      !
      close(72)
      !
    ENDDO
    !
    open(72, file=trim(exportDir)//"/mgrid")
    write(72, '("# Full G-vectors grid")')
    write(72, '("# G-vector index, G-vector(1:3) miller indices. Format: ''(4i10)''")')
    !
    do ink = 1, ngm_g
      write(72, '(4i10)') ink, itmp_g(1:3,ink)
    enddo
    !
    close(72)
    !
    !DEALLOCATE( itmp_g )
    !
    write(50, '("# Spin. Format: ''(i10)''")')
    write(50, '(i10)') nspin
    !
    allocate( nTyp(nsp) )
    nTyp = 0
    do i = 1, nat
      nTyp(ityp(i)) = nTyp(ityp(i)) + 1
    enddo
    !
    DO i = 1, nsp
      !
      call read_upf(upf, grid, ierr, 71, trim(outdir)//'/'//trim(prefix)//'.save/'//trim(psfile(i)))
      !
      if (  upf%typ == 'PAW' ) then
        !
        write(stdout, *) ' PAW type pseudopotential found !'
        !
        write(50, '("# Element")')
        write(50, *) trim(atm(i))
        write(50, '("# Number of Atoms of this type. Format: ''(i10)''")')
        write(50, '(i10)') nTyp(i)
        write(50, '("# Number of projectors. Format: ''(i10)''")')
        write(50, '(i10)') upf%nbeta              ! number of projectors
        !
        write(50, '("# Angular momentum, index of the projectors. Format: ''(2i10)''")')
        ms = 0
        do inb = 1, upf%nbeta
          write(50, '(2i10)') upf%lll(inb), inb
          ms = ms + 2*upf%lll(inb) + 1
        enddo
        !
        write(50, '("# Number of channels. Format: ''(i10)''")')
        write(50, '(i10)') ms
        !
        write(50, '("# Number of radial mesh points. Format: ''(2i10)''")')
        write(50, '(2i10)') upf%mesh, upf%kkbeta ! number of points in the radial mesh, number of point inside the aug sphere
        !
        write(50, '("# Radial grid, Integratable grid. Format: ''(2ES24.15E3)''")')
        do im = 1, upf%mesh
          write(50, '(2ES24.15E3)') upf%r(im), upf%rab(im) ! r(mesh) radial grid, rab(mesh) dr(x)/dx (x=linear grid)
        enddo
        !
        write(50, '("# AE, PS radial wfc for each beta function. Format: ''(2ES24.15E3)''")')
        if ( upf%has_wfc ) then   ! if true, UPF contain AE and PS wfc for each beta
          do inb = 1, upf%nbeta
            do im = 1, upf%mesh
              write(50, '(2ES24.15E3)') upf%aewfc(im, inb), upf%pswfc(im, inb) 
                                        ! wfc(mesh,nbeta) AE wfc, wfc(mesh,nbeta) PS wfc
            enddo
          enddo
        else
          write(50, *) 'UPF does not contain AE and PS wfcs!!'
          stop
        endif
        !
      endif
      !
    enddo
    !
  ENDIF
  !
  DEALLOCATE( rtmp_gg )

!  ! for each k point build and write the global G+k indexes array
!  ALLOCATE( igwk( npwx_g,nkstot ) )
!  !WRITE(0,*) "Writing grids for wfc"
!  !CALL iotk_write_attr (attr,"npwx",npwx_g,first=.true.)
!  !IF(ionode) CALL iotk_write_begin(50,"Wfc_grids",ATTR=attr)
!
!
!  DO ik = 1, nkstot
!    igwk(:,ik) = 0
!    !
!    ALLOCATE( itmp1( npw_g ), STAT= ierr )
!    IF ( ierr/=0 ) CALL errore('pw_export','allocating itmp1', abs(ierr) )
!    itmp1 = 0
!    !
!    IF( ik >= iks .and. ik <= ike ) THEN
!      DO  ig = 1, ngk( ik-iks+1 )
!        itmp1( igk_l2g( ig, ik-iks+1 ) ) = igk_l2g( ig, ik-iks+1 )
!      ENDDO
!    ENDIF
!    !
!    CALL mp_sum( itmp1 )
!    !
!    ngg = 0
!    DO  ig = 1, npw_g
!      IF( itmp1( ig ) == ig ) THEN
!        ngg = ngg + 1
!        igwk( ngg , ik) = ig
!      ENDIF
!    ENDDO
!    IF( ngg /= ngk_g( ik ) ) THEN
!      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
!    ENDIF
!    !
!    DEALLOCATE( itmp1 )
!    !
!  ENDDO
!
!  DEALLOCATE( itmp_g )
!
!
#ifdef __MPI
  CALL poolrecover (et, nbnd, nkstot, nks)
#endif


  WRITE(stdout,*) "Writing Eigenvalues"

  IF( ionode ) THEN
    !
    write(50, '("# Fermi Energy (Hartree). Format: ''(ES24.15E3)''")')
    write(50, '(ES24.15E3)') ef*ryToHartree
    flush(50)
    !
    DO ik = 1, nkstot
      !
      ispin = isk( ik )
      !
      open(72, file=trim(exportDir)//"/eigenvalues"//iotk_index(ik))
      !
      write(72, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
      write(72, '("# Eigenvalues (Hartree), band occupation number. Format: ''(2ES24.15E3)''")')
      !
      do ibnd = 1, nbnd
        if ( wk(ik) == 0.D0 ) then 
          write(72, '(2ES24.15E3)') et(ibnd,ik)*ryToHartree, wg(ibnd,ik)
        else
          write(72, '(2ES24.15E3)') et(ibnd,ik)*ryToHartree, wg(ibnd,ik)/wk(ik)
        endif
      enddo
      !
      close(72)
      !
    ENDDO
    !
  endif
  !
  if ( ionode .and. writeWFC ) WRITE(stdout,*) "Writing Wavefunctions"
  !
  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.
  !
  IF ( nkb > 0 ) THEN
    !
    CALL init_us_1
    CALL init_at_1
    !
    CALL allocate_bec_type (nkb,nbnd, becp)
    !
    DO ik = 1, nkstot
      !
      local_pw = 0
      IF ( (ik >= iks) .and. (ik <= ike) ) THEN
        CALL gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        CALL davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)

        CALL init_us_2(npw, igk, xk(1, ik), vkb)
        local_pw = ngk(ik-iks+1)

        IF ( gamma_only ) THEN
          CALL calbec ( ngk_g(ik), vkb, evc, becp )
          WRITE(0,*) 'Gamma only PW_EXPORT not yet tested'
        ELSE
          CALL calbec ( npw, vkb, evc, becp )
          if ( ionode ) then
            !
            WRITE(stdout,*) "Writing projectors of kpt", ik
            !
!            file_exists = .false.
!            inquire(file =trim(exportDir)//"/projectors"//iotk_index(ik), exist = file_exists)
!            if ( .not. file_exists ) then
!              open(73, file=trim(exportDir)//"/projectors"//iotk_index(ik))
!              write(73, '("# Complex projectors |beta>. Format: ''(2ES24.15E3)''")')
!              write(73,'(2i10)') nkb, ngk_g(ik)
!              do j = 1, nkb
!                do i = 1, ngk_g(ik)
!                  write(73,'(2ES24.15E3)') vkb(i,j)
!                enddo
!              enddo
!              close(73)
!            endif
!            !
!            file_exists = .false.
!            inquire(file =trim(exportDir)//"/evc"//iotk_index(ik), exist = file_exists)
!            if ( .not. file_exists ) then
!              !
!              open(74, file=trim(exportDir)//"/evc"//iotk_index(ik))
!              write(74, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
!              write(74, '("# Complex : wavefunction coefficients (a.u.)^(-3/2). Format: ''(2ES24.15E3)''")')
!              write(74,'(2i10)') nbnd, ngk_g(ik)
!              !
!              do j = 1, nbnd
!                do i = 1, ngk_g(ik)
!                  write(74,'(2ES24.15E3)') evc(i,j) 
!                enddo
!              enddo
!              !
!              close(74)
!              !
!            endif
!            !
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
              !
              close(72)
              !
            endif
          endif
        ENDIF
      ENDIF

      ALLOCATE(l2g_new(local_pw))

      l2g_new = 0
      DO ig = 1, local_pw
        ngg = igk_l2g(ig,ik-iks+1)
        DO ig_ = 1, ngk_g(ik)
          IF(ngg == igwk(ig_,ik)) THEN
            l2g_new(ig) = ig_
            exit
          ENDIF
        ENDDO
      ENDDO
      !
      ispin = isk( ik )
      !
      if ( ionode ) then

        file_exists = .false.
        inquire(file =trim(exportDir)//"/wfc"//iotk_index(ik), exist = file_exists)
        if ( .not. file_exists ) then
          open (72, file=trim(exportDir)//"/wfc"//iotk_index(ik))
          write(72, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
          write(72, '("# Complex : wavefunction coefficients (a.u.)^(-3/2). Format: ''(2ES24.15E3)''")')
          !
          open(73, file=trim(exportDir)//"/projectors"//iotk_index(ik))
          write(73, '("# Complex projectors |beta>. Format: ''(2ES24.15E3)''")')
          write(73,'(2i10)') nkb, ngk_g(ik)
!          WRITE(stdout,*) "Writing Wavefunctions of kpt", ik
!          open(74, file=trim(exportDir)//"/evc"//iotk_index(ik))
!          write(74, '("# Spin : ",i10, " Format: ''(a9, i10)''")') ispin
!          write(74, '("# Complex : wavefunction coefficients (a.u.)^(-3/2). Format: ''(2ES24.15E3)''")')
        endif
      endif
      !
      CALL mp_bcast( file_exists, ionode_id )
      !
      if ( .not. file_exists ) then
        CALL write_restart_wfc(72, exportDir, ik, nkstot, kunit, ispin, nspin, &
                               wfc_scal, evc, twf0, evc, twfm, npw_g, gamma_only, nbnd, &
                               l2g_new(:),local_pw )
        CALL write_restart_wfc(73, exportDir, ik, nkstot, kunit, ispin, nspin, &
                               wfc_scal, vkb, twf0, evc, twfm, npw_g, gamma_only, nkb, &
                               l2g_new(:), local_pw )
      endif
      !
      if ( .not. file_exists .and. ionode ) then
        close(72)
        close(73)
!        close(74)
      endif
      !
      DEALLOCATE(l2g_new)
    ENDDO
    !
    CALL deallocate_bec_type ( becp )
    !
  ENDIF

  DEALLOCATE( igk_l2g )
  DEALLOCATE( igwk )
  DEALLOCATE ( ngk_g )

END SUBROUTINE write_export

END PROGRAM pw_export_for_TME

