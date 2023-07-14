module GVelMod
  
  use constants, only: dp
  use base, only: nSpins
  use errorsAndMPI

  implicit none

  integer, allocatable :: bandL(:), bandR(:)
    !! Indices of left and right bands
  integer :: iBandInit, iBandFinal
    !! Initial and final bands
  integer :: iBaseK
    !! Index of base k-point within displacement pattern
  integer :: nBands
    !! Total number of bands
  integer :: nBaseKPoints
    !! Number of base k-points
  integer :: nDegen
    !! Number of degenerate bands
  integer :: nDispkPerCoord
    !! Number of displaced k-points per coordinate

  real(kind=dp) :: degenTol
    !! Tolerance for testing for degeneracies
  real(kind=dp), allocatable :: eigv(:,:,:,:)
    !! Eigenvalue for each band, displacement 
    !! (negative to positive), direction (x/y/z), 
    !! and spin

  character(len=300) :: exportDir
    !! Directory to be used for export

  contains

!----------------------------------------------------------------------------
  subroutine getBaseKPointIndex(nDispkPerCoord, iBaseK)
    !! Get the index of the base k-point within
    !! the displacement pattern in the `groupedEigenvalues.isp.ik`
    !! files
   
    use miscUtilities, only: hpsort_eps

    implicit none

    ! Input variables:
    integer, intent(in) :: nDispkPerCoord
      !! Number of displaced k-points per coordinate

    ! Output variables:
    integer, intent(out) :: iBaseK
      !! Index of base k-point within displacement pattern

    ! Local variables:
    integer :: idk
      !! Loop index
    integer :: iPattSort(nDispkPerCoord)
      !! Sorted index order for patternArr from 
      !! negative to positive

    real(kind=dp) :: eps8 = 1.0E-8_dp
      !! Double precision zero
    real(kind=dp) :: patternArr(nDispkPerCoord+1)
      !! Displacement pattern for groups of
      !! k-points for group velocity calculations


    if(ionode) then

      patternArr(1) = 0.0_dp

      open(72, file=trim(exportDir)//"/groupedEigenvaluesX.1.1")

      read(72,*)
      read(72,*)
      read(72,*) (patternArr(idk), idk=2,nDispkPerCoord+1)

      close(72)

      do idk = 1, nDispkPerCoord
        iPattSort(idk) = idk
      enddo

      call hpsort_eps(nDispkPerCoord+1, patternArr, iPattSort, eps8)
        ! Get sorted index order for patternArr from negative to positive

      do idk = 1, nDispkPerCoord+1
        if(iPattSort(idk) == 1) then
          iBaseK = idk
          exit
        endif
      enddo

    endif

    call MPI_BCAST(iBaseK, 1, MPI_INTEGER, root, worldComm, ierr)

    return

  end subroutine

!----------------------------------------------------------------------------
  subroutine readGroupedEigenvalues(ikLocal, isp, nBands, nBaseKPoints, nDispkPerCoord, nSpins, eigv)

    use miscUtilities, only: int2str, ignoreNextNLinesFromFile

    implicit none

    ! Input variables:
    integer, intent(in) :: ikLocal
      !! Local k-point index
    integer, intent(in) :: isp
      !! Spin index
    integer, intent(in) :: nBands
      !! Total number of bands
    integer, intent(in) :: nBaseKPoints
      !! Total number of k-points
    integer, intent(in) :: nDispkPerCoord
      !! Number of displaced k-points per coordinate
    integer, intent(in) :: nSpins
      !! Number of spins

    ! Output variables:
    real(kind=dp), intent(out) :: eigv(iBandInit:iBandFinal,nDispkPerCoord+1,3,nSpins)
      !! Eigenvalue for each band, displacement 
      !! (negative to positive), direction (x/y/z), 
      !! and spin


    ! Local variables
    integer :: ikGlobal
      !! Global k-point index
    integer :: ib, idk, ix
      !! Loop indices


    if(indexInPool == 0) then

      ikGlobal = ikLocal+ikStart_pool-1

      do ix = 1, 3

        if(ix == 1) then
          open(72, file=trim(exportDir)//"/groupedEigenvaluesX."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)))
        else if(ix == 2) then
          open(72, file=trim(exportDir)//"/groupedEigenvaluesY."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)))
        else if(ix == 3) then
          open(72, file=trim(exportDir)//"/groupedEigenvaluesZ."//trim(int2str(isp))//"."//trim(int2str(ikGlobal)))
        endif
        
        call ignoreNextNLinesFromFile(72,4)
      
        do ib = 1, nBands

          read(72,*) (eigv(ib,idk,ix,isp), idk=1,nDispkPerCoord+1)

        enddo
      
        close(72)

     enddo
    endif

    call MPI_BCAST(eigv, size(eigv), MPI_DOUBLE_PRECISION, root, intraPoolComm, ierr)
    
    return
  end subroutine readGroupedEigenvalues
 
end module GVelMod
