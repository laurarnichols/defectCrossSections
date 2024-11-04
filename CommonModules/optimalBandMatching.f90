module optimalBandMatching
  
  use miscUtilities, only: int2str
  use constants
  use errorsAndMPI
  use mpi
  
  implicit none
  
  contains

!----------------------------------------------------------------------------
  subroutine findOptimalPairsAndOutput(iBandLBra, iBandHBra, iBandLKet, iBandHKet, ikLocal, isp, nPairs, Ufi, optimalPairsDir)

    use hungarianOptimizationMod, only: hungarianOptimalMatching

    implicit none

    ! Input variables:
    integer, intent(in) :: iBandLBra, iBandHBra, iBandLKet, iBandHKet
      !! Optional band bounds
    integer, intent(in) :: ikLocal
      !! Current local k-point
    integer, intent(in) :: isp
      !! Current spin channel
    integer, intent(in) :: nPairs
      !! Number of pairs of bands to get overlaps for

    complex(kind=dp), intent(in) :: Ufi(nPairs)
      !! All-electron overlap

    character(len=300), intent(in) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file

    ! Local variables:
    integer :: ikGlobal
      !! Current global k-point
    integer :: ip, ibB, ibK
      !! Loop indices
    integer, allocatable :: maxBandPair(:)
      !! Optimal band pair to maximize overlaps
    integer :: nBandsEachSys
      !! Number of bands per system (bra/ket) to line up

    real(kind=dp), allocatable :: norm2Ufi2D(:,:)
      !! 2D version of overlap array (needed for passing
      !! into optimization subroutine)
    real(kind=dp) :: t1, t2


    nBandsEachSys = (iBandHBra - iBandLBra + 1)
      ! Doesn't matter which system you use here because we already
      ! confirmed that they have the same number

    allocate(norm2Ufi2D(nBandsEachSys,nBandsEachSys))

    ip = 0
    do ibB = 1, nBandsEachSys
      do ibK = 1, nBandsEachSys
        ip = ip + 1
        norm2Ufi2D(ibB,ibK) = abs(Ufi(ip))**2
      enddo
    enddo


    allocate(maxBandPair(nBandsEachSys))

    if(ionode) write(*,'("  Beginning optimal matching algorithm to match ",i7," sets of bands.")') nBandsEachSys
    call cpu_time(t1)

    call hungarianOptimalMatching(nBandsEachSys, norm2Ufi2D, .false., maxBandPair)

    call cpu_time(t2)
    if(ionode) write(*,'("  Optimal matching complete! (",f10.2," secs)")') t2-t1


    ikGlobal = ikLocal+ikStart_pool-1

    open(64,file=trim(optimalPairsDir)//'/optimalPairs.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal))//'.out')

    write(64,'("# Number of bands in each system, iBandLBra, iBandHBra, iBandLKet, iBandHKet")')
    write(64,'(5i10)') nBandsEachSys, iBandLBra, iBandHBra, iBandLKet, iBandHKet

    write(64,'("# ibBra, ibKet, Overlap")')

    do ibK = 1, nBandsEachSys ! Loop over second (col) index
      write(64,'(2i10,ES13.4E3)') maxBandPair(ibK)+iBandLBra-1, ibK+iBandLKet-1, norm2Ufi2D(maxBandPair(ibK),ibK) 
    enddo

    close(64)

    return

  end subroutine findOptimalPairsAndOutput

!----------------------------------------------------------------------------
  subroutine readAllOptimalPairs(ikGlobal, nSpins, spin1Skipped, spin2Skipped, optimalPairsDir, iBandLKet, iBandHKet, &
        ibBra_optimal)
    ! The optimal pairs file contains the best match between the bra-system (defect)
    ! and ket-system (perfect crystal) state indices. Because we take the energies 
    ! from the perfect-crystal system, I do not rearrange those states. Instead, 
    ! I find the corresponding match from the defect system. 
    !
    ! This file can be confusing because the indices come from the defect system, 
    ! but the perfect-crystal-system states are not shifted while those from the
    ! defect system are. I tried to think of a way around this, but using the defect
    ! indices just makes sense so that all of the codes can work consistently with
    ! the zeroth- and first-order term, with the latter containing only defect states.

    implicit none

    ! Input variables:
    integer, intent(in) :: ikGlobal
      !! Current global k-point
    integer, intent(in) :: nSpins
      !! Number of spins (tested to be consistent
      !! across all systems)

    logical, intent(in) :: spin1Skipped, spin2Skipped
      !! If spin channels skipped

    character(len=300), intent(in) :: optimalPairsDir
      !! Path to store or read optimalPairs.out file

    ! Output variables:
    integer, intent(out) :: iBandLKet, iBandHKet
      !! Band bounds from optimalPairs file
    integer, allocatable, intent(out) :: ibBra_optimal(:,:)
      !! Optimal index from the bra system corresponding 
      !! to the input index from the ket system

    ! Local variables:
    integer :: iBandLBra, iBandHBra
      !! Band bounds from optimalPairs file
    integer :: isp, ibp, ibKet, ibBra
      !! Loop indices
    integer :: nBandsEachSys
      !! Number of bands per system (bra/ket) to line up


    do isp = 1, nSpins
      if((isp == 1 .and. .not. spin1Skipped) .or. (isp == 2 .and. .not. spin2Skipped)) then

        open(64,file=trim(optimalPairsDir)//'/optimalPairs.'//trim(int2str(isp))//'.'//trim(int2str(ikGlobal))//'.out')

        read(64,*)
        read(64,'(5i10)') nBandsEachSys, iBandLBra, iBandHBra, iBandLKet, iBandHKet

        if(isp == 1 .or. spin1Skipped) allocate(ibBra_optimal(nSpins,iBandLKet:iBandHKet))
          ! Assume that the band bounds are the same for both spin channels, which
          ! should currently be the case unless the user messes with the files


        read(64,*)

        do ibp = 1, nBandsEachSys
          read(64,*) ibBra, ibKet
          ibBra_optimal(isp,ibKet) = ibBra
            ! Storing it like this means it will work even if the 
            ! band pairs get reordered somehow
        enddo

        close(64)
      endif
    enddo

    return

  end subroutine readAllOptimalPairs

end module optimalBandMatching
