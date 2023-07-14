program GVel

  use GVelMod

  implicit none

  integer :: ikLocal, ikGlobal, ib, ix, ibd, iDegen, isp
    !! Loop indices

  real(kind=dp) :: t0, t1, t2
    !! Timers


  call cpu_time(t0)

  call mpiInitialization('GVel')

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools


  ! Get inputs:
  !  * nBaseKPoints (get only number of middle k-points)
  !  * nDispkPerCoord
  !  * nSpins
  !  * iBandInit, iBandFinal
  !  * degenTol
  !  * exportDir

  call distributeItemsInSubgroups(myPoolId, nBaseKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools


  call getBaseKPointIndex(nDispkPerCoord, iBaseK)


  allocate(eigv(iBandInit:iBandFinal,nDispkPerCoord+1,3,nSpins))
    ! Band, displacement (negative to positive), 
    ! direction (x/y/z), and spin
  allocate(bandL(iBandInit:iBandFinal))
  allocate(bandR(iBandInit:iBandFinal))

  do isp = 1, nSpins
    do ikLocal = 1, nkPerPool
      ! Iterate through middle k-points

      call readGroupedEigenvalues(ikLocal, isp, nBands, nBaseKPoints, nDispkPerCoord, nSpins, eigv)
        ! Read eigenvalues for this k-point for all positions and directions

      do ix = 1, 3

        ib = iBandInit
        do while(ib <= iBandFinal)

          ibd = ib + 1

          ! Find any degenerate bands in this group
          do while(ibd <= iBandFinal .and. abs(eigv(ibd,ix,iBaseK,isp) - eigv(ib,ix,iBaseK,isp) < degenTol))
            ibd = ibd + 1
          enddo

          nDegen = ibd - ib

          if(nDegen == 1) then
            ! If there are no degeneracies for this band,
            ! assume that the band goes straight across

            bandL(ib) = ib
            bandR(ib) = ib

  !       * Lock in those with good enough fits
          
          else
            ! If there are degeneracies, assume that bands
            ! form a fan, where lowest bands on one side
            ! correspond to highest bands on the other side
            !
            ! NOTE: I am not sure how bands should be assigned
            ! here so that they properly line up with bands in
            ! other directions. Still waiting on guidance from
            ! Sok and Xiaoguang on how that should work.

            ! Need to handle the case where last band is in 
            ! degeneracy group

        
            do iDegen = 1, nDegen
              bandL(ib+iDegen-1) = ib+iDegen-1
                ! Left band goes lowest to highest
              bandR(ib+iDegen-1) = ib+(nDegen-1)-(iDegen-1)
                ! Right band goes highest to lowest
            enddo

            ! NOTE: As far as I understand, the degeneracy groups
            ! will be symmetric, and we don't care about the sign
            ! of the slope, so the higher side is arbitrary.
            

  !         For each degeneracy:
  !           * Do linear regression on points 
  !           * Lock in those with good enough fits


            ib = ibd - 1
              ! Make the outer band loop skip the degenerate bands 
              ! that have already been handled

          endif
        enddo

  !     * Output degenerate bands info
  !     * Output bands not locked in after first round as potential crossings
  !     While there are bands not locked in or we've done < 5 loops
  !       For each band:
  !         * Check if locked in
  !         * If not locked in, check the fit of this point and the point above (not locked in) with the left or right points swapped 
  !         * Choose the one that makes both fits better
  !         * If the fits now meet the tolerance, lock them in
      enddo
    enddo
  enddo

  deallocate(eigv)
  deallocate(bandL)
  deallocate(bandR)

end program GVel
