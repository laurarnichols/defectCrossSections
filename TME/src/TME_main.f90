program TMEmain
  use TMEmod
  
  implicit none
  
  integer :: ikLocal, ikGlobal, isp
    !! Loop indices

  logical :: bothSpinChannelsExist = .false.
    !! If `allElecOverlap.isp.ik` files exist
    !! for both spin channels at the current 
    !! k-point
  logical :: spin1Read = .false.
    !! If spin channel 1 at the current k-point
    !! was read from a file. Tells us if we need
    !! to read/calculate some things  for the second
    !! that are normally done for the first channel
    !! only

  
  call mpiInitialization('TME')
    !! Initialize MPI

  call getCommandLineArguments()
    !! * Get the number of pools from the command line

  call setUpPools()
    !! * Split up processors between pools and generate MPI
    !!   communicators for pools
  
  if(ionode) call cpu_time(t0)


  if(ionode) write(*, '("Pre-k-loop: [ ] Read inputs  [ ] Read full PW grid  [ ] Set up tables ")')
  call cpu_time(t1)


  call readInput(ispSelect, maxAngMom, maxGIndexGlobal, nKPoints, nGVecsGlobal, realLattVec, recipLattVec, baselineDir, &
        loopSpins, subtractBaseline)
    !! Read input, initialize, check that required variables were set, and
    !! distribute across processes
    !! @todo Figure out if `realLattVec` used anywhere. If not remove. @endtodo


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools

  call distributeItemsInSubgroups(indexInPool, nGVecsGlobal, nProcPerPool, nProcPerPool, nProcPerPool, iGStart_pool, &
          iGEnd_pool, nGVecsLocal)
    !! * Distribute G-vectors across processes in pool


  call cpu_time(t2)
  if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [ ] Read full PW grid  [ ] Set up tables (",f10.2," secs)")') t2-t1
  call cpu_time(t1)

  allocate(mill_local(3,nGVecsLocal))

  call getFullPWGrid(iGStart_pool, nGVecsLocal, nGVecsGlobal, mill_local)

  
  call cpu_time(t2)
  if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [X] Read full PW grid  [ ] Set up tables (",f10.2," secs)")') t2-t1
  call cpu_time(t1)


  allocate(gCart(3,nGVecsLocal))
  allocate(Ylm((maxAngMom+1)**2,nGVecsLocal))

  call setUpTables(maxAngMom, nGVecsLocal, mill_local, numOfTypes, recipLattVec, atoms, gCart, Ylm)

  deallocate(mill_local)

  
  call cpu_time(t2)
  if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [X] Read full PW grid  [X] Set up tables (",f10.2," secs)")') t2-t1
  call cpu_time(t1)
    

  allocate(Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nKPerPool, nSpins))
  Ufi(:,:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
  

  do ikLocal = 1, nkPerPool
    
    if(ionode) write(*,'("Beginning k-point loop ", i4, " of ", i4)') ikLocal, nkPerPool
    
    ikGlobal = ikLocal+ikStart_pool-1
      !! Get the global `ik` index from the local one
    
    if(indexInPool == 0) then

      if(nSpins == 1) then
        bothSpinChannelsExist = overlapFileExists(ikGlobal, 1)
      else if(nSpins == 2) then
        bothSpinChannelsExist = overlapFileExists(ikGlobal, 1) .and. overlapFileExists(ikGlobal, 2)
      endif
      
    endif
    
    call MPI_BCAST(bothSpinChannelsExist, 1, MPI_LOGICAL, root, intraPoolComm, ierr)

    if(.not. bothSpinChannelsExist) then

      !-----------------------------------------------------------------------------------------------
      !> Read projectors

      if(ionode) write(*, '("  Pre-spin-loop: [ ] Read projectors ")') 
      call cpu_time(t1)


      call distributeItemsInSubgroups(indexInPool, npwsPC(ikGlobal), nProcPerPool, nProcPerPool, nProcPerPool, iGkStart_poolPC, &
              iGkEnd_poolPC, nGkVecsLocalPC)

      allocate(betaPC(nGkVecsLocalPC,nProjsPC))

      call readProjectors('PC', iGkStart_poolPC, ikGlobal, nGkVecsLocalPC, nProjsPC, betaPC)


      call distributeItemsInSubgroups(indexInPool, npwsSD(ikGlobal), nProcPerPool, nProcPerPool, nProcPerPool, iGkStart_poolSD, &
              iGkEnd_poolSD, nGkVecsLocalSD)

      allocate(betaSD(nGkVecsLocalSD,nProjsSD))

      call readProjectors('SD', iGkStart_poolSD, ikGlobal, nGkVecsLocalSD, nProjsSD, betaSD)

      call cpu_time(t2)
      if(ionode) write(*, '("  Pre-spin-loop: [X] Read projectors (",f10.2," secs)")') t2-t1

      
      !-----------------------------------------------------------------------------------------------
      !> Allocate arrays that are potentially spin-polarized and start spin loop

      allocate(wfcPC(nGkVecsLocalPC, iBandIinit:iBandIfinal))
      allocate(wfcSD(nGkVecsLocalSD, iBandFinit:iBandFfinal))
      
      allocate(cProjBetaPCPsiSD(nProjsPC, iBandFinit:iBandFfinal))
      allocate(cProjBetaSDPhiPC(nProjsSD, iBandIinit:iBandIfinal))

      allocate(cProjPC(nProjsPC, iBandIinit:iBandIfinal))
      allocate(cProjSD(nProjsSD, iBandFinit:iBandFfinal))

      allocate(paw_PsiPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))
      allocate(paw_SDPhi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))

      allocate(pawKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nGVecsLocal))
      allocate(pawSDK(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nGVecsLocal))

      allocate(paw_id(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))

      
      do isp = 1, nSpins

        if(ionode) write(*,'("  Beginning spin loop ", i2, " of ", i2)') isp, nSpins

        !-----------------------------------------------------------------------------------------------
        !> Check if the `allElecOverlap.isp.ik` file exists

        if(.not. overlapFileExists(ikGlobal, isp)) then

          !-----------------------------------------------------------------------------------------------
          !> Read wave functions and calculate overlap
      
          if(ionode) write(*, '("    Ufi calculation: [ ] Overlap  [ ] Cross projections  [ ] PAW wfc  [ ] PAW k")') 
          call cpu_time(t1)


          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Read) &
            call readWfc('PC', iBandIinit, iBandIfinal, iGkStart_poolPC, ikGlobal, min(isp,nSpinsPC), nGkVecsLocalPC, npwsPC(ikGlobal), wfcPC)
        
          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Read) &
            call readWfc('SD', iBandFinit, iBandFfinal, iGkStart_poolSD, ikGlobal, min(isp,nSpinsSD), nGkVecsLocalSD, npwsSD(ikGlobal), wfcSD)
        
          call calculatePWsOverlap(ikLocal, isp)
        

          call cpu_time(t2)
          if(ionode) write(*, '("    Ufi calculation: [X] Overlap  [ ] Cross projections  [ ] PAW wfc  [ ] PAW k (",f6.2," secs)")') t2-t1
          call cpu_time(t1)


          !-----------------------------------------------------------------------------------------------
          !> Calculate cross projections

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Read) &
            call calculateCrossProjection(iBandFinit, iBandFfinal, nGkVecsLocalPC, nGkVecsLocalSD, nProjsPC, betaPC, wfcSD, cProjBetaPCPsiSD)


          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Read) &
            call calculateCrossProjection(iBandIinit, iBandIfinal, nGkVecsLocalSD, nGkVecsLocalPC, nProjsSD, betaSD, wfcPC, cProjBetaSDPhiPC)
        

          call cpu_time(t2)
          if(ionode) write(*, '("    Ufi calculation: [X] Overlap  [X] Cross projections  [ ] PAW wfc  [ ] PAW k (",f6.2," secs)")') t2-t1
          call cpu_time(t1)


          !-----------------------------------------------------------------------------------------------
          !> Have process 0 in each pool calculate the PAW wave function correction for PC

          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Read) &
            call readProjections('PC', iBandIinit, iBandIfinal, ikGlobal, min(isp,nSpinsPC), nProjsPC, cProjPC)

          if(indexInPool == 0) then

            call pawCorrectionWfc(nIonsPC, TYPNIPC, nProjsPC, numOfTypesPC, cProjPC, cProjBetaPCPsiSD, atomsPC, paw_PsiPC)

          endif


          !-----------------------------------------------------------------------------------------------
          !> Have process 1 in each pool calculate the PAW wave function correction for PC

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Read) &
            call readProjections('SD', iBandFinit, iBandFfinal, ikGlobal, min(isp,nSpinsSD), nProjsSD, cProjSD)

          if(indexInPool == 1) then

            call pawCorrectionWfc(nIonsSD, TYPNISD, nProjsSD, numOfTypes, cProjBetaSDPhiPC, cProjSD, atoms, paw_SDPhi)

          endif

          call cpu_time(t2)
          if(ionode) write(*, '("    Ufi calculation: [X] Overlap  [X] Cross projections  [X] PAW wfc  [ ] PAW k (",f6.2," secs)")') t2-t1
          call cpu_time(t1)
      
      
          !-----------------------------------------------------------------------------------------------
          !> Broadcast wave function corrections to other processes

          call MPI_BCAST(paw_PsiPC, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
          call MPI_BCAST(paw_SDPhi, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)


          !-----------------------------------------------------------------------------------------------
          !> Have all processes calculate the PAW k correction
      
          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Read) &
            call pawCorrectionK('PC', nIonsPC, TYPNIPC, maxAngMom, nGVecsLocal, numOfTypesPC, numOfTypes, posIonPC, gCart, Ylm, atomsPC, atoms, pawKPC)
      

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Read) &
            call pawCorrectionK('SD', nIonsSD, TYPNISD, maxAngMom, nGVecsLocal, numOfTypes, numOfTypes, posIonSD, gCart, Ylm, atoms, atoms, pawSDK)

        
          call cpu_time(t2)
          if(ionode) write(*, '("    Ufi calculation: [X] Overlap  [X] Cross projections  [X] PAW wfc  [X] PAW k (",f8.2," secs)")') t2-t1
          call cpu_time(t1)
      

          !-----------------------------------------------------------------------------------------------
          !> Sum over PAW k corrections

          paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
          do ibi = iBandIinit, iBandIfinal
        
            do ibf = iBandFinit, iBandFfinal   
              paw_id(ibf,ibi) = dot_product(conjg(pawSDK(ibf,ibi,:)),pawKPC(ibf,ibi,:))
            enddo
        
          enddo

          Ufi(:,:,ikLocal,isp) = Ufi(:,:,ikLocal,isp) + paw_id(:,:)*16.0_dp*pi*pi/omega


          call MPI_ALLREDUCE(MPI_IN_PLACE, Ufi(:,:,ikLocal,isp), size(Ufi(:,:,ikLocal,isp)), MPI_DOUBLE_COMPLEX, &
                  MPI_SUM, intraPoolComm, ierr)


          if(indexInPool == 0) then 
          
            Ufi(:,:,ikLocal,isp) = Ufi(:,:,ikLocal,isp) + paw_SDPhi(:,:) + paw_PsiPC(:,:)

            if(order == 1 .and. subtractBaseline) &
              call readAndSubtractBaseline(iBandIinit, iBandIfinal, iBandFinit, iBandFfinal, ikLocal, isp, nSpins, Ufi)
        
            call writeResults(ikLocal,isp, Ufi)
        
          endif
  
      
        endif ! If this spin channel exists


        if(isp == nSpins) then
          deallocate(wfcSD, wfcPC)
          deallocate(betaPC, betaSD)
          deallocate(cProjBetaPCPsiSD, cProjBetaSDPhiPC)
          deallocate(cProjPC, cProjSD)
          deallocate(pawKPC, pawSDK)
          deallocate(paw_id)
          deallocate(paw_PsiPC, paw_SDPhi)
        endif

      enddo ! Spin loop
    endif ! If both spin channels exist
  enddo ! k-point loop


  call MPI_BARRIER(worldComm, ierr)
  if(ionode) write(*,'("Done with k loop!")')

  
  call finalizeCalculation()
  
  call MPI_FINALIZE(ierr)
  
end program TMEmain
