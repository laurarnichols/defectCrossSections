program transitionMatrixElements
  use declarations
  
  implicit none
  
  integer :: ikLocal, ikGlobal, isp
    !! Loop indices

  logical :: bothSpinChannelsExist = .false.
    !! If `allElecOverlap.isp.ik` files exist
    !! for both spin channels at the current 
    !! k-point
  logical :: spin1Exists = .false.
    !! If `allElecOverlap.1.ik` file exists
    !! for spin channel 1 at the current 
    !! k-point
  logical :: thisSpinChannelExists = .false.
    !! If `allElecOverlap.isp.ik` file exists
    !! for single spin channel at the current 
    !! k-point

  character(len=300) :: ikC
    !! Character index

  
  call mpiInitialization()
    !! Initialize MPI
  
  if(ionode) call cpu_time(t0)


  if(ionode) write(*, '("Pre-k-loop: [ ] Read inputs  [ ] Read full PW grid ")')
  call cpu_time(t1)


  call readInput(maxGIndexGlobal, nKPoints, nGVecsGlobal, realLattVec, recipLattVec)
    !! Read input, initialize, check that required variables were set, and
    !! distribute across processes
    !! @todo Figure out if `realLattVec` used anywhere. If not remove. @endtodo


  call distributeItemsInSubgroups(myPoolId, nKPoints, nProcs, nProcPerPool, nPools, ikStart_pool, ikEnd_pool, nkPerPool)
    !! * Distribute k-points in pools

  call distributeItemsInSubgroups(indexInPool, nGVecsGlobal, nProcPerPool, nProcPerPool, nProcPerPool, iGStart_pool, &
          iGEnd_pool, nGVecsLocal)
    !! * Distribute G-vectors across processes in pool


  call cpu_time(t2)
  if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [ ] Read full PW grid (",f10.2," secs)")') t2-t1
  call cpu_time(t1)

  allocate(mill_local(3,nGVecsLocal))

  call getFullPWGrid(iGStart_pool, iGEnd_pool, nGVecsLocal, nGVecsGlobal, mill_local)

  
  call cpu_time(t2)
  if(ionode) write(*, '("Pre-k-loop: [X] Read inputs  [X] Read full PW grid (",f10.2," secs)")') t2-t1
    

  if(indexInPool == 0) allocate(eigvI(iBandIinit:iBandIfinal), eigvF(iBandFinit:iBandFfinal))
  
  allocate(Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nKPerPool, nSpins))
  Ufi(:,:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
  

  do ikLocal = 1, nkPerPool
    
    ikGlobal = ikLocal+ikStart_pool-1
      !! Get the global `ik` index from the local one
    
    if(indexInPool == 0) then

      do isp = 1, nSpins
        call checkIfCalculated(ikGlobal, isp, bothSpinChannelsExist)
      enddo
      
    endif
    
    call MPI_BCAST(bothSpinChannelsExist, 1, MPI_LOGICAL, root, intraPoolComm, ierr)

    if(bothSpinChannelsExist) then

      if(indexInPool == 0) then

        do isp = 1, nSpins

          call readUfis(ikLocal,isp)

        enddo

      endif
    
    else

      !-----------------------------------------------------------------------------------------------
      !> Read projectors

      if(indexInPool == 0) &
        write(*, '("Pre-spin-loop for k-point",i4,": [ ] Read projectors ")') &
              ikGlobal
      call cpu_time(t1)


      call distributeItemsInSubgroups(indexInPool, npwsPC(ikGlobal), nProcPerPool, nProcPerPool, nProcPerPool, iGkStart_poolPC, &
              iGkEnd_poolPC, nGkVecsLocalPC)

      allocate(betaPC(nGkVecsLocalPC,nProjsPC))

      call readProjectors('PC', iGkStart_poolPC, ikGlobal, nGkVecsLocalPC, nProjsPC, npwsPC(ikGlobal), betaPC)


      call distributeItemsInSubgroups(indexInPool, npwsSD(ikGlobal), nProcPerPool, nProcPerPool, nProcPerPool, iGkStart_poolSD, &
              iGkEnd_poolSD, nGkVecsLocalSD)

      allocate(betaSD(nGkVecsLocalSD,nProjsSD))

      call readProjectors('SD', iGkStart_poolSD, ikGlobal, nGkVecsLocalSD, nProjsSD, npwsSD(ikGlobal), betaSD)

      call cpu_time(t2)
      if(indexInPool == 0) &
        write(*, '("Pre-spin-loop for k-point",i4,": [X] Read projectors (",f10.2," secs)")') &
              ikGlobal, t2-t1

      
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

        !-----------------------------------------------------------------------------------------------
        !> Check if the `allElecOverlap.isp.ik` file exists

        if(indexInPool == 0) then

          call checkIfCalculated(ikGlobal, isp, thisSpinChannelExists)
      
        endif
    
        call MPI_BCAST(thisSpinChannelExists, 1, MPI_LOGICAL, root, intraPoolComm, ierr)

        if(thisSpinChannelExists) then

          if(indexInPool == 0) then

            call readUfis(ikLocal,isp)

          endif

          spin1Exists = .true.
    
        else

          !-----------------------------------------------------------------------------------------------
          !> Read wave functions and calculate overlap
      
          if(indexInPool == 0) &
            write(*, '("Ufi calculation for k-point",i4," and spin ",i1,": [ ] Overlap  [ ] Cross projections  [ ] PAW wfc  [ ] PAW k")') &
              ikGlobal, isp
          call cpu_time(t1)


          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call readWfc('PC', iBandIinit, iBandIfinal, iGkStart_poolPC, ikGlobal, min(isp,nSpinsPC), nGkVecsLocalPC, npwsPC(ikGlobal), wfcPC)
        
          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call readWfc('SD', iBandFinit, iBandFfinal, iGkStart_poolSD, ikGlobal, min(isp,nSpinsSD), nGkVecsLocalSD, npwsSD(ikGlobal), wfcSD)
        
          call calculatePWsOverlap(ikLocal, isp)
        

          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("Ufi calculation for k-point",i4," and spin ",i1,": [X] Overlap  [ ] Cross projections  [ ] PAW wfc  [ ] PAW k (",f6.2," secs)")') &
              ikGlobal, isp, t2-t1
          call cpu_time(t1)


          !-----------------------------------------------------------------------------------------------
          !> Calculate cross projections

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call calculateCrossProjection(iBandFinit, iBandFfinal, ikGlobal, nGkVecsLocalPC, nGkVecsLocalSD, nProjsPC, betaPC, wfcSD, cProjBetaPCPsiSD)

          if(isp == nSpins) then
            deallocate(wfcSD)
            deallocate(betaPC)
          endif


          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call calculateCrossProjection(iBandIinit, iBandIfinal, ikGlobal, nGkVecsLocalSD, nGkVecsLocalPC, nProjsSD, betaSD, wfcPC, cProjBetaSDPhiPC)
        
          if(isp == nSpins) then
            deallocate(wfcPC)
            deallocate(betaSD)
          endif


          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("Ufi calculation for k-point",i4," and spin ",i1,": [X] Overlap  [X] Cross projections  [ ] PAW wfc  [ ] PAW k (",f6.2," secs)")') &
              ikGlobal, isp, t2-t1
          call cpu_time(t1)


          !-----------------------------------------------------------------------------------------------
          !> Have process 0 in each pool calculate the PAW wave function correction for PC

          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call readProjections('PC', iBandIinit, iBandIfinal, ikGlobal, min(isp,nSpinsPC), nProjsPC, cProjPC)

          if(indexInPool == 0) then

            call pawCorrectionWfc(nIonsPC, TYPNIPC, nProjsPC, cProjPC, cProjBetaPCPsiSD, atomsPC, paw_PsiPC)

          endif
          
          if(isp == nSpins) deallocate(cProjBetaPCPsiSD)


          !-----------------------------------------------------------------------------------------------
          !> Have process 1 in each pool calculate the PAW wave function correction for PC

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call readProjections('SD', iBandFinit, iBandFfinal, ikGlobal, min(isp,nSpinsSD), nProjsSD, cProjSD)

          if(indexInPool == 1) then

            call pawCorrectionWfc(nIonsSD, TYPNISD, nProjsSD, cProjBetaSDPhiPC, cProjSD, atoms, paw_SDPhi)

          endif

          if(isp == nSpins) deallocate(cProjBetaSDPhiPC)

          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("Ufi calculation for k-point",i4," and spin ",i1,": [X] Overlap  [X] Cross projections  [X] PAW wfc  [ ] PAW k (",f6.2," secs)")') &
              ikGlobal, isp, t2-t1
          call cpu_time(t1)
      
      
          !-----------------------------------------------------------------------------------------------
          !> Broadcast wave function corrections to other processes

          call MPI_BCAST(paw_PsiPC, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
          call MPI_BCAST(paw_SDPhi, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)


          !-----------------------------------------------------------------------------------------------
          !> Have all processes calculate the PAW k correction
      
          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call pawCorrectionK('PC', nIonsPC, TYPNIPC, numOfTypesPC, posIonPC, atomsPC, atoms, pawKPC)

          if(isp == nSpins) then
            deallocate(cProjPC)
          endif
      

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call pawCorrectionK('SD', nIonsSD, TYPNISD, numOfTypes, posIonSD, atoms, atoms, pawSDK)

          if(isp == nSpins) then
            deallocate(cProjSD)
          endif

        
          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("Ufi calculation for k-point",i4," and spin ",i1,": [X] Overlap  [X] Cross projections  [X] PAW wfc  [X] PAW k (",f6.2," secs)")') &
              ikGlobal, isp, t2-t1
          call cpu_time(t1)
      

          !-----------------------------------------------------------------------------------------------
          !> Sum over PAW k corrections

          paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
          do ibi = iBandIinit, iBandIfinal
        
            do ibf = iBandFinit, iBandFfinal   
              paw_id(ibf,ibi) = sum(pawSDK(ibf,ibi,:)*pawKPC(ibf,ibi,:))
            enddo
        
          enddo

          if(isp == nSpins) then
            deallocate(pawKPC)
            deallocate(pawSDK)
          endif

          Ufi(:,:,ikLocal,isp) = Ufi(:,:,ikLocal,isp) + paw_id(:,:)*16.0_dp*pi*pi/omega

          if(isp == nSpins) deallocate(paw_id)

          call MPI_ALLREDUCE(MPI_IN_PLACE, Ufi(:,:,ikLocal,isp), size(Ufi(:,:,ikLocal,isp)), MPI_DOUBLE_COMPLEX, &
                  MPI_SUM, intraPoolComm, ierr)


          if(indexInPool == 0) then 
          
            Ufi(:,:,ikLocal,isp) = Ufi(:,:,ikLocal,isp) + paw_SDPhi(:,:) + paw_PsiPC(:,:)
        
            call writeResults(ikLocal,isp)
        
          endif
  
          if(isp == nSpins) then
            deallocate(paw_PsiPC)
            deallocate(paw_SDPhi)
          endif
      
        endif ! If this spin channel exists
      enddo ! Spin loop
    endif ! If both spin channels exist
  enddo ! k-point loop


  call MPI_BARRIER(worldComm, ierr)
  if(ionode) write(*,'("Done with k loop!")')
    
  if(calculateVfis .and. indexInPool == 0) call calculateVfiElements()

  
  call finalizeCalculation()
  
  call MPI_FINALIZE(ierr)
  
end program transitionMatrixElements
