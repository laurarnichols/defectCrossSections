program transitionMatrixElements
  use declarations
  
  implicit none
  
  integer :: ikLocal, ikGlobal, iType, isp
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
    
  call readInput(maxGIndexGlobal, nKPoints, nGVecsGlobal, realLattVec, recipLattVec)
    !! Read input, initialize, check that required variables were set, and
    !! distribute across processes
    !! @todo Figure out if `realLattVec` used anywhere. If not remove. @endtodo


  call distributeKpointsInPools(nKPoints)
    !! Split the k-points across the pools


  allocate(gIndexGlobalToLocal(nGVecsGlobal), gVecProcId(nGVecsGlobal))

  call getFullPWGrid(nGVecsGlobal, gIndexGlobalToLocal, gVecProcId, mill_local, nGVecsLocal)
    !! Read the full PW grid from `mgrid` and distribute
    !! across processes



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
      !> Read PW grids from `grid.ik` files

      allocate(gKIndexGlobalPC(npwsPC(ikGlobal)))
      allocate(gKIndexGlobalSD(npwsSD(ikGlobal)))

      if(indexInPool == 0) call readGrid('PC', ikGlobal, npwsPC(ikGlobal), gKIndexGlobalPC)
      if(indexInPool == 1) call readGrid('SD', ikGlobal, npwsSD(ikGlobal), gKIndexGlobalSD)

      call MPI_BCAST(gKIndexGlobalPC, size(gKIndexGlobalPC), MPI_INT, 0, intraPoolComm, ierr)
      call MPI_BCAST(gKIndexGlobalSD, size(gKIndexGlobalSD), MPI_INT, 1, intraPoolComm, ierr)

      !-----------------------------------------------------------------------------------------------
      !> Read projectors

      allocate(betaPC(nGVecsLocal,nProjsPC))

      call readProjectors('PC', ikGlobal, nProjsPC, npwsPC(ikGlobal), gKIndexGlobalPC, betaPC)

      allocate(betaSD(nGVecsLocal,nProjsSD))

      call readProjectors('SD', ikGlobal, nProjsSD, npwsSD(ikGlobal), gKIndexGlobalSD, betaSD)

      
      !-----------------------------------------------------------------------------------------------
      !> Allocate arrays that are potentially spin-polarized and start spin loop

      allocate(wfcPC(nGVecsLocal, iBandIinit:iBandIfinal))
      allocate(wfcSD(nGVecsLocal, iBandFinit:iBandFfinal))
      
      allocate(cProjBetaPCPsiSD(nProjsPC, nBands))
      allocate(cProjBetaSDPhiPC(nProjsSD, nBands))

      allocate(cProjPC(nProjsPC, nBands))
      allocate(cProjSD(nProjsSD, nBands))

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
      
          if(indexInPool == 0) write(*, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4, " and spin ", i2)') &
                  ikGlobal, nKPoints, isp
          call cpu_time(t1)


          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call readWfc('PC', iBandIinit, iBandIfinal, ikGlobal, min(isp,nSpinsPC), npwsPC(ikGlobal), gKIndexGlobalPC, wfcPC)
        
          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call readWfc('SD', iBandFinit, iBandFfinal, ikGlobal, min(isp,nSpinsSD), npwsSD(ikGlobal), gKIndexGlobalSD, wfcSD)
        
          call calculatePWsOverlap(ikLocal, isp)
        

          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> for k-point", i4, " and spin ", i1, " done in", f10.2, " secs.")') &
                  ikGlobal, isp, t2-t1
          call cpu_time(t1)


          !-----------------------------------------------------------------------------------------------
          !> Calculate cross projections

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call calculateCrossProjection(iBandFinit, iBandFfinal, ikGlobal, nProjsPC, npwsPC(ikGlobal), betaPC, wfcSD, cProjBetaPCPsiSD)

          if(isp == nSpins) then
            deallocate(wfcSD)
            deallocate(gKIndexGlobalPC)
            deallocate(betaPC)
          endif


          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("      Calculating <betaPC|wfcSD> for k-point", i4, " and spin ", i1, " done in", f10.2, " secs.")') &
                  ikGlobal, isp, t2-t1
          call cpu_time(t1)


          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call calculateCrossProjection(iBandIinit, iBandIfinal, ikGlobal, nProjsSD, npwsSD(ikGlobal), betaSD, wfcPC, cProjBetaSDPhiPC)
        
          if(isp == nSpins) then
            deallocate(wfcPC)
            deallocate(gKIndexGlobalSD)
            deallocate(betaSD)
          endif


          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("      Calculating <betaSD|wfcPC> for k-point", i4, " and spin ", i1, " done in", f10.2, " secs.")') &
                  ikGlobal, isp, t2-t1


          !-----------------------------------------------------------------------------------------------
          !> Have process 0 in each pool calculate the PAW wave function correction for PC

          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call readProjections('PC', ikGlobal, min(isp,nSpinsPC), nProjsPC, cProjPC)

          if(indexInPool == 0) then

            call pawCorrectionWfc(nIonsPC, TYPNIPC, cProjPC, cProjBetaPCPsiSD, atomsPC, paw_PsiPC)

          endif
          
          if(isp == nSpins) deallocate(cProjBetaPCPsiSD)


          !-----------------------------------------------------------------------------------------------
          !> Have process 1 in each pool calculate the PAW wave function correction for PC

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call readProjections('SD', ikGlobal, min(isp,nSpinsSD), nProjsSD, cProjSD)

          if(indexInPool == 1) then

            call pawCorrectionWfc(nIonsSD, TYPNISD, cProjBetaSDPhiPC, cProjSD, atoms, paw_SDPhi)

          endif

          if(isp == nSpins) deallocate(cProjBetaSDPhiPC)
      
      
          !-----------------------------------------------------------------------------------------------
          !> Broadcast wave function corrections to other processes

          call MPI_BCAST(paw_PsiPC, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
          call MPI_BCAST(paw_SDPhi, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)


          !-----------------------------------------------------------------------------------------------
          !> Have all processes calculate the PAW k correction

          call cpu_time(t1)
      
          if(isp == 1 .or. nSpinsPC == 2 .or. spin1Exists) &
            call pawCorrectionK('PC', nIonsPC, TYPNIPC, numOfTypesPC, posIonPC, atomsPC, atoms, pawKPC)

          if(isp == nSpins) deallocate(cProjPC)
      

          call cpu_time(t2)
          if(indexInPool == 0) &
            write(*, '("      <\\vec{k}|PAW_PC> for k-point ", i4, " and spin ", i1, " done in", f10.2, " secs.")') &
                  ikGlobal, isp, t2-t1
          call cpu_time(t1)
      

          if(isp == 1 .or. nSpinsSD == 2 .or. spin1Exists) &
            call pawCorrectionK('SD', nIonsSD, TYPNISD, numOfTypes, posIonSD, atoms, atoms, pawSDK)

          if(isp == nSpins) deallocate(cProjSD)

        
          call cpu_time(t2)
          if(indexInPool == 0) write(*, '("      <PAW_SD|\\vec{k}> for k-point ", i4, " and spin ", i1, " done in", f10.2, " secs.")') &
                  ikGlobal, isp, t2-t1
      

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


  deallocate(npwsPC)
  deallocate(wkPC)
  deallocate(xkPC)
  deallocate(posIonPC)
  deallocate(TYPNIPC)

  do iType = 1, numOfTypesPC
    deallocate(atomsPC(iType)%r)
    deallocate(atomsPC(iType)%lps)
    deallocate(atomsPC(iType)%F)
    deallocate(atomsPC(iType)%F1)
    deallocate(atomsPC(iType)%F2)
    deallocate(atomsPC(iType)%bes_J_qr)
  enddo

  deallocate(atomsPC)

  deallocate(npwsSD)
  deallocate(wk)
  deallocate(xk)
  deallocate(posIonSD)
  deallocate(TYPNISD)

  do iType = 1, numOfTypes
    deallocate(atoms(iType)%lps)
    deallocate(atoms(iType)%r)
    deallocate(atoms(iType)%F)
    deallocate(atoms(iType)%F1)
    deallocate(atoms(iType)%F2)
    deallocate(atoms(iType)%bes_J_qr)
  enddo

  deallocate(atoms)

  deallocate(gIndexGlobalToLocal)
  deallocate(gVecProcId)
  deallocate(mill_local)
  

  if(indexInPool == 0) deallocate(eigvI, eigvF)
  
  ! Calculating Vfi
  
  if (ionode) then
    
    ! Finalize Calculation
     
    call finalizeCalculation()
    
  endif
  
  call MPI_FINALIZE(ierr)
  
end program transitionMatrixElements
