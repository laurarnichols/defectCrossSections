program transitionMatrixElements
  use declarations
  
  implicit none
  
  integer :: ikLocal, ikGlobal, iType
    !! Loop indices

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



  allocate(paw_id(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))
  allocate(paw_PsiPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))
  allocate(paw_SDPhi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))

  if(indexInPool == 0) allocate(eigvI(iBandIinit:iBandIfinal), eigvF(iBandFinit:iBandFfinal))
  
  allocate(Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nKPerPool, nSpins))
  Ufi(:,:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
  

  do ikLocal = 1, nkPerPool
    
    ikGlobal = ikLocal+ikStart_pool-1
      !! Get the global `ik` index from the local one
    
    if(indexInPool == 0) then
      
      tmes_file_exists = .false.
      call checkIfCalculated(ikGlobal, isp, tmes_file_exists)
      
    endif
    
    call MPI_BCAST(tmes_file_exists, 1, MPI_LOGICAL, root, intraPoolComm, ierr)
    
    if(.not. tmes_file_exists) then

      !-----------------------------------------------------------------------------------------------
      !> Read PW grids from `grid.ik` files

      allocate(gKIndexGlobalPC(npwsPC(ikGlobal)))
      allocate(gKIndexGlobalSD(npwsSD(ikGlobal)))

      if(indexInPool == 0) call readGrid('PC', ikGlobal, npwsPC(ikGlobal), gKIndexGlobalPC)
      if(indexInPool == 1) call readGrid('SD', ikGlobal, npwsSD(ikGlobal), gKIndexGlobalSD)

      call MPI_BCAST(gKIndexGlobalPC, size(gKIndexGlobalPC), MPI_INT, 0, intraPoolComm, ierr)
      call MPI_BCAST(gKIndexGlobalSD, size(gKIndexGlobalSD), MPI_INT, 1, intraPoolComm, ierr)

      
      !-----------------------------------------------------------------------------------------------
      !> Read wave functions and calculate overlap
      
      if(indexInPool == 0) write(*, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4)') ikGlobal, nKPoints
      call cpu_time(t1)
      
      allocate(wfcPC(nGVecsLocal, iBandIinit:iBandIfinal), wfcSD(nGVecsLocal, iBandFinit:iBandFfinal))
        
      call calculatePWsOverlap(ikLocal, isp)
        !! Read wave functions and get overlap
        
      call cpu_time(t2)
      if(indexInPool == 0) write(*, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> for k-point", i4, " done in", f10.2, " secs.")') ikGlobal, t2-t1
      call cpu_time(t1)

      !-----------------------------------------------------------------------------------------------
      !> Calculate cross projections
      
      allocate(cProjBetaPCPsiSD(nProjsPC, nBands, nSpins))

      call calculateCrossProjection('PC', iBandFinit, iBandFfinal, ikGlobal, nProjsPC, npwsPC(ikGlobal), gKIndexGlobalPC, &
            wfcSD, cProjBetaPCPsiSD)
        
      deallocate(wfcSD)
      deallocate(gKIndexGlobalPC)

      call cpu_time(t2)
      if(indexInPool == 0) write(*, '("      Calculating <betaPC|wfcSD> for k-point", i4, " done in", f10.2, " secs.")') ikGlobal, t2-t1
      call cpu_time(t1)
      
      allocate(cProjBetaSDPhiPC(nProjsSD, nBands, nSpins))

      call calculateCrossProjection('SD', iBandIinit, iBandIfinal, ikGlobal, nProjsSD, npwsSD(ikGlobal), gKIndexGlobalSD, &
            wfcPC, cProjBetaSDPhiPC)
        
      deallocate(wfcPC)
      deallocate(gKIndexGlobalSD)

      call cpu_time(t2)
      if(indexInPool == 0) write(*, '("      Calculating <betaSD|wfcPC> for k-point", i4, " done in", f10.2, " secs.")') ikGlobal, t2-t1


      !-----------------------------------------------------------------------------------------------
      !> Have process 0 in each pool calculate the PAW wave function correction for PC

      allocate(cProjPC(nProjsPC, nBands, nSpins))

      call readProjections('PC', ikGlobal, nProjsPC, cProjPC)

      if(indexInPool == 0) then

        call pawCorrectionWfc(nIonsPC, TYPNIPC, cProjPC, cProjBetaPCPsiSD, atomsPC, paw_PsiPC)

      endif
          
      deallocate(cProjBetaPCPsiSD)


      !-----------------------------------------------------------------------------------------------
      !> Have process 1 in each pool calculate the PAW wave function correction for PC

      allocate(cProjSD(nProjsSD, nBands, nSpins))

      call readProjections('SD', ikGlobal, nProjsSD, cProjSD)

      if(indexInPool == 1) then

        call pawCorrectionWfc(nIonsSD, TYPNISD, cProjBetaSDPhiPC, cProjSD, atoms, paw_SDPhi)

      endif

      deallocate(cProjBetaSDPhiPC)
      
      
      !-----------------------------------------------------------------------------------------------
      !> Broadcast wave function corrections to other processes

      call MPI_BCAST(paw_PsiPC, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
      call MPI_BCAST(paw_SDPhi, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)


      !-----------------------------------------------------------------------------------------------
      !> Have all processes calculate the PAW k correction

      call cpu_time(t1)

      allocate(pawKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nGVecsLocal))
      
      call pawCorrectionK('PC', nIonsPC, TYPNIPC, numOfTypesPC, posIonPC, atomsPC, atoms, pawKPC)
      

      call cpu_time(t2)
      if(indexInPool == 0) write(*, '("      <\\vec{k}|PAW_PC> for k-point ", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1
      call cpu_time(t1)
      

      allocate(pawSDK(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nGVecsLocal))
      
      call pawCorrectionK('SD', nIonsSD, TYPNISD, numOfTypes, posIonSD, atoms, atoms, pawSDK)

        
      call cpu_time(t2)
      if(indexInPool == 0) write(*, '("      <PAW_SD|\\vec{k}> for k-point ", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1
      

      !-----------------------------------------------------------------------------------------------
      !> Sum over PAW k corrections

      paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
      do ibi = iBandIinit, iBandIfinal
        
        do ibf = iBandFinit, iBandFfinal   
          paw_id(ibf,ibi) = sum(pawSDK(ibf,ibi,:)*pawKPC(ibf,ibi,:))
        enddo
        
      enddo

      Ufi(:,:,ikLocal,isp) = Ufi(:,:,ikLocal,isp) + paw_id(:,:)*16.0_dp*pi*pi/omega

      if(indexInPool == 0) then
        call MPI_REDUCE(MPI_IN_PLACE, Ufi(:,:,ikLocal,isp), size(Ufi(:,:,ikLocal,isp)), MPI_DOUBLE_COMPLEX, &
                MPI_SUM, 0, intraPoolComm, ierr)
      else
        call MPI_REDUCE(Ufi(:,:,ikLocal,isp), Ufi(:,:,ikLocal,isp), size(Ufi(:,:,ikLocal,isp)), MPI_DOUBLE_COMPLEX, &
                MPI_SUM, 0, intraPoolComm, ierr)
      endif


      if(indexInPool == 0) then 
        
        Ufi(:,:,ikLocal,isp) = Ufi(:,:,ikLocal,isp) + paw_SDPhi(:,:) + paw_PsiPC(:,:)
        
        call writeResults(ikLocal,isp)
        
      endif
      
      deallocate(cProjPC, pawKPC)
      deallocate(cProjSD, pawSDK)
      
    else
      
      if(indexInPool == 0) call readUfis(ikLocal,isp)
       
    endif
    
  enddo


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
  
  
  deallocate(paw_id)
  deallocate(paw_PsiPC)
  deallocate(paw_SDPhi)

  if(indexInPool == 0) deallocate(eigvI, eigvF)
  
  ! Calculating Vfi
  
  if (ionode) then
    
    ! Finalize Calculation
     
    call finalizeCalculation()
    
  endif
  
  call MPI_FINALIZE(ierr)
  
end program transitionMatrixElements
