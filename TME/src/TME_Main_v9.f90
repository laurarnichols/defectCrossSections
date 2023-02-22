program transitionMatrixElements
  use declarations
  
  implicit none
  
  real(kind = dp) :: t1, t2
    !! For timing different processes

  integer :: ikLocal, ikGlobal, iType
    !! Loop indices
  
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
  if(ionode) then
    allocate ( paw_SDKKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_fi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( eigvI (iBandIinit:iBandIfinal), eigvF (iBandFinit:iBandFfinal) )
    
    
  endif

  allocate(Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nKPerPool))
  Ufi(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
  
  do ikLocal = 1, nkPerPool
    
    ikGlobal = ikLocal+ikStart_pool-1
      !! Get the global `ik` index from the local one
    
    if(ionode) then
      
      tmes_file_exists = .false.
      call checkIfCalculated(ikGlobal,tmes_file_exists)
      
    endif
    
    call MPI_BCAST(tmes_file_exists, 1, MPI_LOGICAL, root, worldComm, ierr)
    
    if(.not. tmes_file_exists) then
      
      !-----------------------------------------------------------------------------------------------
      !> Read wave functions and calculate overlap
      
      if(indexInPool == 0) write(iostd, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4)') ikGlobal, nKPoints
        
      call cpu_time(t1)
      
      allocate(wfcPC(nGVecsLocal, iBandIinit:iBandIfinal), wfcSD(nGVecsLocal, iBandFinit:iBandFfinal))
        
      call calculatePWsOverlap(ikGlobal)
        !! Read wave functions and get overlap
        
      call cpu_time(t2)
      if(indexInPool == 0) write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> for k-point", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1

      !-----------------------------------------------------------------------------------------------
      !> Read projections
      
      allocate(cProjPC(nProjsPC, nBands, nSpins))

      call readProjections('PC', ikGlobal, nProjsPC, cProjPC)

      allocate(cProjSD(nProjsSD, nBands, nSpins))

      call readProjections('SD', ik, nProjsSD, cProjSD)
        

      !-----------------------------------------------------------------------------------------------
      !> Calculate cross projections

      allocate(cProjBetaPCPsiSD(nProjsPC, nBands, nSpins))

      call calculateCrossProjection('PC', iBandFinit, iBandFfinal, ikGlobal, nProjsPC, wfcSD, cProjBetaPCPsiSD)
        
      deallocate(wfcSD)

      allocate(cProjBetaSDPhiPC(nProjsSD, nBands, nSpins))

      call calculateCrossProjection('SD', iBandIinit, iBandIfinal, ikGlobal, nProjsSD, wfcPC, cProjBetaSDPhiPC)
        
      deallocate(wfcPC)


      !-----------------------------------------------------------------------------------------------
      !> Have process 0 in each pool calculate the PAW wave function correction for PC
      if(indexInPool == 0) then

        call cpu_time(t1)
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> for k-point ", i2, " begun.")') ikGlobal
        flush(iostd)

        call pawCorrectionWfc(nIonsPC, TYPNIPC, cProjPC, cProjBetaPCPsiSD, atomsPC, paw_PsiPC)

        call cpu_time(t2)
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> for k-point ", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1

      endif
          
      deallocate(cProjBetaPCPsiSD)


      !-----------------------------------------------------------------------------------------------
      !> Have process 1 in each pool calculate the PAW wave function correction for PC
      if(indexInPool == 1) then

        call cpu_time(t1)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> for k-point ", i2, " begun.")') ikGlobal
        flush(iostd)

        call pawCorrectionWfc(nIonsSD, TYPNISD, cProjBetaSDPhiPC, cProjSD, atoms, paw_SDPhi)

        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> done in", f10.2, " secs.")') t2-t1
        flush(iostd)

      endif

      deallocate(cProjBetaSDPhiPC)
      
      
      !-----------------------------------------------------------------------------------------------
      !> Broadcast wave function corrections to other processes

      call MPI_BCAST(paw_PsiPC, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 0, intraPoolComm, ierr)
      call MPI_BCAST(paw_SDPhi, size(paw_PsiPC), MPI_DOUBLE_COMPLEX, 1, intraPoolComm, ierr)


      !-----------------------------------------------------------------------------------------------
      !> Have all processes calculate the PAW k correction

      call cpu_time(t1)
      if(indexInPool == 0) write(iostd, '("      <\\vec{k}|PAW_PC> for k-point ", i2, " begun.")') ikGlobal
      

      allocate(pawKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nGVecsLocal))
      
      call pawCorrectionK('PC', nIonsPC, TYPNIPC, numOfTypesPC, posIonPC, atomsPC, pawKPC)
      

      call cpu_time(t2)
      if(indexInPool == 0) write(iostd, '("      <\\vec{k}|PAW_PC> for k-point ", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1
        
      call cpu_time(t1)
      if(indexInPool == 0) write(iostd, '("      <PAW_SD|\\vec{k}> for k-point ", i2, " begun.")') ikGlobal
      

      allocate(pawSDK(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nGVecsLocal))
      
      call pawCorrectionK('SD', nIonsSD, TYPNISD, numOfTypes, posIonSD, atoms, pawSDK)

        
      call cpu_time(t2)
      if(indexInPool == 0) write(iostd, '("      <PAW_SD|\\vec{k}> for k-point ", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1
      

      !-----------------------------------------------------------------------------------------------
      !> Sum over PAW k corrections
        
      call cpu_time(t1)
      if(indexInPool == 0) write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> begun.")')
      

      paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
      do ibi = iBandIinit, iBandIfinal
        
        do ibf = iBandFinit, iBandFfinal   
          paw_id(ibf,ibi) = sum(pawSDK(ibf,ibi,:)*pawKPC(ibf,ibi,:))
        enddo
        
      enddo

      
      if(indexInPool == 0) paw_SDKKPC(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
      call MPI_REDUCE(paw_id, paw_SDKKPC, size(paw_id), MPI_DOUBLE_COMPLEX, MPI_SUM, root, intraPoolComm, ierr)
      

      call cpu_time(t2)
      if(indexInPool == 0) then write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
        
        Ufi(:,:,ik) = Ufi(:,:,ik) + paw_SDPhi(:,:) + paw_PsiPC(:,:) + paw_SDKKPC(:,:)*16.0_dp*pi*pi/omega
        
        call writeResults(ik)
        
      endif
      
      deallocate ( cProjPC, pawKPC )
      deallocate ( cProjSD, pawSDK )
      
    else
      
      if (ionode) call readUfis(ik)
       
    endif
    
  enddo

  deallocate(npwsPC)
  deallocate(wkPC)
  deallocate(xkPC)
  deallocate(posIonPC)
  deallocate(TYPNIPC)

  do iType = 1, numOfTypesPC
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

  do iType = 1, numOfTypesPC
    deallocate(atoms(iType)%lps)
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
  
  ! Calculating Vfi
  
  if (ionode) then
    
    if (calculateVfis ) call calculateVfiElements()
    
    ! Finalize Calculation
     
    call finalizeCalculation()
    
  endif
  
  call MPI_FINALIZE(ierr)
  
end program transitionMatrixElements
