program transitionMatrixElements
  use declarations
  
  implicit none
  
  real(kind = dp) :: t1, t2
    !! For timing different processes

  integer :: ikLocal, ikGlobal
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


  allocate(gIndexLocalToGlobal(nGVecsGlobal), gVecProcId(nGVecsGlobal))

  call getFullPWGrid(nGVecsGlobal, gIndexGlobalToLocal, gVecProcId, mill_local, nGVecsLocal)
    !! Read the full PW grid from `mgrid` and distribute
    !! across processes



  allocate(paw_id(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal))
  if(ionode) then
    allocate ( Ufi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nKPoints) )
    allocate ( paw_SDKKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_PsiPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_SDPhi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( paw_fi(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal) )
    allocate ( eigvI (iBandIinit:iBandIfinal), eigvF (iBandFinit:iBandFfinal) )
    
    
    Ufi(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind = dp)
  endif
  
  do ikLocal = 1, nkPerPool
    
    ikGlobal = ikLocal+ikStart_pool-1
      !! Get the global `ik` index from the local one
    
    if(ionode) then
      
      tmes_file_exists = .false.
      call checkIfCalculated(ikGlobal,tmes_file_exists)
      
    endif
    
    call MPI_BCAST(tmes_file_exists, 1, MPI_LOGICAL, root, worldComm, ierr)
    
    if(.not. tmes_file_exists) then
      
      allocate(cProjPC(nProjsPC, nBands, nSpins))
      allocate(cProjSD(nProjsSD, nBands, nSpins))
      
      write(iostd, '(" Starting Ufi(:,:) calculation for k-point", i4, " of", i4)') ikGlobal, nKPoints
      flush(iostd)
        
      call cpu_time(t1)
      
      allocate(wfcPC(nGVecsLocal, iBandIinit:iBandIfinal), wfcSD(nGVecsLocal, iBandFinit:iBandFfinal))
        
      call calculatePWsOverlap(ikLocal)
        !! Read wave functions and get overlap
        
      call cpu_time(t2)
      write(iostd, '("      <\\tilde{Psi}_f|\\tilde{Phi}_i> for k-point", i2, " done in", f10.2, " secs.")') ikGlobal, t2-t1
      flush(iostd)

      if(ikGlobal == 1) then
        write(100+indexInPool, '(2i10)') iBandIinit, npwsPC(ikGlobal) 
        do ig = 1, npwsPC(ikGlobal)
          write(100+indexInPool,'(2ES24.15E3)') wfcPC(ig,iBandIinit)
        enddo
      endif

      call exitError('main', 'stopping for debugging', 5)
        

      call cpu_time(t1)
      write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> for k-point", i2, " begun.")') ikGlobal
      flush(iostd)
        
        call readProjectionsPC(ik)
        
        allocate ( cProjBetaPCPsiSD(nProjsPC, nBands, nSpins) )
        call projectBetaPCwfcSD(ik)
        
        deallocate ( wfcSD )
        
        call pawCorrectionPsiPC()
        
        deallocate ( cProjBetaPCPsiSD )
        
        call cpu_time(t2)
        write(iostd, '("      <\\tilde{Psi}_f|PAW_PC> done in", f10.2, " secs.")') t2-t1
        call cpu_time(t1)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> begun.")')
        flush(iostd)
        
        call readProjectionsSD(ik)
        
        allocate ( cProjBetaSDPhiPC(nProjsSD, nBands, nSpins) )
        call projectBetaSDwfcPC(ik)
        
        deallocate ( wfcPC )
        
        call pawCorrectionSDPhi()
        deallocate ( cProjBetaSDPhiPC )
        
        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\tilde{Phi}_i> done in", f10.2, " secs.")') t2-t1
        flush(iostd)
        
        call cpu_time(t1)
        write(iostd, '("      <\\vec{k}|PAW_PC> begun.")')
        flush(iostd)
      
      call MPI_BCAST(cProjPC, size(cProjPC), MPI_DOUBLE_COMPLEX, root, worldComm, ierr)
      call MPI_BCAST(cProjSD, size(cProjSD), MPI_DOUBLE_COMPLEX, root, worldComm, ierr)
      
      allocate ( pawKPC(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nPWsI(myid):nPWsF(myid)) )
      
      call pawCorrectionKPC()
      
      if (ionode) then
        call cpu_time(t2)
        write(iostd, '("      <\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
        
        call cpu_time(t1)
        write(iostd, '("      <PAW_SD|\\vec{k}> begun.")')
        flush(iostd)
        
      endif
      
      allocate ( pawSDK(iBandFinit:iBandFfinal, iBandIinit:iBandIfinal, nPWsI(myid):nPWsF(myid) ) )
      
      call pawCorrectionSDK()
      
      if (ionode) then
        
        call cpu_time(t2)
        write(iostd, '("      <PAW_SD|\\vec{k}> done in", f10.2, " secs.")') t2-t1
        
        call cpu_time(t1)
        write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> begun.")')
        flush(iostd)
        
      endif
      
      paw_id(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
      do ibi = iBandIinit, iBandIfinal
        
        do ibf = iBandFinit, iBandFfinal   
          paw_id(ibf,ibi) = sum(pawSDK(ibf,ibi,:)*pawKPC(ibf,ibi,:))
        enddo
        
      enddo
      
      if (ionode) paw_SDKKPC(:,:) = cmplx( 0.0_dp, 0.0_dp, kind = dp )
      
      CALL MPI_REDUCE(paw_id, paw_SDKKPC, size(paw_id), MPI_DOUBLE_COMPLEX, MPI_SUM, root, worldComm, ierr)
      
      if (ionode) then
        
        call cpu_time(t2)
        write(iostd, '("      \\sum_k <PAW_SD|\\vec{k}><\\vec{k}|PAW_PC> done in", f10.2, " secs.")') t2-t1
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

  deallocate(gIndexLocalToGlobal)
  deallocate(gVecProcId)
  deallocate(mill_local)
  
  
  if ( allocated ( paw_id ) ) deallocate ( paw_id )
  if (ionode) then
    if ( allocated ( paw_PsiPC ) ) deallocate ( paw_PsiPC )
    if ( allocated ( paw_SDPhi ) ) deallocate ( paw_SDPhi )
  endif
  
  ! Calculating Vfi
  
  if (ionode) then
    
    if (calculateVfis ) call calculateVfiElements()
    
    ! Finalize Calculation
     
    call finalizeCalculation()
    
  endif
  
  call MPI_FINALIZE(ierr)
  
end program transitionMatrixElements
