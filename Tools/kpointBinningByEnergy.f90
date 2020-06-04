program kpointBinningByEnergy
  !! Take in eigenvalue files for each kpoint
  !! and energy bins from the VfisVsE file, bin
  !! each bands for each kpoint by energy, and
  !! output the results
  !!
  !! All energies are in Hartree
  !!
  !! <h2>Walkthrough</h2>

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision

  real(kind=dp), allocatable :: eBin(:)
    !! Energies for each bin
  real(kind=dp), allocatable :: DEHartree(:,:)
    !! Energy in Hartree

  integer, allocatable :: bandsBinnedByE(:)
    !! Number of kpoints in each bin
  integer :: iBandFinalState
    !! Band index for the final state
  integer :: iBandMax
    !! Max band index for the initial state
  integer :: nBins
    !! Number of energy bins
  integer :: nKpoints
    !! Total number of kpoints

  character(len=256) :: finalStateDir
    !! Path for final state eigenvalues
  character(len=256) :: initialStateDir
    !! Path for initial state eigenvalues

  call readInputs(nKpoints, nBins, finalStateDir, initialStateDir, iBandFinalState, iBandMax)
    !! * Read the input parameters

  allocate(eBin(nBins), DEHartree(iBandFinalState+1:iBandMax,nKpoints), bandsBinnedByE(nBins))

  write(*,*) "Reading eigenvalue files"

  call readEigenvalueFiles(nKpoints, finalStateDir, initialStateDir, iBandFinalState, DEHartree)
    !! * Read eigenvalue file for each kpoint

  write(*,*) "Done reading eigenvalue files"
  write(*,*) "Reading energy bins from VfisVsE file"

  call getEnergyBins(nBins, eBin)
    !! * Read in the energy bins from the VfisVsE file


  write(*,*) "Done reading energy bins"
  write(*,*) "Counting kpoints in each bin"

  call binBands(nBands, nBins, nKpoints, bandsBinnedByE, DEHartree, eBin)
    !! * Assign each band for each kpoint to the energy bins

  write(*,*) "Done binning kpoints"
  write(*,*) "Writing out results"

  call writeResults(nBins, eBin, bandsBinnedByE)
    !! * Write the bin energy and count to `bandsBinnedByE` file

  write(*,*) "Program complete"

  deallocate(eBin)
  deallocate(DEHartree)
  deallocate(bandsBinnedByE)
    ! Deallocate arrays

end program kpointBinningByEnergy

!==============================================================================

subroutine readInputs(nKpoints, nBins, finalStateDir, initialStateDir, iBandFinalState, iBandMax)
  !! Read input parameters

  implicit none

  integer :: iBandFinalState
    !! Band index for the final state
  integer :: iBandMax
    !! Max band index for the initial state
  integer :: nBins
    !! Number of energy bins
  integer :: nKpoints
    !! Total number of kpoints

  character(len=256) :: finalStateDir
    !! Path for final state eigenvalues
  character(len=256) :: initialStateDir
    !! Path for initial state eigenvalues

  namelist /input/ nKpoints, nBins, finalStateDir, initialStateDir, iBandFinalState, iBandMax

  write(*,*) "Reading input parameters"

  read(*,input)

  return
end subroutine readInputs

!------------------------------------------------------------------------------

subroutine readEigenvalueFiles(nKpoints, finalStateDir, initialStateDir, iBandFinalState, DEHartree)
  !! Read eigenvalue file for each kpoint

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision

  integer :: ib, ik
    !! Loop indices
  integer, intent(in) :: iBandFinalState
    !! Band index for the final state
  integer, intent(in) :: iBandMax
    !! Max band index for the initial state
  integer, intent(in) :: nKpoints
    !! Total number of kpoints

  real(kind=dp), intent(out) :: DEHartree(iBandFinalState+1:iBandMax,nKpoints)
    !! Energy difference between initial and
    !! final state in Hartree
  real(kind=dp) :: dummyD
    !! Dummy variable to ignore input
  real(kind=dp) :: eFinalState(nKpoints)
    !! Energy of final band in final state
    !! for each k-point
  real(kind=dp) :: eInitialState
    !! Band energy of initial state

  character(len=256) :: fileName
    !! String to hold eigenvalue file name 
    !! at each k-point
  character(len=256), intent(in) :: finalStateDir
    !! Path for final state eigenvalues
  character(len=256) :: finalStateFile
    !! File name with path for final state
  character(len=256), intent(in) :: initialStateDir
    !! Path for initial state eigenvalues
  character(len=256) :: initialStateFile
    !! File name with path for initial state

  do ik = 1, nKpoints
    ! Read in eigenvalues files for each kpoint

    if(ik < 10) then
      write(fileName, '("eigenvalues.", I1)') ik
    else
      write(fileName, '("eigenvalues.", I2)') ik
    endif

    initialStateFile = trim(trim(initialStateDir) // '/' // fileName)
    finalStateFile = trim(trim(finalStateDir) // '/' // fileName)

    open(unit=18,file=trim(finalStateFile))
      ! Open the eigenvalue file for this kpoint 
      ! for the final state
    
    read(18,*) 
    read(18,*) 
      ! Skip the first 2 blank lines

    do ib = 1, iBandFinalState-1
      ! Ignore all energy bands below final band

      read(18,*) 

    enddo

    read(18,*) eFinalState(ik), dummyD
      ! Read the band energy of the final state

    close(18)
      ! Close the eigenvalue file

    open(unit=18,file=trim(initialStateFile))
      ! Open the eigenvalue file for this kpoint 
      ! for the initial state
    
    read(18,*) 
    read(18,*) 
      ! Skip the first 2 blank lines

    do ib = 1, iBandFinalState
      ! Ignore all energy bands at or below final band

      read(18,*) 

    enddo

    do ib = iBandFinalState+1, iBandMax
      ! Read in all of the energy bands between the 
      ! bounds for this kpoint,

      read(18,*) eInitialState, dummyD

      DEHartree(ib,ik) = eInitialState - eFinalState(ik)

    enddo

    close(18)
      ! Close the eigenvalue file

  enddo

  return
end subroutine readEigenvalueFiles

!------------------------------------------------------------------------------

subroutine getEnergyBins(nBins, eBin)
  !! Read in the energy bins from the VfisVsE file

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision

  integer :: iE
    !! Loop indices
  integer, intent(in) :: nBins
    !! Number of energy bins

  real(kind=dp), intent(out) :: eBin(nBins)
    !! Energies for each bin
  real(kind=dp) :: dummyD
    !! Dummy variable to ignore input

  open(18,file="VfisVsE")
    ! Open VfisVsE file

  read(18,*)
  read(18,*)
  read(18,*)
  read(18,*)
  read(18,*)
  read(18,*)
    ! Ignore comments

  do iE = 1, nBins
    ! Read in the energies for each bin

    read(18,*) eBin(iE), dummyD, dummyD

  enddo

  close(18)

  return
end subroutine getEnergyBins

!------------------------------------------------------------------------------

subroutine binBands(iBandFinalState, iBandMax, nBins, nKpoints, bandsBinnedByE, DEHartree, eBin)
  !! Assign each band for each kpoint to the energy bins

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision

  integer :: ib, iE, ik
    !! Loop indices
  integer, intent(in) :: iBandFinalState
    !! Band index for the final state
  integer, intent(in) :: iBandMax
    !! Max band index for the initial state
  integer :: iEMin
    !! Minimum index for energy loop
  integer, intent(in) :: nBins
    !! Number of energy bins
  integer, intent(in) :: nKpoints
    !! Total number of kpoints

  integer, intent(out) :: bandsBinnedByE(nBins)
    !! Number of kpoints in each bin

  real(kind=dp), intent(in) :: eBin(nBins)
    !! Energies for each bin
  real(kind=dp), intent(in) :: DEHartree(iBandFinalState+1:iBandMax, nKpoints)
    !! Energy in Hartree

  bandsBinnedByE = 0
  iEMin = iBandFinalState + 1

  do ik = 1, nKpoints
    ! Choose a kpoint

    do ib = 2, nBins
      ! Go through each bin

      do iE = iEMin, iBandMax
        ! Cycle through bands starting at last index
        ! This starts at 1 but gets updated as you go 
        ! through bins because the energies are increasing
        ! so you don't need to go through each band
        ! every time

        if(DEHartree(iE,ik) >= eBin(ib-1) .AND. DEHartree(iE,ik) < eBin(ib)) then
          ! Increment the count in the previous bin if the band
          ! energy is less than the current bin energy boundary
        
          bandsBinnedByE(ib-1) = bandsBinnedByE(ib-1) + 1

        else if(ib == nBins) then
          ! If this is the last bin, all remaining energies 
          ! will be greater than the current bin boundary,
          ! so put the rest of the bands in this bin and exit
          
          bandsBinnedByE(ib) = bandsBinnedByE(ib) + nBands - iE + 1
          exit

        else if(DEHartree(iE,ik) >= eBin(ib)) then
          ! If you are not in the last bin and the band energy
          ! goes higher than the current bin boundary, update
          ! the minimum band index and exit to move on to the 
          ! next bin

          iEMin = iE
          exit

        endif

      enddo
    enddo

    iEMin = iBandFinalState + 1
      ! Once you go through all of the bins and are moving
      ! on to another kpoint, restart the minimum band index

  enddo

  return
end subroutine binBands

!------------------------------------------------------------------------------

subroutine writeResults(nBins, eBin, bandsBinnedByE)
  !! Write the bin energy and count to `bandsBinnedByE` file

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision

  integer :: ib
    !! Loop index
  integer, intent(in) :: nBins
    !! Number of energy bins

  integer, intent(out) :: bandsBinnedByE(nBins)
    !! Number of kpoints in each bin

  real(kind=dp), intent(in) :: eBin(nBins)
    !! Energies for each bin

  open(18,file="bandsBinnedByE")

  write(18,*) "# Energy bands for each kpoint binned by energy"
  write(18,*) "# Energy bin (Hartree), Number of bands/kpoints in bin"

  do ib = 1, nBins
    ! Write the energy and count for each bin

    write(18,*) eBin(ib), bandsBinnedByE(ib)

  enddo

  close(18)

  return
end subroutine writeResults
