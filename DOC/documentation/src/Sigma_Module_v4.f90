
module sigma_module
  !
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
  integer, parameter :: int32 = selected_int_kind(5)
  integer, parameter :: iostd = 16
  integer,         parameter :: nOfEnergyBins = 5040
  !
  !
  real(kind = dp), parameter ::         abCM = 0.529177219217e-8_dp
  real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
  real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
  real(kind = dp), parameter ::    maxEnergy = 10.0_dp
  real(kind = dp), parameter ::           pi = 3.1415926535897932_dp
  real(kind = dp), parameter ::        twopi = 2.0_dp*pi
  !
  character(len = 11), parameter :: output = 'sigmaStatus'
  !
  integer(kind = int32) :: ios
  integer :: m
  integer :: nEnergies
  integer :: numOfVfis
  !
  real(kind = dp) :: volume
  real(kind = dp) :: de
  real(kind = dp) :: eifMin
  real(kind = dp) :: DHifMin
  !
  real(kind = dp), allocatable :: E(:)
  real(kind = dp), allocatable :: energy(:)
  real(kind = dp), allocatable :: lorentz(:)
  real(kind = dp), allocatable :: lorentzByPhonon(:)
  real(kind = dp), allocatable :: lsf(:)
  real(kind = dp), allocatable :: lsfVsE(:)
  real(kind = dp), allocatable :: lsfVsEbyPhonon(:)
  real(kind = dp), allocatable :: sigma(:)
  real(kind = dp), allocatable :: sigmaByPhonon(:)
  real(kind = dp), allocatable :: Vfis(:)
  !
  character(len = 256) :: VfisInput
  character(len = 256) :: LSFinput
  character(len = 256) :: crossSectionOutput
  !
  logical :: file_exists
  !
  namelist /elphscat/ VfisInput, LSFinput, crossSectionOutput
  !
  !
contains
  !
  !
  subroutine readInputs()
    !! Read input parameters and read LSF and TME output
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    !> * Check if an output file exists; if it does, delete it
    inquire(file = output, exist = file_exists)
    !
    if ( file_exists ) then
      !
      open (unit = 11, file = output, status = "old")
      !
      close(unit = 11, status = "delete")
      !
    endif
    !
    open (iostd, file = output, status='new')
      !! * Open new output file
    !
    call initialize()
      !! * Set default values of input parameters
    !
    READ (5, elphscat, iostat = ios)
      !! * Read in input parameters
    !
    call checkInputAndUpdateParameters()
      !! * Check if input parameters were updated
      !!   and do some basic checks
    !
    call readLSF()
      !! * Read the LSF output
    !
    call readVfis()
      !! * Read the TME output
    !
    return
    !
  end subroutine readInputs
  !
  !
  subroutine initialize()
    !! Set default values for input parameters  
    !
    implicit none
    !
    VfisInput = ''
    LSFinput = ''
    crossSectionOutput = ''
    !
    de = maxEnergy*eVToHartree/real(nOfEnergyBins, dp)
    !
    return
    !
  end subroutine initialize
  !
  !
  subroutine checkInputAndUpdateParameters()
    !! Check that the input variables don't still have their default
    !! values. If the input files are not defined, end the program.
    !
    implicit none
    !
    if ( VfisInput == '' ) then
      write(iostd, '(" Vfi elements input (input variable VfisInput) is not defined!")')
    else
      inquire(file =trim(VfisInput), exist = file_exists)
      if ( file_exists ) then
        write(iostd, '(" Vfi elements input : ", a)') trim(VfisInput)
      else
        write(iostd, '(" Vfi elements input : ", a, " does not exist!")') trim(VfisInput)
      endif
    endif
    !
    if ( LSFinput == '' ) then
      write(iostd, '(" LSF input (input variable LSFinput) is not defined!")')
    else
      inquire(file =trim(LSFinput), exist = file_exists)
      if ( file_exists ) then
        write(iostd, '(" LSF input : ", a)') trim(LSFinput)
      else
        write(iostd, '(" LSF input : ", a, " does not exist!")') trim(LSFinput)
      endif
    endif
    !
    if ( crossSectionOutput == '' ) then
      write(iostd, '(" crossSectionOutput is not defined! File name : crossSection, will be used.")')
      crossSectionOutput = 'crossSection'
    else
      write(iostd, '(" Cross section output file name : ", a)') trim(crossSectionOutput)
    endif
    !
    if ( ( VfisInput == '' ) .or. ( LSFinput == '' ) ) then
      !
      write(iostd, '(" One or both of the input files is not defined! ")')
      write(iostd, '(" ********************************************** ")')
      write(iostd, '(" *               Program stops!               * ")')
      write(iostd, '(" *       Please check the output file.        * ")')
      write(iostd, '(" ********************************************** ")')
      !
      stop
      !
    endif
    !
    flush(iostd)
    !
    return
    !
  end subroutine checkInputAndUpdateParameters
  !
  !
  subroutine readLSF()
    !! Read LSF output
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    character(len = 1) :: dummyC1
    character(len = 8) :: dummyC8
    character(len = 9) :: dummyC9
      !! @todo Merge these dummy characters @endtodo
    !
    real(kind = dp) :: ee
      !! Energy in eV
    !
    integer :: iE
      !! Loop index over energies
    !
    open(1, file=trim(LSFinput), status='old')
      !! * Open the LSF output file
    !
    read(1,'(a1, i10, a9, i5, a8)') dummyC1, nEnergies, dummyC9, m, dummyC8
      !! * Read in the number of energies and ??
    !
    allocate ( E(-nEnergies:nEnergies), lsfVsE(-nEnergies:nEnergies), lsfVsEbyPhonon(-nEnergies:nEnergies) )
    !
    do iE = -nEnergies, nEnergies
      !! * For each energy, read in the energy value as well as `lsfVsE` and `lsfVsEbyPhonon`
      !
      read(1,'(F16.8,2E18.6e3)') ee, lsfVsE(iE), lsfVsEbyPhonon(iE)
      E(iE) = ee*eVToHartree
      !
    enddo
    !
    close(1)
    !
  end subroutine readLSF
  !
  !
  subroutine readVfis()
    !! Read TME output, get `Vfis` and `lsf` using the same index,
    !! and output the results to two output files
    !!
    !! @note
    !! This subroutine does not line up with what is output from the TME program,
    !! so I'm not 100% sure where the input is coming from or what it's supposed
    !! to be
    !! @endnote
    !!
    !! <h2>Walkthrough</h2>
    !!
    implicit none
    !
    integer :: i, iE0, iE
    real(kind = dp) :: dummyD1, dummyD2, Ee, VfiOfE, VfiOfE0, eBin
    character(len =  1) :: dummyC1
    character(len = 32) :: dummyC32
    character(len = 35) :: dummyC35
      !! @todo Merge dummy variables @endtodo
    !
    open(1, file=trim(VfisInput), status="old")
    !
    !read(1, '(a1, i10, a9, f15.4, a16)') dummyC1, nEVfi, dummyC9, volume, dummyC16
    !
    read(1, *)
      !! * Ignore the first line as it is a comment
    !
    read(1, '(a32, ES24.15E3, a35)') dummyC32, volume, dummyC35
      !! * Read cell volume
    !
    read(1, '(a32, ES24.15E3, a35)') dummyC32, eifMin, dummyC35
      !! * Read minimum transition energy
    !
    read(1, '(a32, ES24.15E3, a35)') dummyC32, DHifMin, dummyC35
      !! * Read \(|\Delta H_{if}|^2\) at min transition energy
    !
    read(1, '(a32, ES24.15E3, a35)') dummyC32, eBin, dummyC35
      !! * Read energy bin size
    !
    read(1, *) 
      !! * Ignore the next line as it is a comment
    !
    read(1, '(i10)') numOfVfis
      !! * Read in `nOfEnergies`? 
    !
    allocate ( Vfis(0:numOfVfis), energy(0:numOfVfis), lsf(0:numOfVfis) )
    !
    Vfis(:) = 0.0_dp
    energy(:) = 0.0_dp
    lsf(:) = 0.0_dp
    !
    read(1, '(3ES24.15E3)' ) Ee, VfiOfE0, dummyD1
    Vfis(1) = VfiOfE0
    energy(1) = Ee
      !! * Read in the initial values of energy and `Vfis`
    !
    iE = int(Ee/de) + 1
      !! * Calculate the energy index
    !
    do i = 2, numOfVfis
      !! * For each energy
      !!    * Read in energy and `VFiOfE`
      !!    * Calculate the indices needed to get `VFis`
      !!      and `lsfVsE` using the same index
      !!    * Average `lsf` over those indices and store 
      !!      in a single index matching that of `Vfis`
      !
      iE0 = iE ! int(energy(i-1)/deltaE) + 1 !  iE
      !
      read(1, '(3ES24.15E3)') Ee, VfiOfE, dummyD2
      !
      energy(i) = Ee
      !
      iE = int(Ee/de) + 1
      !
      !Vfis(iE0:iE) = VfiOfE0
      Vfis(i) = VfiOfE
      !VfiOfE0 = VfiOfE
      !
      lsf(i-1) = sum(lsfVsE(iE0:iE))/(iE-iE0+1)
      !
      write(26,*) E(iE0), Ee, lsf(i) ! sum(lsfVsE(iE0:iE))/(iE-iE0+1)
        !! @todo Figure out where this file is opened @endtodo
      !
    enddo
    !
    close(1)
    close(26)
    !
    do iE = 0, numOfVfis ! -nEnergies, nEnergies
      !! * For each energy, output the energy in eV, `Vfis` and `lsf`
      !
      write(44,*) energy(iE)*HartreeToEv, Vfis(iE), lsf(iE)
        !! @todo Figure out where this file is opened @endtodo
        !
    enddo
    !
    close(44)
    !
    return
    !
  end subroutine readVfis
  !
  !
  subroutine calculateSigma()
    !
    implicit none
    !
    integer :: iE
    real(kind = dp) :: vg, sigma0
    !
    allocate( sigma(numOfVfis) ) ! , sigmaByPhonon(-nEnergies:nEnergies) )
    !allocate( sigma(-nEnergies:nEnergies), sigmaByPhonon(-nEnergies:nEnergies) )
    !
    iE = int(eifMin/de) + 1
    write(6,*) eifMin, eifMin*HartreeToEv, iE
    sigma0 = twoPi*abCM**2*volume*DHifMin*lsfVsE(iE)/sqrt(2.0_dp*E(iE))
    !
    !do iE = 1, numOfVfis ! -nEnergies, nEnergies - 1
    !  if ( (E(iE) < eifMin).and.(E(iE+1) > eifMin) ) sigma0 = twoPi*abCM**2*volume*DHifMin*lsfVsE(iE)/sqrt(2.0_dp*E(iE))
    !enddo
    !
    write(6,*) eifMin*HartreeToEv, sigma0
    !
    sigma(:) = 0.0_dp
    !
    do iE = 0, numOfVfis ! -nEnergies, nEnergies
      vg = 1.0_dp
      if ( energy(iE) > 0.0_dp ) vg = sqrt(2.0_dp*energy(iE))
      !write(6,*) iE, energy(iE), vg, Vfis(iE), lsf(iE)
      sigma(iE)         = twoPi*abCM**2*volume*Vfis(iE)*lsf(iE)/vg
      !sigma(iE)         = twoPi*abCM**2*volume*Vfis(iE)*lsfVsE(iE)/vg
      !sigmaByPhonon(iE) = twoPi*abCM**2*volume*Vfis(iE)*lsfVsEbyPhonon(iE)/vg
    enddo
    !
    return
    !
  end subroutine calculateSigma
  !
  !
  subroutine writeSigma()
    !
    implicit none
    !
    integer :: iE
    !
    open(2, file=trim(crossSectionOutput), status='unknown')
    !
    do iE = 0, numOfVfis ! -nEnergies, nEnergies
      !
      write(2,*) energy(iE)*HartreeToEv, sigma(iE)!, sigmaByPhonon(iE)
      !write(2,*) E(iE), sigma(iE), sigmaByPhonon(iE)
      !
    enddo
    !
    close(2)
    !
    return
    !
  end subroutine writeSigma
  !
  !
end module sigma_module
