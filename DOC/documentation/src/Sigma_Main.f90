program crossSection
  !
  use sigma_module
  !
  implicit none
  !
  ! Reading input, initializing and checking all variables of the calculation.
  !
  call readInputs()
  !
  call calculateSigma()
  !
  call writeSigma()
  !
end program crossSection
