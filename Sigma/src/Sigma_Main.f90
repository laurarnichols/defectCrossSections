program crossSection
  !! Combine the results of TME and LSF
  !! into the cross section
  !
  use sigma_module
  !
  implicit none
  !
  call readInputs()
    !! * Read input, initializing and checking all variables of the calculation
  !
  call calculateSigma()
    !! * Calculate \(\sigma(E)\)
  !
  call writeSigma()
    !! * Output \(\sigma(E)\)
  !
end program crossSection
