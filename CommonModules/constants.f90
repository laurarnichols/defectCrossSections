module constants
  implicit none
  
  integer, parameter :: dp    = selected_real_kind(15, 307)
    !! Used to make reals double precision
  integer, parameter :: iostd = 6
    !! Unit number for output file
  
  real(kind=dp), parameter ::       angToBohr = 1.889725989_dp
  real(kind=dp), parameter ::          angToM = 1e-10_dp
  real(kind=dp), parameter ::     BohrToMeter = 5.29177e-11_dp
  real(kind=dp), parameter ::       elecMToKg = 9.10938e-31_dp
  real(kind=dp), parameter ::     eVToHartree = 1.0_dp/27.21138386_dp
  real(kind=dp), parameter ::           eVToJ = 1.6021766208e-19_dp
  real(kind=dp), parameter ::          eVToRy = 1.0_dp/13.60569301_dp
  real(kind=dp), parameter ::   daltonToElecM = 1.66053907e-27_dp/9.1093837e-31_dp
  real(kind=dp), parameter ::     HartreeToEv = 27.21138386_dp
  real(kind=dp), parameter ::      HartreeToJ = 4.359744650e-18_dp
  real(kind=dp), parameter ::            hbar = 1.0545718e-34_dp
    !! J*s
  real(kind=dp), parameter ::     hbar_atomic = 1.0_dp
  real(kind=dp), parameter ::           kB_SI = 1.38064852e-23_dp
  real(kind=dp), parameter ::       kB_atomic = 3.1668e-6_dp 
    !! Boltzmann constant in different units
  real(kind=dp), parameter ::              pi = atan(1.0_dp)*4.0_dp
  real(kind=dp), parameter ::     ryToHartree = 0.5_dp
  real(kind=dp), parameter ::    THzToHartree = 2.4188843265864e-5_dp
    !! Still frequency, not energy; just Hartree units
  real(kind=dp), parameter ::         THzToHz = 1e12_dp
  real(kind=dp), parameter :: time_atomicToSI = 2.4188843265864e-17_dp
  real(kind=dp), parameter ::           twopi = 2.0_dp*pi

  complex(kind = dp), parameter ::         ii = cmplx(0.0_dp, 1.0_dp, kind=dp)
  
end module constants
