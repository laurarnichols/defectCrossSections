module constants
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
    !! Used to make reals double precision
  integer, parameter :: iostd = 16
    !! Unit number for output file
  !
  real(kind = dp), parameter ::         abCM = 0.529177219217e-8_dp
  real(kind = dp), parameter :: angToBohr = 1.889725989_dp
    !! Conversion factor from Angstrom to Bohr
  real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
    !! Conversion factor from eV to Hartree
  real(kind = dp), parameter :: eVToRy = 1.0_dp/13.60569301_dp
    !! Conversion factor from eV to Rydberg
  real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
    !! Conversion factor from Hartree to eV
  real(kind = dp), parameter ::           pi = 3.1415926535897932_dp
  real(kind = dp), parameter :: ryToHartree = 0.5_dp
    !! Conversion factor from Rydberg to Hartree
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  real(kind = dp), parameter ::        twopi = 2.0_dp*pi
  !
end module constants
