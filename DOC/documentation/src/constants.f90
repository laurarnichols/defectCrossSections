module constants
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
    !! Used to make reals double precision
  integer, parameter :: iostd = 16
    !! Unit number for output file
  !
  real(kind = dp), parameter ::         abCM = 0.529177219217e-8_dp
  real(kind = dp), parameter :: eVToHartree  = 1.0_dp/27.21138386_dp
  real(kind = dp), parameter :: HartreeToEv  = 27.21138386_dp
  real(kind = dp), parameter ::           pi = 3.1415926535897932_dp
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  real(kind = dp), parameter ::        twopi = 2.0_dp*pi
  !
end module constants
