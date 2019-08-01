module constants
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
    !! Used to make reals double precision
  integer, parameter :: iostd = 16
    !! Unit number for output file
  !
  real(kind = dp), parameter :: THzToHartree = 1.0_dp/6579.683920729_dp
  !
end module constants
