module GVelMod
  
  use constants, only: dp
  use base, only: nKPoints
  use errorsAndMPI

  implicit none

  integer :: iBandInit, iBandFinal
    !! Initial and final bands
  integer :: nDegen
    !! Number of degenerate bands

  real(kind=dp) :: degenTol
    !! Tolerance for testing for degeneracies
  real(kind=dp), allocatable :: eigv(:,:,:)
    !! Eigenvalue for each band, direction
    !! (x/y/z), and position (left/middle/right)

end module GVelMod
