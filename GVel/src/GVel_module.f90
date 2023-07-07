module GVelMod
  
  use constants, only: dp
  use base, only: nKPoints
  use errorsAndMPI

  implicit none

  integer :: iBandInit, iBandFinal
    !! Initial and final bands

  real(kind=dp), allocatable :: eigv(:,:,:)
    !! Eigenvalue for each band, position 
    !! (left/middle/right), and direction 
    !! (x/y/z)

end module GVelMod
