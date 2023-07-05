module base

  implicit none

  integer :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !! Energy band bounds for initial and final state
  integer :: nKPoints
    !! Total number of k-points
  integer :: nSpins
    !! Max number of spins for both systems
  integer :: order
    !! Order of matrix element (0 or 1)

end module base
