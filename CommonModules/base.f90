module base

  implicit none

  integer :: iBandIinit, iBandIfinal, iBandFinit, iBandFfinal
    !! Energy band bounds for initial and final state
  integer :: ispSelect
    !! Selection of a single spin channel if input
    !! by the user
  integer :: nKPoints
    !! Total number of k-points
  integer :: nSpins
    !! Max number of spins for both systems
  integer :: order
    !! Order of matrix element (0 or 1)

  logical :: loopSpins
    !! Whether to loop over available spin channels;
    !! otherwise, use selected spin channel

end module base
