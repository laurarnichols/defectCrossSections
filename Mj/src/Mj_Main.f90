program MjME
  !! Calculate some of the main variables in equations 
  !! 42 and 43, calculate new displacements of atoms,
  !! and generate new QE input based on the new positions
  !!
  !! <h2>Walkthrough</h2>
  !!
  use MjModule
  !
  implicit none
  !
  call cpu_time(ti)
    !! @todo Make sure that there is an end timer @endtodo
  !
  call readInputs()
    !! * Initialize variables, read input, and check that 
    !!   all required variables were read and have values 
    !!   that make sense
  !
  call computeGeneralizedDisplacements()
    !! * Calculate \delta q_j^2
  !
  call computeVariables()
    !! * Compute main parts of equations 42 and 43 in paper
    !!   to make whole formula more manageable
  !
  call displaceAtoms()
    !! * Calculate new positions of atoms for each mode 
    !!   based on `maxDisplacement` and `phonD`
  !
  !> * If a base QE input was read, then export the new QE 
  !>   input files with the new postions
  !> * Otherwise, just write out the positions
  if ( readQEInput ) then
    !
    call exportQEInput()
    !
  else
    !
    call writeNewAtomicPositions()
    !
  endif
  !
end program MjME
