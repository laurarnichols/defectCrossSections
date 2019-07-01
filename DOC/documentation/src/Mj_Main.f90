program MjME
  !
  use MjModule
  !
  implicit none
  !
  call cpu_time(ti)
  !
  ! Reading input, initializing and checking all variables of the calculation.
  !
  call readInputs()
  !
  call computeGeneralizedDisplacements()
  !
  call computeVariables()
  !
  call displaceAtoms()
  !
  if ( readQEInput ) then
    call exportQEInput()
  else
    call writeNewAtomicPositions()
  endif
  !
end program MjME
