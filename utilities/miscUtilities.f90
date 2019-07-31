module miscUtilities
  !
  implicit none
  !
  contains
  !
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine int2str(integ, string)
    !! Write a given integer to a string, using only as many digits as needed
    !
    implicit none
    integer :: integ
    character(len = 300) :: string
    !
    if ( integ < 10 ) then
      write(string, '(i1)') integ
    else if ( integ < 100 ) then
      write(string, '(i2)') integ
    else if ( integ < 1000 ) then
      write(string, '(i3)') integ
    else if ( integ < 10000 ) then
      write(string, '(i4)') integ
    endif
    !
    string = trim(string)
    !
    return
    !
  end subroutine int2str
  !
end module miscUtilities
