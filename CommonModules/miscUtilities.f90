module miscUtilities
  !
  implicit none
  !
  integer, parameter :: dp    = selected_real_kind(15, 307)
    !! Used to make reals double precision
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
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine arrangeLargerToSmaller(nModes, x, s2L)
    !! Sort `s2L` based on descending order
    !! of `x`
    !!
    !! @todo Change this to a more efficient algorithm @endtodo
    !
    implicit none
    !
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(inout) :: s2L(nModes)
      !! Indexes phonon information that come in
      !! in order (1, 2, 3, ...) and are rearranged
      !! here so that they are the indices in order
      !! of the magnitude of the argument `x`
    integer :: i, iMode
    !
    real(kind = dp), intent(in) :: x(nModes)
      !! Argument to modified Bessel function
    real(kind = dp), allocatable :: temp(:)
    real(kind = dp) :: tmpr
    integer :: tmpi
    !
    allocate( temp(nModes) )
    !
    temp(:) = 0.0_dp
    temp(:) = x(:)
    !
    do iMode = 1, nModes
      !
      do i = 1, nModes-1
        !
        if ( temp(i) < temp(i+1) ) then
          !
          tmpi = s2L(i)
          s2L(i) = s2L(i+1)
          s2L(i+1) = tmpi
          !
          tmpr = temp(i)
          temp(i) = temp(i+1)
          temp(i+1) = tmpr
          !
        endif
        !
      enddo
      !
    enddo
    !
    deallocate ( temp )
    !
    return
    !
  end subroutine arrangeLargerToSmaller
  !
end module miscUtilities
