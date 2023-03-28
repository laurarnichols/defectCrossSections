module miscUtilities
  
  use constants
  
  implicit none
  
  contains
  
!----------------------------------------------------------------------------
  subroutine int2str(integ, string)
    
    implicit none
    integer :: integ
    character(len = 300) :: string
    
    if ( integ < 10 ) then
      write(string, '(i1)') integ
    else if ( integ < 100 ) then
      write(string, '(i2)') integ
    else if ( integ < 1000 ) then
      write(string, '(i3)') integ
    else if ( integ < 10000 ) then
      write(string, '(i4)') integ
    else if ( integ < 100000 ) then
      write(string, '(i5)') integ
    else if ( integ < 1000000 ) then
      write(string, '(i6)') integ
    else if ( integ < 10000000 ) then
      write(string, '(i7)') integ
    endif
    
    string = trim(string)
    
    return
    
  end subroutine int2str
  
!----------------------------------------------------------------------------
  subroutine int2strLeadZero(integ, length, string)
    
    implicit none
    integer :: integ
    integer :: length
    character(len = 300) :: string
    character(len = 300) :: lengthC
    character(len = 300) :: formatString
    
    call int2str(length, lengthC)
    formatString = '(i'//trim(lengthC)//'.'//trim(lengthC)//')'

    write(string, formatString) integ
    
    string = trim(string)
    
    return
    
  end subroutine int2strLeadZero
  
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine arrangeLargerToSmaller(nModes, x, s2L)
    !! Sort `s2L` based on descending order
    !! of `x`
    !!
    !! @todo Change this to a more efficient algorithm @endtodo
    
    implicit none
    
    integer, intent(in) :: nModes
      !! Number of phonon modes
    integer, intent(inout) :: s2L(nModes)
      !! Indexes phonon information that come in
      !! in order (1, 2, 3, ...) and are rearranged
      !! here so that they are the indices in order
      !! of the magnitude of the argument `x`
    integer :: i, iMode
    
    real(kind = dp), intent(in) :: x(nModes)
      !! Argument to modified Bessel function
    real(kind = dp) :: temp(nModes)
    real(kind = dp) :: tmpr
    integer :: tmpi
    
    temp(:) = 0.0_dp
    temp(:) = x(:)
    
    do iMode = 1, nModes
      
      do i = 1, nModes-1
        
        if ( temp(i) < temp(i+1) ) then
          
          tmpi = s2L(i)
          s2L(i) = s2L(i+1)
          s2L(i+1) = tmpi
          
          tmpr = temp(i)
          temp(i) = temp(i+1)
          temp(i+1) = tmpr
          
        endif
        
      enddo
      
    enddo
    
    return
    
  end subroutine arrangeLargerToSmaller
  
  
  function findloc(stringArr, searchString) result(loc)
    !! Find the first instance of a search string in
    !! an array of strings
    !!
    implicit none
    
    integer :: i
    integer :: loc
    
    character(len=*), intent(in) :: searchString
    character(len=*), intent(in) :: stringArr(:)
    
    loc = -1
    
    do i = 1, size(stringArr)
      
      if (stringArr(i) == searchString) loc = i   
      
    enddo
    
  end function findloc
  
end module miscUtilities
