module miscUtilities
  
  use constants
  
  implicit none
  
  contains
  
!----------------------------------------------------------------------------
  function int2str(integ) result(string)
    
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
    
  end function int2str
  
!----------------------------------------------------------------------------
  function int2strLeadZero(integ, length) result(string)
    
    implicit none
    integer :: integ
    integer :: length
    character(len = 300) :: string
    character(len = 300) :: lengthC
    character(len = 300) :: formatString
    
    formatString = '(i'//trim(int2str(length))//'.'//trim(int2str(length))//')'

    write(string, formatString) integ
    
    string = trim(string)
    
  end function int2strLeadZero

!----------------------------------------------------------------------------
  function getFirstLineWithKeyword(fUnit, keyword) result(line)
    !! Read from file until ecounter keyword

    implicit none

    ! Input variables:
    integer, intent(in) :: fUnit
      !! Unit for file to read from
    
    character(*) :: keyword
      !! Keyword to search for

    ! Output variables:
    character(len=300) :: line
      !! First line with keyword

    ! Local variables:
    logical :: found
      !! If the line has been found


    found = .false.
    do while (.not. found)
      !! * Ignore everything until you get to a
      !!   line with keyword

      read(fUnit, '(A)') line

      if (index(line,keyword) /= 0) found = .true.

    enddo

  end function getFirstLineWithKeyword

!----------------------------------------------------------------------------
  subroutine ignoreNextNLinesFromFile(fUnit, n)
    !! Ignore next n lines from a file

    implicit none

    ! Input variables:
    integer, intent(in) :: fUnit
      !! Unit for file to read from
    integer, intent(in) :: n
      !! Number of lines to ignore

    ! Local variables:
    integer :: iLine
      !! Loop index
    

    do iLine = 1, n

      read(fUnit,*)

    enddo

    return

  end subroutine ignoreNextNLinesFromFile
  
!----------------------------------------------------------------------------
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
