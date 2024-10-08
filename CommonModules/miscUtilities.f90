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

!-----------------------------------------------------------------------------------------------
  subroutine getUniqueInts(arrSize, arr, nUnique, uniqueVals)

    implicit none

    ! Input variables:
    integer, intent(in) :: arrSize
      !! Size of array to get unique values from
    integer, intent(in) :: arr(arrSize)
      !! Array to get unique values from

    ! Output variables:
    integer, intent(out) :: nUnique
      !! Count of unique values in the array
    integer, allocatable, intent(out) :: uniqueVals(:)
      !! Unique values the size of nUnique

    ! Local variables:
    integer :: maxValOverall
      !! Maximum value in the array
    integer :: minValRemain
      !! Min value remaining after removing the values
      !! already found to be unique from consideration
    integer :: uniqueValsOrigSize(arrSize)
      !! Unique values stored in an array the size of
      !! the input array


    if(arrSize > 0) then
      nUnique = 1
      ! If there is at least one value in the array,
      ! there is always at least one unique value

      minValRemain = minval(arr)
      maxValOverall = maxval(arr)

      uniqueValsOrigSize(1) = minValRemain
      
      do while (minValRemain < maxValOverall)
          minValRemain = minval(arr, mask=arr>minValRemain)
            ! Get the minimum value remaining among numbers greater
            ! than the previous minimum value remaining. This will
            ! automatically exclude numbers equal to the previous
            ! minimum value.

          nUnique = nUnique + 1
          uniqueValsOrigSize(nUnique) = minValRemain
            ! Store the new unique value in the next slot of the 
            ! unique-values array
      enddo

      allocate(uniqueVals(nUnique), source=uniqueValsOrigSize(1:nUnique))

    endif

    return

  end subroutine

!----------------------------------------------------------------------------
  subroutine hpsort_eps(n, ra, ind, eps)
    !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
    !! considering two elements equal if their difference is less than `eps`
    !!
    !! n is input, ra is replaced on output by its sorted rearrangement.
    !! Create an index table (ind) by making an exchange in the index array
    !! whenever an exchange is made on the sorted data array (ra).
    !! In case of equal values in the data array (ra) the values in the
    !! index array (ind) are used to order the entries.
    !!
    !! if on input ind(1)  = 0 then indices are initialized in the routine,
    !! if on input ind(1) != 0 then indices are assumed to have been
    !!                initialized before entering the routine and these
    !!                indices are carried around during the sorting process
    !!
    !! To sort other arrays based on the index order in ind, use
    !!      sortedArray(i) = originalArray(ind(i))
    !! If you need to recreate the original array like
    !!      originalArray(i) = sortedArray(indRec(i))
    !! indRec can be generated by sorting a sequential array (routine will
    !! initialize to sequential if not already) based on the order in ind.
    !!      realInd = real(ind)
    !!      hpsort_eps(n, realInd, indRec, eps)
    !!
    !! From QE code, adapted from Numerical Recipes pg. 329 (new edition)
    !!

    implicit none

    ! Input/Output variables:
    real(kind=dp), intent(in) :: eps
    integer, intent(in) :: n

    integer, intent(inout) :: ind(:)
    real(kind=dp), intent(inout) :: ra (:)


    ! Local variables
    integer :: i, ir, j, l, iind
    real(kind=dp) :: rra

    !> Initialize index array, if not already initialized
    if (ind (1) .eq.0) then
      do i = 1, n
        ind(i) = i
      enddo
    endif

    ! nothing to order
    if (n.lt.2) return
    ! initialize indices for hiring and retirement-promotion phase
    l = n / 2 + 1

    ir = n

    sorting: do

      ! still in hiring phase
      if ( l .gt. 1 ) then
        l    = l - 1
        rra  = ra (l)
        iind = ind(l)
        ! in retirement-promotion phase.
      else
        ! clear a space at the end of the array
        rra  = ra (ir)
        !
        iind = ind(ir)
        ! retire the top of the heap into it
        ra (ir) = ra (1)
        !
        ind(ir) = ind(1)
        ! decrease the size of the corporation
        ir = ir - 1
        ! done with the last promotion
        if ( ir .eq. 1 ) then
          ! the least competent worker at all !
          ra (1)  = rra
          !
          ind(1) = iind
          exit sorting
        endif
      endif
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      do while ( j .le. ir )
        if ( j .lt. ir ) then
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then
            if (ra(j).lt.ra(j+1)) j = j + 1
          else
            ! this means ra(j) == ra(j+1) within tolerance
            if (ind(j) .lt.ind(j + 1) ) j = j + 1
          endif
        endif
        ! demote rra
        if ( abs(rra - ra(j)).ge.eps ) then
          if (rra.lt.ra(j)) then
            ra (i) = ra (j)
            ind(i) = ind(j)
            i = j
            j = j + j
          else
            ! set j to terminate do-while loop
            j = ir + 1
          end if
        else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind(j) ) then
            ra (i) = ra (j)
            ind(i) = ind(j)
            i = j
            j = j + j
          else
            ! set j to terminate do-while loop
            j = ir + 1
          endif
        end if
      enddo
      ra (i) = rra
      ind(i) = iind
    end do sorting
    
    return 
  end subroutine hpsort_eps
  
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
