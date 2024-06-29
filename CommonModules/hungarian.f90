module hungarianOptimizationMod
  ! Linear assignment problem (LAP) solver
  ! -----------------------------------------------
  ! The first method (method = 1) is a brute force
  ! solution that computes all the permutations (n!)
  ! using Heap's algorithm - time complexity O(n!)
  ! The second method (method = 2) uses the Munkres
  ! algorithm (also called the Hungarian algo) to
  ! solve the problem - time complexity O(n^3)
  ! -----------------------------------------------
  ! F. Brieuc - March 2017
  !
  ! Implementation of the Munkres algorithm (also referred to as the Hungarian
  ! algorithm). J. Munkres, Journal of the SIAM 5, 1 (1957)

  use errorsAndMPI

  implicit none

  contains

  subroutine hungarianOptimalMatching(nElements, cost_in, minimizeCost, colSolution)

    implicit none
      
    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix

    real(kind=dp), intent(in) :: cost_in(nElements,nElements)
      !! Real cost as input

    logical, intent(in) :: minimizeCost
      !! If cost should be minimized

    ! Output variables:
    integer, allocatable, intent(out) :: colSolution(:)
      !! Optimized row indices for each corresponding column index

    ! Local variables:
    integer :: iRow, iCol
      !! Loop indices
    integer, allocatable :: M(:,:)
      !! Mask matrix to mark special zeros
    integer :: pathRow0, pathCol0  
      !! Starting point for path finding part
    integer, allocatable :: rowCover(:), colCover(:) 
      !! Mark rows and columns covered
    integer :: step
      !! What step we are on

    real(kind=dp), allocatable :: cost(:,:)
      !! Real cost used in algorithm

    logical :: done
      !! If done with algorithm


    if(nElements < 0) call exitError('hungarianOptimalMatching', 'must have at least one element in cost matrix', 1)


    allocate(cost(nElements,nElements))

    if(.not. minimizeCost) then ! maximize sum
      cost = -cost_in
      cost = cost - minval(cost)
    else
      cost = cost_in
    endif


    allocate(colSolution(nElements))
    allocate(M(nElements,nElements),rowCover(nElements),colCover(nElements))

    M = 0
    rowCover = 0
    colCover = 0
    pathRow0 = 0
    pathCol0 = 0
    done = .false.
    step = 1


    do while(.not. done)
      select case(step)
      case(1)
        call step1(nElements, step, cost)
      case(2)
        call step2(nElements, cost, M, rowCover, colCover, step)
      case(3)
        call step3(nElements, M, colCover, step)
      case(4)
        call step4(nElements, cost, M, pathRow0, pathCol0, rowCover, colCover, step)
      case(5)
        call step5(nElements, pathRow0, pathCol0, M, rowCover, colCover, step)
      case(6)
        call step6(nElements, rowCover, colCover, step, cost)
      case default ! done
         do iCol = 1, nElements
            do iRow = 1, nElements
               if(M(iRow,iCol) == 1) colSolution(iCol) = iRow
            enddo
         enddo
         done = .true.
      end select
    enddo

    deallocate(M)
    deallocate(rowCover)
    deallocate(colCover)
    deallocate(cost)

    return

  end subroutine hungarianOptimalMatching

!----------------------------------------------------------------------------
  subroutine step1(nElements, step, cost)
    ! Find the smallest value in each column and subtract it from
    ! all elements in the column. Ensures each column has at least
    ! one zero. Go to step 2.

    implicit none

    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix

    ! Output variables:
    integer, intent(out) :: step
      !! What step we are on

    real(kind=dp), intent(inout) :: cost(nElements,nElements)
      !! Real cost to be maximized

    ! Local variables:
    integer :: iCol
      !! Loop indices


    do iCol = 1, nElements
      cost(:,iCol) = cost(:,iCol) - minval(cost(:,iCol))
    enddo

    step = 2

    return

  end subroutine step1

!----------------------------------------------------------------------------
  subroutine step2(nElements, cost, M, rowCover, colCover, step)
    ! Search for zeros.
    ! Find a zero (Z) in the matrix. If no zeros has been previously starred in
    ! its row and column then star Z. Go to step 3.

    implicit none

    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix

    real(kind=dp), intent(in) :: cost(nElements,nElements)
      !! Real cost to be maximized

    ! Output variables:
    integer, intent(inout) :: M(nElements,nElements)
      !! Mask matrix to mark special zeros
    integer, intent(inout) :: rowCover(nElements), colCover(nElements) 
      !! Mark rows and columns covered
    integer, intent(out) :: step
      !! What step we are on

    ! Local variables:
    integer :: i, j
      !! Loop indices


    do i = 1, nElements
      do j = 1, nElements
        if(cost(j,i) == 0.0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
           M(j,i) = 1
           rowCover(i) = 1
           colCover(j) = 1
        endif
      enddo
    enddo

    ! Uncover the rows and columns
    rowCover = 0
    colCover = 0

    step = 3

    return

  end subroutine step2

!----------------------------------------------------------------------------
  subroutine step3(nElements, M, colCover, step)
    ! cover each column containing a starred zero. If n column are covered
    ! the starred zero describe an optimal assignment and we are done otherwise
    ! go to step 4.
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix
    integer, intent(in) :: M(nElements,nElements)
      !! Mask matrix to mark special zeros

    ! Output variables:
    integer, intent(inout) :: colCover(nElements) 
      !! Mark columns covered
    integer, intent(out) :: step
      !! What step we are on

    ! Local variables:
    integer :: colCount, i, j


    colCount = 0
    do i = 1, nElements
      do j = 1, nElements
        ! if starred and column is uncovered
        if(M(j,i) == 1 .and. colCover(j) == 0) then
          colCover(j) = 1
          colCount = colCount + 1
        endif
      enddo
    enddo

    if(colCount == nElements) then
      step = 0
    else
      step = 4
    endif

    return

  end subroutine step3

!----------------------------------------------------------------------------
  subroutine step4(nElements, cost, M, pathRow0, pathCol0, rowCover, colCover, step)
    ! Find a uncovered zero and prime it. If there is no starred zero in the row
    ! go to step 5. Otherwise, cover the row and uncover the column containing
    ! the starred zero. Continue until no uncovered zeros is left. Go to step 6.

    implicit none

    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix

    real(kind=dp), intent(in) :: cost(nElements,nElements)
      !! Real cost to be maximized

    ! Output variables:
    integer, intent(inout) :: M(nElements,nElements)
      !! Mask matrix to mark special zeros
    integer, intent(out) :: pathRow0, pathCol0  
      !! Starting point for path finding part
    integer, intent(inout) :: rowCover(nElements), colCover(nElements) 
      !! Mark rows and columns covered
    integer, intent(out) :: step
      !! What step we are on

    ! Local variables:
    logical :: done, starInRow
    integer :: i, j, row, col


    done = .false.

    do while (.not. done)
      ! find an uncovered zero
      row = 0; col = 0
      starInRow = .false.
      loop1: do i = 1, nElements
        loop2: do j = 1, nElements
          if(cost(j,i) < 1e-8 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
            row = i
            col = j
            exit loop1
          endif
        enddo loop2
      enddo loop1

      if(row == 0) then !no zero uncoverred left
        done = .true.
        step = 6
      else
        M(col,row) = 2 !primed zero
        ! search if there is a starred zero in the same row
        do j = 1, nElements
          if(M(j,row) == 1) then
            starInRow = .true.
            col = j
          endif
        enddo
        
        if(starInRow) then ! if there is a starred zero in line
          rowCover(row) = 1
          colCover(col) = 0
        else ! if no starred zero in line
          done = .true.
          step = 5
          pathRow0 = row
          pathCol0 = col
        endif
      endif
    enddo

    return

  end subroutine step4

!----------------------------------------------------------------------------
  subroutine step5(nElements, pathRow0, pathCol0, M, rowCover, colCover, step)
    ! Augmenting path algorithm: construct a serie of alternating primed and
    ! starred zeros as follows. Let Z0 be the uncoverd primed zero found in
    ! step 4. Let Z1 be the starred zero in the column of Z0 (if any).
    ! Let Z2 be the primed zero in the row of Z1 (there will always be one).
    ! Continue until the series terminates at a primed zero that has no starred
    ! zero in its column. Then unstar each starred zeros of the series, star
    ! each primed zeros of the series, erase all primes and uncover every line
    ! and columns. Return to step 3.

    implicit none

    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix
    integer, intent(in) :: pathRow0, pathCol0  
      !! Starting point for path finding part

    ! Output variables:
    integer, intent(inout) :: M(nElements,nElements)
      !! Mask matrix to mark special zeros
    integer, intent(inout) :: rowCover(nElements), colCover(nElements) 
      !! Mark rows and columns covered
    integer, intent(out) :: step
      !! What step we are on

    ! Local variables:
    logical :: done
    integer :: i, j
    integer :: row, col
    integer :: pathCount
    integer, allocatable :: path(:,:)


    allocate(path(2*nElements+1,2))

    pathCount = 1

    path(pathCount,1) = pathRow0
    path(pathCount,2) = pathCol0

    done = .false.

    do while(.not. done)
      ! search for a starred zero in column
      row = 0
      col = path(pathCount,2)
      do i = 1, nElements
        if(M(col,i) == 1) row = i
      enddo
      if (row /= 0) then ! update path
        pathCount = pathCount + 1
        path(pathCount,1) = row
        path(pathCount,2) = path(pathCount-1,2)
      else
        done = .true.
      endif
      if(.not. done) then
        ! search for a prime zero in row
        do j = 1, nElements
          if (M(j,row) == 2) col = j
        enddo
        ! update path
        pathCount = pathCount + 1
        path(pathCount,1) = path(pathCount-1,1)
        path(pathCount,2) = col
      endif
    enddo

    ! augment path
    do i = 1, pathCount
      if(M(path(i,2),path(i,1)) == 1) then
        M(path(i,2),path(i,1)) = 0
      else
        M(path(i,2),path(i,1)) = 1
      endif
    enddo

    ! clear covers and erase primes
    do i = 1, nElements
      rowCover(i) = 0
      colCover(i) = 0
      do j = 1, nElements
        if (M(j,i) == 2) M(j,i) = 0
      enddo
    enddo

    step = 3

    deallocate(path)

    return

  end subroutine step5

!----------------------------------------------------------------------------
  subroutine step6(nElements, rowCover, colCover, step, cost)
    ! Search for the smallest uncovered value and add it to covered rows
    ! and substract it from uncovered columns. Return to step 4.
    
    implicit none

    ! Input variables:
    integer, intent(in) :: nElements
      !! Square size of the cost matrix

    ! Output variables:
    integer, intent(inout) :: rowCover(nElements), colCover(nElements) 
      !! Mark rows and columns covered
    integer, intent(out) :: step
      !! What step we are on

    real(kind=dp), intent(inout) :: cost(nElements,nElements)
      !! Real cost to be maximized

    ! Local variables
    integer :: i, j

    real(kind=dp) :: tmp, minVal

    minVal = huge(tmp)

    do i = 1, nElements
      do j = 1, nElements
        if(rowCover(i) == 0 .and. colCover(j) == 0 .and. cost(j,i) < minVal) then
          minVal = cost(j,i)
        endif
      enddo
    enddo
    do i = 1, nElements
      do j = 1, nElements
        if(rowCover(i) == 1) cost(j,i) = cost(j,i) + minVal
        if(colCover(j) == 0) cost(j,i) = cost(j,i) - minVal
      enddo
    enddo

    step = 4

    return

  end subroutine step6

end module hungarianOptimizationMod
