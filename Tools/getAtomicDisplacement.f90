program getAtomicDisplacement
  !! Read in two sets of POSCARS or CONTCARS and output 
  !! the magnitude of the displacement vector for each 
  !! atom and the average displacement vector magnitude
  !!
  !! Will check to make sure that atom counts match
  !!
  !! <h2>Walkthrough</h2>

  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)
    !! Used to set real variables to double precision

  real(kind=dp) :: avgDispMag
    !! Average value of the displacement vector for
    !! all atoms
  real(kind=dp), allocatable :: coord1(:,:)
    !! First set of coordinates
  real(kind=dp), allocatable :: coord2(:,:)
    !! Second set of coordinates
  real(kind=dp) :: coordDiff(3)
    !! Difference between coordinates for a given atom
  real(kind=dp), allocatable :: dispMag(:)
    !! Magnitude of displacement vector for each atom
  real(kind=dp) :: latticeVectors(3,3)
    !! Lattice vectors
  real(kind=dp) :: scalingFactor
    !! Lattice scaling factor

  integer :: nargs
    !! Number of command-line arguments
  integer :: nAtoms
    !! Number of atoms in the system

  character(len=256) :: file1
    !! Coordinate file 1
  character(len=256) :: file2
    !! Coordinate file 2


  write(*,*) "Checking that we have 2 command-line arguments"

  nargs = command_argument_count()

  if(nargs /= 2) call exitError("getAtomicDisplacement", "invalid number of arguments, expected 2", 1)
    !! Make sure that there are two arguments

  write(*,*) "2 command-line arguments found."

  call get_command_argument(1, file1)
  call get_command_argument(2, file2)
    !! Get POSCAR/CONTCAR file names from command line

  write(*,*) "Checking compatibility of ", trim(file1), " and ", trim(file2)

  return 

  call checkCompatibility()
    !! Check that the number of atoms of each type is the same
    !! for both files so that they can be compared
    !! @todo Write `checkCompatibility` subroutine @endtodo

  write(*,*) "Files compatible. Reading coordinates."

  call readCoordinateFile(file1)
  call readCoordinateFile(file2)
    !! Read in POSCAR/CONTCAR files
    !! @todo Write `readCoordinateFile` subroutine @endtodo

  write(*,*) "Getting displacement magnitude per atom and overall average."

  call getDispMag()
    !! Get the displacement vector magnitude for each
    !! atom and the average magnitude
    !! @todo Write `getDispMag` subroutine @endtodo

  write(*,*) "Writing results."

  call writeResults()
    !! Write out results
    !! @todo Write `writeResults` subroutine @endtodo

  write(*,*) "Done."

end program getAtomicDisplacement

!==============================================================================

subroutine exitError(calledFrom, message, ierr)
  !! Output error message and abort if ierr > 0
  !!
  !! Can ensure that error will cause abort by
  !! passing abs(ierror)
  !!
  !! <h2>Walkthrough</h2>
  !!
    
  implicit none

  integer, intent(in) :: ierr
    !! Error code

  character(len=*), intent(in) :: calledFrom
    !! Place where this subroutine was called from
  character(len=*), intent(in) :: message
    !! Error message

  character(len=6) :: cerr
    !! String version of error


  if ( ierr <= 0 ) return
    !! * Do nothing if the error is less than or equal to zero

  write( cerr, fmt = '(I6)' ) ierr
    !! * Write ierr to a string
  write(unit=*, fmt = '(/,1X,78("%"))' )
    !! * Output a dividing line
  write(unit=*, fmt = '(5X,"Error in ",A," (",A,"):")' ) trim(calledFrom), trim(adjustl(cerr))
    !! * Output where the error occurred and the error
  write(unit=*, fmt = '(5X,A)' ) TRIM(message)
    !! * Output the error message
  write(unit=*, fmt = '(1X,78("%"),/)' )
    !! * Output a dividing line

  write( *, '("     stopping ...")' )
  
  call flush( 16 )

  stop 2

  return

end subroutine exitError
