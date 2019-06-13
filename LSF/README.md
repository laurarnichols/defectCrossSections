# Mj Module
## Directory Structure and Files
```markdown
Mj
|	README.md
|_______DOC
|	|	SiVH3.in
|	|	equilibriumAtomicPositions66
|	|	input.in
|	|	phonons_SiVH3newDisp.dat
|_______src
	|	Makefile
	|	Mj.x
	|	Mj_Main.f90
	|	Mj_Main.o
	|	Mj_Module_v1.f90
	|	Mj_Module_v1.o
	|	mjmodule.mod
```

## How to Run
* Change into the `DOC` directory
* Execute `cat input.in | ../src/Mj.x` to send the contents of the `input.in` file into the program as input
* The output messages will be held in the `status` file now in the `DOC` directory
* The program will generate a new QE input file for given mode(s), each in their own directory

## Code Walkthrough

`Mj_main.f90:MjME`
* Pull in `MjModule` that includes the needed (global) variables
* Start a timer
* Call `readInputs()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()`
* Check if output file already exists
* Delete the output file if any exists
* Open a new output file (the name is given in the variable declarations: `output=[file name]`)
* Call `initialize()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()-->initialize()`
* Set default values for all variables to test if they were read from the files
* Return to `readInputs()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()`
* Read input from command line (or from file if you use command given in the [How to Run](#how-to-run) section)

_Note: The input must be in a specific form as the variables are grouped into a [namelist](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc6/index.html)_
* Call `checkAndUpdateInput()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()-->checkAndUpdateInput()`
* Check each variable to see if it still has the default value
* If it does:
	* Output an error message to the output file
	* Tell the program to abort
* Else:
	* Output the new variable value to the output file
* Return to `readInputs()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()`
* Call `readPhonons()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()-->readPhonons()`
* Set up some local variables
* Open the file given as `phononsInput`
* Read the number of q points and the number of atoms from the first line of the file
* Output the read values to the output file
* Calculate the number of modes as `nModes = 3*nAtoms - 3`
* Read in a blank line from file
* Allocate space for the `atomD` and `atomM` arrays
* __What are atomD and atomM arrays representing?__
* Initialize all values in `atomD` and `atomM` arrays to 0.0
* Loop through the next `nAtoms` lines and read in the data from the `atomD` and `atomM` arrays
* Read in a blank line from file
* Allocate space for the `phonQ`, `phonF`, and `phonD` arrays 
* __What are phonQ, phonF, and phonD arrays representing?__
* Initialize all values in `phonQ`, `phonF`, and `phonD` arrays to 0.0
* For each q point:
	* Read in one column of the `phonQ` array
	* For each mode:
		* Read in a blank line from file
		* Read in the frequency in THz
		* Convert the frequency to Hartree units
		* Read in a trash line
		* For each atom:
			* Read in the `phonD` array
* Close the input file
* Use `flush` to make the output file available for other processes
* Return to `readInputs()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()`
* Call `readAtomicPositions()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()-->readAtomicPositions()`
* Open the file given as `equilibriumAtomicPositions`
* Allocate space for the `elements` and `atomPosition` arrays
* Initialize all values in the `atomPosition` array to 0.0
* For each atom:
	* Read in the element, and the 3D position
* Close the input file
* Return to `readInputs()`

`Mj_main.f90:MjME-->Mj_Module_v1.f90:readInputs()`
* Return to `MjME`

`Mj_main.f90:MjME`
* Call `computeGeneralizedDisplacements()`

`Mj_main.f90:MjME-->computeGeneralizedDisplacements()`
_From here, the code gets more scientific. Because there are no comments and the variable names aren't clear,
I can't tell what the code is doing exactly.`
