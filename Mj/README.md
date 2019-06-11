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
* Change into the `Mj/src/` directory
* Execute `cat ../DOC/input.in | ./Mj.x` to send the contents of the `input.in` file into the program as input

## Code Walkthrough

`Mj_main.f90`
* Pull in `MjModule` that includes the needed (global) variables
* Start a timer
* Call `readInputs()`

`Mj_main.f90-->Mj_Module_v1.f90:readInputs()`
* Check if output file already exists
* Delete the output file if any exists
* Open a new output file
* Call `initialize()`

`Mj_main.f90-->Mj_Module_v1.f90:readInputs()-->initialize()`
* Set default values for all variables to test if they were read from the files
* Go back to `readInputs()`

`Mj_main.f90-->Mj_Module_v1.f90:readInputs()`
* Read input from command line (or from file if you use command given in the [How to Run](#how-to-run) section)

_Note: The input must be in a specific form as the variables are grouped into a [namelist](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc6/index.html)_
* Call `checkAndUpdateInput()`

`Mj_main.f90-->Mj_Module_v1.f90:readInputs()-->checkAndUpdateInput()`

