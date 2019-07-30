# Mj
## Directory Structure and Files
```
Mj
|	README.md
|_______DOC
|	|	SiVH3.in
|	|	equilibriumAtomicPositions66
|	|	input.in
|	|	phonons_SiVH3newDisp.dat
|_______src
	|	Makefile
	|	Mj_Main.f90
	|	Mj_Module_v1.f90
```

## How to Run
* Ensure that the executables are up to date by going to the main folder and running `make LSF`
* Change into the `Mj/DOC` directory
* Execute `cat input.in | ../src/Mj.x` to send the contents of the `input.in` file into the program as input
* The output messages will be held in the `status` file now in the `DOC` directory
* The program will generate a new QE input file for given mode(s), each in their own directory

## Inputs
* `input.in`
	* `QEInput` (string) -- directs the program to the QE input file
	* `phononsInput` (string) -- directs the program to the phonons input file
	* `equilibriumAtomicPositions` (string) -- directs the program to a file that has the equilibrium position of each of the atoms
	* `temperature` -- the temperature of the system
	* `maxDisplacement` -- the maximum displacement value for atoms 
	* `modeI` -- initial mode (must be smaller than `modeF`)
	* `modeF` -- final mode
	* `qPoint` -- 
	
_Note: Do not alter the `&lsfInput` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

* `QEInput` file (e.g. SiVH3.in)

* [`phononsInput`](https://github.com/laurarnichols/carrierCrossSections/blob/master/Mj/DOC/phononsInput.md) file (e.g. phonons_SiVH3newDisp.dat)
	* `nOfqPoints` (integer) -- the number of q points
	* `nAtoms` (integer) -- the number of atoms in the system
	* `atomD` (real `3xnAtoms` array) -- 
	* `atomM` (real vector of length `nAtoms`) -- 
	* `phonQ` (real `3xnOfqPoints` array) -- 
	* `freqInTHz` (real) -- the frequency of the mode in THz; this value gets automatically converted to Hartree and put in the `phonF` array
	* `phonD` (real `3xnAtomsxnModesxnOfqPoints` array) -- phonon displacements

## Outputs
* Information for each mode ordered by magnitude of modified Bessel function argument `x`
   * Index of mode based on input
   * Mode frequency in eV
   * Generalized coordinates
   * Argument of modified Bessel function (`x`)
* New postions or new QE input (if give QE input)
