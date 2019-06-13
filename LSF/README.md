# Mj Module
## Directory Structure and Files
```
LSF
|	README.md
|_______DOC
|	|	input.in
|	|	phonons_SiVH3newDisp.dat
|_______src
	|_______linearTerms
	|	|	LSF_linear_Main.f90
	|	|	LSF_linear_Module_v1.f90
	|	|	Makefile
	|_______zerothOrder
		|	LSF_zeroth_Main.f90
		|	LSF_zeroth_Module_v35.f90
		|	Makefile
```

## How to Run
* Ensure that the executables are up to date by going to the main folder and running `make LSF`
* Change into the `LSF/DOC` directory
* Execute `cat input.in | ../src/<order folder>/LSF<0 or 1>.x` to send the contents of the `input.in` file into the program as input
* The output messages will be held in the `status` file now in the `DOC` directory

## Inputs
* `input.in`
	* `phononsInput` (string) -- directs the program to the phonons input file
	* `phononsInputFormat` (string) -- what format the phonons input file has (VASP or QE)
	* `continueLSFfromFile` (string) -- directs the program to continue from a given output file if it exists
	* `maximumNumberOfPhonons` -- the maximum number of phonons to use
	* `temperature` -- the temperature of the system
	* `maxEnergy` -- ??
	* `nMC` -- ??
	
_Note: Do not alter the `&lsfInput` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

* `phononsInput` file (e.g. phonons_SiVH3newDisp.dat)
	* `nOfqPoints` (integer) -- the number of q points
	* `nAtoms` (integer) -- the number of atoms in the system
	* `atomD` (real `3xnAtoms` array) -- ??
	* `atomM` (real vector of length `nAtoms`) -- ??
	* `phonQ` (real `3xnOfqPoints` array) -- ??
	* `freqInTHz` (real) -- the frequency of the mode in THz; this value gets automatically converted to Hartree and put in the `phonF` array
	* `phonD` (real `3xnAtomsxnModesxnOfqPoints` array) -- phonon displacements
