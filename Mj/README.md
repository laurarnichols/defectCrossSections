# Mj Module
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

## Variable Glossary

<table>
	<tr>
		<th>Name</th>
		<th>Type</th>
		<th>Defined In</th>
		<th>Allocated In</th>
		<th>Deallocated In</th>
		<th>Used In</th>
		<th>Meaning</th>
	</tr>
	<tr>
		<td>`abCM`</td>
		<td>`real, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`abortExecution`</td>
		<td>`logical`</td>
		<td>`checkAndUpdateInput()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>`checkAndUpdateInput()`</td>
	</tr>
	<tr>
		<td>`atomD`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td>`readPhonons()`</td>
		<td></td>
		<td>`readPhonons()`</td>
		<td></td>
	</tr>
	<tr>
		<td>`atomM`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td>`readPhonons()`</td>
		<td></td>
		<td>`readPhonons()`</td>
		<td></td>
	</tr>
	<tr>
		<td>`atomPosition`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`besOrderNofModeM`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`coth`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`dp`</td>
		<td>`integer, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Used to set the real numbers to double precision</td>
	</tr>
	<tr>
		<td>`dummyC`</td>
		<td>`character`</td>
		<td>`readPhonons()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>Dummy variable for trash from input file</td>
	</tr>
	<tr>
		<td>`dummyD`</td>
		<td>`real`</td>
		<td>`readPhonons()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>Dummy variable for trash from input file</td>
	</tr>
	<tr>
		<td>`elements`</td>
		<td>`character`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td>Array to hold element names</td>
	</tr>
	<tr>
		<td>`equilibriumAtomicPositions`</td>
		<td>`character`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`
		</td>
		<td>The name of the file that has the equilibrium atomic positions</td>
	</tr>
	<tr>
		<td>`eVToHartree`</td>
		<td>`real, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Conversion factor from eV to Hartree</td>
	</tr>
	<tr>
		<td>`file_exists`</td>
		<td>`logical`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readInputs()`</td>
		<td>Indicate if a file exists</td>
	</tr>
	<tr>
		<td>`freqInTHz`</td>
		<td>`real`</td>
		<td>`readPhonons()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>The phonon frequency in THz</td>
	</tr>
	<tr>
		<td>`genCoord`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`HartreeToEv`</td>
		<td>`real, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Conversion factor from Hartree to eV</td>
	</tr>
	<tr>
		<td>`iAtom`</td>
		<td>`integer`</td>
		<td>`readPhonons()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td>`iMode`</td>
		<td>`integer`</td>
		<td>`readPhonons()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td>`int32`</td>
		<td>`integer, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>To set the kind of integers</td>
	</tr>
	<tr>
		<td>`int64`</td>
		<td>`integer, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>To set the kind of integers</td>
	</tr>
	<tr>
		<td>`ios`</td>
		<td>`integer`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`iostd`</td>
		<td>`integer, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`readInputs()`<br/>
			`checkAndUpdateInput()`<br/>
			`readPhonons()`
		</td>
		<td>The unit of the output file</td>
	</tr>
	<tr>
		<td>`iQ`</td>
		<td>`integer`</td>
		<td>`readPhonons()`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td>`kT`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`checkAndUpdateInput()`</td>
		<td>Boltzmann constant time T then converted to Hartree</td>
	</tr>
	<tr>
		<td>`maxDisplacement`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`
		</td>
		<td>Max atomic displacement in each direction</td>
	</tr>
	<tr>
		<td>`MjInput`</td>
		<td>`namelist`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readInputs()`</td>
		<td>Grouping of several variables for input</td>
	</tr>
	<tr>
		<td>`modeF`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`
		</td>
		<td></td>
	</tr>
	<tr>
		<td>`modeI`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`
		</td>
		<td></td>
	</tr>
	<tr>
		<td>`nAtoms`</td>
		<td>`integer`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>The number of atoms</td>
	</tr>
	<tr>
		<td>`newAtomicPosition`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`newAtomicPositions`</td>
		<td>`character`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`nModes`</td>
		<td>`integer`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>The number of modes</td>
	</tr>
	<tr>
		<td>`nOfqPoints`</td>
		<td>`integer`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>The number of q points</td>
	</tr>
	<tr>
		<td>`output`</td>
		<td>`character, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readInputs()`</td>
		<td>The name of the output file</td>
	</tr>
	<tr>
		<td>`phonD`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td>`readPhonons()`</td>
		<td></td>
		<td>`readPhonons()`</td>
		<td>Phonon displacements</td>
	</tr>
	<tr>
		<td>`phonF`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td>`readPhonons()`</td>
		<td></td>
		<td>`readPhonons()`</td>
		<td>Phonon frequencies</td>
	</tr>
	<tr>
		<td>`phononsInput`</td>
		<td>`character`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`<br/>
			`readPhonons()`
		</td>
		<td>The name of the phonons input file</td>
	</tr>
	<tr>
		<td>`phonQ`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td>`readPhonons()`</td>
		<td></td>
		<td></td>
		<td>`readPhonons()`</td>
	</tr>
	<tr>
		<td>`pi`</td>
		<td>`real, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>The value of pi</td>
	</tr>
	<tr>
		<td>`QEInput`</td>
		<td>`character`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`
		</td>
		<td>The name of the QE input file</td>
	</tr>
	<tr>
		<td>`qPoint`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`readQEInput`</td>
		<td>`logical`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`checkAndUpdateInput()`</td>
		<td>Whether or not QEInput was read</td>
	</tr>
	<tr>
		<td>`s2L`</td>
		<td>`integer, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`Sj`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`t1`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Start time within subprocess</td>
	</tr>
	<tr>
		<td>`t2`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time within subprocess</td>
	</tr>
	<tr>
		<td>`temperature`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			`initialize()`<br/>
			`checkAndUpdateInput()`
		</td>
		<td>Temperature of the system</td>
	</tr>
	<tr>
		<td>`tf`</td>
		<td>`real`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time</td>
	</tr>
	<tr>
		<td>`THzToHartree`</td>
		<td>`real, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>`readPhonons()`</td>
		<td>Conversion factor from THz to Hartree</td>
	</tr>
	<tr>
		<td>`twopi`</td>
		<td>`real, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>2 times pi</td>
	</tr>
	<tr>
		<td>`un`</td>
		<td>`integer, parameter`</td>
		<td>`MjModule`</td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`wby2kT`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td>`x`</td>
		<td>`real, allocatable`</td>
		<td>`MjModule`</td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
</table>
