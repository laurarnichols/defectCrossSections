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
		<td><code>abCM</code></td>
		<td><code>real, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>abortExecution</code></td>
		<td><code>logical</code></td>
		<td><code>checkAndUpdateInput()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td><code>checkAndUpdateInput()</code></td>
	</tr>
	<tr>
		<td><code>atomD</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td><code>readPhonons()</code></td>
		<td></td>
		<td><code>readPhonons()</code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>atomM</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td><code>readPhonons()</code></td>
		<td></td>
		<td><code>readPhonons()</code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>atomPosition</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>besOrderNofModeM</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>coth</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>dp</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Used to set the real numbers to double precision</td>
	</tr>
	<tr>
		<td><code>dummyC</code></td>
		<td><code>character</code></td>
		<td><code>readPhonons()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>Dummy variable for trash from input file</td>
	</tr>
	<tr>
		<td><code>dummyD</code></td>
		<td><code>real</code></td>
		<td><code>readPhonons()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>Dummy variable for trash from input file</td>
	</tr>
	<tr>
		<td><code>elements</code></td>
		<td><code>character</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td>Array to hold element names</td>
	</tr>
	<tr>
		<td><code>equilibriumAtomicPositions</code></td>
		<td><code>character</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()<code>
		</td>
		<td>The name of the file that has the equilibrium atomic positions</td>
	</tr>
	<tr>
		<td><code>eVToHartree</code></td>
		<td><code>real, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Conversion factor from eV to Hartree</td>
	</tr>
	<tr>
		<td><code>file_exists</code></td>
		<td><code>logical</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInputs()</code></td>
		<td>Indicate if a file exists</td>
	</tr>
	<tr>
		<td><code>freqInTHz</code></td>
		<td><code>real</code></td>
		<td><code>readPhonons()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>The phonon frequency in THz</td>
	</tr>
	<tr>
		<td><code>genCoord</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>HartreeToEv</code></td>
		<td><code>real, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Conversion factor from Hartree to eV</td>
	</tr>
	<tr>
		<td><code>iAtom</code></td>
		<td><code>integer</code></td>
		<td><code>readPhonons()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td><code>iMode</code></td>
		<td><code>integer</code></td>
		<td><code>readPhonons()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td><code>int32</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>To set the kind of integers</td>
	</tr>
	<tr>
		<td><code>int64</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>To set the kind of integers</td>
	</tr>
	<tr>
		<td><code>ios</code></td>
		<td><code>integer</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>iostd</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>readInputs()</code><br/>
			<code>checkAndUpdateInput()</code><br/>
			<code>readPhonons()<code>
		</td>
		<td>The unit of the output file</td>
	</tr>
	<tr>
		<td><code>iQ</code></td>
		<td><code>integer</code></td>
		<td><code>readPhonons()</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td><code>kT</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>checkAndUpdateInput()</code></td>
		<td>Boltzmann constant time T then converted to Hartree</td>
	</tr>
	<tr>
		<td><code>maxDisplacement</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()<code>
		</td>
		<td>Max atomic displacement in each direction</td>
	</tr>
	<tr>
		<td><code>MjInput</code></td>
		<td><code>namelist</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInputs()</code></td>
		<td>Grouping of several variables for input</td>
	</tr>
	<tr>
		<td><code>modeF</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()<code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>modeI</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()<code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>nAtoms</code></td>
		<td><code>integer</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>The number of atoms</td>
	</tr>
	<tr>
		<td><code>newAtomicPosition</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>newAtomicPositions</code></td>
		<td><code>character</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nModes</code></td>
		<td><code>integer</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>The number of modes</td>
	</tr>
	<tr>
		<td><code>nOfqPoints</code></td>
		<td><code>integer</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>The number of q points</td>
	</tr>
	<tr>
		<td><code>output</code></td>
		<td><code>character, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInputs()</code></td>
		<td>The name of the output file</td>
	</tr>
	<tr>
		<td><code>phonD</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td><code>readPhonons()</code></td>
		<td></td>
		<td><code>readPhonons()</code></td>
		<td>Phonon displacements</td>
	</tr>
	<tr>
		<td><code>phonF</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td><code>readPhonons()</code></td>
		<td></td>
		<td><code>readPhonons()</code></td>
		<td>Phonon frequencies</td>
	</tr>
	<tr>
		<td><code>phononsInput</code></td>
		<td><code>character</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()</code><br/>
			<code>readPhonons()<code>
		</td>
		<td>The name of the phonons input file</td>
	</tr>
	<tr>
		<td><code>phonQ</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td><code>readPhonons()</code></td>
		<td></td>
		<td></td>
		<td><code>readPhonons()</code></td>
	</tr>
	<tr>
		<td><code>pi</code></td>
		<td><code>real, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>The value of pi</td>
	</tr>
	<tr>
		<td><code>QEInput</code></td>
		<td><code>character</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()<code>
		</td>
		<td>The name of the QE input file</td>
	</tr>
	<tr>
		<td><code>qPoint</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>readQEInput</code></td>
		<td><code>logical</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>checkAndUpdateInput()</code></td>
		<td>Whether or not QEInput was read</td>
	</tr>
	<tr>
		<td><code>s2L</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>Sj</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>t1</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Start time within subprocess</td>
	</tr>
	<tr>
		<td><code>t2</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time within subprocess</td>
	</tr>
	<tr>
		<td><code>temperature</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize()</code><br/>
			<code>checkAndUpdateInput()<code>
		</td>
		<td>Temperature of the system</td>
	</tr>
	<tr>
		<td><code>tf</code></td>
		<td><code>real</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time</td>
	</tr>
	<tr>
		<td><code>THzToHartree</code></td>
		<td><code>real, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readPhonons()</code></td>
		<td>Conversion factor from THz to Hartree</td>
	</tr>
	<tr>
		<td><code>twopi</code></td>
		<td><code>real, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>2 times pi</td>
	</tr>
	<tr>
		<td><code>un</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>MjModule</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>wby2kT</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><code>x</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>MjModule</code></td>
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
