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
		<td><pre>abCM</pre></td>
		<td><pre>real, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>abortExecution</pre></td>
		<td><pre>logical</pre></td>
		<td><pre>checkAndUpdateInput()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td><pre>checkAndUpdateInput()</pre></td>
	</tr>
	<tr>
		<td><pre>atomD</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>atomM</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>atomPosition</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>besOrderNofModeM</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>coth</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>dp</pre></td>
		<td><pre>integer, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Used to set the real numbers to double precision</td>
	</tr>
	<tr>
		<td><pre>dummyC</pre></td>
		<td><pre>character</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>Dummy variable for trash from input file</td>
	</tr>
	<tr>
		<td><pre>dummyD</pre></td>
		<td><pre>real</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>Dummy variable for trash from input file</td>
	</tr>
	<tr>
		<td><pre>elements</pre></td>
		<td><pre>character</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td>Array to hold element names</td>
	</tr>
	<tr>
		<td><pre>equilibriumAtomicPositions</pre></td>
		<td><pre>character</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()<pre>
		</td>
		<td>The name of the file that has the equilibrium atomic positions</td>
	</tr>
	<tr>
		<td><pre>eVToHartree</pre></td>
		<td><pre>real, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Conversion factor from eV to Hartree</td>
	</tr>
	<tr>
		<td><pre>file_exists</pre></td>
		<td><pre>logical</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readInputs()</pre></td>
		<td>Indicate if a file exists</td>
	</tr>
	<tr>
		<td><pre>freqInTHz</pre></td>
		<td><pre>real</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>The phonon frequency in THz</td>
	</tr>
	<tr>
		<td><pre>genCoord</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>HartreeToEv</pre></td>
		<td><pre>real, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Conversion factor from Hartree to eV</td>
	</tr>
	<tr>
		<td><pre>iAtom</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td><pre>iMode</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td><pre>int32</pre></td>
		<td><pre>integer, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>To set the kind of integers</td>
	</tr>
	<tr>
		<td><pre>int64</pre></td>
		<td><pre>integer, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>To set the kind of integers</td>
	</tr>
	<tr>
		<td><pre>ios</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>iostd</pre></td>
		<td><pre>integer, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>readInputs()</pre><br/>
			<pre>checkAndUpdateInput()</pre><br/>
			<pre>readPhonons()<pre>
		</td>
		<td>The unit of the output file</td>
	</tr>
	<tr>
		<td><pre>iQ</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>Integer used for loop</td>
	</tr>
	<tr>
		<td><pre>kT</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>checkAndUpdateInput()</pre></td>
		<td>Boltzmann constant time T then converted to Hartree</td>
	</tr>
	<tr>
		<td><pre>maxDisplacement</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()<pre>
		</td>
		<td>Max atomic displacement in each direction</td>
	</tr>
	<tr>
		<td><pre>MjInput</pre></td>
		<td><pre>namelist</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readInputs()</pre></td>
		<td>Grouping of several variables for input</td>
	</tr>
	<tr>
		<td><pre>modeF</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()<pre>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><pre>modeI</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()<pre>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><pre>nAtoms</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>The number of atoms</td>
	</tr>
	<tr>
		<td><pre>newAtomicPosition</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>newAtomicPositions</pre></td>
		<td><pre>character</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>nModes</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>The number of modes</td>
	</tr>
	<tr>
		<td><pre>nOfqPoints</pre></td>
		<td><pre>integer</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>The number of q points</td>
	</tr>
	<tr>
		<td><pre>output</pre></td>
		<td><pre>character, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readInputs()</pre></td>
		<td>The name of the output file</td>
	</tr>
	<tr>
		<td><pre>phonD</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
		<td><pre>readPhonons()</pre></td>
		<td>Phonon displacements</td>
	</tr>
	<tr>
		<td><pre>phonF</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
		<td><pre>readPhonons()</pre></td>
		<td>Phonon frequencies</td>
	</tr>
	<tr>
		<td><pre>phononsInput</pre></td>
		<td><pre>character</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()</pre><br/>
			<pre>readPhonons()<pre>
		</td>
		<td>The name of the phonons input file</td>
	</tr>
	<tr>
		<td><pre>phonQ</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td><pre>readPhonons()</pre></td>
		<td></td>
		<td></td>
		<td><pre>readPhonons()</pre></td>
	</tr>
	<tr>
		<td><pre>pi</pre></td>
		<td><pre>real, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>The value of pi</td>
	</tr>
	<tr>
		<td><pre>QEInput</pre></td>
		<td><pre>character</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()<pre>
		</td>
		<td>The name of the QE input file</td>
	</tr>
	<tr>
		<td><pre>qPoint</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>readQEInput</pre></td>
		<td><pre>logical</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>checkAndUpdateInput()</pre></td>
		<td>Whether or not QEInput was read</td>
	</tr>
	<tr>
		<td><pre>s2L</pre></td>
		<td><pre>integer, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>Sj</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>t1</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Start time within subprocess</td>
	</tr>
	<tr>
		<td><pre>t2</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time within subprocess</td>
	</tr>
	<tr>
		<td><pre>temperature</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<pre>initialize()</pre><br/>
			<pre>checkAndUpdateInput()<pre>
		</td>
		<td>Temperature of the system</td>
	</tr>
	<tr>
		<td><pre>tf</pre></td>
		<td><pre>real</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time</td>
	</tr>
	<tr>
		<td><pre>THzToHartree</pre></td>
		<td><pre>real, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><pre>readPhonons()</pre></td>
		<td>Conversion factor from THz to Hartree</td>
	</tr>
	<tr>
		<td><pre>twopi</pre></td>
		<td><pre>real, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>2 times pi</td>
	</tr>
	<tr>
		<td><pre>un</pre></td>
		<td><pre>integer, parameter</pre></td>
		<td><pre>MjModule</pre></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>wby2kT</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
		<td></td>
		<td></td>
		<td></td>
		<td></td>
	</tr>
	<tr>
		<td><pre>x</pre></td>
		<td><pre>real, allocatable</pre></td>
		<td><pre>MjModule</pre></td>
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
