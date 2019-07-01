# TME Module
## Directory Structure and Files
```
TME
|	README.md
|_______DOC
|	|	TME_Input.in
|_______src
	|	Makefile
	|	TME_Main_v9.f90
	|	TME_Module_v28.f90
```

## How to Run
* Ensure that the executables are up to date by going to the main folder and running `make TME`
* Make sure that you have already [run `Export_QE-5.3.0.x`](../QE-dependent/QE-5.3.0/Export/README.md)
* Run the program and send the contents of the `TME_Input.in` file in as input (e.g., `./bin/TME.x < ExampleRun/TME/input/TME_Input.in`)

## Inputs
* `TME_Input.in`
	* `exportDirSD` (string) -- directs the program to the output directory from `Export` program
	* `exportDirPC` (string) -- directs the program to the output directory from `Export` program
	* `elementsPath` (string) -- path to store outputs
	* `iBandIinit` (integer) -- initial state initial band
	* `iBandIfinal` (integer) -- initial state final band
	* `iBandFinit` (integer) -- final state initial band
	* `iBandFfinal` (integer) -- final state final band
	* `ki` (integer) -- initial k-point
	* `kf` (integer) -- final k-point
	_Note: if `ki` and `kf` are not set, all k-points will be calculated_

	* `calculateVfis` (boolean) -- set to true if there is an incoming electron that will be captured in the system
	* `eBin` (real) -- ?? in eV
	
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

* `exportDirSD`/`exportDirPC` directories

## Variable Glossary

### A
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
		<td><code>absVfi2</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>abortExecution</code></td>
		<td><code>logical</code></td>
		<td><code>checkInitialization</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>checkInitialization</code></td>
		<td>If the program should abort because the inputs weren't all correct</td>
	</tr>
	<tr>
		<td><code>at</code></td>
		<td><code>real, array</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>atoms</code></td>
		<td><code>atom, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>atomsPC</code></td>
		<td><code>atom, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### B
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
		<td><code>bes_J_qr</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>betaPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>betaSD</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>bg</code></td>
		<td><code>real, array</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### C
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
		<td><code>calculateVfis</code></td>
		<td><code>logical</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>coulomb</code></td>
		<td><code>logical</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>count</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>cProjBetaPCPsiSD</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>cProjBetaSDPhiPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>cProjPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>cProjSD</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### D
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
		<td><code>DE</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>displmnt</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>dp</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>transitionMatrixElements</code></td>
		<td>Set real variables to double precision</td>
	</tr>
</table>

### E
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
		<td><code>eBin</code></td>
		<td><code>real</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>eigvF</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>eigvI</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>elementsPath</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>eVToHartree</code></td>
		<td><code>real, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>checkInitialization</code></td>
		<td>Conversion factor from eV to Hartree</td>
	</tr>
	<tr>
		<td><code>exportDirPC</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td>The output directory from <code>Export</code></td>
	</tr>
	<tr>
		<td><code>exportDirSD</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td>The output directory from <code>Export</code></td>
	</tr>
</table>

### F
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
		<td><code>F</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>F1</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>F2</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>fftxMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>fftxMin</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>fftyMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>fftyMin</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>fftzMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>fftzMin</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>file_exists</code></td>
		<td><code>logical</code></td>
		<td>
			<code>readInput</code><br/>
			<code>checkInitialization</code>
		</td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>readInput</code><br/>
			<code>checkInitialization</code>
		</td>
		<td>If file exists</td>
	</tr>
</table>

### G
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
		<td><code>gamma_only</code></td>
		<td><code>logical</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>groundState</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>gvecs</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>gx</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>gy</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>gz</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### H
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
		<td><code>HartreeToEv</code></td>
		<td><code>real, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Conversion factor from Hartree to eV</td>
	</tr>
</table>

### I
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
		<td><code>i</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>iBandFfinal</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>iBandFinit</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>iBandIfinal</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>iBandIinit</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>initialize</code><br/>
			<code>checkInitialization</code>
		</td>
		<td></td>
	</tr>
	<tr>
		<td><code>ibf</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ibi</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>id</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ierr</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>transitionMatrixElements</code></td>
		<td>Error code returned from MPI</td>
	</tr>
	<tr>
		<td><code>ig</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>igN</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `vec` structure</td>
	</tr>
	<tr>
		<td><code>igM</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `vec` structure</td>
	</tr>
	<tr>
		<td><code>igvs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ii</code></td>
		<td><code>complex, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ik</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ind</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `vec` structure</td>
	</tr>
	<tr>
		<td><code>ind2</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>input</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>inputPC</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ios</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInput</code></td>
		<td>Status returned from I/O commands</td>
	</tr>
	<tr>
		<td><code>iostd</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>readInput</code><br/>
			<code>checkInitialization</code>
		</td>
		<td>Unit number for output file</td>
	</tr>
	<tr>
		<td><code>iPn</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>iRc</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>iqs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>iTypes</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### J
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
		<td><code>j</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>JMAX</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### K
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
		<td><code>kf</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>initialize</code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ki</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>initialize</code></td>
		<td></td>
	</tr>
</table>

### L
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
		<td><code>lMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>lmMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>lps</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
</table>

### M
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
		<td><code>master</code></td>
		<td><code>logical</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>maxL</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>mkdir</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>checkInitialization</code></td>
		<td>The command for creating the elements path directory</td>
	</tr>
	<tr>
		<td><code>myid</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>transitionMatrixElements</code></td>
		<td>ID for each MPI process</td>
	</tr>
</table>

### N
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
		<td><code>n</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>n1</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>n2</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>n3</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>n4</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nBands</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>newVecs</code></td>
		<td><code>vec, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nF</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nGf</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nGi</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nGvsF</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nGvsI</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nI</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nIonsPC</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nIonsSD</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nFs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>ngs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nIs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nKpts</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>initialize</code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nKptsPC</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>np</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nPP</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nProjsPC</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nProjsSD</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>npw</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>npwMf</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>npwMi</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>npwNf</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>npwNi</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nPWsF</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code>transitionMatrixElements</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Number of final plane waves??</td>
	</tr>
	<tr>
		<td><code>nPWsI</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code>transitionMatrixElements</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Number of initial plane waves??</td>
	</tr>
	<tr>
		<td><code>npwsPC</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>npwsSD</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nSpins</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>nSquareProcs</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numOfAtoms</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>numOfGvecs</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numOfPWs</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInput</code></td>
		<td>Number of plane waves??</td>
	</tr>
	<tr>
		<td><code>numOfPWsPC</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInput</code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numOfPWsSD</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInput</code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numOfTypes</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numOfTypesPC</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numOfUsedGvecsPP</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>numprocs</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>transitionMatrixElements</code></td>
		<td>Number of processes in the MPI pool</td>
	</tr>
</table>

### O
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
		<td><code>omega</code></td>
		<td><code>real</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>output</code></td>
		<td><code>character, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInput</code></td>
		<td>Name of output file</td>
	</tr>
</table>

### P
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
		<td><code>paw</code></td>
		<td><code>complex</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>paw2</code></td>
		<td><code>complex</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>paw_fi</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>paw_id</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>paw_PsiPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>paw_SDKKPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>paw_SDPhi</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pawKPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pawPsiPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pawSDK</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pawSDPhi</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pi</code></td>
		<td><code>real, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Pi</td>
	</tr>
	<tr>
		<td><code>posIonPC</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>posIonSD</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>psuedo1</code></td>
		<td><code>complex</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>psuedo2</code></td>
		<td><code>complex</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pwGindPC</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pwGindSD</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pwGs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>pwGvecs</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### R
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
		<td><code>r</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>rab</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>root</code></td>
		<td><code>integer, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>transitionMatrixElements</code></td>
		<td>ID of the root process</td>
	</tr>
</table>

### S
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
		<td><code>sq4pi</code></td>
		<td><code>real, parameter</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>symbol</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure that gives the element name</td>
	</tr>
</table>

### T
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
		<td><code>t0</code></td>
		<td><code>real</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td>
			<code>transitionMatrixElements</code><br/>
			<code>readInput</code>
		</td>
		<td>Start time for program</td>
	</tr>
	<tr>
		<td><code>t1</code></td>
		<td><code>real</code></td>
		<td><code>transitionMatrixElements</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>Start time for??</td>
	</tr>
	<tr>
		<td><code>t2</code></td>
		<td><code>real</code></td>
		<td><code>transitionMatrixElements</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td></td>
		<td>End time for??</td>
	</tr>
	<tr>
		<td><code>textDum</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>tf</code></td>
		<td><code>real</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>End time for the program</td>
	</tr>
	<tr>
		<td><code>threej</code></td>
		<td><code>real</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>TME_Input</code></td>
		<td><code>namelist</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>readInput</code></td>
		<td>Used to gather variables for input</td>
	</tr>
	<tr>
		<td><code>tmes_file_exists</code></td>
		<td><code>logical</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>TYPNIPC</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>TYPNISD</code></td>
		<td><code>integer, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### U
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
		<td><code>Ufi</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>

### V
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
		<td><code>vecs</code></td>
		<td><code>vec, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>VfisOutput</code></td>
		<td><code>character</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code>initialize</code></td>
		<td></td>
	</tr>
</table>

### W
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
		<td><code>wae</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
	<tr>
		<td><code>wfcPC</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>wfcSD</code></td>
		<td><code>complex, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>wk</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>wkPC</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>wps</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
	</tr>
</table>

### X
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
		<td><code>xk</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code>xkPC</code></td>
		<td><code>real, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
	<tr>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td><code></code></td>
		<td></td>
	</tr>
</table>
