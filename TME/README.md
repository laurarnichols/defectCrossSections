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
		<td><code>ind</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `vec` structure</td>
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
		<td><code>newVecs</code></td>
		<td><code>vec, allocatable</code></td>
		<td><code>declarations</code></td>
		<td><code></code></td>
		<td><code></code></td>
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
		<td><code>nMax</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
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
		<td><code>numOfAtoms</code></td>
		<td><code>integer</code></td>
		<td><code>declarations</code></td>
		<td>N/A</td>
		<td>N/A</td>
		<td><code></code></td>
		<td>Part of the `atom` structure</td>
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
