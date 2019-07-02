# TME Module
## Directory Structure and Files
```
TME
|	README.md
|-------DOC
|	|	TME_Input.in
|-------src
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
