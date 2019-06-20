# Export Program
## Directory Structure and Files
```
QE-dependent/QE-5.3.0/Export
|	README.md
|_______DOC
|	|	export.in
|_______src
	|	Export_QE-5.3.0_v3.f90
	|	Makefile
```

## How to Run
* Ensure that the executables are up to date by going to the main folder and running `make Export_QE-5.3.0`
* Make sure you have already done a QE scf calculation (see [RunningQE.md](RunningQE.md))
* Run the program and send the contents of the `export.in` file in as input (e.g., `./bin/Export_QE-5.3.0.x < ExampleRun/Export/input/export.in`)

## Inputs
* `export.in`
	* `prefix` (string) -- the prefix for all output files (must match what was used for QE scf calculation)
	* `outdir` (string) -- the path to the QE output folder
	* `exportDir` (string) -- the path to the output folder for the `Export_QE-5.3.0.x` program
	* `writeWFC` (boolean) -- ??
	
_Note: Do not alter the `&inputpp` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_
