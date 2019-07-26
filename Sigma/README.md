# Sigma
## Directory Structure and Files
```
Sigma
|	README.md
|-------DOC
|	| input.in
|-------src
	|	Makefile
	|	Sigma_Main.f90
	|	Sigma_Module_v4.f90
```

## How to Run
* Ensure that the executables are up to date by going to the main folder and running `make TME`
* Make sure that you have already run [`TME.x`](../TME/README.md) and [`LSF.x`](../LSF/README.md)
* Run the program and send the contents of the `input.in` file in as input (e.g., `./bin/Sigma.x < ExampleRun/Sigma/input/input.in` from the home folder)

## Inputs
* `input.in`
	* `VfisInput` (string) -- directs the program to the output file containing the transition matrix elements (from `TME` program)
	* `LSFInput` (string) -- directs the program to the output file containing the line shape function (from `LSF` program)
	
_Note: Do not alter the `&elphscat` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## Outputs 
Sigma as a function of energy
