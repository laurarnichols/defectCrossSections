# Export From VASP Source

This program uses the VASP source code to generate the projectors and projections to be used by the matrix element code (`TME`). It is currently not integrated with VASP's parallelization scheme, so **you must use 1 MPI process and set `NCORE = 1` and `KPAR = 1`**. Otherwise, the program will not work. 

The following options can be passed in at the command line:
* `-ks` or `-kStart` is the k-point you want to start with
* `-ke` or `-kEnd` is the k-point you want to end with
* `-id` or `-inputDirectory` is the path to the VASP files to be exported
* `-od` or `-outputDirectory` is the path for the output files detailing how the export proceeds (all typical VASP output files)
* `-ed` or `-exportDirectory` is the path for the exported files (`projectors.ik`, `projections.ik`, and `wfc.ik`)

An example of how to run the code is
```
aprun -n 1 path-to-package/bin/Export_FromVASPSrc.x -ks 1 -ke 2 -id '../' -od './trashOutput' -ed './'
```

I believe that the only files this calculation needs from the VASP calculation are `INCAR`, `KPOINTS`, `POSCAR`, `POTCAR`, `WAVECAR`, and `CHGCAR`.

For compiling the code, make a copy of the VASP 5 source directory and update the path in the main `Makefile`. Go to that copy and run `make std`, then compile the `Export_FromVASPSrc` code. The local `Makefile` will copy over the customized files to the VASP-source copy and compile.
