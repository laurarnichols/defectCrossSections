# Export From VASP Source

This program uses the VASP source code to generate the wave function, projectors, and projections to be used by the matrix element code (`TME`). It is currently not integrated with VASP's parallelization scheme, so *you must set `NCORE = number-of-cores` and `KPAR = 1`*. Otherwise, the program will not work. 

Currently, the input files are read from the current directory automatically (like when running VASP) and the output is written to the same directory. I am not 100% sure on what input files are needed. The `OUTCAR` file is overwritten, so make sure to save that file if you need it somewhere else. *You should run the `ExportFromVASPOutput` first because some of the output may change after you run `ExportFromVASPSrc`.*
