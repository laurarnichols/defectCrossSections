# Export From VASP

This program exports the data from VASP output files into a form that can be processed by the matrix element code (`TME`). The program needs to know where the output files are and where to put the exported data.

The input file should look like
```
&inputParams
  VASPDir = 'path-to-VASP-output'
  exportDir = 'path-to-put-exported-files'
  gammaOnly = .false.
/
```

The input variables are:
* `VASPDir`
  * Default: `'./'`
  * The directory where the VASP output files are
* `exportDir`
  * Default: `'./Export'`
  * The directory for where to store the exported files
  * If directory does not exist, it will be created
* `gammaOnly`
  * Default: `.false.`
  * Whether the Gamma-point only version of VASP (`gam`) was used

The output files are
* `eigenvalues.ik`
* `grid.ik`
* `input`
* `mgrid`
* `wfc.ik`

## Running and parallelization

To run the code, use something like 
```
aprun -n num-procs path-to-package/bin/Export_VASP.x -nk num-k-pools < export.in > export.out
```
There must be the same number of processes per pool, so `num-procs` must be evenly divisible by `num-k-pools`. The k-points are split up into `num-k-pools` sequentially (i.e., k-points get assigned to pool 1 first then pool 2, etc.). The plane waves are then split across each process in a given pool. The calculations over plane waves are completely distributed, so they are very fast. The bottleneck of the export is writing out the projectors and wave functions, as they require all of the data to be gathered to a single process per pool that handles the I/O for each. The I/O processors are different for the projectors and wave functions so those files can be output at the same time. Because of this bottleneck, it is most efficient to minimize the number of k-points per pool. 

The only way to solve the bottleneck would be to use binary files and direct access and/or MPI I/O, but that would require updating codes down the line, and we aren't concerned with that right now.
