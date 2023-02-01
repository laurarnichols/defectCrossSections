# Export From VASP Output

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

To run the code, use something like `aprun -n num-procs path-to-package/bin/Export_FromVASPOutput.x < exportFromOutput.in > exportFromOutput.out`. 

Currently, the I/O is only done by the root node (`ionode`), but the processing of the data is split across all processes.

The output files are
* `eigenvalues.ik`
* `grid.ik`
* `input`
* `mgrid`
* `wfc.ik`

