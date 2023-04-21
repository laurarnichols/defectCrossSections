# Export From VASP

The following VASP calculations are needed for a capture calculation:
* relaxed final charge state (plus SCF after)
* relaxed initial charge state (plus SCF after)
* final charge state in initial positions

The results must then be post-processed to be used as input to the other programs in the suite. **The `Export` code is not currently tested on relaxation calculations (the total energy may be extracted incorrectly), but you should always run an SCF calculation after relaxation in VASP to get the charge density correct for the final positions anyways.** The program needs to know where the output files are and where to put the exported data.

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
* `eigenvalues.isp.ik` -- eigenvalues for each band for given spin channel and k-point
* `grid.ik` -- Miller indices for G-vectors such that $G+k <$ cutoff
* `groundState` -- highest occupied band for each spin channel and k-point
* `input` -- main output file with cell, k-point, pseudopotential information, and total energy
* `mgrid` -- full G-vector grid in Miller indices
* `projections.isp.ik` -- $\langle \beta | \psi \rangle$ for given spin channel and k-point
* `projectors.ik` -- $|\beta\rangle$ for given k-point (binary)
* `wfc.isp.ik` -- wave function coefficients for given spin channel and k-point (binary)

## Running and parallelization

To run the code, use something like 
```
aprun -n num-procs path-to-package/bin/Export_VASP.x -nk num-k-pools -nb num-band-groups < export.in > export.out
```
The data is split up as follows: all processors get divided into k-point pools, each pool gets divided into band groups, and the atoms and plane waves get divided over the processes in each band group.

There must be the same number of processes per subgroup, so `num-procs` must be evenly divisible by `num-k-pools` and `num-procs/num-k-pools` must be evenly divisible by `num-band-groups`. The k-points, bands, atoms, and $G+k$ vectors are split up sequentially (i.e., k-points get assigned to pool 1 first then pool 2, etc.). 

Atom parallelization is currently only implemented for outputting the projectors. A single process must handle the I/O for the projectors, so each process in the first band group (only first band group participates here because projectors don't depend on bands) writes out its atoms to its own `projectors.ik.indexInBgrp` file and the root process in the band group merges all of the individual files.

Band parallelization is currently implemented for reading/writing the `wfc.isp.ik` files and the `projections.isp.ik` files. Like with the projectors, the root process in each band group outputs its own file and the root process in the pool merges the individual files. 

The bulk of the Export time is in outputting the `projectors.ik` and `wfc.isp.ik` files. The choice of the number of band groups will need to balance these two calculations. More band groups means the wave function write will go faster, but there will be fewer processes per band group to split up the atoms, so the projector write will go slower. Similarly, the calculations must be repeated for each k-point, so more pools means more loops can be done in parallel, but it also means fewer processors per band group, which will affect the projector writes. Test the different choices to make the best selection for your system.
