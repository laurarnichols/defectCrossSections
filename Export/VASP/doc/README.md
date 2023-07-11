# Export From VASP

The following VASP calculations are needed for a capture calculation:
* relaxed final charge state (plus SCF after)
* relaxed initial charge state (plus SCF after)
* final charge state in initial positions

The results must then be post-processed to be used as input to the other programs in the suite. **The `Export` code is not currently tested on relaxation calculations (the total energy may be extracted incorrectly), but you should always run an SCF calculation after relaxation in VASP to get the charge density correct for the final positions anyways.** The program needs to know where the output files are and where to put the exported data.

The input file should look like
```f90
&inputParams
  VASPDir = 'path-to-VASP-output'               ! default './'
  exportDir = 'path-to-put-exported-files'      ! default './Export'
  gammaOnly = logical                 ! default .false.
  energiesOnly = logical              ! default .false.
  groupForGroupVelocity = logical     ! default .false.
/
```
VASP has a gamma-only version that makes simplifications on the PW grid. The Export code is not fully integrated with those assumptions, but it can still be used to export the energy-related values (`energiesOnly`). The `energiesOnly` option is helpful if getting accurate energies from HSE calculations so that the gamma-only version can be used. It is also helpful for group-velocity calculations where multiple k-points are needed for each original k-point to get the derivative of the energy bands with respect to k in each direction. Both positive and negative displacements are needed in each direction to determine if/where band crossings occur. The `groupForGroupVelocity` option assumes that the input k-points are in the order base, +/-x, +/-y, and +/-z. The `Export` code will output a `groupedEigenvalues.isp.ik` where `ik` is the index of the base k-point. The eigenvalues are output for each band in the same order as the k-points. It is assumed that there are 7 k-points per group.

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
```bash
aprun -n num-procs path-to-package/bin/Export_VASP.x -nk num-k-pools -nb num-band-groups < export.in > export.out
```
The data is split up as follows: all processors get divided into k-point pools, each pool gets divided into band groups, and the atoms and plane waves get divided over the processes in each band group.

There must be the same number of processes per subgroup, so `num-procs` must be evenly divisible by `num-k-pools` and `num-procs/num-k-pools` must be evenly divisible by `num-band-groups`. The k-points, bands, atoms, and $G+k$ vectors are split up sequentially (i.e., k-points get assigned to pool 1 first then pool 2, etc.). 

Atom parallelization is currently only implemented for outputting the projectors. Band parallelization is currently implemented for reading/writing the `wfc.isp.ik` files and the `projections.isp.ik` files. The projectors and wave functions are written out to binary files so that the I/O can be done in parallel.
