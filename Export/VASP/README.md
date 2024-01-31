# Export From VASP

The following VASP calculations are needed for a capture calculation:
* relaxed perfect crystal (plus SCF after)
* relaxed final charge state (plus SCF after)
* relaxed initial charge state (plus SCF after)
* final charge state in initial positions

The results must then be post-processed to be used as input to the other programs in the suite. **The `Export` code is not currently tested on relaxation calculations (the total energy may be extracted incorrectly), but you should always run an SCF calculation after relaxation in VASP to get the charge density correct for the final positions anyways.** The program needs to know where the output files are and where to put the exported data.

The input file should look like
```f90
&inputParams
  ispSelect = integer                                   ! selection of a single spin channel; default unused
  VASPDir = 'path-to-VASP-output'                       ! default './'
  exportDir = 'path-to-put-exported-files'              ! default './Export'
  gammaOnly = logical                                   ! default .false.
  energiesOnly = logical                                ! default .false.
  groupForGroupVelocity = logical                       ! default .false.
  nDispkPerCoord = integer                              ! how many displaced k-points per coordinate if groupForGroupVelocity = .true.
  pattern = 'displacement-pattern-for-group-velocity'   ! e.g., '-0.02 -0.01 0.01 0.02'
/
```
VASP has a gamma-only version that makes simplifications on the PW grid. The Export code is not fully integrated with those assumptions, but it can still be used to export the energy-related values (`energiesOnly`). The `energiesOnly` option is helpful if getting accurate energies from HSE calculations so that the gamma-only version can be used. It is also helpful for group-velocity calculations where multiple k-points are needed for each original k-point to get the derivative of the energy bands with respect to k in each direction. Typically, 5 k-points are used for each k-point and direction (base and two positive and negative) so that you can determine if/where crossings occur and if bands should be fit to lines or parabolas. However, this program allows you to select up to 6 displaced k-points per coordinate (`nDispkPerCoord`). You should also set the displacement pattern used for the k-points. The displacement pattern only applies to the displaced k-points. It is assumed that the base k-point is the first k-point in each group and that the displacement pattern is the same for all coordinate directions. The pattern should line up with the k-point list used in the VASP calculation. The `Export` code will output a `groupedEigenvalues.isp.ikGroup` where `ikGroup` is the index of the k-point group or base k-point. The eigenvalues are output for each band from negative to positive and in x, y, z order.

The output files are
* `eigenvalues.isp.ik` -- eigenvalues for each band for given spin channel and k-point
* `groupedEigenvalues.isp.ik` -- grouped eigenvalues for each base k-point and spin channel if `groupForGroupVelocity = .true.`
* `grid.ik` -- Miller indices for G-vectors such that $G+k < $ cutoff
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
