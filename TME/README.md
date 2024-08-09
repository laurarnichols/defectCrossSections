# Transition matrix element (`TME`)

The `TME` program primarily calculates the all-electron overlap between certain bands in two different systems, `<braWfc|ketWfc>`, based on the PAW method. In the `overlapOnly` mode, only the overlaps are output. Otherwise, the matrix element based on the given `order` is calculated. The matrix elements require the bands to be read from the `energyTable` created by the `EnergyTabulator` code, but the `overlapOnly = .true.` mode allows the specification of bands without using the energy table.

For each calculation, the `braExportDir` and `ketExportDir` must be specified. If the states are to be read from the energy table, `energyTableDir` must be given. For the `overlapOnly` mode, it is also possible to define all of the band pairs explicitly using, e.g.
```f90
nPairs = 5
braBands = '128 129 130 128 129'
ketBands = '128 128 128 129 129'
```
The list of bands must be given in a string, then `nPairs` of bands will be read from the strings. You can also give a range of bands for the bra and ket systems using `iBandLBra`, `iBandHBra`, `iBandLKet`, and `iBandHKet`, where `L` and `H` represent the lower and upper limits on the bands, respectively.

There are a few optional parameters the user can set:
* `capture` -- logical, default `.true.`. For capture, it is assumed that there is a single possible final state and that the energy table was also tabulated with `captured = .true.`. The scattering code has been generalized to allow overlaps between k-points, but the `Export` code would need to be updated because of the way the wave functions are represented in the finite supercell. The overlap between k-points has, therefore, not been validated.
* `intraK` -- logical, default `.false.`. Right now, this only removes k-point parallelization and affects the expected format of the energy table.
* `dqOnly` -- logical, default `.false.`. Allows the option to only divide the first-order overlap by $\delta q_j$ rather than multiplying by the energy. This might be helpful when using the same overlaps to calculate matrix elements for different problems, like equilibrium and nonequilibrium capture. However, the adjustment to the energies can also be done in the `LSF` program.

Output files are 
* `allElecOverlap.isp.ik`/`allElecOverlap.isp` -- for each band in range and a given spin (and k-point for capture or inter-k overlaps): initial and final band, all-electron overlaps, and matrix elements (if `.not. overlapOnly`)
* stdout -- timing information and status updates

The same TME code is used for both the zeroth-order and the first-order matrix elements. For both orders, the overlap part comes from from the `TME` code and the energy difference comes from the [`EnergyTabulator`](../EnergyTabulator) code. The code assumes that the order of the atom types in the two systems matches, with the possible addition of atom types at the end of one of the files.

## Zeroth-order

The zeroth-order term has a single matrix element that is shown in the GaN paper to be $$M_{\text{e}}^{\text{BO}} = \langle \phi_f | \psi_i^0 \rangle (E_f - E_i),$$ where $E_f - E_i$ is the total electronic energy difference of the defect crystal before and after capture and $|\phi_f\rangle$ and $|\psi_i^0\rangle$ are the single-quasiparticle orbitals of the final defect state and the perfect-crystal initial state, respectively. 

The `TME` input file for the zeroth-order should look like
```f90
&TME_Input
  ! Which order of matrix element to calculate
  order = 0

  ! Spin channel selection; default unused
  ispSelect = integer  
  
  ! Systems used for overlaps
  braExportDir    = 'path-to-final-charge-state-initial-positions-export'
  ketExportDir    = 'path-to-perfect-crystal-export'
  energyTableDir = 'path-to-energy-tables'
  
  ! Output info
  outputDir = 'path-to-store-overlap-files' 			! default './TMEs'
/
```
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## First-order

The first-order term has a matrix element for each mode: $$M_j = \frac{\varepsilon_i - \varepsilon_f}{\delta q_j} \langle \phi_f^j | \phi_i\rangle,$$ where $\varepsilon_i - \varepsilon_f$ is the eigenvalue energy difference between the initial and final band states and $|\phi_f^j\rangle$ is the final-state defect orbital after displacing the atoms along the phonon directions by a small $\delta q_j$. The formalism assumes that the phonon modes are the same in the initial and final electronic states, so it doesn't matter which electronic state is used as long as the same charge state is used for both input wave functions and the atoms are in the initial positions. In practice, the ground state is usually simplest and fastest.

The equation above for the matrix element assumes that the wave functions are orthogonal, but the numerical overlaps of the non-displaced wave functions $\langle \phi_{l'}|\phi_l\rangle$ are not actually $\delta_{ll'}$. To increase the numerical accuracy of the matrix elements, it is best to subtract the non-displaced wave function overlaps: $$M_j = \frac{\varepsilon_i - \varepsilon_f}{\delta q_j} \left( \langle \phi_f^j | \phi_i\rangle - \langle \phi_f|\phi_i\rangle \right).$$
These overlaps can be found by inputting the non-displaced system for the bra and ket inputs. These baseline overlaps can then be subtracted from the overlaps from all of the mode-shifted overlaps.

The `TME` input file for the first-order term should look like
```f90
&TME_Input
  ! Which order of matrix element to calculate
  order = 1

  ! Spin channel selection; default unused
  ispSelect = integer                                   
  
  ! Systems used for overlaps
  braExportDir    = 'path-to-displaced-defect-export'
  ketExportDir    = 'path-to-undisplaced-defect-export'
  energyTableDir = 'path-to-energy-tables'

  ! Baseline info:
  subtractBaseline = .true. or .false. ! If "orthogonal" overlap should be subtracted
  baselineDir = 'path-to-non-displaced-overlaps'

  ! Parameters to get dq_j
  dqFName = 'path-to-dq-file-from-shifter'
  phononModeJ = integer           ! phonon-mode index
  
  ! Output info
  outputDir = 'path-to-store-overlap-files' 			! default './TMEs'
/
```
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## Running and parallelization

To run the code, use something like 
```bash
aprun -n num-procs path-to-package/bin/TME.x -nk num-k-pools < TME.in > TME.out
```

All processors get divided into k-point pools, then plane waves get divided over the processes in each pool. There must be the same number of processes per pool, so `num-procs` must be evenly divisible by `num-k-pools`. The k-points, $G$ vectors, and $G+k$ vectors are split up sequentially (i.e., k-points get assigned to pool 1 first then pool 2, etc.). 
